import scipy
import json
from time import perf_counter

import numpy as np

# from functools import singledispatchmethod
from importlib.resources import files

from .units import IndexFormula, Role, Pol

class Material:
    def __init__(self, mat_name, data):
        self.data = data
        self.name = mat_name
        ri_formula = list(data.keys())[0]
        ri_data = data[ri_formula]
        self.formula = IndexFormula.from_str(ri_formula)

        if self.formula is IndexFormula.Const:
            self.bnd = (0, np.inf)
            self._n = lambda wavelength: ri_data['n']
            self._xi = lambda wavelength: ri_data.get('xi', 0.)

        elif self.formula is IndexFormula.Linear_interp:
            self._n = scipy.interpolate.make_interp_spline(ri_data['wavelength'], ri_data['n'], k=1)
            if 'xi' in ri_data:
                self._xi = scipy.interpolate.make_interp_spline(ri_data['wavelength'], ri_data['xi'], k=1)
            else:
                self._xi = lambda wavelength: 0.

        elif self.formula is IndexFormula.Sellmeier:
            coef_A = ri_data['coef_A']
            def Sellmeier(wavelength):
                # wavelength in nm
                wv2 = wavelength ** 2
                pwr2n = (coef_A[0] + wv2 * (coef_A[1] / (wv2 - coef_A[2]) +
                                            coef_A[3] / (wv2 - coef_A[4]) +
                                            coef_A[5] / (wv2 - coef_A[6])))
                return np.sqrt(pwr2n)
            self._n = Sellmeier
            self._xi = lambda wavelength: 0.

        elif self.formula is IndexFormula.Cauchy:
            coef_A = ri_data['coef_A']
            def Cauchy(wavelength):
                # wavelength in nm
                return coef_A[0] + coef_A[1] / wavelength ** 2 + coef_A[2] / wavelength ** 4
            self._n = Cauchy
            self._xi = lambda wavelength: 0.

        else:
            raise ValueError('Uncorrect Index Formula')

    def n(self, wave):
        return self._n(float(wave))

    def xi(self, wave):
        return self._xi(float(wave))


class Design:
    n_a = 1.
    def __init__(self, data):
        class Substrate:
            self._role = Role.Substrate
            def __init__(self, mat, d=None):
                self._material = mat
                self._d = None if d == 0 else d

            @property
            def d(self):
                if self._d is None: return 0
                return self._d

            def n(self, wave):
                return self._material.n(wave)

            def xi(self, wave):
                return self._material.xi(wave)

        class Layer:
            count = 0
            __slots__ = ('_material', '_d_th', '_d_cur', '_witness', '_role', '_number')
            def __init__(self, mat, d, role=None, witness=1):
                Layer.count += 1
                self._material = mat
                self._d_th = d
                self._d_cur = d
                self._witness = witness
                if role is None:
                    self._role = Role.H if Layer.count % 2 else Role.L
                else:
                    self._role = role
                self._number = Layer.count

            def n(self, wave):
                return self._material.n(wave)

            @property
            def d(self):
                return self._d_cur

            @property
            def d_th(self):
                return self._d_th

            @property
            def role(self):
                return self._role

        self._d = None
        self._n = {}
        self._n_s = {}

        self.name = data['name']
        materials = {mat_name: Material(mat_name, data['mat_inf'][mat_name]) for mat_name in data['mat_inf']}

        self.subs = Substrate(materials[data['layers'][0]], data['thicknesses'][0])
        self.layers = [Layer(materials[mat_name], d)
                       for mat_name, d in zip(data['layers'][1:], data['thicknesses'][1:])]


    @property
    def N(self):
        return len(self.layers)

    @property
    def d(self):
        if self._d is None:
            self._d = np.array([layer.d for layer in self.layers])
        return self._d

    def n(self, wave):
        if wave in self._n: return self._n[wave]
        self._n[wave] = np.array([layer.n(wave) for layer in self.layers])
        return self._n[wave]

    def n_s(self, wave):
        if wave in self._n_s: return self._n_s[wave]
        self._n_s[wave] = self.subs.n(wave)
        return self._n_s[wave]

    @classmethod
    def from_name(cls, name):
        fname = files('opticalcoating.resources.Designs').joinpath(f'{name}.json')
        with open(fname, 'r', encoding='utf-8') as file:
            data = json.load(file)
            return cls(data)
        raise FileNotFoundError(f'Design {name} is not found')

    def __call__(self, wave, q_subs=False, backside=False):
        time = []
        time.append(perf_counter())
        def phase(d, n):
            if wave.angle == 0:
                phi = (2.0 * np.pi / wave.wavelength) * d * n
                q = n
            else:
                gamma = np.arcsin(np.sin(wave.angle) * Design.n_a / n)
                phi = (2.0 * np.pi / wave.wavelength) * d * n * gamma
                if wave.polarisation is Pol.S:
                    q = n * np.cos(gamma)
                elif wave.polarisation is Pol.P:
                    q = n / np.cos(gamma)
                else:
                    raise ValueError(f"{wave.polarisation} is unknown polarisation")
            return phi, q

        time.append(perf_counter())
        ph = phase(self.d, self.n(wave))
        time.append(perf_counter())
        M = np.linalg.multi_dot(np.vectorize(lambda phi, q: np.array([[np.cos(phi), 1j * np.sin(phi) / q],
                                                                      [1j * np.sin(phi) * q, np.cos(phi)]],
                                                                     dtype=np.complex128),
                                             signature='(),()->(2,2)')(*ph))
        time.append(perf_counter())
        n_s = self.n_s(wave)
        if wave.angle == 0:
            q_s = n_s
            q_a = Design.n_a
        else:
            gamma_s = np.arcsin(np.sin(wave.angle) * Design.n_a / n_s)
            if wave.polarisation is Pol.S:
                q_s = n_s * np.cos(gamma_s)
                q_a = Design.n_a * np.cos(wave.angle)
            elif wave.polarisation is Pol.P:
                q_s = n_s / np.cos(gamma_s)
                q_a = Design.n_a / np.cos(wave.angle)
            else:
                raise ValueError(f"{wave.polarisation} is unknown polarisation")
        time.append(perf_counter())
        t = 2 * q_a / (M[0][0] * q_a + M[1][1] * q_s + M[0][1] * q_s * q_a + M[1][0])
        r = ((M[0][0] * q_a - M[1][1] * q_s + M[0][1] * q_s * q_a - M[1][0]) /
             (M[0][0] * q_a + M[1][1] * q_s + M[0][1] * q_s * q_a + M[1][0]))
        time.append(perf_counter())
        T_1 = abs(t)**2 * q_s / q_a
        R_1 = abs(r)**2
        if q_subs:
            T_2 = 4.0 * q_a * q_s / (q_a + q_s) ** 2
            R_2 = ((q_a - q_s) / (q_a + q_s)) ** 2
            if backside:
                T_1, T_2 = T_2, T_1
                R_1, R_2 = R_2, R_1
            # формула для коэффициента поглощения верна для нормального падения
            sgm = np.exp(-4 * np.pi * self.subs.xi(wave) * self.subs.d / wave.wavelength)
            T = sgm * T_1 * T_2 / (1.0 - R_1 * R_2 * sgm ** 2)
            R = (R_1 + (sgm ** 2 * R_2 * T_1 * T_1) / (1.0 - R_1 * R_2 * sgm ** 2))
        else:
            T, R = T_1, R_1
        time.append(perf_counter())
        print(*[end - start for end, start in zip(time[1:], time[:-1])])
        return T, R
