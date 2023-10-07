import scipy
import json

import numpy as np

from importlib.resources import files

from .units import IndexFormula, Role

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

            @property
            def d(self):
                return self._d_cur

            @property
            def d_th(self):
                return self._d_th

            @property
            def role(self):
                return self._role

        self.name = data['name']
        materials = {mat_name: Material(mat_name, data['mat_inf'][mat_name]) for mat_name in data['mat_inf']}

        self.subs = Substrate(materials[data['layers'][0]], data['thicknesses'][0])
        self.layers = [Layer(materials[mat_name], d)
                       for mat_name, d in zip(data['layers'][1:], data['thicknesses'][1:])]


    @classmethod
    def from_name(cls, name):
        fname = files('opticalcoating.resources.Designs').joinpath(f'{name}.json')
        with open(fname, 'r', encoding='utf-8') as file:
            data = json.load(file)
            return cls(data)
        raise FileNotFoundError(f'Design {name} is not found')