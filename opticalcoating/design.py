from .units import IndexFormula
import numpy as np
import scipy

class Material:
    def __init__(self, data):
        self.data = data
        self.name = list(data.keys())[0]
        ri_formula = list(data[self.name].keys())[0]
        ri_data = data[self.name][ri_formula]
        self.formula = IndexFormula(ri_formula)

        if self.formula is IndexFormula.Const:
            self.bnd = (0, np.inf)
            self._n = lambda wavelength: ri_data['n']
            self._xi = lambda wavelength: ri_data.get('xi', 0.)

        elif self.formula is IndexFormula.Linear_interp:
            self._n = scipy.interpolate.make_interp_spline(ri_data['wavelength'], ri_data['n'], k=1)
            self._xi = scipy.interpolate.make_interp_spline(ri_data['wavelength'], ri_data.get('xi', 0.), k=1)

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


class Layer:
    def __init__(self):
        self._material
        self._thickness
        self._witness

class Design:
    def __init__(self):
        pass