from functools import singledispatchmethod, total_ordering
from enum import Enum
from math import pi

class Role(Enum):
    # Material (layer) roles
    Substrate, S, Subs, subs  = 1, 1, 1, 1
    H, L, A, B, C = 2, 3, 4, 5, 6

class IndexFormula(Enum):
    Const = 1
    Linear_interp = 2
    Sellmeier = 3
    Cauchy = 4

    @staticmethod
    def from_str(name):
        if name in 'Const': return IndexFormula.Const
        if name in 'Table': return IndexFormula.Linear_interp
        if name in 'Sellmeier': return IndexFormula.Sellmeier
        if name in 'Cauchy': return IndexFormula.Cauchy
        raise ValueError(f'{name} is absent in IndexFormula')

class TermCase(Enum):
    OnTime = 1
    LateStop = 2
    EmergencyStop = 3

class WidthForm(Enum):
    rect, gauss = 1, 2
    Rect, Gauss = 1, 2
    r, g = 1, 2
    R, G = 1, 2

class Pol(Enum):
    S, P = 1, 2
    s, p = 1, 2

    @staticmethod
    def from_str(s):
        if s in ('S', 's'): return Pol.S
        if s in ('P', 'p'): return Pol.P

@total_ordering
class Wave:
    __slots__ = ('_wavelength', '_polarisation', '_angle')

    @singledispatchmethod
    def __init__(self, wavelength, polarisation=Pol.S, angle=0.):
        """Incident wave class
        :param wavelength: wavelength in nm
        :param polarisation: 'S' or 'P' polarisation
        :param angle: angle in radians
        """
        self._wavelength = wavelength
        self._polarisation = polarisation
        self._angle = angle

    @__init__.register(list)
    @__init__.register(tuple)
    def _from_list_tuple(self, wave_tuple):
        self._wavelength = wave_tuple[0]
        self._polarisation = wave_tuple[1] if len(wave_tuple) > 1 and wave_tuple[1] is not None else Pol.S
        self._angle = wave_tuple[2] if len(wave_tuple) > 2 and wave_tuple[2] is not None else 0

    @property
    def wavelength(self):
        return self._wavelength

    @property
    def polarisation(self):
        return self._polarisation

    @property
    def angle(self):
        return self._angle

    def __call__(self, wavelength):
        return Wave(wavelength, self.polarisation, self.angle)

    def __float__(self):
        return float(self.wavelength)

    def __repr__(self):
        return f"Wave({self.wavelength}, '{self.polarisation}', {self.angle})"

    def __eq__(self, other):
        if isinstance(other, Wave):
            return (self.wavelength == other.wavelength and
                    self.polarisation is other.polarisation and
                    self.angle == other.angle)
        return NotImplemented

    def __hash__(self):
        return hash((self.wavelength, self.polarisation, self.angle))

    def __lt__(self, other):
        if isinstance(other, Wave):
            if self.angle > other.angle: return True
            if self.polarisation is not other.polarisation and other.polarisation is Pol.P: return True
            return self.wavelength < other.wavelength
        return NotImplemented

    def __add__(self, other):
        if isinstance(1., (int, float)):
            return Wave(self.wavelength + other, self.polarisation, self.angle)
        elif isinstance(other, Wave):
            if (self.polarisation is other.polarisation and self.angle == other.angle):
                return Wave(self.wavelength + other.wavelength, self.polarisation, self.angle)
            else:
                return NotImplemented
        else:
            return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(1., (int, float)):
            return Wave(self.wavelength - other, self.polarisation, self.angle)
        elif isinstance(other, Wave):
            if (self.polarisation is other.polarisation and self.angle == other.angle):
                return Wave(self.wavelength - other.wavelength, self.polarisation, self.angle)
            else:
                return NotImplemented
        else:
            return NotImplemented

    def __rsub__(self, other):
        if isinstance(1., (int, float)):
            return Wave(other - self.wavelength, self.polarisation, self.angle)
        elif isinstance(other, Wave):
            if (self.polarisation is other.polarisation and self.angle == other.angle):
                return Wave(other.wavelength - self.wavelength, self.polarisation, self.angle)
            else:
                return NotImplemented
        else:
            return NotImplemented

    @staticmethod
    def deg2rad(deg): return deg * pi / 180.

    @staticmethod
    def rad2deg(rad): return rad * 180. / pi
