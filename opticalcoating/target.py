import json
from importlib.resources import files
from datetime import datetime
from .calc_flux import Wave

class Target:
    def __init__(self, name, waves, flux_target, creation_time=None, q_TR='T', dTa=None):
        """Target file
        :param waves: tuple of Wave class objects
        :param flux_target: absolute value, i.e. range 0.0 - 1.0
        """
        self.name = name
        self.waves = waves
        self.q_TR = q_TR
        self.flux_target = flux_target
        self.creation_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S") if creation_time is None else creation_time
        self.dTa = dTa


    @classmethod
    def from_json(cls, name):
        fname = files(f'opticalcoating.resources.Targets').joinpath(f'Target_{name}.json')
        with open(fname, 'r', encoding='utf-8') as file:
            dct = json.load(file)
        return cls(name=name,
                   creation_time = dct['creation_time'],
                   waves=[Wave.from_tuple(wv_tuple)
                          for wv_tuple in zip(dct['wv_list'], dct.get('pol_list'), dct.get('angle_list'))],
                   q_TR = dct['q_TR'],
                   flux_target = dct['flux_target'],
                   dTa = dct.get('dTa'))


    def save(self):
        dct = {'name': self.name,
               'creation_time': self.creation_time,
               'q_TR': self.q_TR,
               'wv_list': [wave.wavelength for wave in self.waves],
               'pol_list': [wave.polarisation for wave in self.waves],
               'angle_list': [wave.angle for wave in self.waves],
               'flux_target': self.flux_target}

        with open(f"../../opticalcoating/resources/Targets/Target_{self.name}.json",
                  'w', encoding='utf-8') as file:
            json.dump(dct, file, indent=3)