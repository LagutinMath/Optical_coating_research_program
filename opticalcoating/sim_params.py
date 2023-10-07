from time import time
import random
from .units import Wave
from .design_class import Design
from .deposition_simulation import SetUpParameters


class SimParams:
    def __init__(self, info):
        self.des = Design(name=info['design'])
        self.trg = info.get('target')
        self.term_algs = SimParams.expand(info['term_alg'], self.des.N)
        self.start_seed = random.randint(1, 10**16)

        wv_list = SimParams.expand(info['wave'], self.des.N)
        polarisation_list = SimParams.expand(info['polarisation'], self.des.N)
        angle_list = SimParams.expand(info['angle'], self.des.N)

        waves = [Wave(wv, pol, ang) for wv, pol, ang in zip(wv_list, polarisation_list, angle_list)]
        rates = SimParams.expand(info['rate'], self.des.N)

        if not isinstance(info['std rate'], list):
            rates_sigmas = [None] + [0.01 * info['std rate'] * rate for rate in rates[1:]]
        else:
            temp = [None] + (self.des.N // len(info['std rate']) + 1) * info['std rate']
            rates_sigmas = [None] + [0.01 * std_rate * rate for rate, std_rate in zip(rates[1:], temp[1:])]

        if not isinstance(info['std meas'], list):
            meas_sigmas = [None] + [0.01 * info['std meas'] for _ in range(self.des.N)]
        else:
            temp = [None] + (self.des.N // len(info['std meas']) + 1) * info['std meas']
            meas_sigmas = [None] + [0.01 * std_meas for std_meas in temp[1:]]

        self.set_up_pars = SetUpParameters(N=self.des.N,
                                           waves=waves,
                                           rates=rates,
                                           q_TR=info['T/R'],
                                           backside=info['backside'],
                                           meas_sigmas=meas_sigmas,
                                           rates_sigmas=rates_sigmas,
                                           width=info.get('width', None))


    def get_tasks(self, N_sim):
        tasks = [{'des_th': self.des,
                  'term_algs': self.term_algs,
                  'set_up_pars': self.set_up_pars,
                  'rnd_seed': rnd_seed}
                 for rnd_seed in range(self.start_seed, self.start_seed + N_sim)]
        return tasks


    @classmethod
    def from_dict(cls, dct_params, default=None):
        if default is None:
            default = {}
            default['term_alg'] = 'QS'
            default['rate'] = (0.5, 0.8)  # nm / sec
            default['std rate'] = (5, 10)  # in percent
            default['std meas'] = 0.1  # in percent
            default['T/R'] = 'T'
            default['polarisation'] = 'S'
            default['angle'] = 0
            default['backside'] = False
        dct = default | dct_params
        return cls(dct)


    @staticmethod
    def expand(x, N):
        if not isinstance(x, list): return [None] + N * [x]
        temp = [None] + (N // len(x) + 1) * x
        return temp[:N + 1]