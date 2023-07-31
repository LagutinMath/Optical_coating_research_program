import json
import numpy as np
from copy import deepcopy
from importlib.resources import files
from .calc_flux import calc_flux
from .statistics_info import StatInfo


def merit(des, target, error=None, MF_d_th=0.0):
    if error is not None:
        des = deepcopy(des)
        des.d = [des.d[0]] + [des.d[i] + error[i - 1] for i in range(1, des.N + 1)]

    waves, T_target = np.array(target[0]), np.array(target[1])
    dTa = np.array(target[2]) if len(target) > 2 else None
    T_des = np.vectorize(lambda x: calc_flux(des, x, q_subs=False, q_percent=False, q_TR='T'))(waves)

    if dTa is None: delta_T = T_des - T_target
    else: delta_T = (T_des - T_target) / dTa
    return np.sqrt(np.dot(delta_T, delta_T) / len(waves)) - MF_d_th


class ProcessedStatistics:
    def __init__(self, info):
        """:param: M - число анализируемых симуляций"""
        self.statistic_num = info['statistic_num']
        self.c_value = info['c_value']
        self.mean_delta_MF_rnd = info['mean_delta_MF_rnd']
        self.c_array = np.array(info['c_array'])
        self.MF_d_th = info['MF_d_th']


    @classmethod
    def calc(cls, statistic_num, des, target, M=None):
        info = {}
        info['statistic_num'] = statistic_num
        info['MF_d_th'] = merit(des, target)
        errors = np.array(StatInfo.load(statistic_num).error_list)
        if M is None or M > len(errors): M = len(errors)
        errors = errors[:M]

        median_val = np.median(np.linalg.norm(errors, axis=1))
        rnd_errors = np.random.normal(0., median_val / np.sqrt(des.N), size=(M, des.N))
        delta_MF_rnd = np.vectorize(lambda x: merit(des, target, error=x, MF_d_th=info['MF_d_th']),
                                    signature='(n)->()')(rnd_errors)
        info['mean_delta_MF_rnd'] = np.mean(delta_MF_rnd)

        delta_MF = np.vectorize(lambda x: merit(des, target, error=x, MF_d_th=info['MF_d_th']),
                                signature='(n)->()')(errors)
        info['c_array'] = delta_MF / info['mean_delta_MF_rnd']
        info['c_value'] = np.mean(info['c_array'])
        return cls(info)


    @classmethod
    def load(cls, statistic_num):
        fname = files('opticalcoating.resources.c_values').joinpath(
            f'c_value{str(statistic_num).zfill(3)}.json')
        with open(fname, 'r') as file: info = json.load(file)
        info['statistic_num'] = statistic_num
        return cls(info)


    def save(self):
        info = {'c_value': self.c_value,
                'MF_d_th': self.MF_d_th,
                'mean_delta_MF_rnd': self.mean_delta_MF_rnd,
                'c_array': self.c_array.tolist()}

        fname = files('opticalcoating.resources.c_values').joinpath(
            f'c_value{str(self.statistic_num).zfill(3)}.json')
        with open(fname, 'w') as file:
            json.dump(info, file, indent=3)
        print(f'"c_value{str(self.statistic_num).zfill(3)}.json" is successfully saved')