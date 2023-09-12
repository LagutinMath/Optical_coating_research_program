import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import seaborn as sns
from datetime import datetime
from importlib.resources import files
from .save_data import find_file_name
import opticalcoating.visualisation as vis


class StatInfo:
    def __init__(self, stat_dict):
        self.des_name = stat_dict['design']
        self.trg_name = stat_dict.get('target', stat_dict['design'])
        self.creation_time = stat_dict['creation_time']
        self.start_rnd_seed = stat_dict['start_rnd_seed']
        self.waves = stat_dict['waves']
        self.term_algs = [None] + stat_dict['term_algs']
        self.rates = [None] + stat_dict['rates']
        self.rates_sigmas = [None] + stat_dict['rates_sigmas']
        self.meas_sigmas = [None] + stat_dict['meas_sigmas']
        self.error_list = stat_dict['error list']


    @classmethod
    def legacy(cls, des_th, term_algs, set_up_pars, error_list, start_rnd_seed, target_name=None):
        stat_dict = {'design': des_th.name,
                     'target': target_name if target_name is not None else des_th.name,
                     'creation_time': datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                     'start_rnd_seed': start_rnd_seed,
                     'waves': [wave.wavelength for wave in set_up_pars.waves[1:]],
                     'term_algs': term_algs[1:],
                     'rates': set_up_pars.rates[1:],
                     'rates_sigmas': set_up_pars.rates_sigmas[1:],
                     'meas_sigmas': set_up_pars.meas_sigmas[1:],
                     'error list': error_list}
        return cls(stat_dict)


    @classmethod
    def from_sim_param(cls, sim_param, err_list):
        stat_dict = {'design': sim_param.des.name,
                     'target': sim_param.trg if sim_param.trg is not None else sim_param.des.name,
                     'creation_time': datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
                     'start_rnd_seed': sim_param.start_seed,
                     'waves': [wave.wavelength for wave in sim_param.set_up_pars.waves[1:]],
                     'term_algs': sim_param.term_algs[1:],
                     'rates': sim_param.set_up_pars.rates[1:],
                     'rates_sigmas': sim_param.set_up_pars.rates_sigmas[1:],
                     'meas_sigmas': sim_param.set_up_pars.meas_sigmas[1:],
                     'error list': err_list}
        return cls(stat_dict)


    @classmethod
    def load(cls, statistic_num):
        """Загружает данные проведенных симуляций как словарь"""
        fname = files(f'opticalcoating.resources.Statistics').joinpath(
            f'Statistic{str(statistic_num).zfill(3)}.json')
        with open(fname, 'r') as file: stat_dict = json.load(file)
        return cls(stat_dict)


    def make_dict(self):
        stat_dict = {'design': self.des_name,
                     'target': self.trg_name,
                     'creation_time': self.creation_time,
                     'start_rnd_seed': self.start_rnd_seed,
                     'waves': self.waves,
                     'term_algs': self.term_algs[1:],
                     'rates': self.rates[1:],
                     'rates_sigmas': self.rates_sigmas[1:],
                     'meas_sigmas': self.meas_sigmas[1:],
                     'error list': self.error_list}
        return stat_dict


    def save(self):
        file_name = find_file_name('Statistic')
        with open(file_name, 'w') as file:
            json.dump(self.make_dict(), file, indent=3)
        print(f'"{file_name}" is successfully saved')


    def save_plain_txt(self):
        file_name = find_file_name('Statistic', ext='.txt')

        n = len(self.error_list)
        m = len(self.error_list[0])

        with open(file_name, 'w') as file:
            for i in range(n):
                for j in range(m):
                    print(self.error_list[i][j], end='\t', file=file)
                print('', file=file)


    def mean_error_norm(self):
        errors = pd.DataFrame(self.error_list)
        M, N = errors.shape
        errors_norm = pd.Series([np.linalg.norm(errors.iloc[i, :]) for i in range(M)])
        return errors_norm.mean()


    def error_rms(self):
        errors = pd.DataFrame(self.error_list)
        M, N = errors.shape
        errors_rms = pd.Series([np.linalg.norm(errors.iloc[:, j]) / np.sqrt(M) for j in range(N)])
        return errors_rms

    def delta(self):
        errors = pd.DataFrame(self.error_list)
        _, N = errors.shape
        return np.linalg.norm(self.error_rms()) / np.sqrt(N)

    # Visualisation
    def rms_bar(self, ymax=None, lang='en', pic_ext=None, **kwargs):
        vis.rms_bar(self, ymax, lang, pic_ext, **kwargs)



def mean_error_norm(num):
    errors = pd.DataFrame(StatInfo.load(num).error_list)
    M, N = errors.shape
    errors_norm = pd.Series([np.linalg.norm(errors.iloc[i, :]) for i in range(M)])
    return errors_norm.mean()


def error_norm_hist(num, *, xmax=None):
    errors = pd.DataFrame(StatInfo.load(num).error_list)
    M, N = errors.shape
    errors_norm = pd.Series([np.linalg.norm(errors.iloc[i, :]) for i in range(M)])

    font_properties = {'size': 22,
                       'family': 'Times New Roman'}
    rc('font', **font_properties)

    plt.figure(figsize=(16, 9))
    if xmax is None:
        sns.histplot(data=errors_norm, bins=40)
        xmax = errors_norm.max()
    else:
        sns.histplot(data=errors_norm, binwidth=1.05 * xmax / 40)
    plt.xlim(0., 1.05 * xmax)
    plt.xlabel('Значение нормы вектора ошибок')
    plt.ylabel('Число симуляций')





