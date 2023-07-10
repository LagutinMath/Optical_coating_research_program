import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import seaborn as sns
from datetime import datetime
from importlib.resources import files
from .save_data import find_file_name


class StatInfo:
    def __init__(self, stat_dict):
        self.des_name = stat_dict['design']
        self.creation_time = stat_dict['creation_time']
        self.start_rnd_seed = stat_dict['start_rnd_seed']
        self.waves = stat_dict['waves']
        self.term_algs = [None] + stat_dict['term_algs']
        self.rates = [None] + stat_dict['rates']
        self.rates_sigmas = [None] + stat_dict['rates_sigmas']
        self.meas_sigmas = [None] + stat_dict['meas_sigmas']
        self.error_list = stat_dict['error list']


    @classmethod
    def legacy(cls, des_th, term_algs, set_up_pars, error_list, start_rnd_seed):
        stat_dict = {'design': des_th.name,
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
    def load(cls, statistic_num):
        """Загружает данные проведенных симуляций как словарь"""
        fname = files(f'opticalcoating.resources.Statistics').joinpath(
            f'Statistic{str(statistic_num).zfill(3)}.json')
        with open(fname, 'r') as file: stat_dict = json.load(file)
        return cls(stat_dict)


    def make_dict(self):
        stat_dict = {'design': self.des_name,
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
            file.close()


    def save_plain_txt(self):
        file_name = find_file_name('Statistic', ext='.txt')

        n = len(self.error_list)
        m = len(self.error_list[0])

        with open(file_name, 'w') as file:
            for i in range(n):
                for j in range(m):
                    print(self.error_list[i][j], end='\t', file=file)
                print('', file=file)
            file.close()


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


def error_rms_bar(num, *, ymax=None, colored=True, special_layers=None):
    errors = pd.DataFrame(StatInfo.load(num).error_list)
    M, N = errors.shape
    errors_rms = pd.Series([np.linalg.norm(errors.iloc[:, j])/np.sqrt(M) for j in range(N)])

    font_properties = {'size': 22,
                       'family': 'Times New Roman'}
    rc('font', **font_properties)

    fig = plt.figure(figsize=(16, 9))
    if colored:
        ax = fig.add_subplot()
        ax.bar(range(1, N + 1, 2), errors_rms[0:N:2], color='b')
        ax.bar(range(2, N + 1, 2), errors_rms[1:N:2], color='r')
        if special_layers is not None:
            sp_errors = [errors_rms[i - 1] for i in special_layers]
            ax.bar(special_layers, sp_errors, color='m')
    else:
        plt.bar(x=range(1, N + 1), height=errors_rms)

    plt.xlim(1 - 0.5, N + 0.5)
    if ymax is None:
        ymax = max(errors_rms)
    plt.ylim(0., 1.05 * ymax)
    plt.xlabel('Номер слоя')
    plt.ylabel('Среднеквадратичная ошибка на слое, нм')


