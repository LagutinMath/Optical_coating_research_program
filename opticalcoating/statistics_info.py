import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from opticalcoating.save_data import find_file_name
import numpy as np
import seaborn as sns


class StatInfo:
    def __init__(self, des_th, term_algs, set_up_pars, err_list, start_rnd_seed):
        self.err_list = err_list
        self.start_rnd_seed = start_rnd_seed
        self.des = des_th
        self.set_up_pars = set_up_pars
        self.term_algs = term_algs

        self.N_sim = len(err_list)

    def make_dict(self):
        stat_dict = {'design': self.des.name, 'start_rnd_seed': self.start_rnd_seed, 'error list': self.err_list,
                     'term_algs': self.term_algs}
        return stat_dict

    def save(self):
        file_name = find_file_name('Statistic')

        with open(file_name, 'w') as file:
            json.dump(self.make_dict(), file, indent=3)
            file.close()

    def save_plain_txt(self):
        file_name = find_file_name('Statistic', ext='.txt')

        n = len(self.err_list)
        m = len(self.err_list[0])

        with open(file_name, 'w') as file:
            for i in range(n):
                for j in range(m):
                    print(self.err_list[i][j], end='\t', file=file)
                print('', file=file)
            file.close()


def load_dict(num):
    """Загружает данные проведенных симуляций как словарь"""
    fname = 'Statistics/Statistic' + str(num).zfill(3) + '.json'
    with open(fname, 'r') as file:
        return json.load(file)


def mean_error_norm(num):
    info = load_dict(num)
    errors = pd.DataFrame(info['error list'])
    M, N = errors.shape
    errors_norm = pd.Series([np.linalg.norm(errors.iloc[i, :]) for i in range(M)])

    return errors_norm.mean()


def error_norm_hist(num, *, show=False, xmax=None):
    info = load_dict(num)
    errors = pd.DataFrame(info['error list'])
    M, N = errors.shape
    errors_norm = pd.Series([np.linalg.norm(errors.iloc[i, :]) for i in range(M)])

    font_properties = {'size': 22,
                       'family': 'Times New Roman'}
    rc('font', **font_properties)

    plt.figure(figsize=(16, 9))
    sns.histplot(data=errors_norm, bins=40)
    if xmax is None:
        xmax = errors_norm.max()
    plt.xlim(0., 1.05 * xmax)
    plt.xlabel('Значение нормы вектора ошибок')
    plt.ylabel('Число симуляций')

    if show:
        plt.show()
    else:
        plt.savefig(find_file_name('Picture', '.png'))


def error_rms_bar(num, *, show=False, ymax=None):
    info = load_dict(num)
    errors = pd.DataFrame(info['error list'])
    M, N = errors.shape
    errors_rms = pd.Series([np.linalg.norm(errors.iloc[:, j])/np.sqrt(M) for j in range(N)])

    font_properties = {'size': 22,
                       'family': 'Times New Roman'}
    rc('font', **font_properties)

    plt.figure(figsize=(16, 9))
    plt.bar(x=range(1, N + 1), height=errors_rms)
    plt.xlim(1 - 0.5, N + 0.5)
    if ymax is None:
        ymax = max(errors_rms)
    plt.ylim(0., 1.05 * ymax)
    plt.xlabel('Номер слоя')
    plt.ylabel('Среднеквадратичная ошибка на слое, нм')

    if show:
        plt.show()
    else:
        plt.savefig(find_file_name('Picture', '.png'))

