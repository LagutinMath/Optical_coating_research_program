from opticalcoating.selfcompensation import *
import pandas as pd
from opticalcoating.statistics_info import load_dict
import json

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc
from opticalcoating.save_data import find_file_name


class ProcessedStatistics:
    def __init__(self, *, des=None, target=None, statistic_num=None):
        init1 = (des is not None and target is not None and statistic_num is not None)
        init2 = (des is None and target is None and statistic_num is not None)
        if init1:
            waves, T_target = target
            self.statistic_num = statistic_num
            self.MF_d_th = MF_calc(des, waves, T_target)
            sim_error_list = load_dict(statistic_num)
            sim_errors = pd.DataFrame(sim_error_list['error list'])
            # Число симуляций:
            self.M = len(sim_errors)

            sim_errors_norm = [np.linalg.norm(sim_errors.iloc[i,:]) for i in range(self.M)]
            median_val = np.median(sim_errors_norm)
            self.rnd_errors = pd.DataFrame(rnd_generation(self.M, des.N, median_val))

            delta_MF_rnd = [delta_MF_calc(des, waves, T_target, self.MF_d_th, self.rnd_errors.iloc[i, :]) for i in range(self.M)]
            self.mean_delta_MF_rnd = np.mean(delta_MF_rnd)

            delta_MF_sim = np.array([delta_MF_calc(des, waves, T_target, self.MF_d_th, sim_errors.iloc[i, :]) for i in
                                 range(self.M)])
            self.c_array = delta_MF_sim / self.mean_delta_MF_rnd
            self.c_value = self.c_array.mean()
        elif init2:
            file_name = 'c_values/c_value' + str(statistic_num).zfill(3) + '.json'
            with open(file_name, 'r') as file:
                info_dict = json.load(file)
                file.close()

            self.statistic_num = statistic_num
            self.MF_d_th = info_dict['MF_d_th']
            self.rnd_errors = pd.DataFrame(info_dict['rnd_error_list'])
            self.mean_delta_MF_rnd = info_dict['mean_delta_MF_rnd']
            self.c_array = np.array(info_dict['c_array'])
            self.c_value = info_dict['c_value']
        else:
            raise NameError('Wrong Initialization')


    def make_dict(self):
        sim_dict = {'MF_d_th': self.MF_d_th, 'rnd_error_list': self.rnd_errors.values.tolist(),
                    'mean_delta_MF_rnd': self.mean_delta_MF_rnd, 'c_array': self.c_array.tolist(), 'c_value': self.c_value}
        return sim_dict


    def save(self):
        file_name = 'c_values/c_value' + str(self.statistic_num).zfill(3) + '.json'

        with open(file_name, 'w') as file:
            json.dump(self.make_dict(), file, indent=3)
            file.close()


    def c_hist(self, *, show=False, xmax=None):
        font_properties = {'size': 22,
                           'family': 'Times New Roman'}
        rc('font', **font_properties)

        lolik = pd.Series(self.c_array)

        plt.figure(figsize=(16, 9))
        sns.histplot(data=lolik, bins=40)
        if xmax is None:
            xmax = lolik.max()
        plt.xlim(0., 1.05 * xmax)
        plt.xlabel('Значение c')
        plt.ylabel('Число симуляций')
        plt.savefig(find_file_name('Picture', '.png'))

        if show:
            plt.show()

