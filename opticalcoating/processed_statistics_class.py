from opticalcoating.selfcompensation import *
import pandas as pd
from opticalcoating.statistics_info import load_dict
import json

class ProcessedStatistics:
    def __init__(self, des, target, statistic_num):
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


    def make_dict(self):
        sim_dict = {'MF_d_th': self.MF_d_th, 'rnd_error_list': self.rnd_errors.values.tolist(),
                    'mean_delta_MF_rnd': self.mean_delta_MF_rnd, 'c_array': self.c_array.tolist(), 'c_value': self.c_value}
        return sim_dict


    def save(self):
        file_name = 'c_values/c_value' + str(self.statistic_num).zfill(3) + '.json'

        with open(file_name, 'w') as file:
            json.dump(self.make_dict(), file, indent=3)
            file.close()

