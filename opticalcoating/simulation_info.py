import json
from .tools import find_file_name


class SimInfo:
    def __init__(self, des_th_d, des_act_d, time_list_res, flux_meas_res, term_cond_case, wavelength,
                 rnd_seed=None, d_j_act_t=None, set_up_pars=None, des=None, nonloccoef=None):
        self.N_layers = len(des_th_d) - 1
        self.rnd_seed = rnd_seed
        self.d_th = des_th_d
        self.d_act = des_act_d
        self.d_j_act_t = d_j_act_t
        self.time_list = time_list_res
        self.flux_meas = flux_meas_res
        self.term_cond_case = term_cond_case
        self.errors_d = [d_act - d_th for (d_act, d_th) in list(zip(*[des_act_d, des_th_d]))]
        self.wavelength = wavelength
        # Сложные класс, которые в json не записываем
        self.set_up_pars = set_up_pars
        self.des = des
        self.nonloccoef = nonloccoef

    def d_act_t(self, j, i):
        if self.d_j_act_t is not None:
            res = [0.0 for _ in range(self.N_layers + 1)]
            res[0] = self.d_th[0]
            for layer in range(1, j):
                res[layer] = self.d_j_act_t[layer][-1]
            res[j] = self.d_j_act_t[j][i]
            return res

    def make_dict(self):
        sim_dict = {'time_list': self.time_list,
                    'flux_meas': self.flux_meas,
                    'wavelength': self.wavelength,
                    'actual thicnesses': self.d_act}
        if self.rnd_seed is not None:
            sim_dict['rnd_seed'] = self.rnd_seed
        if self.d_j_act_t is not None:
            sim_dict['d_j_act_t'] = self.d_j_act_t
        return sim_dict

    def save(self):
        file_name = find_file_name('Simulation')

        with open(file_name, 'w') as file:
            json.dump(self.make_dict(), file, indent=3)
            file.close()
