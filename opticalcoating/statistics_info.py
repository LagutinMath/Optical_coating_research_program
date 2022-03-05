from opticalcoating.save_data import find_file_name
import json

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