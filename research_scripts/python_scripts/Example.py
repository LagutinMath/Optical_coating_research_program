import sys
sys.path.append('../../')
import opticalcoating as oc
import multiprocessing as mp
from os import listdir
from importlib.resources import files


def stat_nums():
    names = listdir(files(f'opticalcoating.resources.Statistics'))
    return {int(name[9:12]) for name in names if name.startswith('Statistic')}

def simulation_mult(kwarg):
    # kwarg = {'des_th', 'term_algs', 'set_up_pars', 'rnd_seed'}
    return oc.simulation(**kwarg).errors_d[1:]

if __name__ == '__main__':
    N_proc = mp.cpu_count() - 1
    N_sim = 100

    default = {}
    default['term_alg'] = 'QS'
    default['rate'] = [0.5, 0.8] # nm / sec
    default['std rate'] = [5, 10] # in percent
    default['std meas'] = 0.1 # in percent
    default['T/R'] = 'T'
    default['polarisation'] = 'S'
    default['angle'] = 0
    default['backside'] = False

    sim_params = ({'design': "AR4_1",
                   'target': "AR4",
                   'wave': 400,
                   'T/R': 'R'}, )

    for sim_param in [oc.SimParams.from_dict(param, default) for param in sim_params]:
        with mp.Pool(processes=N_proc) as p:
            err_list = oc.timer()(p.map)(simulation_mult, sim_param.get_tasks(N_sim))
        oc.StatInfo.from_sim_param(sim_param, err_list).save()

    for stat_num in stat_nums():
        oc.timer()(oc.ProcessedStatistics.calc)(stat_num).save()

    oc.process_statistics()






