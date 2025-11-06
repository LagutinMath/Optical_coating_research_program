# import sys
# sys.path.append('../..')
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
    # N_proc = mp.cpu_count() - 1
    # N_sim = 100

    default = {}
    default['term_alg'] = 'QS'
    default['rate'] = [0.5, 0.8] # nm / sec
    default['std rate'] = [5, 10] # in percent
    default['std meas'] = 0.1 # in percent
    default['T/R'] = 'R'
    default['polarisation'] = 'S'
    default['angle'] = 0
    default['backside'] = False

    sim_params = {'design': "QWM20",
                   'target': "QWM20",
                   'wave': 510,
                   'T/R': 'R'}
    
    sim_params.update(default)
    
    # print(oc.simulation(**sim_params).errors_d[1:])

    set_up_pars = oc.SetUpParameters(N=20, waves=510, q_TR='R', tau=2., rates=0.5, mono_width=0., meas_sigmas=None, meas_syst=0.,
                 rates_sigmas=None, rates_syst=0., delay_time=0., delay_time_sigma=0., r_index_syst=None,
                 r_index_sigmas=None, backside=False, width=None)

    print(oc.simulation(des_th="QWM20", term_algs="QS", set_up_pars=set_up_pars, rnd_seed=None, q_subs=True, save_M=False))

    print()

    # for sim_param in [oc.SimParams.from_dict(param, default) for param in sim_params]:
    #     with mp.Pool(processes=N_proc) as p:
    #         err_list = oc.timer()(p.map)(simulation_mult, sim_param.get_tasks(N_sim))
    #     oc.StatInfo.from_sim_param(sim_param, err_list).save()

    # for stat_num in stat_nums():
    #     oc.timer()(oc.ProcessedStatistics.calc)(stat_num).save()

    # oc.process_statistics()
