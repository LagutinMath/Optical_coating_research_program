# Статистика номер 6
sim_num = 6
# С новым дизайном WDM302

from time import perf_counter
import multiprocessing as mp

from opticalcoating.deposition_simulation import simulation, SetUpParameters
from opticalcoating.design_class import Design
from opticalcoating.calc_flux import Wave
from opticalcoating.statistics_info import StatInfo, mean_error_norm, error_norm_hist, error_rms_bar

from opticalcoating.processed_statistics_class import ProcessedStatistics


def simulation_mult(kwarg):
    # {'des_th', 'term_algs', 'set_up_pars', 'rnd_seed'}
    return simulation(**kwarg).errors_d[1:]


if __name__ == '__main__':
    N_proc = mp.cpu_count()
    N_sim = 2000

    test_des = Design(name='WDM302')
    term_algs = (test_des.N + 1) * ['Quasiswing']

    waves = [None] + test_des.N * [Wave(1550)]

    r_H, r_L = 0.24, 0.18
    rates = [None] + int(test_des.N / 2) * [r_H, r_L]
    rates_sigmas = [None] + int(test_des.N / 2) * [0.05 * r_H, 0.1 * r_L]

    set_up_pars = SetUpParameters(N=test_des.N, waves=waves, rates=rates, q_TR='T', meas_sigmas=0.001,
                                  rates_sigmas=rates_sigmas)

    tasks = [{'des_th': test_des, 'term_algs': term_algs, 'set_up_pars': set_up_pars, 'rnd_seed': rnd_seed} for rnd_seed in range(10000000, 10000000 + N_sim)]
    p = mp.Pool(processes=N_proc)


    start = perf_counter()
    err_list = p.map(simulation_mult, tasks)
    p.close()
    stop = perf_counter()
    print(f'mean sim time = {(stop - start) / N_sim} second')
    print(f'total time of calc = {(stop - start)} second')


    res = StatInfo(test_des, term_algs, set_up_pars, err_list, 10000000)
    res.save()

    # -----------------------------------------------

    wv_list = [1546.6075, 1547.9627, 1548.1323, 1548.3019, 1548.5015, 1548.6084, 1548.7154, 1548.8223, 1548.9293,
               1549.0363, 1549.1433, 1549.2504, 1549.3574, 1549.4645, 1549.5715, 1549.6786, 1549.7857, 1549.8929,
               1550, 1550.1072, 1550.2143, 1550.3215, 1550.4287, 1550.5359, 1550.6431, 1550.7504, 1550.8576, 1550.9649,
               1551.0722, 1551.1795, 1551.2868, 1551.3941, 1551.5015, 1551.7019, 1551.8723, 1552.0427, 1553.4075]
    wv = [Wave(x) for x in wv_list]
    T_target = 3 * [0.0] + 31 * [1.0] + 3 * [0.0]

    start = perf_counter()
    c_vals = ProcessedStatistics(des=test_des, target=(wv, T_target), statistic_num=sim_num)
    stop = perf_counter()
    print(f'time to calc c = {stop - start} second')

    c_vals.save()

# на 2000 случаях
# mean sim time = 15.6 second
# total time of calc = 31248 second
# time to calc c = 141 second

