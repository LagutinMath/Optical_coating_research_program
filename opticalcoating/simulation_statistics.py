from opticalcoating.statistics_info import StatInfo
from opticalcoating.deposition_simulation import simulation

def collect_statistics(des_th, term_algs, set_up_pars, N_sim, start_rnd_seed=10000000):

    err_list = N_sim * [[]]

    for x in range(N_sim):
        rnd_seed = start_rnd_seed + x
        res = simulation(des_th, term_algs, set_up_pars, rnd_seed)
        err_list[x] = res.errors_d[1:des_th.N + 1]

    return StatInfo(des_th, term_algs, set_up_pars, err_list, start_rnd_seed)