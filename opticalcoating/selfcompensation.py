import numpy as np
from opticalcoating.calc_flux import calc_flux


def MF_calc(des, wv, T_target, q_percent=False):
    # T_target not in persent too
    L = len(wv)
    T_array = np.zeros((L, 1))
    A_array = np.zeros((L, 1))
    for i in range(L):
        T_array[i] = calc_flux(des, wv[i], q_subs=False, q_percent=q_percent, q_TR='T')
        A_array[i] = (T_array[i] - T_target[i])**2

    MF = np.sqrt(np.sum(A_array) / L)
    return MF


def delta_MF_calc(des, wv, T_target, MF_d_th=0, delta_error=None):
    # delta_d = [sum(n) for n in zip_longest(des.d, delta_error, fillvalue=0)]
    # delta_error - 36 layers
    d_th = des.d
    if delta_error is not None:
        des.d = [des.d[0]] + [des.d[i] + delta_error[i - 1] for i in range(1, des.N + 1)]
    # реализовать если процентах T
    delta_MF = MF_calc(des, wv, T_target) - MF_d_th
    des.d = d_th
    return delta_MF


def mean_delta_MF_rnd_calc(des, wv, T_target, MF_d_th, delta_random):
    #M number of simulations
    M, m = delta_random.shape
    delta_MF_rnd = np.zeros((M, 1))
    for i in range(M):
        delta_MF_rnd[i] = delta_MF_calc(des, wv, T_target, MF_d_th, delta_random[i,:])
    mean_value = np.mean(delta_MF_rnd)
    return mean_value


def c_calc(des, wv, T_target, MF_d_th, delta_sim, mean_delta_MF_rnd):
    c = delta_MF_calc(des, wv, T_target, MF_d_th, delta_sim) / mean_delta_MF_rnd
    return c


def rnd_generation(M, N, inv):
    # additive errors case
    return np.random.normal(0., inv/np.sqrt(N), size=(M, N))

