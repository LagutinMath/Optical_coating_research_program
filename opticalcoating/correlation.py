import numpy as np


def correlation(A):
    """Вычисление числа beta --- коэф.корреляции"""
    # M --- number of simulations
    # m --- number of layers
    M, m = A.shape
    mu = (1/M) * A.T @ A
    vals, vect = np.linalg.eig(mu)
    sgm_list = np.sqrt(vals)
    mse = np.square(sgm_list).mean()
    rmse = np.sqrt(mse)
    sgm_av = rmse
    p_array = np.zeros(m)
    for i in range(m):
        p_array[i] = sgm_av/sgm_list[i]
    return np.prod(p_array)**(1/m)


def std_values(A):
    # M = len(A[:, 1])
    # m = len(A[1, :])
    M, m = A.shape
    # print('sim=', M)
    # print('layers=', m)

    mu = (1/M) * A.T @ A
    # print('mu=', mu)
    vals, vect = np.linalg.eig(mu)

    # print('eig_value=', vals)
    sgm_list = np.sqrt(vals)
    # print('sqm=', sgm_list)
    mse = np.square(sgm_list).mean()
    rmse = np.sqrt(mse)
    # print('rmse =',rmse)
    sgm_av = rmse
    sgm_list[::-1].sort()
    return sgm_list

# temp
def eig_vec_0(A):
    # M = len(A[:, 1])
    # m = len(A[1, :])
    M, m = A.shape
    # print('sim=', M)
    # print('layers=', m)

    mu = (1/M) * A.T @ A
    # print('mu=', mu)
    vals, vect = np.linalg.eig(mu)

    # print('eig_value=', vals)
    sgm_list = np.sqrt(vals)
    return vect[0]