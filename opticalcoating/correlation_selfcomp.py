from opticalcoating.deposition_simulation import *


def correlation(A):
    M = len(A[:, 1])
    m = len(A[1, :])
    print('sim=',M)
    print('layers=',m)

    mu = (1/M)*A.T@A
    print('mu=',mu)
    vals, vect = np.linalg.eig(mu)
    print('eig_value=',vals)
    #vals=np.real(vals) ???
    sgm_list = np.sqrt(vals)
    print('sqm=', sgm_list)
    mse=np.square(sgm_list).mean()
    rmse=math.sqrt(mse)
    print('rmse =',rmse)
    sgm_av=rmse
    sgm_list[::-1].sort()
    print('sorted eig value=',sgm_list)

    p_array=np.zeros((m,1))
    for i in range(m):
        p_array[i] = sgm_av/sgm_list[i]
    print('p_array =',p_array)

    beta=np.prod(p_array)**(1/m)
    return beta


def MF_calc(des, wv, T_target):
    #T_target not in persent too
    L=len(wv)
    T_array=np.zeros((L,1))
    A_array = np.zeros((L, 1))
    for i in range(L):
        T_array[i]=calc_flux(des, wv[i], q_subs=False, q_percent=False, n_a=1, q_TR='T')
        A_array[i] = (T_array[i]-T_target[i])**2

    MF = math.sqrt(np.sum(A_array)/L)
    return MF

def delta_MF_calc(des, wv, T_target, MF_d_th=0, delta_error=None):
    #delta_d = [sum(n) for n in zip_longest(des.d, delta_error, fillvalue=0)]
    #delta_error - 36 layers
    if delta_error is not None:
        des.d = [des.d[0]]+[des.d[i]+delta_error[i-1] for i in range(1,des.N+1)]


    #реализовать если процентах T
    delta_MF = MF_calc(des, wv, T_target) - MF_d_th
    return delta_MF

def mean_delta_MF_rnd_calc(M, des, wv, T_target, MF_d_th, delta_random):
    #M number of simulations
    delta_MF_rnd = np.zeros((M, 1))
    for i in range(M):
        delta_MF_rnd[i] = delta_MF_calc(des, wv, T_target, MF_d_th, delta_random[i,:])
    mean_value = np.mean(delta_MF_rnd)
    return mean_value

def c_calc(des, wv, T_target, MF_d_th, delta_sim, mean_delta_MF_rnd):
    c = delta_MF_calc(des, wv, T_target, MF_d_th, delta_sim)/mean_delta_MF_rnd
    return c

def rnd_generation(M,N, inv):
    #additive errors case
    rnd_matrix = np.zeros((M,N))
    for i in range(M):
        for j in range(N):
            rnd_matrix[i,j] = np.random.normal(0, inv/math.sqrt(N))

    return rnd_matrix






# des = Design('Z36')
# a=np.array([1,2,3])
# b=np.array([4,5])
# result = [sum(n) for n in zip_longest(des.d, b, fillvalue=0)]
# des.d=result
# print(result)
# print(des.d)


# # Загрузка данных корреляции
# d_act_mat = scipy.io.loadmat("target_list.mat")
# print(d_act_mat.keys())
# print(d_act_mat)



















# -------test of reverse R and R_b------------------------
# def R_test(R, R_b, sgm):
#     R_value=R_b+(sgm**2 * R * (1-R_b)**2)/(1-sgm**2 * R * R_b)
#     return R_value
#
#
# print(R_test(0.1,0.2,0),'is equal',R_test(0.2,0.1,0))
# print(R_test(0.1,0.2,0.8) - R_test(0.2,0.1,0.8))