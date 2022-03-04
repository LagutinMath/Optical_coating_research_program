import numpy as np

from correlation_selfcomp import *

#-------correlation-----------------
A = np.array([[1,2,3],[6,5,9]])
B = np.loadtxt('Statistics/Statistic003.txt')
#print(correlation(B))


#------selfcompensation-------------
delta_sim=np.loadtxt('Statistics/Statistic003.txt')
M = len(delta_sim[:, 1]) #number of sim
m = len(delta_sim[1, :]) #number of layers
des = Design('Z36')
# wv=np.loadtxt('lambda_list_test.txt')
wv=[Wave(x) for x in range(400,642,2)]+[Wave(x) for x in range(680,802,2)]
T_target=121*[1,]+61*[0,]
# print(len(wv))
# print(len(T_target))

# k=0.09186683
# a=MF_calc(des, wv, T_target)
# print(a*100-k)
# print(a)



MF_d_th=0.0009186697517661457
sim_norm=np.zeros((M,1))
for i in range(M):
    sim_norm[i]=np.linalg.norm(delta_sim[i,:])
inv=np.median(sim_norm)

delta_random=rnd_generation(M,m,inv)
mean_delta_MF_rnd = mean_delta_MF_rnd_calc(M, des, wv, T_target, MF_d_th, delta_random)
print(mean_delta_MF_rnd)

c_array = np.zeros((M,1))
for i in range(M):
    c_array[i] = c_calc(des, wv, T_target, MF_d_th, delta_sim[i,:], mean_delta_MF_rnd)
mean_value_c = np.mean(c_array)

print(mean_value_c)
print(c_array)

#print(delta_MF_calc(des, wv, T_target, MF_d_th, delta_error))



#not right c value, app 2, should be<1 !!!!!!!!!!!!!!!!!