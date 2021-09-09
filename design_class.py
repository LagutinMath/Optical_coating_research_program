from scipy.interpolate import interp1d
import math
import numpy as np
import json
# from numba import njit
import matplotlib.pyplot as plt


class Design:
    def __init__(self, des_name=None, thicknesses=None, n_const=None, des_json=None):
        self.q_n_const = False

        if des_json is not None:
            self.materials = des_json
        elif des_name is not None:
            with open('Designs/' + des_name + '.json', 'r') as file:
                self.materials = json.load(file)
        elif n_const is not None:
            self.q_n_const = True

        if thicknesses is None:
            self.d = self.materials["thicknesses"]
        else:
            self.d = thicknesses

        self.N = len(self.d) - 1
        # self.n_fix: 
        # первый столбец --- посчитанные ранее показатели преломления 
        # второй столбец --- для какой длины волны они посчитаны
        self.n_fix = np.zeros((2, self.N + 1), dtype=float)
        if self.q_n_const:
            self.n_fix[0] = n_const

    def increase_layer_thickness(self, j, delta_d):
        self.d[j] += delta_d

    def n(self, layer_num, wv):
        if self.n_fix[1][layer_num] != wv.wavelength:
            # при следующем вызове с той же длиной волны n не будет пересчитываться
            self.n_fix[1][layer_num] = wv.wavelength
            # В случае истинности self.q_n_const ничего делать не надо. Значения инициализированы конструктором
            if not self.q_n_const:
                layer_material = self.materials["mat_inf"][self.materials["layers"][layer_num]]
                if "Cauchy" in layer_material:
                    self.n_fix[0][layer_num] = layer_material["Cauchy"][0] + layer_material["Cauchy"][
                        1] / wv.wavelength ** 2 + layer_material["Cauchy"][2] / wv.wavelength ** 4
                else:
                    f_interp = interp1d(layer_material["Table"]["wavelength"], layer_material["Table"]["n"])
                    self.n_fix[0][layer_num] = float(f_interp(wv.wavelength))
        return self.n_fix[0][layer_num]

    def thickness_plot(self):
        plt.bar(range(1, self.N + 1), self.d[1:self.N + 1])
        plt.show()


class Wave:
    def __init__(self, wavelength, polarisation='S', angle=0):
        self.wavelength = wavelength
        self.polarisation = polarisation
        self.angle = angle


def sp_mat_mult(A, B):  # special matrix multiplication
    m00 = A[0][0] * B[0][0] - A[0][1] * B[1][0]
    m01 = A[0][0] * B[0][1] + A[0][1] * B[1][1]
    m10 = A[1][0] * B[0][0] + A[1][1] * B[1][0]
    m11 = A[1][1] * B[1][1] - A[1][0] * B[0][1]
    return [[m00, m01], [m10, m11]]


def calc_flux(des, wv, q_subs=True, q_percent=False, n_a=1, q_TR='both'):
    M = [[1.0, 0.0], [0.0, 1.0]]
    for j in range(des.N + 1):
        if wv.angle == 0:
            q_j = n_j = des.n(j, wv)
            phi_j = 2.0 * np.pi * n_j * des.d[j] / wv.wavelength
        else:
            n_j = des.n(j, wv)
            gamma_j = np.arcsin(np.sin(wv.angle) * n_a / n_j)
            q_j = n_j * np.cos(gamma_j) if wv.polarisation == 'S' else n_j / np.cos(gamma_j)
            phi_j = 2.0 * np.pi * n_j * des.d[j] * np.cos(gamma_j) / wv.wavelength

        cos_phi = math.cos(phi_j)
        sin_phi = math.sin(phi_j)

        M_j_00 = cos_phi
        M_j_01 = sin_phi / q_j
        M_j_10 = sin_phi * q_j
        M_j_11 = cos_phi

        M = sp_mat_mult([[M_j_00, M_j_01], [M_j_10, M_j_11]], M)

    n_s = des.n(0, wv)
    if wv.angle == 0:
        q_s = n_s
        q_a = n_a
    else:
        gamma_s = np.arcsin(np.sin(wv.angle) * n_a / n_s)
        q_s = n_s * np.cos(gamma_s) if wv.polarisation == 'S' else n_s / np.cos(gamma_s)
        q_a = n_a * np.cos(wv.angle) if wv.polarisation == 'S' else n_a / np.cos(wv.angle)

    T_1 = 4.0 * q_a * q_s / ((M[0][0] * q_a + M[1][1] * q_s) ** 2 + (M[1][0] + M[0][1] * q_a * q_s) ** 2)
    R_1 = (((M[0][0] * q_a - M[1][1] * q_s) ** 2 + (M[0][1] * q_a * q_s - M[1][0]) ** 2) /
           ((M[0][0] * q_a + M[1][1] * q_s) ** 2 + (M[0][1] * q_a * q_s + M[1][0]) ** 2))
    if q_subs:
        T_2 = 4.0 * q_a * q_s / (q_a + q_s) ** 2
        R_2 = ((q_a - q_s) / (q_a + q_s)) ** 2
        T = T_1 * T_2 / (1.0 - R_1 * R_2)
        R = (R_1 + (R_2 * T_1 * T_1) / (1.0 - R_1 * R_2))
    else:
        T, R = T_1, R_1

    if q_TR == 'both':
        return (100.0 * T, 100.0 * R) if q_percent else (T, R)
    elif q_TR == 'T':
        return 100.0 * T if q_percent else T
    elif q_TR == 'R':
        return 100.0 * R if q_percent else R

