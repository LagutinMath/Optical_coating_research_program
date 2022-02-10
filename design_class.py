from scipy.interpolate import interp1d
import math
import numpy as np
import json
import os.path
# from numba import njit
import matplotlib.pyplot as plt

class Design:
    def __init__(self, des_name=None, thicknesses=None, n_const=None, des_json=None):
        self.q_n_const = False
        self.name = ''
        self.n_const = []

        if des_json is not None:
            self.materials = des_json
            self.name = des_json['name']
        elif des_name is not None:
            fname = 'Designs/' + des_name + '.json'
            if os.path.isfile(fname):
                with open(fname, 'r') as file:
                    self.materials = json.load(file)
                    self.name = self.materials['name']
                    if 'n_const' in self.materials:
                        self.q_n_const = True
                        n_const = self.materials['n_const']
            else:
                self.name = des_name

        if n_const is not None:
            self.q_n_const = True
            self.n_const = n_const

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

    def create_simple_json(self):
        # подразумевается, что есть только имя и вектора d и n
        fname = 'Designs/' + self.name + '.json'
        if not os.path.isfile(fname):
            with open(fname, 'w') as file:
                inf = json.dumps({'name': self.name, 'thicknesses': self.d, 'n_const': self.n_const}, indent=4)
                print(inf, file=file)

    def create_json(self):
        print(1)
        fname = 'Designs/' + self.name + '.json'
        print(fname)
        if not os.path.isfile(fname):
            with open(fname, 'w') as file:
                inf = json.dumps(self.materials, indent=4)
                print(inf, file=file)

    def increase_layer_thickness(self, j, delta_d):
        self.d[j] += delta_d

    def n(self, layer_num, wv):
        if self.n_fix[1][layer_num] != wv.wavelength:
            # при следующем вызове с той же длиной волны n не будет пересчитываться
            self.n_fix[1][layer_num] = wv.wavelength
            # В случае истинности self.q_n_const ничего делать не надо. Значения инициализированы конструктором
            if not self.q_n_const:
                layer_material = self.materials["mat_inf"][self.materials["layers"][layer_num]]
                if "Sellmeier" in layer_material:
                    self.n_fix[0][layer_num] = Sellmeier_n(layer_material["Sellmeier"], wv.wavelength)
                elif "Cauchy" in layer_material:
                    self.n_fix[0][layer_num] = layer_material["Cauchy"][0] + layer_material["Cauchy"][
                        1] / wv.wavelength ** 2 + layer_material["Cauchy"][2] / wv.wavelength ** 4
                else:
                    f_interp = interp1d(layer_material["Table"]["wavelength"], layer_material["Table"]["n"])
                    self.n_fix[0][layer_num] = float(f_interp(wv.wavelength))
        return self.n_fix[0][layer_num]

    def thickness_plot(self):
        # Сделать картинку во весь экран для 27 дюймового монитора
        # diag = math.sqrt(16**2 + 9**2)
        # ipseg = 27/diag  # inches per segment
        # fig = plt.figure('Thicknesses', figsize=(16 * ipseg, 9 * ipseg))

        fig = plt.figure('Thicknesses', figsize=(16, 9))
        ax = fig.add_subplot()
        ax.bar(range(1, self.N + 1, 2), self.d[1:self.N + 1:2], color='b')
        ax.bar(range(2, self.N + 1, 2), self.d[2:self.N + 1:2], color='r')
        plt.xlabel('Layer number')
        plt.ylabel('Physical thickness, nm')
        # plt.title('Design physical thicknesses')
        plt.show()

    def spectral_plot(self, *, q_TR='T', wv_range=[380, 760]):
        fig = plt.figure('Spectral plot', figsize=(16, 9))
        ax = fig.add_subplot()
        N_pts = 1000
        x_range = np.linspace(wv_range[0], wv_range[1], N_pts)
        y_range = N_pts * [0.0]
        for i in range(N_pts):
            y_range[i] = calc_flux(self, Wave(x_range[i]), q_percent=True, q_TR=q_TR)
        ax.plot(x_range, y_range)
        plt.xlabel('Wavelength, nm')
        plt.ylabel('T, %')
        # plt.title('Design physical thicknesses')
        plt.show()


def Sellmeier_n(coef_A, wvlen):
    # wvlen in nm
    wv2 = wvlen ** 2
    pwr2n = (coef_A[0] + wv2 * (coef_A[1] / (wv2 - coef_A[2]) +
                                coef_A[3] / (wv2 - coef_A[4]) +
                                coef_A[5] / (wv2 - coef_A[6])))
    return np.sqrt(pwr2n)


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


def calc_flux(des, wv, *, q_subs=True, q_percent=False, n_a=1, q_TR='both'):
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

