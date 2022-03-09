from scipy.interpolate import interp1d
import numpy as np
import json
import os.path
import copy
import matplotlib.pyplot as plt
import opticalcoating.calc_flux as cf
from opticalcoating.calc_flux import Wave
from math import inf


class Design:
    def __init__(self, *, name=None, d=None, n_const=None, info=None, witness_layers=None):

        # ПРОВЕРКА ПРАВИЛЬНОСТИ ИНИЦИАЛИЗАЦИИ
        # Инициализация по имени дизайна (имя файла)
        init_1 = (name is not None)
        # Инициализация толщинами и коэфициентами преломления без дисперсии
        init_2 = ((d is not None) and (n_const is not None))
        # Инициализация словарём с параметрами
        init_3 = (info is not None)
        if not (init_1 ^ init_2 ^ init_3):
            raise NameError('Wrong initialisation Design.')

        if init_1:
            fname = 'Designs/' + name + '.json'
            with open(fname, 'r') as file:
                self.info = json.load(file)
            self.name = self.info['name']
            self.d = self.info["thicknesses"]
            if 'n_const' in self.info:
                self.q_n_const = True
                self.n_const = self.info['n_const']
            else:
                self.q_n_const = False
        elif init_2:
            self.d = d
            self.q_n_const = True
            self.n_const = n_const
        elif init_3:
            self.info = info
            self.name = info['name']
            self.d = self.info["thicknesses"]
            if 'n_const' in self.info:
                self.q_n_const = True
                self.n_const = self.info['n_const']
            else:
                self.q_n_const = False
        else:
            raise NameError('Impossible. Wrong initialisation Design.')

        self.N = len(self.d) - 1
        if witness_layers is None:
            self.witnesses = [list(range(1,self.N + 1))]
        else:
            self.witnesses = witness_layers
        # текущий наблюдаемый свидетель
        self.w_cur = 0
        self.n_memory = {}
        self.xi_memory = {}


    def witness_num(self, layer):
        """Возвращает номер свидетеля на который будет напыляться слой layer
        Нумерация свидетелей начинается начинается с нуля"""
        for i, whitness in enumerate(self.witnesses):
            if layer in whitness:
                return i


    def witness_layers(self, layer):
        """Возвращает список слоев на свидетеле со слоем layer до слоя layer не включительно"""
        all_layers_on_witness = self.witnesses[self.witness_num(layer)]
        ans = []
        for j in all_layers_on_witness:
            if j < layer:
                ans.append(j)
            else:
                break
        return ans


    def create_simple_json(self, *, name):
        # подразумевается, что есть только имя и вектора d и n
        fname = 'Designs/' + name + '.json'
        if not os.path.isfile(fname):
            with open(fname, 'w') as file:
                inf = json.dumps({'name': name, 'thicknesses': self.d, 'n_const': self.n_const}, indent=4)
                print(inf, file=file)


    def create_json(self):
        print(1)
        fname = 'Designs/' + self.name + '.json'
        print(fname)
        if not os.path.isfile(fname):
            with open(fname, 'w') as file:
                inf = json.dumps(self.info, indent=4)
                print(inf, file=file)


    def increase_layer_thickness(self, j, delta_d):
        self.d[j] += delta_d


    def n(self, layer_num, wv):
        if self.q_n_const:
            ans = self.n_const[layer_num]
        else:
            layer_material = self.info["layers"][layer_num]
            layer_material_info = self.info["mat_inf"][layer_material]
            if (layer_material, wv) in self.n_memory:
                ans = self.n_memory[(layer_material, wv)]
            else:
                if "Sellmeier" in layer_material_info:
                    ans = Sellmeier_n(layer_material_info["Sellmeier"], wv.wavelength)
                elif "Cauchy" in layer_material_info:
                    ans = Cauchy_n(layer_material_info["Cauchy"], wv.wavelength)
                else:
                    f_interp = interp1d(layer_material_info["Table"]["wavelength"], layer_material_info["Table"]["n"])
                    ans = float(f_interp(wv.wavelength))
                self.n_memory[(layer_material, wv)] = ans
        return ans


    def xi(self, layer_num, wv):
        if self.q_n_const:
            ans = 0.0
        else:
            layer_material = self.info["layers"][layer_num]
            layer_material_info = self.info["mat_inf"][layer_material]
            if (layer_material, wv) in self.xi_memory:
                ans = self.xi_memory[(layer_material, wv)]
            else:
                if "xi" in layer_material_info["Table"]:
                    f_interp = interp1d(layer_material_info["Table"]["wavelength"], layer_material_info["Table"]["xi"])
                    ans = float(f_interp(wv.wavelength))
                else:
                    ans = 0.0
            self.xi_memory[(layer_material, wv)] = ans
        return ans


    def wv_bnd(self, layer_num):
        # return min and max correct approximation values of wavelength for n and xi
        if self.q_n_const:
            return [0., inf]
        else:
            layer_material = self.info["mat_inf"][self.info["layers"][layer_num]]
            if "Sellmeier" in layer_material or "Cauchy" in layer_material:
                return [0., inf]
                # There is not restrict in Sellmeier case (you cant use wv = sqrt(A[2]), sqrt(A[4]), sqrt(A[6]))
            else:
                return [layer_material["Table"]["wavelength"][0], layer_material["Table"]["wavelength"][-1]]


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
        # plt.title('Design physical d')
        plt.show()


    def calc_flux(self, wv, *, q_subs=True, backside=False, q_percent=False, n_a=1, q_TR='R', layer=None):
        return cf.calc_flux(self, wv, q_subs=q_subs, backside=backside, q_percent=q_percent, n_a=n_a, q_TR=q_TR,
                            layer=layer)


    def spectral_plot(self, *, q_TR='T', wv_range=[380, 760], N_pts=1000):
        fig = plt.figure('Spectral plot', figsize=(16, 9))
        ax = fig.add_subplot()
        x_range = np.linspace(wv_range[0], wv_range[1], N_pts)
        y_range = N_pts * [0.0]
        for i in range(N_pts):
            y_range[i] = self.calc_flux(Wave(x_range[i]), q_percent=True, q_TR=q_TR)
        ax.plot(x_range, y_range)
        plt.xlabel('Wavelength, nm')
        plt.ylabel('T, %')
        # plt.title('Design physical d')
        plt.show()


def Sellmeier_n(coef_A, wvlen):
    # wvlen in nm
    wv2 = wvlen ** 2
    pwr2n = (coef_A[0] + wv2 * (coef_A[1] / (wv2 - coef_A[2]) +
                                coef_A[3] / (wv2 - coef_A[4]) +
                                coef_A[5] / (wv2 - coef_A[6])))
    return np.sqrt(pwr2n)


def Cauchy_n(coef_A, wvlen):
    # wvlen in nm
    return coef_A[0] + coef_A[1] / wvlen ** 2 + coef_A[2] / wvlen ** 4
