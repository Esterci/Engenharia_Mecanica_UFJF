import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
from sqlalchemy import false
import time


class Goodman:
    def __init__(
        self, Su_vida_inf, f_plast_prev, F_min, F_max, tipo="normal", f_cor=0.8
    ) -> None:
        self.tipo = tipo
        self.Su_vida_inf = Su_vida_inf
        self.f_cor = f_cor
        self.f_plast_prev = f_plast_prev
        self.F_min = F_min
        self.F_max = F_max
        pass

    # plotagem diagrama de goodman

    def show(self):
        if self.tipo == "modificado":

            fig = plt.figure(figsize=[16, 9])
            fig.suptitle("Diagrama de Goodman modificado", fontsize=16)

            linha_sobrecarga_1 = np.array(
                ([0, self.Su_vida_inf], [self.f_cor, self.f_cor])
            )

            linha_sobrecarga_2 = np.array(
                ([0, self.f_plast_prev], [self.f_cor, self.f_plast_prev])
            )

            linha_carga = np.array(
                ([0, 0], [self.f_cor * (self.F_min / self.F_max), self.f_cor])
            )

            # Plotando 2D

            ax = fig.add_subplot(1, 1, 1)

            ax.plot(
                linha_sobrecarga_1[:, 0], linha_sobrecarga_1[:, 1], "r", linewidth=2
            )
            ax.plot(
                linha_sobrecarga_2[:, 0], linha_sobrecarga_2[:, 1], "r", linewidth=2
            )
            ax.plot(linha_carga[:, 0], linha_carga[:, 1], "b", linewidth=2)

            ax.set_xlim(0)
            ax.set_ylim(0)
            ax.set_ylabel("$\sigma_{min}/S_u$", fontsize=16)
            ax.set_xlabel("$\sigma_{max}/S_u$", fontsize=16)
            ax.locator_params(axis="y", nbins=20)
            ax.grid()

            plt.show()


def C_converge(C_t, d, tau, F_max, verbose=False):

    k = 0

    erro_antigo = 10000

    while k < 100:

        Ks = (4 * C_t - 1) / (4 * C_t - 4) + 0.675 / C_t

        C_p = tau * np.pi * d ** 2 / (8 * F_max * Ks)

        error = ((C_p - C_t) ** 2) ** 0.5

        C_t = (C_p + C_t) / 2

        if error <= C_t * 0.01:
            return C_t

        elif error > erro_antigo:
            return C_t

        erro_antigo = error

        if verbose:
            print(C_t)
            print(C_p)
            print(error)
            print("*" * 20)
            time.sleep(0.25)


def deg_to_rad(deg):
    return deg * np.pi / (180)


def embreagem_hiperespaco(r_pinhao, r_i_lim ):

    r_i_min, r_i_max = r_i_lim
    r_i_eval = np.linspace(r_i_min, r_i_max, num=15)

    r_i_list = []
    r_o_list = []
    N_list = []
    P_list = []

    for r_i in r_i_eval:

        # constantes embragem

        Pot = 7.5e3
        w = 900 / 60
        torque = Pot / (w * 2 * np.pi)
        f = 0.45
        p_max = 2070e3
        b = 3e-3
        rho_emb = 79

        # calculo de embreagem

        N = 2

        r_o = r_i * (3) ** 0.5

        while r_o > 2 * r_pinhao:

            r_o = (torque / (np.pi * p_max * f * r_i * N) + r_i ** 2) ** 0.5

            N += 1

        F = 2 * np.pi * p_max * r_i * (r_o - r_i)

        P = (np.pi * r_o ** 2 - np.pi * r_i ** 2) * b * rho_emb * N * 2 - 1

        r_i_list.append(r_i)
        r_o_list.append(r_o)
        N_list.append(N)

    df = pd.DataFrame(
        np.vstack((r_i_list, r_o_list, N_list,)).T, columns=["r_i", "r_o", "N",]
    )

    return df


# n_molas_min, n_molas_max = n_molas_lim
#     d_min,d_max = d_lim

#     n_molas_eval = np.linspace(n_molas_min,n_molas_max,num=15,dtype=int)
#     d_mm_list = np.linspace(d_min,d_max,num=15,dtype=int)
#  for n_molas in n_molas_eval:

#         ## constantes molas

#         Su = 1250e6
#         F_min = F / n_molas
#         F_max = F * 1.2
#         Fs = 1.1 * F_max
#         delta_p = 1e-3
#         G = 79e9
#         E = 207e9
#         rho = 7700

#         flambagem_seed = np.array((
#             [0,10000],
#             [5.2,0.75],
#             [10,0.1375],
#             [9,0.1625],
#             [6,0.3667],
#             [8,0.2],
#             [7,0.252]
#         ))


#         f2 = interp1d(
#             flambagem_seed[:,0],
#             flambagem_seed[:,1],
#             kind='cubic'
#         )

#         # calculo molas

#         # encontrando tensão limite

#         tensao_max_teorico = 0.65 * Su

#         f_seguranca = 0.1

#         f_mola_solida = 0.1

#         tensao_max = tensao_max_teorico/(1+f_seguranca)/(1+f_mola_solida)

#         # estabelecendo um diãmetro alvo

#         C_t = 5

#         Kw = (4*C_t - 1)/(4*C_t - 4) + 0.615/C_t

#         d_alvo = (8*F_max*C_t*Kw/(np.pi*tensao_max))**0.5

#         # determinando rigidez

#         k = (F_max-F_min)/delta_p


#         C_list = []
#         D_list = []
#         Nt_list = []
#         Lf_D_list = []
#         delta_Lf_list = []
#         V_list = []
#         P_list = []
#         fn_list = []
#         d_list = []
#         r_i_list = []
#         r_o_list = []
#         N_list = []
#         n_molas_list = []

#         for i,d_mm in enumerate(d_mm_list):

#             d = d_mm*1e-3

#             C = C_converge(C_t,d,tensao_max,F_max)

#             D = C * d

#             D_lim = 2 * np.pi * r_i / (n_molas + np.pi) - d

#             N = d * G/(8 * C**3 * k)

#             N_t = N + 2

#             Ls = N_t * d

#             delta = Fs/k

#             L_f = Ls + delta

#             Lf_D = L_f/D

#             ponto_curva_a = f2(Lf_D)

#             delta_Lf = delta/L_f

#             V = (np.pi * d **2 / 4) * 2*np.pi *D/2

#             P += V * rho * n_molas

#             fn = 353e3*d/(N*D**2) * 60

#             if D > D_lim:
#                 C_list.append("Circ")


#             elif ponto_curva_a < delta_Lf:
#                 C_list.append("Flm")
#                 break

#             else:
#                 C_list.append(C)

#             D_list.append(D)
#             Nt_list.append(N_t)
#             Lf_D_list.append(Lf_D)
#             delta_Lf_list.append(delta_Lf)
#             V_list.append(V)
#             P_list.append(P)
#             fn_list.append(fn)
#             d_list.append(d)
#             r_i_list.append(r_i)
#             r_o_list.append(r_o)
#             N_list.append(N)
#             n_molas_list.append(n_molas)

#             print(i)

# df = pd.DataFrame(np.vstack((
#     d_list,
#     C_list,
#     D_list,
#     Nt_list,
#     Lf_D_list,
#     delta_Lf_list,
#     V_list,
#     P_list,
#     fn_list,
#     r_i_list,
#     r_o_list,
#     N_list,
#     n_molas_list
# )).T,columns=[
#         'd',
#         'C',
#         'D',
#         'Nt',
#         'Lf/D',
#         'delta/Lf',
#         'V',
#         'P',
#         'fn',
#         'r_i',
#         'r_o',
#         'N',
#         'n_molas'
#     ])

# return df
