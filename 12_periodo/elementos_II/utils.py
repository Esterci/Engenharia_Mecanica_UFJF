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


class embreagem:
    def __init__(self, r_pinhao, n_molas_lim = (4,8)) -> None:
        self.r_pinhao = r_pinhao
        self.n_molas_min, self.n_molas_max = n_molas_lim
        pass

    def calc_embreagem(self):

        r_i_eval = np.linspace(10, self.r_pinhao, num=10, dtype=int )
        r_i_list = []
        r_o_list = []
        N_list = []
        P_list = []
        F_emb_list = []
        torque_list = []
        Perim_list = []

        for r_i_mm in r_i_eval:

            r_i = r_i_mm * 1e-3

            # constantes embragem

            Pot = 7.5e3
            w = 900 / 60
            torque_proj = Pot / (w * 2 * np.pi)
            f = 0.45
            p_max = 2070e3
            b = 3e-3
            rho_emb = 7070
            r_pinhao = 75
            torque_proj

            # calculo de embreagem

            N = 2

            r_o = r_i * (3) ** 0.5

            torque = np.pi * p_max * f * r_i * (r_o ** 2 - r_i ** 2) * N

            while torque > torque_proj:

                torque = np.pi * p_max * f * r_i * (r_o ** 2 - r_i ** 2) * N

                N += 1

                if N > 10:
                    break

            F_emb = 2 * np.pi * p_max * r_i * (r_o - r_i)

            P = (np.pi * r_o ** 2 - np.pi * r_i ** 2) * b * rho_emb * (N * 2 - 1)

            Perim_emb = 2 * np.pi * ((r_o - r_i) / 2 + r_i)

            r_i_list.append(r_i)
            r_o_list.append(r_o)
            N_list.append(N)
            P_list.append(P)
            F_emb_list.append(F_emb)
            Perim_list.append(2*np.pi*((r_o-r_i)/2 + r_i))

        self.df_emb = pd.DataFrame(
            np.vstack((
                r_i_list, 
                r_o_list, 
                N_list,
                P_list, 
                F_emb_list,
                torque_list,
                Perim_list,
                )).T, 
                columns=["r_i", "r_o", "N","P","F","torque","Perim"]
        )

    def plot_peso_embreagem(self):
        
        fig = plt.figure(figsize=[16, 9])
        fig.suptitle('Peso da embreagem em função do raio interno', fontsize=16)

        ax = fig.add_subplot(1,1,1)

        # Plotando 2D

        ax.plot(self.df_emb.r_i,
            self.df_emb.P,
            'o', 
            self.df_emb.r_i, 
            self.df_emb.P, 
            '-',
            color='b',
        )

        ax.locator_params(axis='y', nbins=30)
        ax.locator_params(axis='x', nbins=20)
        ax.set_ylabel('$Peso (kg)$', fontsize=16)
        ax.set_xlabel('$r_i (m)$', fontsize=16)
        ax.grid()

        plt.savefig("plots/Peso_da_embreagem_vs_ri.png")

    def calc_mola(self):

        # constantes
        Su = 1250e6
        Su_lim = 0.80
        Su_vida_inf = 0.62
        Sy_lim = 0.65
        delta_p = 1e-3
        G = 79e9
        E = 207e9
        rho = 7700

        tensao_max_teorico = 0.65 * Su

        f_seguranca = 0.1

        f_mola_solida = 0.1

        tensao_max = tensao_max_teorico/(1+f_seguranca)/(1+f_mola_solida)

        n_molas_list = np.linspace(self.n_molas_min,self.n_molas_max,num=5,dtype=int)
        d_mm_list = np.linspace(5,18,num=6,dtype=float)


        for i in range(len(self.df_emb)):

            pd_serie = self.df_emb.loc[i]
            
            for j,n_molas in enumerate(n_molas_list):

                F = pd_serie.F / n_molas

                F_max = pd_serie.F * 1.2

                Fs = 1.1 * F_max
                
                if i ==0 and j == 0:
                    self.plot_goodman(self, pd_serie.F, F_max, Su_vida_inf,Su_lim,Sy_lim)

                for i,d_mm in enumerate(d_mm_list):

                    d = d_mm*1e-3

                    C = C_converge(C_t,d,tensao_max,F_max)

                    D = C * d

                    N = d * G/(8 * C**3 * k)

                    N_t = N + 2

                    Ls = N_t * d

                    delta = Fs/k

                    L_f = Ls + delta

                    Lf_D = L_f/D


                    delta_Lf = delta/L_f


                    V = (np.pi * d **2 / 4) * 2*np.pi *D/2

                    P += V * rho * n_molas

                    fn = 353e3*d/(N*D**2) * 60

                    C_list.append(C)
                    D_list.append(D)
                    Nt_list.append(N_t)
                    Lf_D_list.append(Lf_D)
                    delta_Lf_list.append(delta_Lf)
                    V_list.append(V)
                    P_list.append(P)
                    fn_list.append(fn)
                    d_list.append(d)
                    N_list.append(n_molas)
                    Diam_list.append(D*n_molas)

            df = pd.DataFrame(np.vstack((
                    d_list,
                    C_list,
                    N_list,
                    Diam_list,
                    D_list,
                    Nt_list,
                    Lf_D_list,
                    delta_Lf_list,
                    P_list,
                    fn_list,
                )).T,columns=[
                        'd',
                        'C',
                        'n_molas',
                        'Diam',
                        'D',
                        'Nt',
                        'Lf/D',
                        'delta/Lf',
                        'P',
                        'fn',
                    ], dtype=float)

            df

    def plot_goodman(self, F_min, F_max, Su_vida_inf,Su_lim,Sy_lim):
        
        # plotagem diagrama de goodman

        fig = plt.figure(figsize=[16, 25])

        linha_sobrecarga_1 = np.array(
            ([0,Su_vida_inf],
            [Su_lim,Su_lim])
        ) 

        linha_sobrecarga_2 = np.array(
            ([0,Sy_lim],
            [Su_lim,Sy_lim])
        ) 

        linha_carga = np.array(
            ([0,0],
            [Su_lim*(F_min/F_max),Su_lim])
        )


        # Plotando 2D

        ax = fig.add_subplot(3, 2, (i+1) )

        ax.set_title('Diagrama de Goodman', fontsize=16)

        ax.plot(linha_sobrecarga_1[:,0], linha_sobrecarga_1[:,1], 'r', linewidth=2)
        ax.plot(linha_sobrecarga_2[:,0], linha_sobrecarga_2[:,1], 'r', linewidth=2)
        ax.plot(linha_carga[:,0], linha_carga[:,1], 'b', linewidth=2)

        ax.set_xlim(0)
        ax.set_ylim(0)
        ax.locator_params(axis='y', nbins=20)
        ax.grid()

        plt.savefig("plots/Diagrama_de_Goodman.png")

    def estudo (self):

        self.calc_embreagem()

        self.plot_peso_embreagem()


