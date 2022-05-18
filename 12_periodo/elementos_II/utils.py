import numpy as np
import matplotlib.pyplot as plt


class Goodman:
    def __init__(
        self, Su_vida_inf, f_plast_prev, F_min, F_max, tipo="normal",f_cor=0.8
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

            linha_sobrecarga_2 = np.array(([0, self.f_plast_prev], [self.f_cor, self.f_plast_prev]))

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
            ax.set_ylabel('$\sigma_{min}/S_u$', fontsize=16)
            ax.set_xlabel('$\sigma_{max}/S_u$', fontsize=16)
            ax.locator_params(axis="y", nbins=20)
            ax.grid()

            plt.show()


def C_converge(C_t, d, tau, F_max):

    k = 0

    while k < 100:

        Ks = (4 * C_t - 1) / (4 * C_t - 4) + 0.675 / C_t

        C_p = tau * np.pi * d ** 2 / (8 * F_max * Ks)

        error = ((C_p - C_t) ** 2) ** 0.5

        C_t = (C_p + C_t) / 2

        if error <= C_t * 0.01:
            return C_t

def deg_to_rad(deg):
    return deg*np.pi/(180)