#### Definição do Frame Work ####

from types import ClassMethodDescriptorType
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from CoolProp.Plots import PropertyPlot


class Experimento:

    # definindo funções de conversão
    @classmethod
    def celcius_to_kelvin(cls, T):
        return T + 273.15

    @classmethod
    def psi_to_pascal(cls, P, manometro=True):
        if manometro:
            return P * 6894 + 101325
        return P * 6894

    @classmethod
    def bar_to_pascal(cls, P, manometro=True):
        if manometro:
            return P * 100000 + 101325
        return P * 100000

    @classmethod
    def kelvin_to_celcius(cls, T):
        return T - 273.15
    
    @classmethod
    def pascal_to_psi(cls, P):
        return P / 6894

    @classmethod
    def pascal_to_bar(cls, P):
        return P / 100000

    @classmethod
    def vasao_massica(cls, V, rho):
        return V * rho

    def __init__(self, Temp_dict, Press_dict, I, V, fluido, cool=True, correcao=False):
        # passando valores de temperatura
        self.Ts1 = self.celcius_to_kelvin(Temp_dict["Ts1"])
        self.Ts2 = self.celcius_to_kelvin(Temp_dict["Ts2"])
        self.Ts3 = self.celcius_to_kelvin(Temp_dict["Ts3"])
        self.Ts4 = self.celcius_to_kelvin(Temp_dict["Ts4"])

        # passando valores de pressão
        self.Ps1 = self.bar_to_pascal(Press_dict["Ps1"])
        self.Ps2 = self.bar_to_pascal(Press_dict["Ps2"])
        self.Ps3 = self.bar_to_pascal(Press_dict["Ps3"])
        self.Ps4 = self.bar_to_pascal(Press_dict["Ps4"])

        # passando demais valores
        self.cool = cool
        self.I = I
        self.V = V
        self.fluido = fluido

        # determinando entalpias
        if cool:
            self.Hs1 = PropsSI("H", "P", self.Ps1, "T", self.Ts1, self.fluido)
            self.Hs3 = PropsSI("H", "P", self.Ps3, "T", self.Ts3, self.fluido)
            self.Hs4 = self.Hs3

            if correcao:
                S1 = PropsSI("S", "P", self.Ps1, "T", self.Ts1, self.fluido)
                self.Hs2 = PropsSI("H", "P", self.Ps2, "S", S1, self.fluido)
                self.Ts2 = PropsSI("T", "P", self.Ps2, "S", S1, self.fluido)
                self.Ps4 = self.Ps1

            else:
                self.Hs2 = PropsSI("H", "P", self.Ps2, "T", self.Ts2, self.fluido)

            # determinando densidades

            self.rho1 = PropsSI("D", "P", self.Ps1, "H", self.Hs1, self.fluido)
            self.rho2 = PropsSI("D", "P", self.Ps2, "H", self.Hs2, self.fluido)
            self.rho3 = PropsSI("D", "P", self.Ps3, "H", self.Hs3, self.fluido)
            self.rho4 = PropsSI("D", "P", self.Ps4, "H", self.Hs4, self.fluido)

            # determinando calor especifico

            self.c1 = PropsSI("C", "P", self.Ps1, "T", self.Ts1, self.fluido)
            self.c2 = PropsSI("C", "P", self.Ps2, "T", self.Ts2, self.fluido)
            self.c3 = PropsSI("C", "P", self.Ps3, "T", self.Ts3, self.fluido)
            self.c4 = PropsSI("C", "P", self.Ps4, "T", self.Ts4, self.fluido)

            # determinando vasao massica

            self.dot_m1 = self.vasao_massica(V,self.rho1)
            self.dot_m2 = self.vasao_massica(V,self.rho2)
            self.dot_m3 = self.vasao_massica(V,self.rho3)
            self.dot_m4 = self.vasao_massica(V,self.rho4)

            input_array = np.array(
                (
                    [self.Ts1, self.Ps1, self.Hs1, self.rho1, self.c1, self.dot_m1, 1],
                    [self.Ts2, self.Ps2, self.Hs2, self.rho2, self.c2, self.dot_m2, 2],
                    [self.Ts3, self.Ps3, self.Hs3, self.rho3, self.c3, self.dot_m3, 3],
                    [self.Ts4, self.Ps4, self.Hs4, self.rho4, self.c4, self.dot_m4, 4],
                )
            )

            self.input_df = pd.DataFrame(
                input_array,
                columns=[
                    "Temperatura",
                    "Pressao",
                    "Entalpia",
                    "Densidade",
                    "Calor Esp.",
                    "Vasao Massica",
                    "Ponto na bancada",
                ],
            )

            conv_input_array = np.array(
                (
                    [
                        np.around(self.kelvin_to_celcius(self.Ts1), 1),
                        np.around(self.pascal_to_bar(self.Ps1), 1),
                        np.around(self.Hs1 / 1e3, 1),
                        np.around(self.rho1, 1),
                        np.around(self.c1, 1),
                        1,
                        self.dot_m1,
                    ],
                    [
                        np.around(self.kelvin_to_celcius(self.Ts2), 1),
                        np.around(self.pascal_to_bar(self.Ps2), 1),
                        np.around(self.Hs2 / 1e3, 1),
                        np.around(self.rho2, 1),
                        np.around(self.c2, 1),
                        np.around(self.kelvin_to_celcius(PropsSI("T", "P", self.Ps2, "Q", 1, self.fluido))),
                        self.dot_m2,
                    ],
                    [
                        np.around(self.kelvin_to_celcius(self.Ts3), 1),
                        np.around(self.pascal_to_bar(self.Ps3), 1),
                        np.around(self.Hs3 / 1e3, 1),
                        np.around(self.rho3, 1),
                        np.around(self.c3, 1),
                        1,
                        self.dot_m3,
                    ],
                    [
                        np.around(self.kelvin_to_celcius(self.Ts4), 1),
                        np.around(self.pascal_to_bar(self.Ps4), 1),
                        np.around(self.Hs4 / 1e3, 1),
                        np.around(self.rho4, 1),
                        np.around(self.c4, 1),
                        1,
                        self.dot_m4,
                    ],
                )
            )

            self.conv_input_df = pd.DataFrame(
                conv_input_array,
                columns=[
                    "Temperatura",
                    "Pressao",
                    "Entalpia",
                    "Densidade",
                    "Calor Esp.",
                    "Temp. Sat.",
                    "Vasao massica"
                ],
            )

        else:
            self.Hs2 = PropsSI("H", "P", self.Ps2, "T", self.Ts2, self.fluido)
            self.Hs4 = PropsSI("H", "P", self.Ps4, "T", self.Ts4, self.fluido)
            self.Hs3 = self.Hs4

            if correcao:
                S2 = PropsSI("S", "P", self.Ps2, "T", self.Ts2, self.fluido)
                self.Hs1 = PropsSI("S", "P", self.Ps1, "S", S2, self.fluido)
                self.Ps4 = self.Ps1

            else:
                self.Hs1 = PropsSI("H", "P", self.Ps1, "T", self.Ts1, self.fluido)
                self.Ps4 = self.bar_to_pascal(Press_dict["Ps4"])

            input_array = np.array(
                (
                    [self.Ts2, self.Ps2, self.Hs2, 2],
                    [self.Ts1, self.Ps1, self.Hs1, 1],
                    [self.Ts4, self.Ps4, self.Hs4, 4],
                    [self.Ts3, self.Ps3, self.Hs3, 3],
                )
            )

            self.input_df = pd.DataFrame(
                input_array,
                columns=["Temperatura", "Pressao", "Entalpia", "Ponto na bancada"],
            )

            conv_input_array = np.array(
                (
                    [
                        np.around(self.kelvin_to_celcius(self.Ts2), 1),
                        np.around(self.pascal_to_bar(self.Ps2), 1),
                        np.around(self.Hs2 / 1e3, 1),
                        2,
                    ],
                    [
                        np.around(self.kelvin_to_celcius(self.Ts1), 1),
                        np.around(self.pascal_to_bar(self.Ps1), 1),
                        np.around(self.Hs1 / 1e3, 1),
                        1,
                    ],
                    [
                        np.around(self.kelvin_to_celcius(self.Ts4), 1),
                        np.around(self.pascal_to_bar(self.Ps4), 1),
                        np.around(self.Hs4 / 1e3, 1),
                        4,
                    ],
                    [
                        np.around(self.kelvin_to_celcius(self.Ts3), 1),
                        np.around(self.pascal_to_bar(self.Ps3), 1),
                        np.around(self.Hs3 / 1e3, 1),
                        3,
                    ],
                )
            )

            self.conv_input_df = pd.DataFrame(
                conv_input_array,
                columns=["Temperatura", "Pressao", "Entalpia", "Ponto no bancada"],
            )

    def variacoes(self):

        if self.cool:
            self.delta_T_eva = self.Ts1 - self.Ts4
            self.delta_T_cond = self.Ts3 - self.Ts2

            self.delta_P_eva = self.Ps1 - self.Ps4
            self.delta_P_cond = self.Ps3 - self.Ps2

        else:
            self.delta_T_eva = self.Ts2 - self.Ts3
            self.delta_T_cond = self.Ts4 - self.Ts1

            self.delta_P_eva = self.Ps2 - self.Ps3
            self.delta_P_cond = self.Ps4 - self.Ps1

    def calc_COP(self):

        if self.cool:
            self.COP = (self.Hs1 - self.Hs4) / (self.Hs2 - self.Hs1)

        else:
            self.COP = (self.Hs4 - self.Hs1) / (self.Hs2 - self.Hs1)

    def calc_ef_isentropica(self):
        self.Ss1 = PropsSI("S", "T", self.Ts1, "P", self.Ps1, self.fluido)
        self.Hs2s = PropsSI("H", "S", self.Ss1, "P", self.Ps2, self.fluido)

        self.ef_isentropica = (self.Hs2s - self.Hs1) / (self.Hs2 - self.Hs1)

    def calc_trab_compressao(self):
        if self.cool:
            self.trab_compressao = (self.dot_m2 + self.dot_m1)/2 * (self.Hs2 - self.Hs1)

        else:
            self.trab_compressao = (self.dot_m2 + self.dot_m1)/2 * (self.Hs1 - self.Hs2)

    def calc_calor_evap(self):
        if self.cool:
            self.calor_evap = (self.dot_m1 + self.dot_m4)/2 * (self.Hs1 - self.Hs4)

    def calc_potencia_compressor(self):
        self.pot_comp = self.I * 127

    def calc_sub_resfriamento(self):
        if self.cool:
            T_liq_sat = PropsSI("T", "P", self.Ps3, "Q", 0, self.fluido)
            self.sub_ref = T_liq_sat - self.Ts3

        else:
            print("not implemented yet.")

    def calc_super_aquecimento(self):
        if self.cool:
            T_vap_sat = PropsSI("T", "P", self.Ps1, "Q", 1, self.fluido)
            self.sup_aqc = self.Ts1 - T_vap_sat

        else:
            print("not implemented yet.")

    def analise_defaut(self):

        self.variacoes()

        self.calc_COP()

        self.calc_ef_isentropica()

        self.calc_potencia_compressor()

        prop1 = {"name": "T", "value": self.Ts3}

        prop2 = {"name": "P", "value": self.Ps3}

        self.calc_trab_compressao()

        self.calc_sub_resfriamento()

        self.calc_super_aquecimento()

        if self.cool:
            output_array = np.array(
                (
                    [
                        np.around(self.delta_T_eva, 1),
                        np.around(self.delta_T_cond, 1),
                        np.around(self.pascal_to_bar(self.delta_P_eva), 1),
                        np.around(self.pascal_to_bar(self.delta_P_cond), 1),
                        np.around(self.COP, 1),
                        np.around(self.pot_comp, 1),
                        np.around(self.sub_ref, 1),
                        np.around(self.sup_aqc, 1),
                        np.around(self.ef_isentropica, 1),
                        np.around(self.trab_compressao, 4),
                    ]
                )
            )

            variables_list = np.array(
                (
                    [
                        "delta_T_eva",
                        "delta_T_cond",
                        "delta_P_eva",
                        "delta_P_cond",
                        "COP",
                        "pot_comp",
                        "sub_resfriamento",
                        "sup_aquecimento",
                        "ef_isentropica",
                        "trab_compressao",
                    ]
                )
            )

            self.output_df = pd.DataFrame(
                np.hstack(
                    (
                        variables_list.reshape(len(variables_list), 1),
                        output_array.reshape(len(output_array), 1),
                    )
                ),
                columns=["Propriedade", "Valor"],
            )

        else:
            output_array = np.array(
                (
                    [
                        self.delta_T_eva,
                        self.delta_T_cond,
                        self.pascal_to_bar(self.delta_P_eva),
                        self.pascal_to_bar(self.delta_P_cond),
                        self.COP,
                        self.ef_isentropica,
                        self.trab_compressao,
                        self.pot_comp,
                    ]
                )
            )

            self.output_df = pd.DataFrame(
                output_array.reshape(1, len(output_array)),
                columns=[
                    "delta_T_eva",
                    "delta_T_cond",
                    "delta_P_eva",
                    "delta_P_cond",
                    "COP",
                    "ef_isentropica",
                    "trab_compressao",
                    "pot_comp",
                ],
            )

    def plot_diagrama_ph(self):
        ph_plot = PropertyPlot(self.fluido, "Ph")
        ph_plot.calc_isolines()
        fig = ph_plot.figure
        fig.set_size_inches(16, 9)
        ax = ph_plot.axis

        # definindo listas para fechar o ciclo
        plot_entalpia = self.input_df.Entalpia.to_list()
        plot_pressao = self.input_df.Pressao.to_list()

        plot_entalpia.append(self.input_df.Entalpia[0])
        plot_pressao.append(self.input_df.Pressao[0])

        plot_entalpia = [h / 10 ** 3 for h in plot_entalpia]
        plot_pressao = [h / 10 ** 3 for h in plot_pressao]

        ax.plot(
            plot_entalpia,
            plot_pressao,
            "o",
            plot_entalpia,
            plot_pressao,
            "-",
            color="b",
        )

        ax.tick_params(axis="x", labelsize=24)
        ax.tick_params(axis="y", labelsize=24)
        ax.set_xlabel("Entalpia específica(kJ/kg)", fontsize=24)
        ax.set_ylabel("Pressão (kPa)", fontsize=24)

        # anotando os pontos no diagrama

        if self.cool:
            ax.annotate(
                "S1",
                (
                    self.input_df.Entalpia[0] / 1e3 + 5,
                    self.input_df.Pressao[0] / 1e3 + 5,
                ),
                size=24,
            )
            ax.annotate(
                "S2",
                (
                    self.input_df.Entalpia[1] / 1e3 + 5,
                    self.input_df.Pressao[1] / 1e3 + 5,
                ),
                size=24,
            )
            ax.annotate(
                "S3",
                (
                    self.input_df.Entalpia[2] / 1e3 + 5,
                    self.input_df.Pressao[2] / 1e3 + 5,
                ),
                size=24,
            )
            ax.annotate(
                "S4",
                (
                    self.input_df.Entalpia[3] / 1e3 + 5,
                    self.input_df.Pressao[3] / 1e3 + 5,
                ),
                size=24,
            )
            ph_plot.show()

        else:
            ax.annotate(
                "S2",
                (
                    self.input_df.Entalpia[0] / 1e3 + 5,
                    self.input_df.Pressao[0] / 1e3 + 5,
                ),
                size=24,
            )
            ax.annotate(
                "S1",
                (
                    self.input_df.Entalpia[1] / 1e3 + 5,
                    self.input_df.Pressao[1] / 1e3 + 5,
                ),
                size=24,
            )
            ax.annotate(
                "S4",
                (
                    self.input_df.Entalpia[2] / 1e3 + 5,
                    self.input_df.Pressao[2] / 1e3 + 5,
                ),
                size=24,
            )
            ax.annotate(
                "S3",
                (
                    self.input_df.Entalpia[3] / 1e3 + 5,
                    self.input_df.Pressao[3] / 1e3 + 5,
                ),
                size=24,
            )
            ph_plot.show()

    def plot_P_saida_crompressor_vs_COP(self, min, max):
        if self.cool:
            P_list = np.linspace(min, max, 250)
            original_Hs2 = self.Hs2
            S_ref = PropsSI("S", "T", self.Ts2, "P", self.Ps2, self.fluido)
            results_list = []
            for P in P_list:
                Ts2 = PropsSI("T", "S", S_ref, "P", P, self.fluido)
                self.Hs2 = PropsSI("H", "T", Ts2, "P", P, self.fluido)

                self.calc_COP()

                results_list.append(self.COP)

            # retornando valores originais

            self.Hs2 = original_Hs2
            self.calc_COP()

            # plotando resultados
            fig = plt.figure(figsize=[16, 9])

            # Plotando 2D

            ax = fig.add_subplot(1, 1, 1)

            ax.set_xlabel("$Ps2~[Pa]$", fontsize=16)
            ax.set_ylabel("$COP$", fontsize=16)

            ax.plot(P_list, results_list, "r", linewidth=2)

            ax.grid()

            plt.show()

        else:
            P_list = np.linspace(min, max, 250)
            original_Hs1 = self.Hs1
            S_ref = PropsSI("S", "T", self.Ts1, "P", self.Ps1, self.fluido)
            results_list = []
            for P in P_list:
                Ts1 = PropsSI("T", "S", S_ref, "P", P, self.fluido)
                self.Hs1 = PropsSI("H", "T", Ts1, "P", P, self.fluido)

                self.calc_COP()

                results_list.append(self.COP)

            # retornando valores originais

            self.Hs1 = original_Hs1
            self.calc_COP()

            # plotando resultados
            fig = plt.figure(figsize=[16, 9])

            # Plotando 2D

            ax = fig.add_subplot(1, 1, 1)

            ax.set_xlabel("$Ps1~[Pa]$", fontsize=16)
            ax.set_ylabel("$COP$", fontsize=16)

            ax.plot(P_list, results_list, "r", linewidth=2)

            ax.grid()

            plt.show()

    def plot_sub_resfriamento_vs_COP(self, min, max):
        if self.cool:
            sub_res_list = np.linspace(min, max, 250)
            original_Hs2 = self.Hs2
            results_list = []
            T_liq_sat = PropsSI("T", "P", self.Ps3, "Q", 0, self.fluido)
            
            for sub_res in sub_res_list:
                Ts2 = T_liq_sat + sub_res
                self.Hs2 = PropsSI("H", "T", Ts2, "P", self.Ps2, self.fluido)

                self.calc_COP()

                results_list.append(self.COP)

            # retornando valores originais

            self.Hs2 = original_Hs2
            self.calc_COP()

            # plotando resultados
            fig = plt.figure(figsize=[16, 9])
            fig.suptitle('COP em função do subresfriamento', fontsize=24)

            # Plotando 2D

            ax = fig.add_subplot(1, 1, 1)

            ax.tick_params(axis="x", labelsize=24)
            ax.tick_params(axis="y", labelsize=24)

            ax.set_xlabel("$Subresfriamento~[^oC]$", fontsize=24)
            ax.set_ylabel("$COP$", fontsize=24)

            ax.plot(sub_res_list, results_list, "r", linewidth=2)

            ax.grid()

            plt.show()

        else:
            print("Not implemented yet")

    def plot_calor_evap_vs_P(self, min, max):
        if self.cool:
            P_list = np.linspace(min, max, 250)
            original_Hs1 = self.Hs1
            original_Hs4 = self.Hs4
            results_list = []
            
            for P in P_list:
                self.Hs1 = PropsSI("H", "T", self.Ts1, "P", P, self.fluido)
                self.Hs4 = PropsSI("H", "T", self.Ts4, "P", P, self.fluido)

                self.analise_defaut()

                self.calc_calor_evap()

                results_list.append(self.calor_evap)

            # retornando valores originais

            self.Hs1 = original_Hs1
            self.Hs4 = original_Hs4
            self.analise_defaut()
            self.calc_calor_evap()

            # plotando resultados
            fig = plt.figure(figsize=[16, 9])
            fig.suptitle('Calor de evaporação em função da pressão do evaporador', fontsize=24)

            # Plotando 2D

            ax = fig.add_subplot(1, 1, 1)

            ax.tick_params(axis="x", labelsize=24)
            ax.tick_params(axis="y", labelsize=24)

            ax.set_xlabel("$Q (J)$", fontsize=24)
            ax.set_ylabel("$COP$", fontsize=24)

            ax.plot(P_list, results_list, "r", linewidth=2)

            ax.grid()

            plt.show()

        else:
            print("Not implemented yet")

    def plot_P_saida_crompressor_vs_trab_compressao(self, min, max):
        if self.cool:
            P_list = np.linspace(min, max, 250)
            original_Hs2 = self.Hs2
            S_ref = PropsSI("S", "T", self.Ts2, "P", self.Ps2, self.fluido)
            results_list = []
            for P in P_list:
                Ts2 = PropsSI("T", "S", S_ref, "P", P, self.fluido)
                self.Hs2 = PropsSI("H", "T", Ts2, "P", P, self.fluido)

                self.calc_trab_compressao()

                results_list.append(self.trab_compressao)

            # retornando valores originais

            self.Hs2 = original_Hs2
            self.calc_trab_compressao()

            # plotando resultados
            fig = plt.figure(figsize=[16, 9])

            # Plotando 2D

            ax = fig.add_subplot(1, 1, 1)

            ax.set_xlabel("$Ps2~[Pa]$", fontsize=16)
            ax.set_ylabel("$W_{comp}~[J/s]$", fontsize=16)

            ax.plot(P_list, results_list, "r", linewidth=2)

            ax.grid()

            plt.show()

        else:
            print("not implemented")
        return

    def plot_P_saida_crompressor_vs_ef_isentropica(self, min, max):
        if self.cool:
            P_list = np.linspace(min, max, 250)
            original_Hs2 = self.Hs2
            S_ref = PropsSI("S", "T", self.Ts2, "P", self.Ps2, self.fluido)
            results_list = []
            for P in P_list:
                Ts2 = PropsSI("T", "S", S_ref, "P", P, self.fluido)
                self.Hs2 = PropsSI("H", "T", Ts2, "P", P, self.fluido)

                self.calc_ef_isentropica()

                results_list.append(self.ef_isentropica)

            # retornando valores originais

            self.Hs2 = original_Hs2
            self.calc_ef_isentropica()

            # plotando resultados
            fig = plt.figure(figsize=[16, 9])

            # Plotando 2D

            ax = fig.add_subplot(1, 1, 1)

            ax.set_xlabel("$Ps2(Pa)$", fontsize=16)
            ax.set_ylabel("$\eta$", fontsize=16)

            ax.plot(P_list, results_list, "r", linewidth=2)

            ax.grid()

            plt.show()

        else:
            print("not implemented")
        return
