#### Definição do Frame Work ####

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
    def ft3_h_to_m3_s(cls, V):
        return V * 7.8658e-6

    @classmethod
    def calc_prop_error(cls, prop_res, prop_1, prop_2, fluido="R134a"):

        nome_1, valor_1, erro_1 = prop_1
        nome_2, valor_2, erro_2 = prop_2

        res_1 = PropsSI(
            prop_res, nome_1, valor_1 + erro_1, nome_2, valor_2 + erro_2, fluido
        )
        res_2 = PropsSI(
            prop_res, nome_1, valor_1 + erro_1, nome_2, valor_2 - erro_2, fluido
        )
        res_3 = PropsSI(
            prop_res, nome_1, valor_1 - erro_1, nome_2, valor_2 + erro_2, fluido
        )
        res_4 = PropsSI(
            prop_res, nome_1, valor_1 - erro_1, nome_2, valor_2 - erro_2, fluido
        )

        aux = [res_1, res_2, res_3, res_4]

        res_min = np.min(aux)
        res_max = np.max(aux)

        med = PropsSI(prop_res, nome_1, valor_1, nome_2, valor_2, fluido)

        erro_min = ((med - res_min) ** 2) ** 0.5
        erro_max = ((med - res_max) ** 2) ** 0.5

        if np.isclose(erro_min, erro_max):
            return {"med": med, "erro": erro_min}

        erro = np.min((erro_max, erro_min))

        return {"med": med, "erro": erro}

    @classmethod
    def format_error(cls, prop, pow=0, decimal=1):
        valor = prop["med"]
        erro = prop["erro"]
        valor = np.around(valor * 10 ** pow, decimal)
        erro = np.around(erro * 10 ** pow, decimal)
        return str(valor) + " +/- " + str(erro)

    def __init__(
        self, Temp_dict, Press_dict, I, V, fluido="R134a", cool=True, correcao=True
    ):
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
        self.V = self.ft3_h_to_m3_s(V)
        self.fluido = fluido

        # definindo erros fixos

        self.Temp_error = 0.5
        self.Press_error = 0.5

        # determinando entalpias
        if cool:
            self.Hs1 = self.calc_prop_error(
                "H", ("P", self.Ps1, self.Press_error), ("T", self.Ts1, self.Temp_error)
            )
            self.Hs3 = self.calc_prop_error(
                "H", ("P", self.Ps3, self.Press_error), ("T", self.Ts3, self.Temp_error)
            )
            self.Hs4 = self.Hs3

            if correcao:
                S1 = self.calc_prop_error(
                    "S",
                    ("P", self.Ps1, self.Press_error),
                    ("T", self.Ts1, self.Temp_error),
                )
                self.Hs2 = self.calc_prop_error(
                    "H", ("P", self.Ps2, self.Press_error), ("S", S1["med"], S1["erro"])
                )
                self.Ts2 = PropsSI("T", "P", self.Ps2, "S", S1["med"], self.fluido)

                self.Ps4 = self.Ps1

            else:
                self.Hs2 = self.calc_prop_error(
                    "H",
                    ("P", self.Ps2, self.Press_error),
                    ("T", self.Ts2, self.Temp_error),
                )

            # determinando densidades

            self.rho1 = self.calc_prop_error(
                "D",
                ("P", self.Ps1, self.Press_error),
                ("H", self.Hs1["med"], self.Hs1["erro"]),
            )
            self.rho2 = self.calc_prop_error(
                "D",
                ("P", self.Ps2, self.Press_error),
                ("H", self.Hs2["med"], self.Hs2["erro"]),
            )
            self.rho3 = self.calc_prop_error(
                "D",
                ("P", self.Ps3, self.Press_error),
                ("H", self.Hs3["med"], self.Hs3["erro"]),
            )
            self.rho4 = self.calc_prop_error(
                "D",
                ("P", self.Ps4, self.Press_error),
                ("H", self.Hs4["med"], self.Hs4["erro"]),
            )

            # determinando temperaturas de saturação

            self.Ts1_sat = self.calc_prop_error(
                "T", ("P", self.Ps1, self.Press_error), ("Q", 1, 0)
            )
            self.Ts2_sat = self.calc_prop_error(
                "T", ("P", self.Ps2, self.Press_error), ("Q", 1, 0)
            )
            self.Ts3_sat = self.calc_prop_error(
                "T", ("P", self.Ps3, self.Press_error), ("Q", 1, 0)
            )
            self.Ts4_sat = self.calc_prop_error(
                "T", ("P", self.Ps4, self.Press_error), ("Q", 1, 0)
            )

            # determinando calor especifico

            self.c1 = self.calc_prop_error(
                "C", ("P", self.Ps1, self.Press_error), ("T", self.Ts1, self.Temp_error)
            )
            self.c2 = self.calc_prop_error(
                "C", ("P", self.Ps2, self.Press_error), ("T", self.Ts2, self.Temp_error)
            )
            self.c3 = self.calc_prop_error(
                "C", ("P", self.Ps3, self.Press_error), ("T", self.Ts3, self.Temp_error)
            )
            self.c4 = self.calc_prop_error(
                "C", ("P", self.Ps4, self.Press_error), ("T", self.Ts4, self.Temp_error)
            )

            input_array = np.array(
                (
                    [self.Ts1, self.Ps1, self.Hs1["med"], self.rho1["med"], self.c1["med"], 1],
                    [self.Ts2, self.Ps2, self.Hs2["med"], self.rho2["med"], self.c2["med"], 2],
                    [self.Ts3, self.Ps3, self.Hs3["med"], self.rho3["med"], self.c3["med"], 3],
                    [self.Ts4, self.Ps4, self.Hs4["med"], self.rho4["med"], self.c4["med"], 4],
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
                    "Ponto na bancada",
                ],
            )

            conv_input_array = np.array(
                (
                    [
                        self.format_error(
                            {"med": self.kelvin_to_celcius(self.Ts1), "erro": 0.5,}
                        ),
                        self.format_error(
                            {"med": self.pascal_to_bar(self.Ps1), "erro": 0.5,}
                        ),
                        self.format_error(self.Hs1, -3),
                        self.format_error(self.rho1),
                        self.format_error(self.c1),
                        self.format_error(self.Ts1_sat, decimal=5),
                    ],
                    [
                        self.format_error(
                            {"med": self.kelvin_to_celcius(self.Ts2), "erro": 0.5,}
                        ),
                        self.format_error(
                            {"med": self.pascal_to_bar(self.Ps2), "erro": 0.5,}
                        ),
                        self.format_error(self.Hs2, -3),
                        self.format_error(self.rho2),
                        self.format_error(self.c2),
                        self.format_error(self.Ts2_sat, decimal=5),
                    ],
                    [
                        self.format_error(
                            {"med": self.kelvin_to_celcius(self.Ts3), "erro": 0.5,}
                        ),
                        self.format_error(
                            {"med": self.pascal_to_bar(self.Ps3), "erro": 0.5,}
                        ),
                        self.format_error(self.Hs3, -3),
                        self.format_error(self.rho3),
                        self.format_error(self.c3),
                        self.format_error(self.Ts3_sat, decimal=5),
                    ],
                    [
                        self.format_error(
                            {"med": self.kelvin_to_celcius(self.Ts4), "erro": 0.5,}
                        ),
                        self.format_error(
                            {"med": self.pascal_to_bar(self.Ps4), "erro": 0.5,}
                        ),
                        self.format_error(self.Hs4, -3),
                        self.format_error(self.rho4),
                        self.format_error(self.c4),
                        self.format_error(self.Ts4_sat, decimal=5),
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
                        self.format_error(self.kelvin_to_celcius(self.Ts2), 1),
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
            self.delat_P_val = self.Ps4 - self.Ps3

        else:
            self.delta_T_eva = self.Ts2 - self.Ts3
            self.delta_T_cond = self.Ts4 - self.Ts1

            self.delta_P_eva = self.Ps2 - self.Ps3
            self.delta_P_cond = self.Ps4 - self.Ps1

    def calc_COP(self):

        if self.cool:

            COP = (self.Hs1["med"] - self.Hs4["med"]) / (
                self.Hs2["med"] - self.Hs1["med"]
            )
            erro = (
                1
                / (self.Hs2["med"] - self.Hs1["med"]) ** 2
                * (
                    (self.Hs2["med"] - self.Hs1["med"]) ** 2
                    * ((self.Hs1["erro"] + self.Hs4["erro"]))
                    + (self.Hs1["med"] - self.Hs4["med"]) ** 2
                    * ((self.Hs1["erro"] + self.Hs2["erro"]))
                )
                ** 0.5
            )
            self.COP = {
                "med": COP,
                "erro": erro,
            }

        else:
            self.COP = (self.Hs4 - self.Hs1) / (self.Hs2 - self.Hs1)

    def calc_ef_isentropica(self):
        self.Ss1 = PropsSI("S", "T", self.Ts1, "P", self.Ps1, self.fluido)
        self.Hs2s = PropsSI("H", "S", self.Ss1, "P", self.Ps2, self.fluido)

        self.ef_isentropica = (self.Hs2s - self.Hs1) / (self.Hs2 - self.Hs1)

    def calc_vasao_massica(self):

        res_1 = self.V * (self.rho3["med"] + self.rho3["erro"])
        res_2 = self.V * (self.rho3["med"] - self.rho3["erro"])

        aux = [res_1, res_2]

        res_min = np.min(aux)
        res_max = np.max(aux)

        med = self.V * self.rho3["med"]

        erro_min = ((med - res_min) ** 2) ** 0.5
        erro_max = ((med - res_max) ** 2) ** 0.5

        if np.isclose(erro_min, erro_max):
            self.dot_m = {"med": med, "erro": erro_min}

        else:
            erro = np.min((erro_max, erro_min))

            self.dot_m = {"med": med, "erro": erro}

    def calc_trab_compressao(self):
        if self.cool:
            self.trab_compressao = ()

        else:
            self.trab_compressao = (
                (self.dot_m2 + self.dot_m1) / 2 * (self.Hs1 - self.Hs2)
            )

    def calc_calor_evap(self):
        if self.cool:
            calor_evap = self.dot_m["med"] * (self.Hs1["med"] - self.Hs4["med"])
            erro = (
                self.dot_m["med"] ** 2 * (self.Hs1["erro"] ** 2 + self.Hs4["erro"] ** 2)
                + (self.Hs1["med"] - self.Hs4["med"]) ** 2 * self.dot_m["erro"] ** 2
            ) ** 0.5
            self.calor_evap = {
                "med": calor_evap,
                "erro": erro,
            }

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
        self.calc_vasao_massica()

        # self.calc_ef_isentropica()

        # self.calc_potencia_compressor()

        # prop1 = {"name": "T", "value": self.Ts3}

        # prop2 = {"name": "P", "value": self.Ps3}

        # self.calc_trab_compressao()

        self.calc_calor_evap()

        self.calc_sub_resfriamento()

        self.calc_super_aquecimento()

        if self.cool:
            output_array = np.array(
                (
                    [
                        self.format_error({"med": self.delta_T_eva, "erro": 0.7}),
                        self.format_error({"med": self.delta_T_cond, "erro": 0.7}),
                        self.format_error(
                            {"med": self.pascal_to_bar(self.delta_P_eva), "erro": 0.7}
                        ),
                        self.format_error(
                            {"med": self.pascal_to_bar(self.delta_P_cond), "erro": 0.7}
                        ),
                        self.format_error(self.COP, decimal=4),
                        self.format_error(self.dot_m, decimal=5),
                        self.format_error(self.calor_evap, pow=-3, decimal=4),
                        self.format_error({"med": self.sub_ref, "erro": 0.7}),
                        self.format_error({"med": self.sup_aqc, "erro": 0.7}),
                        # self.format_error(self.pot_comp),
                        # self.format_error(self.ef_isentropica),
                        # self.format_error(self.trab_compressao, 4),
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
                        "vazao_massica",
                        "calor_evap",
                        "sub_resfriamento",
                        "sup_aquecimento",
                        # "pot_comp",
                        # "ef_isentropica",
                        # "trab_compressao",
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
            fig.suptitle("COP em função do subresfriamento", fontsize=24)

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
            fig.suptitle(
                "Calor de evaporação em função da pressão do evaporador", fontsize=24
            )

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
