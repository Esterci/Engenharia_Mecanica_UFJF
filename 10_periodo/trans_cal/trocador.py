PROCEDURE temp_it (T_f_ent_1;T_f_sai_chut_5;T_f_sai_chut_6;T_q_ent_2;T_q_sai_chut_4;T_q_sai_chut_5;T_q_sai_chut_6;m_dot_f;m_dot_q;P_q;P_f : error_q_4 ;error_q_5 ; error_q_6 ; error_q_med ; T_it_4 ; T_it_5 ; T_it_6)

    REPEAT
        "############### Volume de Controle IV ##############"

        "parametros do fluido frio"

        T_f_ent_4 = T_f_ent_1
        T_f_sai_4 = 373,15
        T_med_f_4 = (T_f_ent_4 +T_f_sai_4) / 2
        c_p_f_4 = Cp (Water;T=T_med_f_4;P=P_f)
        I_f_ent_4 = Enthalpy(Water;T=T_f_ent_4;P=P_f)
        I_f_sai_4 = Enthalpy(Water;T=T_f_sai_4;P=P_f)


        "parametros do fluido quente"

        T_q_ent_4 = T_q_ent_2
        T_med_q_4 = (T_q_ent_4 + T_q_sai_chut_4)/2
        c_p_q_4 =Cp(Air_ha;T=T_med_q_4;P=P_q)

        "Balanço de massas e energias"
        
        T_q_sai_4 := T_q_ent_4 - m_dot_f*(I_f_sai_4-I_f_ent_4)/(m_dot_q * c_p_q_4)

        "############### Volume de Controle V ##############"

        "parametros do fluido frio"

        T_f_ent_5 = T_f_sai_4
        T_med_f_5 = (T_f_ent_5 + T_f_sai_chut_5) / 2

        " tirando média dos parâmetros termo-físicos"

        c_p_f_a = Cp (Water;T=T_med_f_5;P=P_f)
        I_f_ent_a = Enthalpy_vaporization(Water;T=T_f_ent_4)
        I_f_sai_a = Enthalpy_vaporization(Water;T=T_f_sai_4)

        c_p_f_v = Cp (Steam;T=T_med_f_5;P=P_f)
        I_f_ent_v = Enthalpy_vaporization(Steam;T=T_f_ent_4)
        I_f_sai_v = Enthalpy_vaporization(Steam;T=T_f_sai_4)

        c_p_f_5 = (c_p_f_v + c_p_f_a)/2
        I_f_ent_5 = (I_f_ent_v + I_f_ent_a)/2
        I_f_sai_5 = (I_f_sai_v + I_f_sai_a)/2

        "parametros do fluido quente"

        T_q_ent_5 = T_q_sai_4
        T_med_q_5 = (T_q_ent_5 + T_q_sai_chut_5)/2
        c_p_q_5 =Cp(Air_ha;T=T_med_q_5;P=P_q)

        "Balanço de massas e energias"

        T_q_sai_5 := T_q_ent_5 - m_dot_f*(I_f_sai_5-I_f_ent_5)/(m_dot_q * c_p_q_5)


        "############### Volume de Controle VI ##############"

        "parametros do fluido frio"

        T_f_ent_6 = T_f_sai_chut_5
        T_med_f_6 = (T_f_ent_6 + T_f_sai_chut_6) / 2
        c_p_f_6 = Cp (Steam;T=T_med_f_6;P=P_f)
        I_f_ent_6 = Enthalpy(Steam;T=T_f_ent_4;P=P_f)
        I_f_sai_6 = Enthalpy(Steam;T=T_f_sai_4;P=P_f)

        "parametros do fluido quente"

        T_q_ent_6 = T_q_sai_5
        T_med_q_6 = abs((T_q_ent_6 + T_q_sai_chut_6)/2)
        c_p_q_6 =Cp(Air_ha;T=T_med_q_6;P=P_q)

        "Balanço de massas e energias"

        T_q_sai_6 := T_q_ent_6 - m_dot_f*(I_f_sai_6-I_f_ent_6)/(m_dot_q * c_p_q_6)

        "############### parametros do loop ##############"

        error_q_4 := ((T_q_sai_chut_4 - T_q_sai_4)^2)^0,5
        error_q_5 := ((T_q_sai_chut_5- T_q_sai_5)^2)^0,5
        error_q_6 := ((T_q_sai_chut_6- T_q_sai_6)^2)^0,5
        error_q_med := (error_q_4 + error_q_5 + error_q_6)/3

        T_q_sai_chut_4 := (T_q_sai_chut_4 + T_q_sai_4)/2
        T_q_sai_chut_5 := (T_q_sai_chut_5 + T_q_sai_5)/2
        T_q_sai_chut_6 := (T_q_sai_chut_6 + T_q_sai_6)/2

    UNTIL(error_q_4<1 AND error_q_5<1 AND error_q_6<1)

    T_it_4 = T_q_sai_chut_4
    T_it_5 = T_q_sai_chut_5
    T_it_6 = T_q_sai_chut_6
END

"############### parametros do trocador ##############"

L = 2
D = 0,076
d_i = 10^(-2)
d_e = 12^(-2)
n_tubos = 19
A_i = pi# * d_i * L * n_tubos
P_i_m = pi# * d_i
A_e = pi# * d_e * L * n_tubos
P_e_m = pi# * d_e
A_i_t = pi# * d_i^2/4
A_c_t = pi# * D^2/4
k_aco = 504

"############### chute inicial de temp. ##############"

T_f_ent_1 = 28 + 273,15
T_f_sai_chut_1 =  600 + 273,15
T_f_sai_chut_5 = (T_f_sai_chut_1 - 373,15 [K])/2 + 373,15 [K]
T_f_sai_chut_6 = T_f_sai_chut_1

T_q_ent_2 = 665 + 273,15
T_q_sai_chut_2 =  50 + 273,15
T_q_sai_chut_4 = (T_q_sai_chut_2 - T_q_ent_2)/3 + T_q_ent_2
T_q_sai_chut_5 = (T_q_sai_chut_2 - T_q_ent_2)/3 * 2 + T_q_ent_2
T_q_sai_chut_6 = T_q_sai_chut_2


"############### Volumes de Controle IV,V,VI ##############"

m_dot_f = 0,16
m_dot_q = 0,102115
P_q = 101,325
P_f = 101,325

"Diâmetro hidraulico"

D_h =  A_i_t/P_i_m

"Defininindo temperaturas intermediarias pro iteração"

CALL temp_it (L;D;d_i;d_e;n_tubos;A_i;P_i_m;A_e;P_e_m;A_i_t;A_c_t;k_aco;T_f_ent_1;T_f_sai_chut_1;T_f_sai_chut_5;T_f_sai_chut_6;T_q_ent_2;T_q_sai_chut_2;T_q_sai_chut_4;T_q_sai_chut_5;T_q_sai_chut_6;m_dot_f;m_dot_q;P_q;P_f;D_h : error_q_4 ;error_q_5 ; error_q_6 ; error_q_med ; T_it_4 ; T_it_5 ; T_it_6)

"############### Volume de Controle IV ##############"

"parametros do fluido frio"

T_f_ent_4 = T_f_ent_1
T_f_sai_4 = 373,15
T_med_f_4 = (T_f_ent_4 +T_f_sai_4) / 2
c_p_f_4 = Cp (Water;T=T_med_f_4;P=P_f)
rho_f_4 = Density(Water;T=T_med_f_4;P=P_f)
mu_f_4 = Viscosity(Water;T=T_med_f_4;P=P_f)
k_f_4 = Conductivity(Water;T=T_med_f_4;P=P_f)
Pr_f_4 = Prandtl(Water;T=T_med_f_4;P=P_f)
I_f_ent_4 = Enthalpy(Water;T=T_f_ent_4;P=P_f)
I_f_sai_4 = Enthalpy(Water;T=T_f_sai_4;P=P_f)


"parametros do fluido quente"

T_q_ent_4 = T_q_ent_2
T_med_q_4 = (T_q_ent_4 + T_it_4)/2
c_p_q_4 =Cp(Air_ha;T=T_med_q_4;P=P_q)
rho_q_4 = Density(Air_ha;T=T_med_q_4;P=P_f)
mu_q_4 =Viscosity(Air_ha;T=T_med_q_4;P=P_q)
k_q_4 =Conductivity(Air_ha;T=T_med_q_4;P=P_q)
Pr_q_4 = Prandtl(Air_ha;T=T_med_q_4;P=P_q)

"Balanço de massas e energias"

m_dot_q * c_p_q_4*(T_q_ent_4-T_q_sai_4) = m_dot_f*(I_f_sai_4-I_f_ent_4)

"Numero de Reynolds"

Re_D_i_4 = (m_dot_f*D_h)/(A_i_t*mu_f_4)

"############### Volume de Controle V ##############"

"parametros do fluido frio"

T_f_ent_5 = T_f_sai_4
T_med_f_5 = (T_f_ent_5 + T_f_sai_chut_5) / 2

" tirando média dos parâmetros termo-físicos"

c_p_f_a = Cp (Water;T=T_med_f_5;P=P_f)
rho_f_a = Density(Water;T=T_med_f_5;P=P_f)
mu_f_a = Viscosity(Water;T=T_med_f_5;P=P_f)
k_f_a = Conductivity(Water;T=T_med_f_5;P=P_f)
Pr_f_a = Prandtl(Water;T=T_med_f_5;P=P_f)
I_f_ent_a = Enthalpy_vaporization(Water;T=T_f_ent_4)
I_f_sai_a = Enthalpy_vaporization(Water;T=T_f_sai_4)

c_p_f_v = Cp (Steam;T=T_med_f_5;P=P_f)
rho_f_v = Density(Steam;T=T_med_f_5;P=P_f)
mu_f_v = Viscosity(Steam;T=T_med_f_5;P=P_f)
k_f_v = Conductivity(Steam;T=T_med_f_5;P=P_f)
Pr_f_v = Prandtl(Steam;T=T_med_f_5;P=P_f)
I_f_ent_v = Enthalpy_vaporization(Steam;T=T_f_ent_4)
I_f_sai_v = Enthalpy_vaporization(Steam;T=T_f_sai_4)

c_p_f_5 = (c_p_f_v + c_p_f_a)/2
rho_f_5 = (rho_f_v + rho_f_a)/2
mu_f_5 = (mu_f_v + mu_f_a)/2
k_f_5 = (k_f_v + k_f_a)/2
Pr_f_5 = (Pr_f_v + Pr_f_a)/2
I_f_ent_5 = (I_f_ent_v + I_f_ent_a)/2
I_f_sai_5 = (I_f_sai_v + I_f_sai_a)/2

"parametros do fluido quente"

T_q_ent_5 = T_q_sai_4
T_med_q_5 = (T_q_ent_5 + T_it_5)/2
c_p_q_5 =Cp(Air_ha;T=T_med_q_5;P=P_q)
rho_q_5 = Density(Air_ha;T=T_med_q_5;P=P_q)
mu_q_5 =Viscosity(Air_ha;T=T_med_q_5;P=P_q)
k_q_5 =Conductivity(Air_ha;T=T_med_q_5;P=P_q)
Pr_q_5 = Prandtl(Air_ha;T=T_med_q_5;P=P_q)

"Balanço de massas e energias"

m_dot_q * c_p_q_5*(T_q_ent_5-T_q_sai_5) = m_dot_f*(I_f_sai_5-I_f_ent_5)

"Numero de Reynolds"

Re_D_i_5 = (m_dot_f*D_h)/(A_i_t*mu_f_5)

"############### Volume de Controle VI ##############"

"parametros do fluido frio"

T_f_ent_6 = T_f_sai_chut_5
T_med_f_6 = (T_f_ent_6 + T_f_sai_chut_6) / 2
c_p_f_6 = Cp (Steam;T=T_med_f_6;P=P_f)
rho_f_6 = Density(Steam;T=T_med_f_6;P=P_f)
mu_f_6 = Viscosity(Steam;T=T_med_f_6;P=P_f)
k_f_6 = Conductivity(Steam;T=T_med_f_6;P=P_f)
Pr_f_6 = Prandtl(Steam;T=T_med_f_6;P=P_f)
I_f_ent_6 = Enthalpy(Steam;T=T_f_ent_4;P=P_f)
I_f_sai_6 = Enthalpy(Steam;T=T_f_sai_4;P=P_f)

"parametros do fluido quente"

T_q_ent_6 = T_q_sai_5
T_med_q_6 = abs((T_q_ent_6 + T_it_6)/2)
c_p_q_6 =Cp(Air_ha;T=T_med_q_6;P=P_q)
rho_q_6 = Density(Air_ha;T=T_med_q_6;P=P_q)
mu_q_6 =Viscosity(Air_ha;T=T_med_q_6;P=P_q)
k_q_6 =Conductivity(Air_ha;T=T_med_q_6;P=P_q)
Pr_q_6 = Prandtl(Air_ha;T=T_med_q_6;P=P_q)

"Balanço de massas e energias"

m_dot_q * c_p_q_6*(T_q_ent_6-T_q_sai_6) = m_dot_f*(I_f_sai_6-I_f_ent_6)

"Numero de Reynolds"

Re_D_i_6 = (m_dot_f*D_h)/(A_i_t*mu_f_6)