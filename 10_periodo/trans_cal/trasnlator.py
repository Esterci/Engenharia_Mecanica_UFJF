from CoolProp.CoolProp import PropsSI
import numpy as np


"############### chute inicial de temp. ##############"

T_f_ent_1 = 28 + 273.15
T_f_sai_chut_1 =  600 
T_f_sai_chut_5 = 399.23628103678413
T_f_sai_chut_6 = 540

T_q_ent_2 = 665 + 273.15
T_q_sai_chut_2 =  473.46962790388955
T_q_sai_chut_4 = 665.8282281750314
T_q_sai_chut_5 = 579.7925504755126
T_q_sai_chut_6 = T_q_sai_chut_2

"############### chute inicial de comprimentos. ##############"

L_chut_5 = 0.6

#m_dot_f_eval = np.linspace(0.009,0.013,num=5)
#d_i_eval = np.linspace(8**(-3),10**(-3),num=5)
m_dot_f_eval = [0.102115]
d_i_eval = [10e-3]


for d_i in d_i_eval:

    d_e = d_i + 1e-3

    for m_dot_f in m_dot_f_eval:



        "############### parametros do trocador ##############"

        L = 2
        D = 0.076
        n_tubos = 19
        P_T = 15 * 10**(-3)
        B_chi = L / 8
        A_i = np.pi * d_i * L
        P_i_m = np.pi * d_i
        A_e = np.pi * d_e * L
        P_e_m = np.pi * d_e
        A_i_t = np.pi * d_i**2/4
        A_c_t = np.pi * D**2/4
        k_aco = 401

        "############### Volumes de Controle IV.V.VI ##############"

        m_dot_q = 0.102115
        P_q = 101325
        P_f = 101325

        "Diâmetro hidraulico"

        D_h =  A_i_t/P_i_m


        "parametros do fluido quente"

        T_q_ent_kern = T_q_ent_2
        T_med_q_kern = (T_q_ent_kern + T_q_sai_chut_5)/2
        mu_q_kern = PropsSI('viscosity','P',P_q,'T',T_med_q_kern,'air')
        Pr_q_kern = PropsSI('Prandtl','P',P_q,'T',T_med_q_kern,'air')
        k_q_kern = PropsSI('conductivity','P',P_q,'T',T_med_q_kern,'air')

        D_e_kern = (4*((P_T**2 * 3**(1/2))/4 - (np.pi * d_e**2)/8))/(np.pi * d_e/2)

        C_kern = P_T - d_e

        A_c = D_e_kern/P_T * C_kern * B_chi

        Re_c = m_dot_q/A_c * D_e_kern/mu_q_kern

        h_e = (0.36 * k_q_kern)/D_e_kern * Re_c**0.55 * Pr_q_kern**(1/3)

        "############### Volume de Controle IV ##############"

        "parametros do fluido frio"

        T_f_ent_4 = T_f_ent_1
        T_f_sai_4 = 373.15
        T_med_f_4 = (T_f_ent_4 +T_f_sai_4) / 2
        c_p_f_4 = PropsSI('CPMASS','P',P_f,'T',T_med_f_4,'water')
        I_f_ent_4 = PropsSI('H','P',P_f,'T',T_f_ent_1,'water')
        I_f_sai_4 = PropsSI('H','P',P_f,'T',T_f_sai_4,'water')
        mu_f_4 = PropsSI('viscosity','P',P_f,'T',T_med_f_4,'water')
        Pr_f_4 = PropsSI('Prandtl','P',P_f,'T',T_med_f_4,'water')
        k_f_4 = PropsSI('conductivity','P',P_f,'T',T_med_f_4,'water')


        "parametros do fluido quente"

        T_q_ent_4 = T_q_ent_2
        T_med_q_4 = (T_q_ent_4 + T_q_sai_chut_4)/2
        c_p_q_4 = PropsSI('CPMASS','P',P_q,'T',T_med_q_4,'air')

        "parametros adicionais"

        mu_f_s_4 = PropsSI('viscosity','P',P_q,'T',T_med_q_4,'air')


        "Balanço de massas e energias"

        T_q_sai_4 = T_q_ent_4 - m_dot_f*(I_f_sai_4-I_f_ent_4)/(m_dot_q * c_p_q_4)

        q_4 = m_dot_q * c_p_q_4 * (T_q_ent_4 - T_q_sai_4)

        "calculo da efetividade"

        c_q_4 = m_dot_q * c_p_q_4
        c_f_4 = m_dot_f * c_p_f_4

        if c_q_4<c_f_4 :
            c_min_4 = c_q_4 
            c_max_4 = c_f_4

        else:
            c_min_4 = c_f_4
            c_max_4 = c_q_4

        c_r_4 = c_min_4/c_max_4

        q_max_4 = c_min_4 * (T_q_ent_4 - T_f_ent_4)

        epsilon_4 = q_4 / q_max_4

        e_c_4 = (2/epsilon_4 - 1 + c_r_4)/(1 + c_r_4**2)**(1/2)

        NUT_4 = -(1 + c_r_4**2)**(-1/2) * np.log((e_c_4 - 1)/(e_c_4 + 1))

        "Numero de Reynolds"

        Re_D_i_4 = (m_dot_f*D_h)/(A_i_t*mu_f_4)

        "numero de Nusselt"

        Nu_D_4 = 3.66

        "coeficiente de convecção interna"

        h_i_4 = (Nu_D_4*k_f_4)/(d_i)

        "coeficientes globais"

        R_d_i_4 = 0.0002
        R_d_e = 0.0009

        U_i_4 = (1/h_i_4 + R_d_i_4 + d_i * (np.log(d_e/d_i))/(2 * k_aco) + (d_i/d_e) * R_d_e + d_i/(d_e * h_e))**(-1)

        "encontrando comprimento do volume de controle"

        L_4 = (NUT_4 * c_min_4)/(n_tubos * U_i_4 * d_i * np.pi)

        "############### Volume de Controle V ##############"

        "parametros do fluido frio"

        T_f_ent_5 = T_f_sai_4
        T_med_f_5 = (T_f_ent_5 + T_f_sai_chut_5) / 2
        I_f_ent_5 = PropsSI('H','Q',0,'T',T_f_ent_5,'water')
        I_f_sai_5 = PropsSI('H','Q',1,'T',T_f_sai_chut_5,'water')

        mu_f_a = PropsSI('viscosity','Q',0,'T',T_med_f_5,'water')
        Pr_f_a = PropsSI('Prandtl','Q',0,'T',T_med_f_5,'water')
        k_f_a = PropsSI('conductivity','Q',0,'T',T_med_f_5,'water')
        c_p_f_a = PropsSI('CPMASS','Q',0,'T',T_med_f_5,'water')
        I_f_a = PropsSI('H','Q',0,'T',T_med_f_5,'water')
        rho_f_a = PropsSI('D','Q',0,'T',T_med_f_5,'water')

        mu_f_v = PropsSI('viscosity','Q',1,'T',T_med_f_5,'water')
        Pr_f_v = PropsSI('Prandtl','Q',1,'T',T_med_f_5,'water')
        k_f_v = PropsSI('conductivity','Q',1,'T',T_med_f_5,'water')
        c_p_f_v = PropsSI('CPMASS','Q',1,'T',T_med_f_5,'water')
        I_f_v = PropsSI('H','Q',1,'T',T_med_f_5,'water')
        rho_f_v = PropsSI('D','Q',0,'T',T_med_f_5,'water')

        "parametros do fluido quente"

        T_q_ent_5 = T_q_sai_4
        T_med_q_5 = (T_q_ent_5 + T_q_sai_chut_5)/2
        c_p_q_5 = PropsSI('CPMASS','P',P_q,'T',T_med_f_5,'air')

        "parametros adicionais"

        mu_f_s_v = PropsSI('viscosity','Q',0,'T',647.096,'water')


        "Balanço de massas e energias"

        T_q_sai_5 = T_q_ent_5 - m_dot_f*(I_f_sai_5 - I_f_ent_5)/(m_dot_q * c_p_q_5)

        q_5 = m_dot_q * c_p_q_5 * (T_q_ent_5 - T_q_sai_5)

        "calculo da efetividade"

        c_q_5 = m_dot_q * c_p_q_5

        c_min_5 = c_q_5 

        c_r_5 = 0

        q_max_5 = c_min_5 * (T_q_ent_5 - T_f_ent_5)

        epsilon_5 = q_5 / q_max_5

        NUT_5 = -np.log(1-epsilon_5)

        "Numero de Reynolds"

        Re_D_i_5_v = (m_dot_f*D_h)/(A_i_t*mu_f_v)

        Re_D_i_5_a = (m_dot_f*D_h)/(A_i_t*mu_f_v)

        "Calculo do coeficiente convectivo"

        q_s_flux = q_5/( n_tubos * np.pi * d_i * L_chut_5)

        X_bar = (q_s_flux * np.pi * d_i * L_chut_5)/(m_dot_f*(I_f_v-I_f_a))

        G_s_f = 1

        Fr = ((m_dot_f/A_i_t) / rho_f_a)/(2*10*d_i)

        if Fr >= 0.04:
            f_Fr =1

        else:
            f_Fr = 2.63*Fr**0.3

        "numero de Nusselt"

        Nu_D_v = 0.27 * Re_D_i_5_v**(4/5) * Pr_f_v**(1/3) * (mu_f_v * mu_f_s_v)**(0.14)

        "coeficiente de convecção interna"

        h_mf = (Nu_D_v*k_f_v)/(d_i)

        if (X_bar > 0) and (X_bar <= 0.8):

            aux = (1.136*(rho_f_a/rho_f_v)**0.45 * X_bar**0.72 * (1-X_bar)**0.08 * f_Fr +
            667*(q_s_flux/((m_dot_f/A_i_t)*(I_f_v-I_f_a)))**0.7 * (1-X_bar)**0.8 * G_s_f)

        else:

            aux = (0.6683*(rho_f_a/rho_f_v)**0.1 * X_bar**0.16 * (1-X_bar)**0.64 * f_Fr + 
            1058*(q_s_flux/((m_dot_f/A_i_t)*(I_f_v-I_f_a)))**0.7 * (1-X_bar)**0.8 * G_s_f)

        h_i_5 = aux * h_mf

        "coeficientes globais"

        R_d_i_5 = 0.0002

        U_i_5 = (1/h_i_5 + R_d_i_5 + d_i * (np.log(d_e/d_i))/(2 * k_aco) + (d_i/d_e) * R_d_e + d_i/(d_e * h_e))**(-1)

        "encontrando comprimento do volume de controle"

        L_5 = (NUT_5 * c_min_5)/(n_tubos * U_i_5 * d_i * np.pi)

        "lei de resfriamento de Newton"

        T_f_sai_5 = T_f_ent_5 + q_5/h_i_5
        
        print('Re_D_i_5_v: {}\nepsilon_5: {}'.format(Re_D_i_5_v,epsilon_5))
        print('T_q_sai_5 = {}'.format(T_q_sai_5))
        print('NUT_5 = {}'.format(NUT_5))
        print('X_bar = {}'.format(X_bar))
        print('G_s_f = {}'.format(G_s_f))
        print('Fr = {}'.format(Fr))
        print('h_mf = {}'.format(h_mf))
        print('aux ={}'.format(aux))
        print('h_i_5 = {}'.format(h_i_5))
        print('U_i_5 = {}'.format(U_i_5))
        print('L_5 = {}'.format(L_5))
        print('T_f_ent_5 = {}'.format(T_f_ent_5))
        print('T_f_sai_5 = {}'.format(T_f_sai_5))
        print('T_q_ent_5 = {}'.format(T_q_ent_5))
        print('T_q_sai_5 = {}'.format(T_q_sai_5))
        print('  ')

        "############### Volume de Controle VI ##############"

        "parametros do fluido frio"

        T_f_ent_6 = T_f_sai_5
        T_med_f_6 = (T_f_ent_6 + T_f_sai_chut_6) / 2
        c_p_f_6 = PropsSI('CPMASS','P',P_f,'T',T_med_f_6,'water')
        I_f_ent_6 = PropsSI('H','P',P_f,'T',T_f_ent_6,'water')
        I_f_sai_6 = PropsSI('H','P',P_f,'T',T_f_sai_chut_6,'water')
        mu_f_6 = PropsSI('viscosity','P',P_f,'T',T_med_f_6,'water')
        Pr_f_6 = PropsSI('Prandtl','P',P_f,'T',T_med_f_6,'water')
        k_f_6 = PropsSI('conductivity','P',P_f,'T',T_med_f_6,'water')


        "parametros do fluido quente"

        T_q_ent_6 = T_q_sai_5
        T_med_q_6 = (T_q_ent_6 + T_q_sai_chut_6)/2
        c_p_q_6 = PropsSI('CPMASS','P',P_q,'T',T_med_q_6,'air')


        "parametros adicionais"

        mu_f_s_6 = PropsSI('viscosity','P',P_f,'T',T_med_f_6,'water')


        "Balanço de massas e energias"

        T_q_sai_6 = T_q_ent_6 - m_dot_f*(I_f_sai_6-I_f_ent_6)/(m_dot_q * c_p_q_6)

        q_6 = m_dot_q * c_p_q_6 * (T_q_ent_6 - T_q_sai_6)

        "calculo da efetividade"

        c_q_6 = m_dot_q * c_p_q_6
        c_f_6 = n_tubos * m_dot_f * c_p_f_6

        if (c_q_6<c_f_6) : 
            c_min_6 = c_q_6 
            c_max_6 = c_f_6

        else: 
            c_min_6 = c_f_6
            c_max_6 = c_q_6


        c_r_6 = c_min_6/c_max_6

        q_max_6 = c_min_6 * (T_q_ent_6 - T_f_ent_6)

        epsilon_6 = q_6 / q_max_6

        e_c_6 = (2/epsilon_6 - 1 + c_r_6)/(1 + c_r_6**2)**(1/2)

        NUT_6 = -(1 + c_r_6**2)**(-1/2) * np.log((e_c_6 - 1)/(e_c_6 + 1))

        "Numero de Reynolds"

        Re_D_i_6 = (m_dot_f*D_h)/(A_i_t*mu_f_6)

        "numero de Nusselt"

        Nu_D_6 = 0.27 * Re_D_i_6**(4/5) * Pr_f_6**(1/3) * (mu_f_6 * mu_f_s_6)**(0.14)

        "coeficiente de convecção interna"

        h_i_6 = (Nu_D_6*k_f_6)/(d_i)

        "coeficientes globais"

        R_d_i_6 = 0.0002

        U_i_6 = (1/h_i_6 + R_d_i_6 + d_i * (np.log(d_e/d_i))/(2 * k_aco) + (d_i/d_e) * R_d_e + d_i/(d_e * h_e))**(-1)

        "encontrando comprimento do volume de controle"

        L_6 = (NUT_6 * c_min_6)/(n_tubos * U_i_6 * d_i * np.pi)

        "lei de resfriamento de Newton"

        T_f_sai_6 = T_f_ent_6 + q_6/h_i_6
        
        print('Volume 6\n')
        
        print('Re_D_i_6_v: {}\nepsilon_5: {}'.format(Re_D_i_6,epsilon_6))
        print('T_q_sai_6 = {}'.format(T_q_sai_6))
        print('NUT_6 = {}'.format(NUT_6))
        print('h_i_6 = {}'.format(h_i_6))
        print('U_i_6 = {}'.format(U_i_6))
        print('L_6 = {}'.format(L_6))
        print('T_f_ent_6 = {}'.format(T_f_ent_6))
        print('T_f_sai_6 = {}'.format(T_f_sai_6))
        print('T_q_ent_6 = {}'.format(T_q_ent_6))
        print('T_q_sai_6 = {}'.format(T_q_sai_6))
        print('  ')
        print(L_4+L_5+L_6)
        print('  ')
        
        
        "############### parametros do loop ##############"

        error_q_4 = ((T_q_sai_chut_4 - T_q_sai_4)**2)**0.5
        error_q_5 = ((T_q_sai_chut_5- T_q_sai_5)**2)**0.5
        error_q_6 = ((T_q_sai_chut_6- T_q_sai_6)**2)**0.5

        error_f_5 = ((T_f_sai_chut_5- T_f_sai_5)**2)**0.5
        error_f_6 = ((T_f_sai_chut_6- T_f_sai_6)**2)**0.5

        erro_L_5 = ((L_chut_5 - L_5)**2)**0.5

        print('  ')

        T_q_sai_chut_4 = (T_q_sai_chut_4 + T_q_sai_4)/2
        T_q_sai_chut_5 = (T_q_sai_chut_5 + T_q_sai_5)/2
        T_q_sai_chut_6 = (T_q_sai_chut_6 + T_q_sai_6)/2

        T_f_sai_chut_5 = (T_f_sai_chut_5 + T_f_sai_5)/2
        T_f_sai_chut_6 = (T_f_sai_chut_6 + T_f_sai_6)/2

        L_chut_5 = (L_chut_5 + L_5)/2
        
        print(error_q_4)
        print(error_q_6)
        print(error_q_6)
        print(' ')
        print(error_f_5)
        print(error_f_6)
        print(' ')
        print(erro_L_5)
        print(' ')
        print(T_q_sai_chut_4)
        print(T_q_sai_chut_5)
        print(T_q_sai_chut_6)
        print(' ')
        print(T_f_sai_chut_5)
        print(T_f_sai_chut_6)
        print(' ')
        print(L_chut_5)  
        print('____________')


