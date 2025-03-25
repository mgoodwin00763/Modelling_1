

import matplotlib.pyplot as plt
import numpy as np


def f(h, T, r, alpha, b, zeta1):
    h_grad = -(r * h) - (alpha * b * T) - (alpha * zeta1)
    return h_grad


def g(h, T, R, gamma, en, b, zeta1, zeta2):
    T_grad = (R * T) + (gamma * h) - en * (h + (b * T))**3 + (gamma * zeta1) + zeta2
    return T_grad

def ENSO_Euler_new(h, T, h_star, T_star, b, alpha, r, dt, R, gamma, en, zeta1, zeta2):
    
        h_new = h + dt * f(h_star, T_star, r, alpha, b, zeta1)
        T_new = T + dt * g(h_star, T_star, R, gamma, zeta1, zeta2)
        
        return h_new, T_new

def ENSO_Matsuno_Heun(h, T, h_star, T_star, b, alpha, r, dt, R, gamma, en, zeta1, zeta2, beta):
    
        h_new = h + dt * (beta * f(h_star, T_star, r, alpha, b, zeta1) + (1 - beta) * f(h, T, r, alpha, b, zeta1) )
        T_new = T + dt * (beta * g(h_star, T_star, R, gamma, en, b, zeta1, zeta2)    + (1 - beta) * g(h, T, R, gamma, en, b, zeta1, zeta2) )
        
        return h_new, T_new
    
""" DEFINE SCALING AND INITIAL CONDITIONS """
    
h_scale = 150
T_scale = 7.5
t_scale = 2

h0 = 0 / h_scale
T0 = 1.125 / T_scale

def main(mu0, total_time, en, mu_ann, tau, f_ann=0, f_ran=0, tau_cor=0, dt = 0.1, h0 = h0, T0 = T0):
    
    """ DEFINE VARIABLES """
    b0 = 2.5 #high end value of coupling parameter
    b = b0 * mu0 #thermocline slope
    gamma = 0.75 #feedback of thermocline gradient on SST gradient
    c = 1 #damping rte of SST anomalies
    R = gamma * b - c #Bjerknes positive feedback process
    r = 0.25 #damping of upper ocean heat content
    alpha = 0.125 #relates enchanced easterly wind stress to recharge of ocean heat content
    zeta1 = 0 #random wind stress forcing
    zeta2 = 0 #random heating added to system

    

    """DEFINE INITAL CONDITIONS"""
    h = np.zeros(total_time)
    T = np.zeros(total_time)
    t = np.zeros(total_time)

    h[0] = h0
    
    T[0] = T0

    for i in range(total_time):
        t[i] = i*dt
        
    W = np.random.uniform(-1, 1, size = total_time)
    
    
        
    """ RUN SIMULATION AND PLOT GRAPHS """
    
    for n in range(total_time - 1):
        
        
        
        zeta1 = f_ann * np.cos(2 * np.pi * t[n] * t_scale / tau) + f_ran * W[n] * tau_cor/dt 
        
        mu = mu0 * (1 + mu_ann * np.cos( (2 * np.pi * t[n] * t_scale / tau) - (5 * np.pi / 6)))
        
        b = b0 * mu      
        
        R = (gamma * b) - c
        
        h_star, T_star = ENSO_Matsuno_Heun(h[n], T[n], h[n], T[n], b, alpha, r, dt, R, gamma, en, zeta1, zeta2, beta = 0) 
        
        h[n+1], T[n+1] = ENSO_Matsuno_Heun(h[n], T[n], h_star, T_star, b, alpha, r, dt, R, gamma, en, zeta1, zeta2, beta = 0.5) 
       
    return h, T, t, h_scale, T_scale, t_scale




def plotting(h, T, t, h_scale, T_scale, t_scale, mu0, title):

    fig3, (ax3, ax4) = plt.subplots(1, 2, figsize=(10,7)) # create a new figure and axes
    ax3.plot(t * t_scale, T * T_scale, label = 'SST') 
    ax3.plot(t * t_scale, h * h_scale, label = 'Height anomaly')
    ax3.legend(bbox_to_anchor=(0.5,-0.1))# create a scatter plot on our new axes
    ax3.set_xlabel("Time (months)")
    ax3.set_ylabel("T, h")
    ax3.set_title("Time series plot")
    
    mu = round(mu0, 2)
    fig3.suptitle(f' {title} \n Coupling parameter: mu = {mu} \n Non linearity = {en}')

    ax4.plot(h * h_scale, T * T_scale)   # create a scatter plot on our new axes
    ax4.set_xlabel("Thermocline depth (h) (metres)")
    ax4.set_ylabel("SST anomaly (T) (deg C)")
    ax4.set_title("Phase plot")
    
    fig3.tight_layout() 
    fig3.savefig("huen2.svg")
    
def plotting_twice(h, T, t, h_scale, T_scale, t_scale, mu0, h2, T2, mu02, Title, label1, label2):

    
    
    fig3, (ax3, ax4) = plt.subplots(1, 2, figsize=(10,7)) # create a new figure and axes
    ax3.plot(t * t_scale, T * T_scale, label = f'SST, \n {label1}', color = 'tab:blue') 
    ax3.plot(t * t_scale, h * h_scale, label = f'Height anomaly,\n {label1}', color = 'tab:orange')
    ax3.plot(t * t_scale, T2 * T_scale, label = f'SST, \n {label2}', color="tab:pink", dashes=(4,2)) 
    ax3.plot(t * t_scale, h2 * h_scale, label = f'Height anomaly,  \n {label2} ', color = 'tab:green', dashes=(4,2))
    ax3.legend(bbox_to_anchor=(0.5,-0.2))# create a scatter plot on our new axes
    ax3.set_xlabel("Time (months)")
    ax3.set_ylabel("T, h")
    ax3.set_title("Time series plot")
    
    mu = round(mu0, 2)
    fig3.suptitle(f' {Title} \n mu = {mu}'  )

    ax4.plot(h * h_scale, T * T_scale, label = f'{label1}', color = 'tab:red')
    ax4.plot(h2 * h_scale, T2 * T_scale, label = f'{label2}', color = 'tab:green')  # create a scatter plot on our new axes
    ax4.set_xlabel("Thermocline depth (h) (metres)")
    ax4.set_ylabel("SST anomaly (T) (deg C)")
    ax4.legend(bbox_to_anchor=(0.5,-0.2))
    ax4.set_title("Phase plot")
    
    fig3.tight_layout() 
    fig3.savefig("graph.svg")
    
    

    
""" TASK A""" 

#SETUP INITIAL PARAMETERS

en = 0 #degree of non-linearity of ROM
frequency = np.sqrt(3/32)
period = 2 * np.pi / frequency
total_time =1000 # 5 * round(period / dt)
dt = 0.1
mu0 = 2/3
mu_ann = 0
tau = 1
title = "Task A - Critical value of mu"

h_a, T_a, t_a, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_a, T_a, t_a, h_scale, T_scale, t_scale, mu0, title) 

""" TASK B """
"""TESTING AROUND SUB CRITICAL AND SUPER CRITCAL VALUES OF MU"""

title = "Task B - subcritical coupling parameter"

mu0 = 0.63
total_time = 1000
h_b1, T_b1, t_b1, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_b1, T_b1, t_b1, h_scale, T_scale, t_scale, mu0, title)     

title = "Task B - supercritical coupling parameter"

mu0 = 0.70
total_time = 1000
h_b2, T_b2, t_b2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_b2, T_b2, t_b2, h_scale, T_scale, t_scale, mu0, title) 

""" TASK C """
""" TURN ON NON - LINEARITY """

mu0 = 2/3
en = 0.1
total_time = 1000
Title = 'Task C: Effects of Non-linearity'
label1 = 'Without non-linearity'
label2 = f'en = {en} '
h_c1, T_c1, t_c1, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting_twice(h_a, T_a, t_a, h_scale, T_scale, t_scale, mu0, h_c1, T_c1, mu0, Title, label1, label2)

mu0 = 0.70
en = 0.1
Title = 'Task C: Effects of Non-linearity'
label1 = 'Without non-linearity'
label2 = f'en = {en} '
h_c2, T_c2, t_c2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting_twice(h_b2, T_b2, t_b2, h_scale, T_scale, t_scale, mu0, h_c2, T_c2, mu0, Title, label1, label2)

mu0 = 0.63
en = 0.1
Title = 'Task C: Effects of Non-linearity'
label1 = 'Without non-linearity'
label2 = f'en = {en} '
h_c3, T_c3, t_c3, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting_twice(h_b1, T_b1, t_b1, h_scale, T_scale, t_scale, mu0, h_c3, T_c3, mu0, Title, label1, label2)

mu0 = 0.72
total_time = 1000
en = 0
h_b5, T_b5, t_b5, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
en = 0.1
h_c5, T_c5, t_c5, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting_twice(h_b5, T_b5, t_b5, h_scale, T_scale, t_scale, mu0, h_c5, T_c5, mu0, Title, label1, label2)

mu0 = 0.75
total_time = 1000
en = 0
h_b4, T_b4, t_b4, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
en = 0.1 
h_c4, T_c4, t_c4, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting_twice(h_b4, T_b4, t_b4, h_scale, T_scale, t_scale, mu0, h_c4, T_c4, mu0, Title, label1, label2)

total_time = 3000
title = 'Task C: Effects of Non-linearity'
mu0 = 0.66
en = 0.1
h_b2, T_b2, t_b2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_b2, T_b2, t_b2, h_scale, T_scale, t_scale, mu0, title) 

title = 'Task C: Effects of Non-linearity'
mu0 = 0.67
en = 0.1
h_b2, T_b2, t_b2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_b2, T_b2, t_b2, h_scale, T_scale, t_scale, mu0, title) 

title = 'Task C: Effects of Non-linearity'
mu0 = 0.672
en = 0.1
h_b2, T_b2, t_b2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_b2, T_b2, t_b2, h_scale, T_scale, t_scale, mu0, title) 

title = 'Task C: Effects of Non-linearity'
mu0 = 0.68
en = 0.1
h_b2, T_b2, t_b2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_b2, T_b2, t_b2, h_scale, T_scale, t_scale, mu0, title) 

title = 'Task C: Effects of Non-linearity'
mu0 = 0.69
en = 0.1
h_b2, T_b2, t_b2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_b2, T_b2, t_b2, h_scale, T_scale, t_scale, mu0, title) 

title = 'Task C: Effects of Non-linearity'
mu0 = 0.70
en = 0.1
h_b2, T_b2, t_b2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_b2, T_b2, t_b2, h_scale, T_scale, t_scale, mu0, title) 



title = 'Task C: Effects of Non-linearity'
mu0 = 0.75
en = 0.1
h_b2, T_b2, t_b2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting(h_b2, T_b2, t_b2, h_scale, T_scale, t_scale, mu0, title) 



""" TASK D """
""" SELF - EXCITATION HYPOTHESIS """

total_time = 1000
mu0 = 0.75
en = 0.1
h_d0, T_d0, t_d0, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

mu0 = 2/3
h_d02, T_d02, t_d02, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

mu0 = 0.80
h_d03, T_d03, t_d03, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

mu0 = 0.68
h_d04, T_d04, t_d04, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

total_time = 1000
mu0 = 0.75
mu_ann = 0.2
en = 0.1
t_scale = 2
#tau = 12 / (t_scale * dt)
tau = 24 / t_scale
h_d1, T_d1, t_d1, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

Title = 'Task D: Self-excitation'
label1 = 'Non-varying coupling parameter'
label2 = 'Annual cycle for coupling parameter'
plotting_twice(h_d0, T_d0, t_d0, h_scale, T_scale, t_scale, mu0, h_d1, T_d1, mu0, Title, label1, label2)

mu0 = 2/3
h_d2, T_d2, t_d2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting_twice(h_d02, T_d02, t_d02, h_scale, T_scale, t_scale, mu0, h_d2, T_d2, mu0, Title, label1, label2)

mu0 = 0.68
h_d4, T_d4, t_d4, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting_twice(h_d04, T_d04, t_d04, h_scale, T_scale, t_scale, mu0, h_d4, T_d4, mu0, Title, label1, label2)

mu0 = 0.80
h_d3, T_d3, t_d3, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)
plotting_twice(h_d03, T_d03, t_d03, h_scale, T_scale, t_scale, mu0, h_d3, T_d3, mu0, Title, label1, label2)

moo = np.zeros(total_time)

#PLOT VARIATION OF MU

for n in range(total_time):
    moo[n] = mu0 * (1 + mu_ann * np.cos( (2 * np.pi * n * dt * t_scale / tau) - (5 * np.pi / 6)))
   

fig10, ax10 = plt.subplots() # create a new figure and axes
ax10.plot(t_d1 * t_scale, moo) 
ax10.set_xlabel("Time (months)")
ax10.set_ylabel("Coupling parameter")
ax10.set_title("Time series plot")


""" TASK E """
""" STOCHASTIC INITATION HYPOTHESIS"""

en = 0
mu_ann = 0.2
dt = tau/365 #timestep of one day

#PLOT WITHOUT WIND FORCING
mu0 = 2/3
h_e01, T_e01, t_e01, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

mu0 = 0.67
h_e05, T_e05, t_e05, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

mu0 = 0.68
h_e06, T_e06, t_e06, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

mu0 = 0.69
h_e07, T_e07, t_e07, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

mu0 = 0.70
h_e02, T_e02, t_e02, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

mu0 = 0.75
h_e03, T_e03, t_e03, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)

mu0 = 0.65
h_e04, T_e04, t_e04, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau)



Title = 'TASK E: Noisy wind forcing (No non-linearity)'
label1 = 'Without noise'
label2 = 'With stochastic wind noise'


#PLOTS WITH WIND FORCING 

en = 0
mu_ann = 0.2
f_ann = 0.02
f_ran = 0.2
tau_cor = tau/365

mu0 = 0.66
h_e1, T_e1, t_e1, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau, f_ann = f_ann, f_ran = f_ran, tau_cor = tau_cor)
plotting_twice(h_e01, T_e01, t_e01, h_scale, T_scale, t_scale, mu0, h_e1, T_e1, mu0, Title, label1, label2)

mu0 = 0.67
h_e5, T_e5, t_e5, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau, f_ann = f_ann, f_ran = f_ran, tau_cor = tau_cor)
plotting_twice(h_e05, T_e05, t_e05, h_scale, T_scale, t_scale, mu0, h_e5, T_e5, mu0, Title, label1, label2)

mu0 = 0.68
h_e6, T_e6, t_e6, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau, f_ann = f_ann, f_ran = f_ran, tau_cor = tau_cor)
plotting_twice(h_e06, T_e06, t_e06, h_scale, T_scale, t_scale, mu0, h_e6, T_e6, mu0, Title, label1, label2)

mu0 = 0.69
h_e7, T_e7, t_e7, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau, f_ann = f_ann, f_ran = f_ran, tau_cor = tau_cor)
plotting_twice(h_e07, T_e07, t_e07, h_scale, T_scale, t_scale, mu0, h_e7, T_e7, mu0, Title, label1, label2)

mu0 = 0.70
h_e2, T_e2, t_e2, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau, f_ann = f_ann, f_ran = f_ran, tau_cor = tau_cor)
plotting_twice(h_e02, T_e02, t_e02, h_scale, T_scale, t_scale, mu0, h_e2, T_e2, mu0, Title, label1, label2)

mu0 = 0.75
h_e3, T_e3, t_e3, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau, f_ann = f_ann, f_ran = f_ran, tau_cor = tau_cor)
plotting_twice(h_e03, T_e03, t_e03, h_scale, T_scale, t_scale, mu0, h_e3, T_e3, mu0, Title, label1, label2)

mu0 = 0.65
h_e4, T_e4, t_e4, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau, f_ann = f_ann, f_ran = f_ran, tau_cor = tau_cor)
plotting_twice(h_e04, T_e04, t_e04, h_scale, T_scale, t_scale, mu0, h_e4, T_e4, mu0, Title, label1, label2)





zeeta = np.zeros(total_time)

        #PLOT VARIATION OF zeta1
W = np.random.uniform(-1, 1, size = total_time)


for n in range(total_time):
    zeeta[n] = f_ann * np.cos(2 * np.pi * n * dt * t_scale / tau) + f_ran * W[n] * tau_cor/dt 
            
fig11, ax11 = plt.subplots() # create a new figure and axes
ax11.plot(t_e2 * t_scale, zeeta) 
ax11.set_xlabel("Time (months)")
ax11.set_ylabel("Zeta_1")
ax11.set_title("Wind stress forcing")





""" TASK E BONUS - CHANGING TIMESTEP OR TAU_COR """

""" TASK F """
""" NON LINEARITY + STOCHASTIC FORCING TOGETHER """

en = 0.1
mu0 = 0.75
mu_ann = 0.2
f_ann = 0.02
f_ran = 0.2
tau_cor = tau/365

"""
Title = 'TASK F: Non linearity + Stochastic forcing'
label1 = 'Linear model \n with stochastic wind noise'
label2 = 'Non-linear model \n with stochastic wind noise'
"""

h_f, T_f, t_f, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau, f_ann = f_ann, f_ran = f_ran, tau_cor = tau_cor)

"""
plotting_twice(h_e3, T_e3, t_e3, h_scale, T_scale, t_scale, mu0, h_f, T_f, mu0, Title, label1, label2)
"""

Title = 'TASK F: Non linearity + Stochastic forcing'
label1 = 'Non-Linear model \n without stochastic wind noise'
label2 = 'Non-linear model \n with stochastic wind noise'
plotting_twice(h_c2, T_c2, t_c2, h_scale, T_scale, t_scale, mu0, h_f, T_f, mu0, Title, label1, label2)

""" TASK G """
""" ENSEMBLE MEMBERS, ASSEMBLE! """


en = 0.1
mu0 = 0.75
mu_ann = 0.2
f_ann = 0.02
f_ran = 0.2
tau_cor = tau/365

def taskG(f_ran, h_scale, T_scale, title, iu):
        
    Ensemble_members = 3

    fig9, ax9 = plt.subplots( figsize=(10,7)) # create a new figure and axes
    fig10, ax10 = plt.subplots( figsize=(10,7))
    fig11, ax11 = plt.subplots( figsize=(10,7))

    for i in range(Ensemble_members):
        
        h0 = 0 / h_scale
        T0 = 1.125/ T_scale
        
        if iu == 1:
        
        #Perturbing initial T and h used to start each forecast
            h0 = (0 + np.random.uniform(-150, 150) ) / h_scale
            
            T0 = (1.125 + np.random.uniform(-15, 15) ) / T_scale
        
        
        
        h_g, T_g, t_g, h_scale, T_scale, t_scale = main(mu0, total_time, en, mu_ann, tau, f_ann = f_ann, f_ran = f_ran, tau_cor = tau_cor, h0 = h0, T0 = T0)
        
        ax9.plot(t_g * t_scale, T_g * T_scale, label = f'SST Ensemble member {i}, Initial SST = {round(T0 * T_scale, 3)} deg C') 
        
        ax11.plot(h_g * h_scale, T_g * T_scale, label = f'Ensemble member {i}')
        
        
        ax10.plot(t_g * t_scale, h_g * h_scale, label = f'SST Ensemble member {i}, Initial height anomaly = {round(h0 * h_scale, 3)} m')


    ax9.legend(bbox_to_anchor=(0.5,-0.2))
    ax9.set_xlabel("Time (months)")
    ax9.set_ylabel("SST anomaly")
    ax9.set_title(f'TASK G: Plume diagram of Ensemble forecast \n {title}')
    
    ax10.legend(bbox_to_anchor=(0.5,-0.2))
    ax10.set_xlabel("Time (months)")
    ax10.set_ylabel("Height anomaly")
    ax10.set_title(f'TASK G: Plume diagram of Ensemble forecast \n {title}')
    
    ax11.set_xlabel("Thermocline depth (h) (metres)")
    ax11.set_ylabel("SST anomaly (T) (deg C)")
    ax11.legend(bbox_to_anchor=(0.5,-0.2))
    ax11.set_title("Phase plot")

f_ran = 0.2
iu = 0
title = "Random wind forcing only"
taskG(f_ran, h_scale, T_scale, title, iu)

f_ran = 0
iu = 1
title = "Initial condition uncertainty only"
taskG(f_ran, h_scale, T_scale, title, iu)

f_ran = 0.2
iu = 1
title = "Initial condition uncertainty \n plus random wind forcing"
taskG(f_ran, h_scale, T_scale, title, iu)


    


