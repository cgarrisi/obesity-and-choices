#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 21:19:28 2022

@author: charlesgarrisi
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

def deriv(y, t, beta, gamma, sigma):
    """
    y : liste contenant les 3 fonctions inconnus 
    t : le temps 
    beta, gamma : les deux facteurs du modèle
    """
    S,I,R,V = y 
 
    # Description des 3 equations differentielles 
    dSdt = -S * I  * beta 
    dIdt = S * I  * beta - gamma * I + sigma * beta * V * I
    dRdt = gamma * I 
    dVdt = gamma * I - sigma * beta * I * V

    return dSdt, dIdt, dRdt, dVdt

# Au temps t0,  70% sains, 30% infecté, 0 guéri, 0 vacciné
y0 = 0.7, 0.2, 0.1, 0

# Evolution sur 365 jours 
t = np.linspace(0, 28)

# Paramètres du modèle 
beta = 0.5
gamma =  0.1
sigma = 1.2 #Choix du sigma difficile 



# Resolution des équations differentielles 
ret = odeint(deriv, y0, t, args = (beta, gamma, sigma))
S,I,R,V = ret.T

plt.plot(t, S, label="Sains")
plt.plot(t, V, label="Ex-obèses")
plt.plot(t, R, label="Vaccinés")
plt.plot(t, I,label="Infectés")

plt.xlabel("Temps")
plt.ylabel("Nombre d'individus")
plt.legend()
plt.title(f"Proportion des individus durant une épidémie modélisé par SIRV avec β = {beta} et γ = {gamma}")
