"""
This program calculates the reverse rate coefficients of the reactions in a chemical mechanism.
The mechanism typically has a forward reaction rate coefficient defined along with equilibrium constant.
This also includes Ea, Beta and A, i.e. activation energy, temperature exponent and Arrhenius coefficient for the
forward reaction, while this info is not readily available for the backward reaction. These quantitites are
calculated for the backward reaction here by simulating a 1D flamelet thereby having the reactions at different
temperatures and then fitting a polynomial to obtain the constants.

"""
import matplotlib.pyplot as plt
import cantera as ct
import pickle
import math
import numpy as np
import os
def premix_freeflame(phi, chem_mech):
    # Simulation parameters
    p = ct.one_atm  # pressure [Pa]
    Tin = 291.0+100.0  # unburned gas temperature [K]
    #phi = 1.0
    a = 2.0/phi
    width = 0.03 # m
    loglevel = 1  # amount of diagnostic output (0 to 8)

    # IdealGasMix object used to compute mixture properties, set to the state of the
    # upstream fuel-air mixture
    gas = ct.Solution(chem_mech)
    xch4 = 0.06
    reactants = {'CH4': 0.16, 'O2': (1+1 - 0.16) * 0.23, 'N2': (1+1 - 0.16) * 0.77}
    #reactants = {'CH4': xch4, 'O2': (1 - xch4) * 0.23, 'N2': (1 - xch4) * 0.77}  # premixed gas composition
    gas.TPX = Tin, p, reactants

    # Set up flame object
    f = ct.FreeFlame(gas, width=width)
    #f=ct.CounterflowPremixedFlame(gas,width=width)
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    #f.show_solution()

    # Solve with mixture-averaged transport model
    f.transport_model = 'Mix'
    f.solve(loglevel=loglevel, auto=True)
    return f

def FlameInp(phi, chem_mech):
    f=premix_freeflame(phi,chem_mech)
    points=f.flame.n_points
    T=f.T
    gas=f.gas
    K_dict ={}
    R = ct.Reaction.listFromFile(chem_mech)
    for i in range(points):
        f.set_gas_state(i)
        Temp=T[i]
        equilib_consts=gas.equilibrium_constants
        fwd_const = gas.forward_rate_constants
        for j in range(len(R)):
            y = Temp * math.log(fwd_const[j] / equilib_consts[j])
            try:
                K_dict[j].append(y)
            except:
                K_dict[j] = [y]
    return K_dict, T

def coeff_extract(path):
    with open(path + "/Kdict.pkl", 'rb') as f:
        K_dict = pickle.load(f)
    with open(path + "/T_Kdict.pkl", 'rb') as f:
        T = pickle.load(f)
    krev_coeff={}
    Ru = 8314.46261815324  # J/K/kmol
    for i in range(len(K_dict.keys())):
        p = np.polyfit(T, K_dict[i],1)
        A = math.exp(p[0])
        Ea = -1.0*Ru*p[1]
        krev_coeff[i] = [A, 0.0, Ea]

    with open(path + "/krev_coeff.pkl", 'wb') as f:
        pickle.dump(krev_coeff, f, pickle.HIGHEST_PROTOCOL)

    return 0

def main(chem_mech):
    #chem_mech = "gri30.cti"
    path = os.getcwd()
    phi = [0.639]  # , 0.82, 1.278]#[0.639, 0.684, 0.73, 0.82, 1.278]
    try:
        with open(path + "/Kdict.pkl", 'rb') as f:
            K_dict = pickle.load(f)
        with open(path + "/T_Kdict.pkl", 'rb') as f:
            T = pickle.load(f)
    except:
        K_collect = {}
        for p in phi:
            K_dict, T = FlameInp(p, chem_mech)
        with open(path + "/Kdict.pkl", 'wb') as f:
            pickle.dump(K_dict, f, pickle.HIGHEST_PROTOCOL)
        with open(path + "/T_Kdict.pkl", 'wb') as f:
            pickle.dump(T, f, pickle.HIGHEST_PROTOCOL)
    coeff_extract(path)

if __name__=='__main__':
    chem_mech = "drm19.cti"
    path=os.getcwd()
    phi = [0.639]#, 0.82, 1.278]#[0.639, 0.684, 0.73, 0.82, 1.278]
    try:
        with open(path+"/Kdict.pkl", 'rb') as f:
            K_dict=pickle.load(f)
        with open(path+"/T_Kdict.pkl", 'rb') as f:
            T=pickle.load(f)
    except:
        K_collect = {}
        for p in phi:
            K_dict,T = FlameInp(p, chem_mech)
        with open(path + "/Kdict.pkl", 'wb') as f:
            pickle.dump(K_dict,f,pickle.HIGHEST_PROTOCOL)
        with open(path + "/T_Kdict.pkl", 'wb') as f:
            pickle.dump(T, f, pickle.HIGHEST_PROTOCOL)
    y = K_dict.keys()
    coeff_extract(path)
    for i in range((len(y))):
        plt.scatter(T,K_dict[i], s=10)

    plt.legend(y)
    plt.xlabel('Temperature(K)')
    plt.ylabel('Intensity Ratio')
    plt.savefig('Ratio_Tabulation', dpi=300, bbox_inches='tight')
    plt.show()
