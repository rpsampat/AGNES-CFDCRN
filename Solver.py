import SaveReactors as sr
import cantera as ct
import pickle
import time
from numpy import array, add, linalg, subtract, amax, absolute, zeros, matrix, transpose, matmul, double, divide, \
    multiply, ones
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve, gmres, cg
from scipy.integrate import BDF
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import math
import os
import rev_rate_coeff_calc as rrc
from joblib import Parallel, delayed
import multiprocessing
import sys
from functools import partial
#from GoverningEquations import GoverningEquations


class CRN_Solver:

    def __init__(self, path, PSR, chem_mech, energy, mfc_gen, boundary, massimb, process):
        self.error_log = []  # logging errors
        self.error_rate_log = []  # logging relative rate of change of error
        self.error_spec_log = []  # logging species associated with epsilon
        self.solver_type = [0]  # 0: local solver, 1: Newton global solver, 2: Time-stepping global solver
        self.errorlog_co2 = []
        self.errorlog_no = []
        self.errorlog_no2 = []
        self.network = []

        self.f_vect = []
        self.Cw = []
        self.fvect_Cw = []
        self.RatesofProduction =[]
        self.J_diffu = []
        self.process = process

        self.dt_fac = 10  # int(tmax/tmin)
        self.g_old = 1.0
        self.epsilon = 1.0  # error
        self.epsilon_omega = 1.0  # rate of error
        self.epsilon_g_omega = 1.0  # rate of change of conserved property.
        self.epsilon_g_omega_rel = 0.0  # rate of change of conserved property.
        self.g_omega = 1.0  # total conserved property
        self.newton_count = 0
        self.RHS_Jac0 = 0.0
        self.alpha = 0.01  # 1e-6
        self.alpha_default = 0.01
        self.max_obj = 0
        self.max_obj_rate = 0

        self.PSR = PSR
        self.path = path
        self.chem_mech = chem_mech
        self.energy = energy
        self.mfc_gen = mfc_gen
        self.boundary = boundary
        self.massimb = massimb
        #self.ge = GoverningEquations(process, PSR, chem_mech, energy, mfc_gen,self)
        self.tol = 1e-15
        self.T_norm = 1.0  # (2500.0 - 0.0) * 4000
        self.Y_high = []
        self.Y_low = []
        self.track_jacobian = "on"  # always keep "on" as required for Jacobian sparse matrix generation
        self.track_jacobian_compo = "off"  # on/off : to track Js and Jw individually
        if self.process == "PIVCRN":
            self.subvect_size = self.mfc_gen.vectsize

        else:
            self.subvect_size = self.mfc_gen.vectsize
        self.constraints()
        rrc.main(self.chem_mech)
        path1 = os.getcwd()
        with open(path1 + "/krev_coeff.pkl", 'rb') as f:
            self.krev = pickle.load(f)



    def constraints(self):
        p = len(self.PSR)
        len_species = len(self.PSR[0].thermo.species_names)
        if self.process == "PIVCRN":
            self.Y_high = ones(p * self.subvect_size)
            self.Y_low = zeros(p * self.subvect_size)
            for i in range(p):
                self.Y_high[i * self.subvect_size + len_species] = 4500
                self.Y_low[i * self.subvect_size + len_species] = 300

        elif self.process == "CFDCRN":
            if self.energy == "on":
                self.Y_high = ones(p * self.subvect_size)
                self.Y_low = zeros(p * self.subvect_size)
                for i in range(p):
                    self.Y_high[i * self.subvect_size + len_species] = 4500
                    self.Y_low[i * self.subvect_size + len_species] = 300
            else:
                self.Y_high = ones(p * self.subvect_size)
                self.Y_low = zeros(p * self.subvect_size)

    def PSR_insert_gas(self, Y, species):
        """
        Inserting new gas state in PSR object
        :param Y:
        :param PSR:
        :param vectsize:
        :param species:
        :return:
        """
        for p in range(len(self.PSR)):
            begin = p * self.subvect_size
            y0 = []
            for i in range(len(species)):
                conc = abs(Y[begin + i])
                y0.append(conc)
            # mass compatability
            mass_total = sum(y0)
            y0 = [y / mass_total for y in y0]
            g = self.PSR[p].thermo
            if self.energy=="on":
                Temp = Y[begin + len(species)]
            else:
                Temp = g.T
            #print "Temperature Calc=",Temp_calc
            if self.process == "CFDCRN":
                Pr = g.P
                g.TPY = Temp, Pr, y0
            elif self.process == "PIVCRN":
                d = self.mfc_gen.rho_calc[p]
                R = 8314.46261815324  # J/K/kmol
                MMW = self.PSR[p].thermo.mean_molecular_weight
                Temp_calc = ct.one_atm / (d * R / MMW)
                try:
                    g.TDY = Temp, d, y0
                except:
                    print Temp
                    print d
                """if abs(Temp-Temp_calc)>=0.0:
                    #g.TPY = Temp, ct.one_atm, y0
                    g.TDY = Temp, d, y0
                else:
                    g.TPY = Temp_calc, ct.one_atm, y0"""

            self.PSR[p].insert(g)
        Yfinal = self.PSR_gas_conc(Y)

        return Yfinal

    def PSR_gas_conc(self, Y):
        Yfinal = []
        for p in range(len(self.PSR)):
            if self.process == "CFDCRN":
                if self.energy == 'on':
                    Yfinal += list(list(self.PSR[p].thermo.Y) + [self.PSR[p].thermo.T])
                else:
                    Yfinal += list(self.PSR[p].thermo.Y)
            elif self.process == "PIVCRN":
                if self.energy == 'on':
                    Yfinal += list(
                        list(self.PSR[p].thermo.Y) + [self.PSR[p].thermo.T])
                else:
                    Yfinal += list(
                        list(self.PSR[p].thermo.Y))
        Yfinal = np.array(Yfinal)
        return Yfinal

    def objective_func(self, t, Y):
        """
        Calculates the value of the function, net production rates of species in this case,
        for current state of all reactors. This is used as the RHS in the Newton's method.
        :param Y:
        :param species:
        :return:
        """
        # Todo: This function has already been converted to numpy arrays
        # print "Objective function Calculation
        vectorsize = self.subvect_size * len(self.PSR)
        species = self.PSR[0].thermo.species_names
        self.Cw = np.zeros(vectorsize)
        self.RatesofProduction = np.zeros(vectorsize)
        self.J_diffu = np.zeros(vectorsize)
        Y_temp = self.PSR_insert_gas(Y, species)  # required to update species in reactors
        colptr = self.mfc_gen.coeffmat.indptr
        rowind = self.mfc_gen.coeffmat.indices
        # Obtaining rates of production
        # print "Obtaining production rates"
        for r_ind in range(len(self.PSR)):
            # TODO: This part could be parallelised
            try:
                r = self.PSR[r_ind]
                # print "Temperature=", r.thermo.T
                # print "Pressure=",r.thermo.P
                # print "Volume=", r.volume
                rates = np.array(list(r.thermo.net_production_rates))  # kmol/m^3/s
                # print "Max rate=", max(abs(rates))
                # print "O2 rate of prod=", rates[r.thermo.species_index("O2")]
                mw = np.array(list(r.thermo.molecular_weights))  # kg/kmol
                # print "O2 mw weight=", mw[r.thermo.species_index("O2")]
                #u_k = np.array(list(r.thermo.partial_molar_int_energies))  # J/kmol

                # ToDo: this is wrong. Need to find quantity for standard enthalpy of formation of species!!
                #  This is causing temperature to become unstable and rise nonlinearly

                #u_k = np.array(r.thermo.partial_molar_enthalpies)# - r.thermo.partial_molar_cp * r.thermo.T)
                u_k = np.array(r.thermo.standard_enthalpies_RT-(r.thermo.standard_cp_R*r.T))

                #u_k = np.array(list(r.thermo.standard_int_energies_RT))
                energy_production = 0.0
                mass_curr = r.mass
                # print "Mass current =",mass_curr
                rowrange = range(colptr[r_ind], colptr[r_ind + 1])
                mflux = 0.0
                for ind in range(self.subvect_size):
                    if ind < len(species):
                        # TODO: include turbulence diffusion effects
                        # Species Mass conservation Equation: production (omega)
                        self.RatesofProduction[
                            r_ind * self.subvect_size + ind] = rates[ind] * mw[ind] * r.volume  # /mass_curr
                        energy_production -= rates[ind] * r.volume * u_k[ind]  # /mass_curr
                        for pt in rowrange:
                            coeff = self.mfc_gen.coeffmat.data[pt]
                            q = rowind[pt]
                            # Species Mass conservation Equation: advection (Cw_1)
                            if self.process == "PIVCRN":
                                # massflowreactor array of velocities for PIVCRN
                                mflux_y = -1.0 * coeff * self.massimb.massflowreactor[r_ind] * self.PSR[r_ind].thermo.Y[
                                    ind]
                                # Mass Conservation Equation
                                #self.Cw[q * self.subvect_size + (
                                 #       self.subvect_size - 1)] += mflux_y  # conservation of mass enforced ; mflux
                                self.Cw[q * self.subvect_size + ind] += mflux_y
                            else:  # CFDCRN
                                # Species Mass conservation Equation: advection (Cw_1)
                                # coeff * massflow[r_ind] as coeffmat generated by normalising columns i.e. total outflow
                                # from a reactor distributed over various other reactors
                                self.Cw[q * self.subvect_size + ind] += -1.0 * coeff * \
                                                                   self.massimb.massflowreactor[r_ind] * \
                                                                   self.PSR[r_ind].thermo.Y[
                                                                       ind]  # /mass_curr
                            """if r_ind * self.mfc_gen.vectsize + ind ==0:
                                print self.Cw[0]"""
                            """if r_ind ==112 and ind == r.thermo.species_index("O2"):
                                #print self.massimb.massflowreactor[pt]
                                print Cw[r_ind * self.mfc_gen.vectsize + ind]"""
                            # Species Mass conservation Equation: diffusion (Cw2)
                            try:
                                # ToDO: check diffusion equation.
                                S_pq = self.mfc_gen.wall_map[r_ind][q].area
                                if S_pq < 0.0:
                                    print "Negative area!!"
                                dx = (self.PSR[r_ind].volume ** 0.333 + self.PSR[q].volume ** 0.333) / 2.0
                                dx1 = self.PSR[q].volume ** 0.333
                                dx2 = self.PSR[r_ind].volume ** 0.333
                                D_spec = self.PSR[q].thermo.mix_diff_coeffs_mass[ind] # m^2/s
                                rho = self.mfc_gen.rho_calc[q]
                                delta_Y = self.PSR[q].thermo.Y[ind] - self.PSR[r_ind].thermo.Y[ind]
                                self.J_diffu[q * self.subvect_size + ind] += -1.0 * rho * D_spec * S_pq * (
                                        delta_Y / dx)  # /mass_curr
                                # diffusion gradient in opposite direction
                                self.J_diffu[r_ind * self.subvect_size + ind] += self.mfc_gen.rho_calc[r_ind] * D_spec * S_pq * (
                                        delta_Y / dx)  # /mass_curr
                            except:
                                pass

                    else:
                        # Energy Conservation Equation: advection (Cw_1)
                        # TODO: include thermal diffusion effects
                        self.RatesofProduction[r_ind * self.subvect_size + ind] = energy_production  # energy production
                        # print "energy_production=", energy_production
                        # TODO: include thermal diffusion effects
                        for pt in rowrange:
                            coeff = self.mfc_gen.coeffmat.data[pt]
                            q = rowind[pt]
                            # energy advection
                            if self.process == "PIVCRN":
                                self.Cw[q * self.subvect_size + ind] += -1 * coeff * \
                                                                   self.massimb.massflowreactor[r_ind] * \
                                                                   self.PSR[r_ind].thermo.cp_mass * self.PSR[
                                                                       r_ind].T
                            else:
                                self.Cw[q * self.subvect_size + ind] += -1 * coeff * \
                                                                   self.massimb.massflowreactor[r_ind] * \
                                                                   self.PSR[r_ind].thermo.cp_mass * self.PSR[
                                                                       r_ind].T
                            try:
                                # wall_map: dictionary of wall objects mapped from reactor to neighbouring reactors.
                                S_pq = self.mfc_gen.wall_map[r_ind][q].area
                                dx = (self.PSR[r_ind].volume ** 0.333 + self.PSR[q].volume ** 0.333) / 2.0
                                K_react = self.PSR[q].thermo.thermal_conductivity
                                rho = self.mfc_gen.rho_calc[q]
                                delta_T = self.PSR[q].T - self.PSR[r_ind].T
                                Js_diffu = -1.0 * rho * K_react * S_pq * (delta_T / dx)
                                self.J_diffu[q * self.subvect_size + ind] += Js_diffu
                                # diffusion gradient in opposite direction
                                self.J_diffu[r_ind * self.subvect_size + ind] += -Js_diffu

                            except:
                                pass

            except Warning:
                print "Failed at = ", r_ind
                continue

        rHS_Jac = add(self.fvect_Cw, array(self.RatesofProduction))
        rHS_Jac1 = add(array(rHS_Jac), self.J_diffu)
        RHS_Jac = add(array(rHS_Jac1), array(self.Cw))
        #RHS_Jac = add(self.fvect_Cw, array(self.Cw))

        return RHS_Jac

    def source_analytic(self, gas, vol, mass):
        """
        Analytic Jacobian of chemical source term
        :param gas:
        :param vol:
        :param mass:
        :return:
        """
        #print "Evaluating Source jac"
        y = list(gas.concentrations)
        mw = list(gas.molecular_weights)  # kg/kmol
        mwu_k = list(gas.partial_molar_int_energies)  # J/kmol
        R = ct.Reaction.listFromFile(self.chem_mech)  # list of reactions in specified mechanism
        R_rev = gas.reverse_rate_constants
        Jw = zeros((len(y), len(y)))
        Jw_energy = zeros((len(y), len(y)))
        if self.energy == 'on':
            Jw = zeros((len(y) + 1, len(y) + 1))
            Jw_energy = zeros((len(y) + 1, len(y) + 1))
        forward = gas.forward_rates_of_progress  # kmol/m^3/s
        # print "Forward = ", forward
        backward = gas.reverse_rates_of_progress  # kmol/m^3/s
        #print "Backward = ", backward
        production_rate = gas.net_production_rates  # kmol/m^3/s

        test_rate_grad = 0
        test2_rate_grad = 0
        test3_rate_grad = 0
        specific_heat = gas.partial_molar_cp
        #hf = list(gas.partial_molar_int_energies)#gas.partial_molar_enthalpies #- gas.partial_molar_cp * gas.T
        #hf = list(gas.partial_molar_enthalpies)
        hf = list(gas.standard_enthalpies_RT-(gas.standard_cp_R*gas.T))
        rates = list(gas.net_production_rates) # kmol/m^3/s

        for r in range(len(R)):
            # dictionary of reactant stoichiometric coeff
            react = R[r].reactants
            prod = R[r].products
            rate_grad = {}
            rate_grad_energy = {}
            Gas_constant = 8314.46261815324  # J/K/kmol
            try:
                E_a = R[r].rate.activation_energy  # J/kmol
                beta = R[r].rate.temperature_exponent
                E_a2 = self.krev[r][2]
                beta2 = self.krev[r][1]
            except:
                # in case of falloff reactions, i.e. reactions with different coefficient of [M](third body)
                # at different pressures.
                # set to low_pressure by default as current application for low pressure systems.
                E_a = R[r].low_rate.activation_energy  # J/kmol
                beta = R[r].low_rate.temperature_exponent


            # Calculating Reaction rate gradients
            for s in react.keys():
                ind = gas.species_index(s)
                try:
                    rate_grad[ind] += react[s] * forward[r] / ((y[ind]) + 1e-15)
                except:
                    rate_grad[ind] = react[s] * forward[r] / ((y[ind]) + 1e-15)
                if self.energy == 'on':
                    try:
                        rate_grad_energy[ind] += react[s] * ((beta / gas.T + E_a / (Gas_constant * gas.T ** 2)) * \
                                                 forward[r] - 1.0 * (beta2 / gas.T + E_a2 / (Gas_constant * gas.T ** 2)) * \
                                                 backward[r])
                    except:
                        rate_grad_energy[ind] = react[s] * ((beta / gas.T + E_a / (Gas_constant * gas.T ** 2)) * \
                                                forward[r] - 1.0 * (beta2 / gas.T + E_a2 / (Gas_constant * gas.T ** 2)) * \
                                                 backward[r])
                        # r] * vol * mwu_k[ind] * specific_heat * gas.T
                    test2_rate_grad += forward[r] * vol * mwu_k[ind] / gas.T
                    test3_rate_grad += -1 * react[s] * forward[r] * vol * mwu_k[ind] / gas.T

            for s in prod.keys():
                ind = gas.species_index(s)
                try:
                    rate_grad[ind] += -1.0 * prod[s] * backward[r] / ((y[ind]) + 1e-15)
                except:
                    rate_grad[ind] = -1.0 * prod[s] * backward[r] / ((y[ind]) + 1e-15)
                if self.energy == 'on':
                    # TODO: include wall thermal diffusion
                    try:
                        rate_grad_energy[ind] += -1.0 * prod[s] * ((beta / gas.T + E_a / (Gas_constant * gas.T ** 2)) * \
                                                forward[r] - 1.0 * (beta2 / gas.T + E_a2 / (Gas_constant * gas.T ** 2)) * \
                                                 backward[r])
                    except:
                        rate_grad_energy[ind] = -1.0 * prod[s] * ((beta / gas.T + E_a / (Gas_constant * gas.T ** 2)) * \
                                                forward[r] - 1.0 * (beta2 / gas.T + E_a2 / (Gas_constant * gas.T ** 2)) * \
                                                 backward[r])
                        # r] * vol * mwu_k[ind] * specific_heat * gas.T
                    test2_rate_grad -= backward[r] * vol * mwu_k[ind] / gas.T
                    test3_rate_grad -= prod[s] * backward[r] * vol * mwu_k[ind] / gas.T

            # Filling in the Jacobian
            # reactant stoichiometric coeff
            for s in react.keys():
                ind = gas.species_index(s)
                react_coeff = -1.0 * react[s]
                for spec in rate_grad:
                    gamma = mass / (vol * mw[spec])
                    # gradient of rate of production of ind wrt spec
                    Jw[ind][spec] += react_coeff * rate_grad[spec] * gamma
                    if self.energy == 'on':
                        # Energy  Conservation Equation: derivative of source term wrt species mass fraction
                        Jw_energy[ind][spec] += react_coeff * rate_grad[spec]* gamma * vol * hf[ind] *  mw[ind]
                if self.energy == 'on':
                    # Energy  Conservation Equation: derivative of source term wrt temperature
                    Jw_energy[ind][len(y)] += rate_grad_energy[ind] * gamma * vol * hf[ind] * mw[ind] # + \
                                              # react_coeff * forward[r] * vol * specific_heat[ind]

                    # Species  Conservation Equation: derivative of source term wrt temperature
                    Jw[ind][len(y)] += rate_grad_energy[ind]

            # product stoichiometric coeff
            for s in prod.keys():
                ind = gas.species_index(s)
                prod_coeff = prod[s]
                for spec in rate_grad:
                    gamma = mass / (vol * mw[spec])
                    # gradient of rate of production of ind wrt spec
                    Jw[ind][spec] += prod_coeff * rate_grad[spec] * gamma
                    if self.energy == 'on':
                        # Energy  Conservation Equation: derivative of source term wrt species mass fraction
                        Jw_energy[ind][spec] += prod_coeff * rate_grad[spec]* gamma * vol * hf[ind] * mw[ind]
                if self.energy == 'on':
                    # Energy  Conservation Equation: derivative of source term wrt temperature
                    Jw_energy[ind][len(y)] += rate_grad_energy[ind]* gamma * vol * hf[ind] * mw[ind]  # + \
                                              # prod_coeff * -1.0 * backward[r] * vol * specific_heat[ind]

                    # Species  Conservation Equation: derivative of source term wrt temperature
                    Jw[ind][len(y)] += rate_grad_energy[ind]

        if self.energy == 'on':
            # Energy  Conservation Equation: derivative of source term wrt temperature
            for i in range(len(gas.species_names)):
                Jw[i][len(y)] += rates[i]*specific_heat[i]*vol*hf[i]/abs(hf[i])
            Jw[len(y)][len(y)] = Jw[len(y)][len(y)] / self.T_norm
        return Jw, Jw_energy

    def jacobian_store(self, J, Jac, p, q, subvect_size, spec, spec0, outer_loop, row_Jac, col_Jac):
        """
        Stores the jacobian data with respective row and column ids used to form a sparse matrix later
        :param J:
        :param Jac:
        :param p:
        :param q:
        :param subvect_size:
        :param spec:
        :param spec0:
        :param outer_loop:
        :param row_Jac:
        :param col_Jac:
        :return:
        """
        if J != 0.0:
            # print "jw=", Jw_add
            # Jac = np.append(Jac, J)
            Jac.append(J)
            # Row in jacobian has constant numerator species
            if outer_loop == 1:
                row = p * subvect_size + spec
                col = q * subvect_size + spec
            else:
                row = p * subvect_size + spec0
                col = q * subvect_size + spec

            # print "rowj=", spec
            # row_Jac = np.append(row_Jac, row)
            # col_Jac = np.append(col_Jac, col)
            row_Jac.append(row)
            col_Jac.append(col)

    def jacobian_print(self, t0, Y0, species):
        #PSR_Jac, PSR_Jac_s, PSR_Jac_w = self.Jacobian_sparse(t0, Y0)
        PSR_Jac = self.Jacobian_sparse(t0, Y0)
        img = csc_matrix.todense(PSR_Jac)
        # img = csc_matrix.todense(self.mfc_gen.coeffmat)
        print "Jacobian temperature=", img[len(species), len(species)]
        # img = img[0:48,0:48]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.imshow(img[0:len(species)*200, 0:len(species)*200], vmin=-1e-7, vmax=1e-7, interpolation='none', cmap=plt.cm.gray, aspect='auto')
        fig.colorbar(im)
        # imgplot = plt.imshow(img)
        # plt.clim(-1, 1)
        plt.show()
        """det_Jac = linalg.det(img)
        print "Jacobian determinant = ", det_Jac
        print "Jacobian trace = ", matrix.trace(img)
        # FD Jacobian
        gas = self.PSR[0]
        # Jw = self.source_analytic(gas.thermo, gas.volume, gas.mass)
        fd_Jac = self.source_Jac_finitediff(gas.thermo, gas.volume, 0, self.chem_mech)
        print fd_Jac
        print "Finite diff Jac det = ", linalg.det(fd_Jac)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.imshow(fd_Jac, vmin=-1e-7, vmax=1e-7, interpolation='none', cmap=plt.cm.gray, aspect='auto')
        fig.colorbar(im)
        # imgplot = plt.imshow(img)
        # plt.clim(-1, 1)
        plt.show()
        samp_jac = [[1e-10, 0], [0, 1e-10]]
        print linalg.det(samp_jac)"""


    def Jacobian_PIVCRN(self, coeff, Jac, Jac_s, Jac_w, row_Jac, col_Jac, row_Jac_s, col_Jac_s, row_Jac_w, col_Jac_w, \
                        subvect_size, q, p, species, Jw, Jw_energy, mw, reactor_q, reactor_p):
        drho_dT = -1.0*(reactor_q.thermo.density / reactor_q.T)
        dT_drho = -1.0 * (reactor_q.T / reactor_q.thermo.density)
        mflux = coeff * self.massimb.massflowreactor[q]  # converting velocity flux to mass flux
        y = reactor_q.Y
        if p == q:
            outer_loop = subvect_size
            # print "coeff=", coeff
        else:
            outer_loop = 1
        for spec0 in range(outer_loop):
            # spec0 iss row index in Jacobian
            for spec in range(subvect_size):
                # spec is column index in Jacobian eand also row index for p!=q
                Js = 0.0
                Jw_add = 0.0
                # Flux terms
                if spec == spec0 or outer_loop == 1:
                    if spec0 == len(species) or (spec == len(species) and outer_loop == 1):
                        if spec == len(species):
                            # Energy  Conservation Equation: derivative of advective flux wrt temperature
                            Js += - 1.0 * mflux * reactor_q.thermo.cp_mass
                            # Energy  Conservation Equation: derivative of diffusive flux wrt temperature
                            try:
                                # wall_map: dictionary of wall objects mapped from reactor to neighbouring reactors.
                                S_pq = self.mfc_gen.wall_map[p][q].area
                                dx = (reactor_p.volume ** 0.333 + reactor_q.volume ** 0.333) / 2.0
                                K_react = (reactor_p.thermo.thermal_conductivity + reactor_q.thermo.thermal_conductivity)/2.0
                                rho = (self.mfc_gen.rho_calc[p] + self.mfc_gen.rho_calc[q])/2.0
                                Js_diffu = -1.0 * (coeff / abs(coeff)) * rho * K_react * S_pq / dx
                                Js += Js_diffu
                                # ToDo: Recording diffusion gradient in incoming
                                if self.track_jacobian == "on":
                                    self.jacobian_store(-Js_diffu, Jac, q, p, subvect_size, spec, spec0, outer_loop, row_Jac,
                                                    col_Jac)
                                if self.track_jacobian_compo == "on":
                                    self.jacobian_store(-Js_diffu, Jac_s, q, p, subvect_size, spec, spec0, outer_loop,
                                                        row_Jac_s, col_Jac_s)
                                # TODO: include turbulence diff
                            except:
                                pass
                        elif spec == len(species) + 1:
                            # Energy  Conservation Equation: derivative of advective flux wrt density
                            Js += -1.0 * reactor_q.thermo.cp_mass * reactor_q.T * coeff \
                            #      - reactor_q.thermo.density_mass * reactor_q.thermo.cp_mass * dT_drho * coeff
                            # Energy  Conservation Equation: derivative of diffusive flux wrt density
                            try:
                                # wall_map: dictionary of wall objects mapped from reactor to neighbouring reactors.
                                S_pq = self.mfc_gen.wall_map[p][q].area
                                dx = (reactor_p.volume ** 0.333 + reactor_q.volume ** 0.333) / 2.0
                                K_react =  (reactor_p.thermo.thermal_conductivity + reactor_q.thermo.thermal_conductivity)/2.0
                                rho = (self.mfc_gen.rho_calc[p] + self.mfc_gen.rho_calc[q])/2.0
                                Js_diffu = -1.0 * (coeff / abs(coeff)) * rho * K_react * S_pq / dx
                                Js += Js_diffu * dT_drho
                                # diffusion gradient in opposite direction
                                if self.track_jacobian == "on":
                                    self.jacobian_store(-Js_diffu * dT_drho, Jac, q, p, subvect_size, spec, spec0, outer_loop, row_Jac,
                                                    col_Jac)
                                if self.track_jacobian_compo == "on":
                                    self.jacobian_store(-Js_diffu * dT_drho, Jac_s, q, p, subvect_size, spec, spec0, outer_loop,
                                                        row_Jac_s, col_Jac_s)
                                # TODO: include turbulence diff
                            except:
                                pass
                    elif spec0 == len(species) + 1 or (spec == len(species) + 1 and outer_loop == 1):
                        if spec == len(species):
                            # Mass Conservation Equation: derivative of advective flux wrt temperature
                            Js += 0.0
                        elif spec == len(species) + 1:
                            # Mass Conservation Equation: derivative of advective flux wrt density
                            Js += -1.0 * coeff
                    else:
                        if spec == len(species):
                            # Species Conservation Equation: derivative of advective flux wrt temperature
                            Js += 0.0
                        elif spec == len(species) + 1:
                            # Species Conservation Equation: derivative of advective flux wrt density
                            Js += (-y[spec] * coeff)
                        else:
                            # Species Conservation Equation: derivative of advective flux wrt species mass fraction
                            Js += -1.0 * mflux
                            try:
                                # Species Conservation Equation: derivative of diffusive flux wrt species mass fraction
                                # wall_map: dictionary of wall objects mapped from reactor to neighbouring reactors.
                                S_pq = self.mfc_gen.wall_map[p][q].area
                                dx = (reactor_p.volume ** 0.333 + reactor_q.volume ** 0.333) / 2.0
                                D_spec = (reactor_p.thermo.mix_diff_coeffs_mass[spec0] + reactor_q.thermo.mix_diff_coeffs_mass[spec0])/2.0
                                rho = (self.mfc_gen.rho_calc[p] + self.mfc_gen.rho_calc[p])/2.0
                                Js_diffu = -1.0 * (coeff / abs(coeff)) * rho * D_spec * S_pq / dx
                                Js += Js_diffu
                                # diffusion gradient in opposite direction
                                if self.track_jacobian == "on":
                                    self.jacobian_store(-Js_diffu, Jac, q, p, subvect_size, spec, spec0, outer_loop, row_Jac,
                                                    col_Jac)
                                if self.track_jacobian_compo == "on":
                                    self.jacobian_store(-Js_diffu, Jac_s, q, p, subvect_size, spec, spec0, outer_loop,
                                                        row_Jac_s, col_Jac_s)
                                # TODO: include turbulence diff
                            except:
                                pass
                # Source terms
                if spec0 == len(species) and p == q:
                    if spec == len(species):
                        # Energy  Conservation Equation: derivative of source term wrt temperature
                        Jw_add -= np.sum([Jw_energy[a, spec] for a in range(len(species))])
                    elif spec == len(species) + 1:
                        # Energy  Conservation Equation: derivative of source term wrt density
                        Jw_add -= np.sum([Jw_energy[a, spec - 1] * dT_drho for a in range(len(species))])
                    else:
                        # Energy  Conservation Equation: derivative of source term wrt species mass fraction
                        Jw_add -= np.sum([Jw_energy[a, spec] for a in range(len(species))])
                if spec0 < len(species) and p == q:
                    if spec == len(species):
                        # Species Conservation Equation: derivative of source term wrt temperature
                        Jw_add += Jw[spec0, spec] * mw[spec0] * reactor_p.volume
                    elif spec == len(species) + 1:
                        # Species Conservation Equation: derivative of source term wrt density
                        Jw_add += Jw[spec0, spec - 1] * mw[spec0] * reactor_p.volume * dT_drho
                    else:
                        # Species Conservation Equation: derivative of source term wrt species mass fraction
                        # Column in species source term jacobian has constant denominator species
                        Jw_add += Jw[spec0, spec] * mw[spec0] * reactor_p.volume

                # Jacobian assembly
                J = (Js + Jw_add)  # /(abs(divisor[spec0+(self.mfc_gen.vectsize*p)])+self.tol)
                if self.track_jacobian == "on":
                    self.jacobian_store(J, Jac, p, q, subvect_size, spec, spec0, outer_loop, row_Jac, col_Jac)
                if self.track_jacobian_compo == "on":
                    self.jacobian_store(Js, Jac_s, p, q, subvect_size, spec, spec0, outer_loop, row_Jac_s, col_Jac_s)
                    self.jacobian_store(Jw_add, Jac_w, p, q, subvect_size, spec, spec0, outer_loop, row_Jac_w,
                                        col_Jac_w)

    def Jacobian_CFDCRN(self, coeff, Jac, Jac_s, Jac_w, row_Jac, col_Jac, row_Jac_s, col_Jac_s, row_Jac_w, col_Jac_w, \
                        subvect_size, q, p, species, Jw, Jw_energy, mw):
        drho_dT = -1.0 * (self.PSR[q].thermo.density / self.PSR[q].T)
        dT_drho = -1.0 * (self.PSR[q].T / self.PSR[q].thermo.density)
        mflux = coeff * self.massimb.massflowreactor[q]
        y = self.PSR[q].Y
        if p == q:
            outer_loop = subvect_size
            # print "coeff=", coeff
        else:
            outer_loop = 1
        for spec0 in range(outer_loop):
            # spec0 iss row index in Jacobian
            for spec in range(subvect_size):
                # spec is column index in Jacobian eand also row index for p!=q
                Js = 0.0
                Jw_add = 0.0
                # Flux terms
                if spec == spec0 or outer_loop == 1:
                    if spec0 == len(species) or (spec == len(species) and outer_loop == 1):
                        if spec == len(species):
                            # Energy  Conservation Equation: derivative of advective flux wrt temperature
                            Js += - 1.0 * mflux * self.PSR[q].thermo.cp_mass
                            # Energy  Conservation Equation: derivative of diffusive flux wrt temperature
                            try:
                                # wall_map: dictionary of wall objects mapped from reactor to neighbouring reactors.
                                S_pq = self.mfc_gen.wall_map[p][q].area
                                dx = (self.PSR[p].volume ** 0.333 + self.PSR[q].volume ** 0.333) / 2.0
                                K_react = (self.PSR[p].thermo.thermal_conductivity + self.PSR[
                                    q].thermo.thermal_conductivity) / 2.0
                                rho = (self.mfc_gen.rho_calc[p] + self.mfc_gen.rho_calc[q]) / 2.0
                                Js_diffu = -1.0 * (coeff / abs(coeff)) * rho * K_react * S_pq / dx
                                Js += Js_diffu
                                # ToDo: Recording diffusion gradient in incoming
                                if self.track_jacobian == "on":
                                    self.jacobian_store(-Js_diffu, Jac, q, p, subvect_size, spec, spec0, outer_loop,
                                                    row_Jac, col_Jac)
                                if self.track_jacobian_compo == "on":
                                    self.jacobian_store(-Js_diffu, Jac_s, q, p, subvect_size, spec, spec0, outer_loop,
                                                        row_Jac_s, col_Jac_s)
                                # TODO: include turbulence diff
                            except:
                                pass
                        elif spec == len(species) + 1:
                            # Energy  Conservation Equation: derivative of advective flux wrt density
                            Js += -1.0 * self.PSR[q].thermo.cp_mass * self.PSR[q].T * coeff \
                                #      - self.PSR[q].thermo.density_mass * self.PSR[q].thermo.cp_mass * dT_drho * coeff
                            # Energy  Conservation Equation: derivative of diffusive flux wrt density
                            try:
                                # wall_map: dictionary of wall objects mapped from reactor to neighbouring reactors.
                                S_pq = self.mfc_gen.wall_map[p][q].area
                                dx = (self.PSR[p].volume ** 0.333 + self.PSR[q].volume ** 0.333) / 2.0
                                K_react = (self.PSR[p].thermo.thermal_conductivity + self.PSR[
                                    q].thermo.thermal_conductivity) / 2.0
                                rho = (self.mfc_gen.rho_calc[p] + self.mfc_gen.rho_calc[q]) / 2.0
                                Js_diffu = -1.0 * (coeff / abs(coeff)) * rho * K_react * S_pq / dx
                                Js += Js_diffu * dT_drho
                                # diffusion gradient in opposite direction
                                if self.track_jacobian == "on":
                                    self.jacobian_store(-Js_diffu * dT_drho, Jac, q, p, subvect_size, spec, spec0,
                                                    outer_loop, row_Jac, col_Jac)
                                if self.track_jacobian_compo == "on":
                                    self.jacobian_store(-Js_diffu * dT_drho, Jac_s, q, p, subvect_size, spec, spec0,
                                                        outer_loop,
                                                        row_Jac_s, col_Jac_s)
                                # TODO: include turbulence diff
                            except:
                                pass
                    elif spec0 == len(species) + 1 or (spec == len(species) + 1 and outer_loop == 1):
                        if spec == len(species):
                            # Mass Conservation Equation: derivative of advective flux wrt temperature
                            Js += 0.0
                        elif spec == len(species) + 1:
                            # Mass Conservation Equation: derivative of advective flux wrt density
                            Js += -1.0 * coeff
                    else:
                        if spec == len(species):
                            # Species Conservation Equation: derivative of advective flux wrt temperature
                            Js += 0.0
                        elif spec == len(species) + 1:
                            # Species Conservation Equation: derivative of advective flux wrt density
                            Js += (-y[spec] * coeff)
                        else:
                            # Species Conservation Equation: derivative of advective flux wrt species mass fraction
                            Js += -1.0 * mflux
                            try:
                                # Species Conservation Equation: derivative of diffusive flux wrt species mass fraction
                                # wall_map: dictionary of wall objects mapped from reactor to neighbouring reactors.
                                S_pq = self.mfc_gen.wall_map[p][q].area
                                dx = (self.PSR[p].volume ** 0.333 + self.PSR[q].volume ** 0.333) / 2.0
                                D_spec = (self.PSR[p].thermo.mix_diff_coeffs_mass[spec0] +
                                          self.PSR[q].thermo.mix_diff_coeffs_mass[spec0]) / 2.0
                                rho = (self.mfc_gen.rho_calc[p] + self.mfc_gen.rho_calc[p]) / 2.0
                                Js_diffu = -1.0 * (coeff / abs(coeff)) * rho * D_spec * S_pq / dx
                                Js += Js_diffu
                                # diffusion gradient in opposite direction
                                if self.track_jacobian == "on":
                                    self.jacobian_store(-Js_diffu, Jac, q, p, subvect_size, spec, spec0, outer_loop,
                                                    row_Jac, col_Jac)
                                if self.track_jacobian_compo == "on":
                                    self.jacobian_store(-Js_diffu, Jac_s, q, p, subvect_size, spec, spec0, outer_loop,
                                                        row_Jac_s, col_Jac_s)
                                # TODO: include turbulence diff
                            except:
                                pass
                # Source terms
                if spec0 == len(species) and p == q:
                    if spec == len(species):
                        # Energy  Conservation Equation: derivative of source term wrt temperature
                        Jw_add -= np.sum([Jw_energy[a, spec] for a in range(len(species))])
                    elif spec == len(species) + 1:
                        # Energy  Conservation Equation: derivative of source term wrt density
                        Jw_add -= np.sum([Jw_energy[a, spec - 1] * dT_drho for a in range(len(species))])
                    else:
                        # Energy  Conservation Equation: derivative of source term wrt species mass fraction
                        Jw_add -= np.sum([Jw_energy[a, spec] for a in range(len(species))])
                if spec0 < len(species) and p == q:
                    if spec == len(species):
                        # Species Conservation Equation: derivative of source term wrt temperature
                        Jw_add += Jw[spec0, spec] * mw[spec0] * self.PSR[p].volume
                    elif spec == len(species) + 1:
                        # Species Conservation Equation: derivative of source term wrt density
                        Jw_add += Jw[spec0, spec - 1] * mw[spec0] * self.PSR[p].volume * dT_drho
                    else:
                        # Species Conservation Equation: derivative of source term wrt species mass fraction
                        # Column in species source term jacobian has constant denominator species
                        Jw_add += Jw[spec0, spec] * mw[spec0] * self.PSR[p].volume

                # Jacobian assembly
                J = (Js + Jw_add)  # /(abs(divisor[spec0+(self.mfc_gen.vectsize*p)])+self.tol)
                if self.track_jacobian == "on":
                    self.jacobian_store(J, Jac, p, q, subvect_size, spec, spec0, outer_loop, row_Jac, col_Jac)
                if self.track_jacobian_compo == "on":
                    self.jacobian_store(Js, Jac_s, p, q, subvect_size, spec, spec0, outer_loop, row_Jac_s, col_Jac_s)
                    self.jacobian_store(Jw_add, Jac_w, p, q, subvect_size, spec, spec0, outer_loop, row_Jac_w,
                                        col_Jac_w)

    def jacobian_parallel(self,p,colptr,rowind, Jac, Jac_s, Jac_w, row_Jac, col_Jac, row_Jac_s, col_Jac_s, row_Jac_w,
                                     col_Jac_w, subvect_size,species):
        gas = self.PSR[p]
        Jw, Jw_energy = self.source_analytic(gas.thermo, gas.volume, gas.mass)
        mw = list(gas.thermo.molecular_weights)
        rowrange = range(colptr[p], colptr[p + 1])
        for quo in rowrange:
            # elem = array([..list of elements..], dtype=...)
            coeff = self.mfc_gen.coeffmat.data[quo]
            q = rowind[quo]
            # flow from reactor p to reactor q
            if self.process == "PIVCRN":
                self.Jacobian_PIVCRN(coeff, Jac, Jac_s, Jac_w, row_Jac, col_Jac, row_Jac_s, col_Jac_s, row_Jac_w,
                                     col_Jac_w, subvect_size, p, q, species, Jw, Jw_energy, mw, self.PSR[p],self.PSR[q])
            elif self.process == "CFDCRN":
                self.Jacobian_CFDCRN(coeff, Jac, Jac_s, Jac_w, row_Jac, col_Jac, row_Jac_s, col_Jac_s, row_Jac_w,
                                     col_Jac_w, subvect_size, p, q, species, Jw, Jw_energy, mw)


    def Jacobian_sparse(self, t, Y):
        """

        :param t:
        :param Y:
        :param species:
        :return:
        """
        # Todo: This function has already been converted to numpy arrays
        """Jac = np.array([])
        row_Jac = np.array([])
        col_Jac = np.array([])"""
        Jac = []
        Jac_s = []
        Jac_w = []
        row_Jac = []
        col_Jac = []
        row_Jac_s = []
        col_Jac_s = []
        row_Jac_w = []
        col_Jac_w = []

        species = self.PSR[0].thermo.species_names
        Y_temp = self.PSR_insert_gas(Y, species)  # required to update species in reactors
        colptr = self.mfc_gen.coeffmat.indptr
        rowind = self.mfc_gen.coeffmat.indices
        if self.process == "PIVCRN":
            vectorsize = (self.mfc_gen.vectsize) * len(self.PSR)  # extra term for mass conservation equation
            subvect_size = self.mfc_gen.vectsize

        else:
            vectorsize = (self.mfc_gen.vectsize) * len(self.PSR)
            subvect_size = self.mfc_gen.vectsize
        print "Jacobian Calculation"
        # divisor = self.RHS_Jac0
        num_cores = multiprocessing.cpu_count()
        #print "Parallel cores=",num_cores
        """Parallel(n_jobs=num_cores-1, backend="threading")(
            delayed(self.jacobian_parallel)(p,colptr,rowind, Jac, Jac_s, Jac_w, row_Jac, col_Jac, row_Jac_s, col_Jac_s, row_Jac_w,
                                     col_Jac_w, subvect_size,species) for p in range(len(self.PSR)))"""
        for p in range(len(self.PSR)):
            # ToDo: Parallelise this operation
            gas = self.PSR[p]
            Jw, Jw_energy = self.source_analytic(gas.thermo, gas.volume, gas.mass)
            mw = list(gas.thermo.molecular_weights)
            rowrange = range(colptr[p], colptr[p + 1])
            for quo in rowrange:
                # elem = array([..list of elements..], dtype=...)
                coeff = self.mfc_gen.coeffmat.data[quo]
                q = rowind[quo]
                # flow from reactor p to reactor q
                if self.process == "PIVCRN":
                    self.Jacobian_PIVCRN(coeff, Jac, Jac_s, Jac_w, row_Jac, col_Jac, row_Jac_s, col_Jac_s, row_Jac_w,
                                         col_Jac_w, subvect_size, p, q, species, Jw, Jw_energy, mw, self.PSR[p],self.PSR[q])
                elif self.process == "CFDCRN":
                    self.Jacobian_CFDCRN(coeff, Jac, Jac_s, Jac_w, row_Jac, col_Jac, row_Jac_s, col_Jac_s, row_Jac_w,
                                         col_Jac_w, subvect_size, p, q, species, Jw, Jw_energy, mw)

        print "Length row = ", len(row_Jac)
        print "Length column = ", len(col_Jac)
        print "Length Jac = ", len(Jac)
        PSR_Jac = csc_matrix((Jac, (row_Jac, col_Jac)), shape=(vectorsize, vectorsize), dtype=double)
        PSR_Jac_s = []
        PSR_Jac_w = []
        if self.track_jacobian_compo == "on":
            PSR_Jac_s = csc_matrix((Jac_s, (row_Jac_s, col_Jac_s)), shape=(vectorsize, vectorsize), dtype=double)
            PSR_Jac_w = csc_matrix((Jac_w, (row_Jac_w, col_Jac_w)), shape=(vectorsize, vectorsize), dtype=double)
        return PSR_Jac #, PSR_Jac_s, PSR_Jac_w


    def localsolver_section(self, PSRInd_updated, Y0, error_cum, tstart_local, dt_fac, dt):
        """

        :param PSRInd_updated:
        :param tres:
        :param error_cum:
        :param tstart_local:
        :param dt_fac:
        :param dt:
        :return:
        """
        # Solving individual reactors
        if self.energy == 'on':
            local_length = len(self.PSR) * 2
        else:
            local_length = len(self.PSR)

        tres = self.massimb.tres
        tol = 1e-15
        for i in range(1):
            error, traverse, PSRInd_updated = self.localsolver(PSRInd_updated, error_cum, tstart_local, dt_fac)
            Y1 = self.PSR_gas_ini()
            t0= dt
            D_Y =subtract(Y1,Y0)
            err = np.divide(np.abs(D_Y), (np.array(Y0) + self.tol))
            max_err = np.ndarray.max(err)
            max_ind = np.argmax(err)
            spec_num = max_ind % self.subvect_size
            self.RHS_Jac0 = self.objective_func(t0, Y1)
            self.g_omega = np.ndarray.max(np.abs(self.RHS_Jac0))  # np.linalg.norm(RHS_Jac1)/len(RHS_Jac1)
            print "Max rate = ", self.g_omega
            self.epsilon_omega = (max_err - self.epsilon) / (self.epsilon + tol)
            self.epsilon = max_err
            self.epsilon_g_omega = (abs(self.g_omega) - abs(self.g_old)) / (abs(self.g_omega) + tol)
            self.error_rate_log.append(self.g_omega)
            self.error_log.append(max_err)
            self.error_spec_log.append(spec_num)
            self.solver_type.append(0)
            print "error local=", error
            # Checking valve mass flows
            m_valve = 0
            m_v_max = 0.0

            for v in self.mfc_gen.ValveRecord:
                mfc_v = self.mfc_gen.ValveRecord[v]
                rel_massflow_valve = v.mdot(0) / mfc_v.mdot(0)
                if rel_massflow_valve > m_v_max:
                    m_v_max = rel_massflow_valve
                    m_valve = v
            print "Max valve flow=", m_v_max

            # To check valve operation
            """try:
                mfc_max = self.mfc_gen.ValveRecord[m_valve]
                reactor_max = self.mfc_gen.mfc_rec[mfc_max][0]
                reactor_max_to = self.mfc_gen.mfc_rec[mfc_max][1]
                print "Reactor max to=", reactor_max_to
                if (reactor_max in self.mfc_gen.mfc_per) or (reactor_max_to in self.mfc_gen.mfc_per):
                    print "In periodic"
                print "Max reactor=", reactor_max
                print "Max reactor pressure=", self.PSR[reactor_max].thermo.P
                print "Max reactor to=", reactor_max_to
                print "Max reactor to pressure=", self.PSR[reactor_max_to].thermo.P
            except:
                print "Inlet Reactor Valve" """

            print "dt_fac=", dt_fac
            if traverse != len(self.PSR):
                dt = dt / 10
            tmin = min(tres)
            tmin_dtfac = dt_fac * tmin
            sum_tres = sum(tres)
            if error < 0.1 and tmin_dtfac < sum_tres:
                print "res min fac=", tmin_dtfac
                print "sum max=", sum_tres
                dt_fac += 1
            if error < 1e-04 or traverse == len(self.PSR):
                break

        if (error < 1e-06) or dt_fac >= 0.85 * len(self.PSR):
            # Correcting Valve Coefficients based on related MFCs
            for v in self.mfc_gen.ValveRecord:
                v.set_valve_coeff(0)

        return Y1, m_v_max, error

    def localsolver(self, PSRInd_updated, error_cum, tstart_local, dt_fac):
        """

        :param PSRInd_updated:
        :param error_cum:
        :param tstart_local:
        :param dt_fac:
        :return:
        """

        traverse = 0
        e_max = self.tol

        for r_ind in range(len(PSRInd_updated)):
            r = self.PSR[r_ind]
            psr_index = r_ind
            net = self.network[r_ind]
            try:
                xin = list(r.thermo.Y)
                if self.energy == 'on':
                    xin = list(list(r.thermo.Y) + [r.thermo.T])
                dt_loc = self.massimb.tres[psr_index]
                net.advance(dt_loc * dt_fac)
                xfin = list(r.thermo.Y)
                if self.energy == 'on':
                    xfin = list(list(r.thermo.Y) + [r.thermo.T])

                # updating reactor residence time
                if self.process == "PIVCRN":
                    self.massimb.tres[psr_index] = r.mass / (self.massimb.massflowreactor[psr_index] * \
                                                             self.PSR[psr_index].thermo.density)
                else:
                    self.massimb.tres[psr_index] = r.mass / self.massimb.massflowreactor[psr_index]
                # func = [(xfin[n] - xin[n]) for n in range(len(xin))]
                func = r.thermo.net_production_rates
                e = []
                e_div = max(xin)

                for ind in range(len(func)):
                    error_cum += (abs(func[ind]) / e_div) ** 2
                    e.append(abs(func[ind]))  # / e_div)
                e_rel = max(e)

                if e_rel > e_max:
                    e_max = e_rel
                tstart_local[psr_index] += dt_loc
                traverse += 1

                for num in self.mfc_gen.Res_dict[r_ind]:
                    g = r.thermo
                    num[0].insert(g)

            except RuntimeError:
                #print "marker"
                """div_list.append(r_ind)
                res = Res_dict[r_ind][0]
                g = Reserv[res][0].thermo
                Reserv[res][1].insert(g)"""
                continue

        error = e_max

        # Residual analysis
        if traverse != len(self.PSR):
            divg = len(self.PSR) - traverse
            print "Diverging reactors=", divg
            error = 1.0
        print "Cummu error=", error_cum

        return error, traverse, PSRInd_updated

    def massflow_const(self, species):
        # Todo: This function has already been converted to numpy arrays
        print "Evaluating constant massflows"
        if self.process == "PIVCRN":
            vectorsize = (self.mfc_gen.vectsize) * len(self.PSR)  # extra term for mass conservation equation
            subvect_size = self.mfc_gen.vectsize

        else:
            vectorsize = (self.mfc_gen.vectsize) * len(self.PSR)
            subvect_size = self.mfc_gen.vectsize
        self.f_vect = np.zeros(subvect_size * len(self.PSR))
        for r_ind in range(len(self.PSR)):
            # TODO: This part could be parallelised
            r = self.PSR[r_ind]
            mass_curr = r.mass
            mass_sum = 0
            for ind in range(self.mfc_gen.vectsize):
                """if self.process == "PIVCRN":
                    speciesflux = self.mfc_gen.fvect[r_ind * self.mfc_gen.vectsize + ind] * self.PSR[
                        r_ind].thermo.density
                else:
                    speciesflux = self.mfc_gen.fvect[r_ind * self.mfc_gen.vectsize + ind]"""
                speciesflux = self.mfc_gen.fvect[r_ind * self.mfc_gen.vectsize + ind]
                self.f_vect[r_ind * subvect_size + ind] = speciesflux
                """if ind == len(species):
                    # energy equation
                    self.f_vect[r_ind * subvect_size + ind] = self.f_vect[r_ind * subvect_size + ind] / self.T_norm"""
                if ind < len(species):
                    mass_sum += speciesflux

        self.fvect_Cw = array(self.f_vect)
        # print self.f_vect
        return 0

    def contraint_check(self, X, X0, len_species):
        """
        Checks constraints of 0<y<1, Tmin<T<Tmax and Sum(Y)=1
        :param X: X0+dx
        :param X0: initial value of vector
        :param len_species: number of species being considered
        :return: X2
        """
        sum = {}
        X2 = array(X)
        for i in range(len(X)):
            p = int(math.floor(i / self.subvect_size))
            spec = i % self.subvect_size
            if spec == 0:
                sum[p] = 0
            if spec == len_species and self.process == "PIVCRN":
                P = ct.one_atm
                Ru = 8314.46261815324  # J/K/kmol
                MW = self.PSR[p].thermo.mean_molecular_weight
                d = self.mfc_gen.rho_calc[p]
                T_calc = P / (d*Ru / MW)
                calc_wt = 0.6
                if abs(T_calc-X2[i])> 500:
                    X2[i] = X2[i] * (1.0 - calc_wt) + T_calc * calc_wt

            if X[i] > self.Y_high[i]:
                X2[i] = X0[i] + 0.5 * (self.Y_high[i] - X0[i])

            if X[i] < self.Y_low[i]:
                X2[i] = X0[i] - 0.5 * (X0[i] - self.Y_low[i])
            if spec < len_species:
                sum[p] += X2[i]

            else:
                pass

        # Sum of mass fractions in reactor = 1
        for j in range(len(X2) / self.subvect_size):
            for k in range(len_species):
                X2[j * self.subvect_size + k] = X2[j * self.subvect_size + k] / sum[j]

        return X2

    def Armijo_Goldstein_rule(self, xin, dx, fx_in, J, t0, beta, len_species):
        """
        Calculates the effective step size based on Armijo-Goldstein rule
        :param xin: initial value of x
        :param dx: descent direction
        :param fx_in: objective function at initial value of x
        :param J: Jacobian at xin
        :param t0:
        :param beta: factor used to calculate alpha in Armijo-Goldstein rule
        :return: x, f(x) or False, False
        """
        x1 = self.contraint_check(np.add(xin,dx), xin, len_species)  # checking for constraints
        dx1 = np.subtract(x1,xin)
        #f1 = self.objective_func(t0, x1)
        #print f1
        J_transp = J.transpose()
        descent = J_transp.dot(np.array(dx1))
        print descent.shape
        #print descent#descent = descent.todense()
        sigma = 0.5 # (0 < sigma < 0.5) factor in Armijo condition
        mu = 1.0 # (0.5 < mu < 1) factor in Goldstein condition
        alpha = np.ones(len(xin))
        rmax = 10
        for i in range(rmax):
            x_test0 = add(xin,multiply(alpha,dx1))
            x_test = self.contraint_check(x_test0, xin, len_species)  # checking for constraints
            f_test = self.objective_func(t0, x_test)
            f_condn_upper = fx_in + sigma * multiply(alpha,descent)
            f_condn_lower = fx_in + mu * multiply(alpha,descent)
            armijo_condn = subtract(f_test,f_condn_upper)  # inv condn; Armijo when f_test<=f_condn
            # goldstein_condn = np.less(f_test,f_condn_lower)  # inv condn; Goldstein when f_test>=f_condn
            count = 0
            for j in range(len(armijo_condn)):
               # print "Armijo\n"
                #print armijo_condn
                #print x_test
                # as function is not positive definite and Armijo rule is for function minimisation,
                # we check the condition such that we minimise the absolute value
                if (armijo_condn[j] > 0.0 and fx_in[j]>0.0) or (armijo_condn[j] < 0.0 and fx_in[j]<0.0) :#descent[j]<0.0:# or any(goldstein_condn):
                    #print descent[j]
                    count += 1
                    alpha[j] = beta ** i + 1
                    """if i == rmax-2:
                        alpha[j] = 0.0
                    else:
                        alpha[j] = beta ** i + 1"""

            #print count
            if count == 0:
                print "iter armijo = ", i
                #print "ftest=", f_test
                return x_test, f_test

            if i == rmax-1:
                print "iter armijo = ", i
                x_test0 = add(xin, (beta ** 1.0)*multiply(alpha, dx1))
                x_test = self.contraint_check(x_test0, xin, len_species)  # checking for constraints
                f_test = self.objective_func(t0, x_test)
                print "Max diff ", np.max(np.abs(np.subtract(x_test0,xin)))
                return x_test, f_test

        return [0], [0]

    def Newton_sparse(self, g, J, x0, t0, alpha, len_species):
        """
        Newton's method on sparse matrix
        :param g:
        :param J:
        :param x0:
        :param alpha:
        :return:
        """
        print "Solving "
        f = -1.0*np.array(g)  # divide(array(g), add(absolute(g),self.tol))
        # J_norm = J/array(g)
        delta_x0 = spsolve(A=J, b=f, use_umfpack=True)
        #val = np.allclose(J.dot(delta_x0), f)
        #print "Closure=", val
        exitcode = 0.0
        if exitcode != 0.0:
            print "Exit code=", exitcode
        #print array(delta_x0).shape
        #print array(delta_x0)
        beta = alpha  # factor used to calculate alpha in Armijo-Goldstein rule
        for i in range(5):
            xfinal, obj_func_final = self.Armijo_Goldstein_rule(x0, array(delta_x0), g, J, t0, beta, len_species)
            if len(xfinal) == 1:
                beta = beta / 2.0
                print "Beta=", beta
            else:
                break

        if len(xfinal) == 1:
            xfinal = x0
            obj_func_final = self.objective_func(t0, x0)
        spec_count = 0
        spec_sum = 0
        """for i in range(len(xfinal)):
            spec_sum += xfinal[i]
            spec_count += 1
            if spec_count == len_species:
                spec_count = 0
                xfinal[i - len_species:i] = xfinal[i - len_species:i] / spec_sum
                spec_sum = 0"""

        #xfinal =np.add(x0,beta*delta_x0)
        #obj_func_final = self.objective_func(t0, xfinal)
        delta_xf = np.subtract(xfinal, x0)
        print "Max diff newton=", np.max(np.abs(delta_xf))
        # delta_x = delta_xf
        print "delta x=", list(delta_xf)
        print "delta x temp = ", delta_xf[len_species]
        #print obj_func_final
        # xfinal = divide(add(xfinal, x0), 2.0)
        return delta_xf, xfinal, obj_func_final, exitcode

    def globalsolver(self, t0, Y0, dt, species, g_old):
        """
        Global solver that currently calls a global newton's method on a sparse system of equations
        :param t0:
        :param Y0:
        :param dt:
        :param species:
        :param g_old:
        :return:
        """
        tol = 1e-15
        print "Newton Solver"
        PSR_Jac = self.Jacobian_sparse(t0, Y0)
        #a = PSR_Jac.data
        #print a
        # RHS_Jac0 = self.objective_func(t0, Y0, species)
        Yinitial = array(Y0)
        #self.jacobian_print(t0, Y0, species)
        # self.alpha = 0.1
        D_Y, Yfinal, RHS_Jac1, exitcode = self.Newton_sparse(self.RHS_Jac0, PSR_Jac, Yinitial, t0, self.alpha, len(species))
        delta_RHS = np.subtract(RHS_Jac1,self.RHS_Jac0)
        self.RHS_Jac0 = RHS_Jac1
        psr_max = 0
        psr_max_epsi = 0
        max_ind_epsi = 0
        #print D_Y
        err = np.divide(np.abs(D_Y), (np.array(Yinitial) + tol))
        max_err = np.max(err)
        max_ind = np.argmax(err)
        spec_num = max_ind % self.subvect_size

        self.max_obj_rate = psr_max * self.mfc_gen.vectsize + max_ind
        print "Reactor num = ", math.floor(max_ind/self.subvect_size)
        try:
            print "Max Species = ", species[spec_num]
        except:
            print "Max Species = Temperature!!"
        print "Number Species = ", len(species)
        print "Vectsize =", self.mfc_gen.vectsize

        tstep = t0 + dt
        self.g_omega = np.linalg.norm(RHS_Jac1)/len(RHS_Jac1)#np.ndarray.max(np.abs(RHS_Jac1))#
        print "Max rate = ", self.g_omega
        max_react = np.argmax(np.abs(RHS_Jac1))
        print "Max rate react= ", max_react
        print "Temperature max react= ", Yfinal[max_react]
        print "Energy rate: BC=",self.fvect_Cw[max_react]
        print "Energy rate: Cw=",self.Cw[max_react]
        print "Energy rate: Production=", self.RatesofProduction[max_react]
        print "Energy rate: Diffusion=", self.J_diffu[max_react]
        positive_rates = np.where(delta_RHS>0.0)
        #if len(positive_rates) == 0:
        self.epsilon_omega = (max_err - self.epsilon) / (self.epsilon + tol)
        self.epsilon = max_err
        self.epsilon_g_omega = (abs(self.g_omega) - abs(self.g_old))/ (abs(self.g_omega)+tol)
        self.error_rate_log.append(self.g_omega)
        self.error_log.append(max_err)
        self.error_spec_log.append(spec_num)
        self.max_obj = psr_max_epsi * self.mfc_gen.vectsize + max_ind_epsi
        self.solver_type.append(1)

        return array(Yfinal), tstep, max_err

    def globalsolver_odeint(self, t0, Y0, tmin, tmax, species):
        """
        Global time integrator: needs to be improved as not working properly currently
        :param t0:
        :param Y0:
        :param tmax:
        :param PSR:
        :param species:
        :param coeffmat:
        :param massflowreactor:
        :param fvect:
        :param vectsize:
        :param wall_map:
        :param chem_mech:
        :param energy:
        :param error_log:
        :param solver_type:
        :return:
        """

        print "Time-stepping ODE Integration"
        print "Time=", t0
        # print "Vol before bdf = ", self.PSR[0].volume
        Yinitial = np.array(Y0)

        func = lambda t, y: self.objective_func(t, y)  # , add(absolute(self.RHS_Jac0),self.tol))
        Jac = lambda t, y: self.Jacobian_sparse(t, y)
        #Jac, Jac_s, Jac_w = self.Jacobian_sparse(t0, Yinitial)
        # Jac_in = self.Jacobian_sparse(t0, Yinitial, species)
        # objfunc_in = self.objective_func(t0, Yinitial, species)
        # step_initial = 1/csc_matrix.max(Jac_in)
        # print "Initial time step =",step_initial
        r = BDF(func, t0=t0, y0=Yinitial, t_bound=t0+1e-03, rtol=1e-12,
                atol=1e-15, jac=Jac)

        tol = 1e-15
        e_max = 0
        Yfinal = None
        # initial mass fractions and temperatures
        xin = []
        xfin = []

        for p in self.PSR:
            if self.energy == 'on':
                xin += list(list(p.thermo.Y) + [p.thermo.T])
            else:
                xin += list(p.thermo.Y)
        count = 0
        t_lim = t0 + tmax
        while r.status == 'running' and count<200:# and r.t < t_lim:  # and self.g_omega > 1e-12 and count<1:  # abs(e_prev - e_max) > 1e-8:
            #print "Count = ", count
            count += 1
            xin = []

            for p in self.PSR:
                if self.energy == 'on':
                    xin += list(list(p.thermo.Y) + [p.thermo.T])
                else:
                    xin += list(p.thermo.Y)
            msg = r.step()
            print "Message=", msg
            print "Status=",r.status
            print "Step-size=",r.step_size
            print "Time=", r.t
            time_end = r.t
            Yfinal = self.contraint_check(r.y, xin, len(species))
            #self.jacobian_print(t0, Yfinal, species)
            r.y = Yfinal
            # print r.y
            xfin = Yfinal
            """for p in self.PSR:
                if self.energy == 'on':
                    xfin += list(list(p.thermo.Y) + [p.thermo.T])
                else:
                    xfin += list(p.thermo.Y)"""
            func = [(xfin[n] - xin[n]) for n in range(len(xin))]  # change in state
            objfunc0 = r.fun(r.t, r.y)
            objfunc = (objfunc0[len(objfunc0) - len(xin):len(objfunc0)])
            self.RHS_Jac0 = objfunc
            objfunc = abs(objfunc)
            err = np.divide(np.abs(func),(np.array(xin) + tol))
            max_err = np.ndarray.max(err)
            max_ind = np.argmax(err)
            spec_num = max_ind % self.subvect_size

            # print "Max species= ", species[max_ind]

            # Logging iteration convergence error
            self.g_omega = np.linalg.norm(objfunc)/len(objfunc)
            self.epsilon_g_omega = (abs(self.g_omega) - abs(self.g_old))
            self.epsilon_g_omega_rel = (abs(self.g_omega) - abs(self.g_old)) / abs(self.g_omega)
            self.epsilon_omega = (max_err - self.epsilon) / (self.epsilon + tol)
            self.epsilon = max_err  # maximum magnitude of change in state
            self.error_rate_log.append(self.g_omega)
            self.error_log.append(self.epsilon)
            self.error_spec_log.append(spec_num)
            e_max = self.g_omega
            self.solver_type.append(2)
            print "Epsilon=",self.epsilon
            print "G omega=",self.g_omega

        print 'End global'

        return Yfinal, time_end, e_max

    def PSR_gas_ini(self):
        Yfinal = []
        for p in range(len(self.PSR)):
            if self.process == "CFDCRN":
                if self.energy == 'on':
                    Yfinal += list(list(self.PSR[p].thermo.Y) + [self.PSR[p].thermo.T])
                else:
                    Yfinal += list(self.PSR[p].thermo.Y)
            elif self.process == "PIVCRN":
                if self.energy == 'on':
                    Yfinal += list(
                        list(self.PSR[p].thermo.Y) + [self.PSR[p].thermo.T])
                else:
                    Yfinal += list(
                        list(self.PSR[p].thermo.Y))
        Yfinal = np.array(Yfinal)
        return Yfinal

    # noinspection PyUnreachableCode
    def solver(self, header, data, data_num, symm_axis):
        # print self.mfc_gen.coeffmat
        # print self.mfc_gen.fvect
        species = self.PSR[0].thermo.species_names
        print "Species=", len(species)
        print species
        print "Network size=", len(self.PSR)
        print "Exit massflow=", self.massimb.m_exit  # should be equal to the inlet massflow
        print "Massflow in outlet reactors", self.massimb.mass_out  # should be double the total inlet massflow
        print "Number of outlet mfc=", self.massimb.count
        MassImb = [0] * len(self.PSR)
        for i in range(len(self.PSR)):
            mfc_in = self.PSR[i].inlets
            mfc_out = self.PSR[i].outlets
            m_in = 0
            m_out = 0
            for j in mfc_in:
                if 'set_mass_flow_rate' in dir(j):
                    m_in += j.mdot(0)
            for j in mfc_out:
                if 'set_mass_flow_rate' in dir(j):
                    m_out += j.mdot(0)
            MassImb[i] += m_in - m_out
            self.network.append(ct.ReactorNet([self.PSR[i]]))
        print "Saving initial state"

        sr.save_data({1: self.mfc_gen.cells_react}, 'ReactorCells', self.path)
        dt = min(self.massimb.tres) / 2  # 0.001
        print "Tres min=", dt
        t0 = time.time()
        t1 = time.clock()

        sr.save(self.PSR, header, data, self.massimb.MassImb2, self.boundary.OutletReact, self.boundary.InletReact,
                self.massimb.tres, self.path, 0, data_num, symm_axis)
        """sr.save(self.PSR, header, data, MassImb, self.boundary.OutletReact, self.boundary.InletReact,
                self.massimb.tres, self.path, 1, data_num, symm_axis)"""

        tstart = 0.0
        changespec = []
        error = 1.0
        error_local = [1] * len(self.PSR)
        tstart_local = [0] * len(self.PSR)

        for i in range(len(species)):
            changespec.append(i)

        # Initialising constant mass flow terms of objective function
        self.massflow_const(species)
        # self.f_vect = self.mfc_gen.fvect  # [0.0] * (self.mfc_gen.vectsize * len(self.PSR))
        # self.massflow_const(species)
        print "Length BCs=", len(self.f_vect)
        print "Sum BCs=", sum(self.f_vect)
        # print self.f_vect

        fact = 1.0
        # Cantera Network
        """net = ct.ReactorNet(self.PSR)
        t0 = 0.0
        dt = 1e-03
        for i in range(100000):
            net.advance(t0+dt)
            t0 = t0+dt"""
        # net.advance_to_steady_state()
        print self.PSR[2].mass
        print self.PSR[2].volume
        print "Max rate cantera = ", max(abs(self.PSR[2].thermo.net_production_rates))

        Y0 = self.PSR_gas_ini()
        #self.ge = GoverningEquations(self.process, self.PSR, self.chem_mech, self.energy, self.mfc_gen, self)
        self.RHS_Jac0 = self.objective_func(t0, Y0)
        print "Initial objective func:", self.RHS_Jac0
        print "Obj func sum =", sum(self.RHS_Jac0)
        print "Max rate = ", max(abs(self.RHS_Jac0))
        epsi_count = 0
        save_count = 0
        t0 = 0.0
        dt = 1e-03
        epsilon2 =1.0
        Y0 = self.PSR_gas_ini()
        # Solver iterations
        for iter_g in range(30000):
            """try:
                sys.stdin.read()
            except KeyboardInterrupt:
                inp = raw_input("Enter val")
                if inp == "end":
                    break"""
            tol = 1e-15
            print "Iteration ", iter_g
            tmin = min(self.massimb.tres)
            tmax = max(self.massimb.tres)
            dt = tmax
            print "dt=", dt
            error_cum = 0
            dt = tmin
            if iter_g < 5 and (abs(self.epsilon_omega) > 1e-05 or abs(self.epsilon) > 1e-04):
                # Local solver
                Ynew, m_v_max, epsilon2 = self.localsolver_section(self.PSR, Y0, error_cum, tstart_local, self.dt_fac, dt)
                Y0 = Ynew
               # print "Initial objective func:", self.RHS_Jac0
                print "Epsilon = ", self.epsilon
                print "Epsilon omega=", self.epsilon_omega
                t0 = t0+dt

            # elif self.epsilon < 1e-02 and self.epsilon>1e-09:
            else:
                # Global Solver
                if self.solver_type[len(self.solver_type) - 1] == 0:
                    print "Initialising local solver params"
                    Y0 = self.PSR_gas_ini()
                    self.RHS_Jac0 = self.objective_func(t0, Y0)
                    # print "Initial objective func:", self.RHS_Jac0
                    print "Obj func sum =", sum(self.RHS_Jac0)
                    g_max = linalg.norm(self.RHS_Jac0) / self.RHS_Jac0.size
                    self.g_old = g_max
                for loop in range(1):
                    print "Global Solver"
                    # TODO: decide epsilon_omega cutoff value
                    # decide Newton/Time stepping based on rate of decrease of residuals. If rate high then Newtom
                    if (abs(self.epsilon<1.0) or (abs(self.g_omega) > 1.0)) and epsi_count<0 and iter_g>1:# and self.alpha > 1.e-14and epsi_count < 5 and abs(self.epsilon_omega<1)\
                            #and abs(self.epsilon<1) and abs(self.epsilon_g_omega<1)):
                        # Newton Raphson Global solver(Faster stepping)
                        if epsi_count < 5:
                            self.solver_type.append(1)
                            Yout, tnew, epsilon2 = self.globalsolver(t0, Y0, dt, species, self.g_old)
                            self.newton_count += 1
                            self.g_old = self.g_omega
                            print "Epsilon = ", self.epsilon
                            print "Epsilon omega=", self.epsilon_omega
                            if self.epsilon_omega >= 0.0 or self.epsilon_g_omega > 0.0:
                                epsi_count += 1
                                print "Epsi count=", epsi_count
                                self.alpha = max(self.alpha / 1.1, 1e-7)
                            #if self.epsilon < 0.0 or self.epsilon_g_omega < 0.0:
                             #   self.alpha = min(self.alpha * 1.3, 1.0)
                            print "Alpha=", self.alpha
                        else:
                            epsi_count = 0
                            self.alpha = max(self.alpha / 3, 1e-7)
                            if self.alpha<=1e-7:
                                self.alpha = self.alpha_default
                            print self.alpha

                    else:
                        # ODE time-stepping  (Slower but more stable and accurate)
                        epsi_count = 0
                        self.alpha = self.alpha_default
                        print "Tmin=", tmin
                        print "Time start=", t0
                        self.newton_count = 0
                        Yout, tnew, epsilon2 = self.globalsolver_odeint(t0, Y0, tmin, tmax, species)
                        self.g_old = self.g_omega

                    # self.epsilon_omega = (epsilon2 - self.epsilon) / (self.epsilon + tol)
                    # self.epsilon = epsilon2
                    t0 = tnew
                    Y0 = np.array(Yout)
                    #print "O2 amount=", Y0[(self.mfc_gen.vectsize + 1) * 2 + self.PSR[2].thermo.species_index("O2")]
                    print "Omega = ", self.g_omega
                    print "Epsilon=", self.epsilon
                    print "Epsilon omega = ", self.epsilon_omega
                    #print "Epsilon g omega = ", self.epsilon_g_omega

            # ToDO: error log issue in local solver
            if (iter_g == 0 or iter_g == 1 or iter_g == 2 or iter_g % 5 == 0):
                if save_count > 5:
                    save_count = 0
                sr.save(self.PSR, header, data, MassImb, self.boundary.OutletReact, self.boundary.InletReact,
                        self.massimb.tres, self.path, save_count, data_num, symm_axis)
                save_dict = {1: self.error_log, 2: self.error_spec_log, 3: self.solver_type, 4: self.error_rate_log}
                sr.save_data(save_dict, "ErrorLog" + str(save_count), self.path)
                sr.save_data({1: self.errorlog_co2}, "ErrorLogCO2" + str(save_count), self.path)
                sr.save_data({1: self.errorlog_no}, "ErrorLogNO" + str(save_count), self.path)
                sr.save_data({1: self.errorlog_no2}, "ErrorLogNO2" + str(save_count), self.path)
                save_count += 1

            # Residual analysis
            if abs(self.epsilon) < 1e-08 and abs(self.g_omega) < 1e-10:  # or abs(self.epsilon_omega) < 1e-12:
                print "Reaction rate=", self.g_omega
                print "tstart=", tstart
                print "error=", self.epsilon
                print "max tres=", max(self.massimb.tres)

                break

        print "Solution completed"
        t01 = time.time() - t0
        t11 = time.clock() - t1
        print "time to Solve0=", t01
        print "time to Solve1=", t11
        sr.save(self.PSR, header, data, MassImb, self.boundary.OutletReact, self.boundary.InletReact,
                self.massimb.tres, self.path, 10, data_num, symm_axis)
        save_dict = {1: self.error_log, 2: self.error_spec_log, 3: self.solver_type, 4: self.error_rate_log}
        print "Lenght rate error", len(self.error_rate_log)
        print "Length error", len(self.error_log)
        print "Length solver type", len(self.solver_type)
        sr.save_data(save_dict, "ErrorLog", self.path)
        sr.save_data({1: self.errorlog_co2}, "ErrorLogCO2", self.path)
        sr.save_data({1: self.errorlog_no}, "ErrorLogNO", self.path)
        sr.save_data({1: self.errorlog_no2}, "ErrorLogNO2", self.path)
        with open(self.path + '/errorlog.pkl', 'wb') as file:
            pickle.dump(self.error_log, file, pickle.HIGHEST_PROTOCOL)
        """with open(path+'/dtlog.pkl','wb') as file:
            pickle.dump(dt_log,file,pickle.HIGHEST_PROTOCOL)"""
        """with open(path + '/readme.txt', 'a') as file:
            file.write("Number of reactors=%s\n" % (len(PSR)))
            file.write("Energy=%s\n" % energy)
            file.write("Time to solve CRN=%f\n" % t11)
            file.write("Total inlet flow=%f\n" % mflow_tot)
            file.write("Total inlet air flow=%f\n" % mflow_tot_air)
            file.write("Total inlet fuel flow=%f\n" % mflow_tot_fuel)
            file.write("Total outlet flow=%f\n" % mflow_out)
            file.write("Total periodic flow=%f\n" % mflow_per)
            file.write("Inlet faces=%i\n" % inlet_faces)
            file.write("Outlet faces=%i\n" % outlet_faces)
            file.write("Periodic faces=%i\n" % per_face)
            file.write("Walls=%i\n" % walls)
            file.write("Inlet area total=%f\n" % f_area_tot)
            file.write("Total Interior faces=%i\n" % interior_faces)
            file.write("MFC total flow=%f\n" % mfc_flow)
            file.write("Total Volume=%f\n" % volume_total)
            file.write("Number of Reservoirs=%f\n" % res_num)
            file.write("Heat loss=%f\n" % heat)"""
