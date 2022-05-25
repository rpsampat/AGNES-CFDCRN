import SaveReactors as sr
import cantera as ct
import pickle
import time
from numpy import array, add, linalg, subtract, amax, absolute, zeros
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.integrate import BDF


def localsolver_section(energy, PSR, PSRInd_updated, network, tres, massflowreactor, error_cum,
                        error_local, tstart_local, Res_dict, dt_fac, mfc_gen, error_log, error, dt, solver_type):
    # Solving individual reactors
    if energy == 'on':
        local_length = len(PSR) * 2
    else:
        local_length = len(PSR)
    for i in range(1):  # local_length):
        error_bef = error
        error, traverse, PSRInd_updated = localsolver(PSRInd_updated, PSR, network, tres, massflowreactor,
                                                      error_cum, error_local, tstart_local, Res_dict,
                                                      dt_fac, energy)
        error_log.append(error)
        solver_type.append(0)
        print "error local=", error
        # Checking valve mass flows
        m_valve = 0
        m_v_max = 0.0

        for v in mfc_gen.ValveRecord:
            mfc_v = mfc_gen.ValveRecord[v]
            rel_massflow_valve = v.mdot(0) / mfc_v.mdot(0)
            if rel_massflow_valve > m_v_max:
                m_v_max = rel_massflow_valve
                m_valve = v
        print "Max valve flow=", m_v_max

        try:
            mfc_max = mfc_gen.ValveRecord[m_valve]
            reactor_max = mfc_gen.mfc_rec[mfc_max][0]
            reactor_max_to = mfc_gen.mfc_rec[mfc_max][1]
            print "Reactor max to=", reactor_max_to
            if (reactor_max in mfc_gen.mfc_per) or (reactor_max_to in mfc_gen.mfc_per):
                print "In periodic"
            print "Max reactor=", reactor_max
            print "Max reactor pressure=", PSR[reactor_max].thermo.P
            print "Max reactor to=", reactor_max_to
            print "Max reactor to pressure=", PSR[reactor_max_to].thermo.P
        except:
            print "Inlet Reactor Valve"

        print "dt_fac=", dt_fac
        if traverse != len(PSR):
            dt = dt / 10
        tmin = min(tres)
        tmin_dtfac = dt_fac * tmin
        sum_tres = sum(tres)
        # if error < 1.0 and error_bef < 1.0:
        # error_log.append(error)
        if error < 0.1 and tmin_dtfac < sum_tres:
            print "res min fac=", tmin_dtfac
            print "sum max=", sum_tres
            dt_fac += 1
        if error < 1e-09 and traverse == len(PSR) and error_bef < 1e-09:
            break

    if (error < 1e-06 and error_bef < 1e-06) or dt_fac >= 0.85 * len(PSR):
        # Correcting Valve Coefficients based on related MFCs
        for v in mfc_gen.ValveRecord:
            v.set_valve_coeff(0)

    return m_v_max, error


def localsolver(PSRInd_updated, PSR, network, tres, massflowreactor, error_cum, error_local, tstart_local,
                Res_dict, dt_fac, energy):
    """
    Local solver
    :param PSRInd_updated:
    :param PSR:
    :param network:
    :param tres:
    :param massflowreactor:
    :param error_cum:
    :param error_local:
    :param tstart_local:
    :param Res_dict:
    :param Reserv:
    :param dt_fac:
    :return:
    """

    traverse = 0
    e_max = 1e-20
    div_list = []

    for r_ind in PSRInd_updated:
        r = PSR[r_ind]
        psr_index = r_ind
        net = network[r_ind]
        try:
            xin = list(r.thermo.Y)
            if energy == 'on':
                xin = list(list(r.thermo.Y) + [r.thermo.T])
            dt_loc = tres[psr_index]
            net.advance(dt_loc * dt_fac)
            xfin = list(r.thermo.Y)
            if energy == 'on':
                xfin = list(list(r.thermo.Y) + [r.thermo.T])

            # updating reactor residence time
            tres[psr_index] = r.mass / massflowreactor[psr_index]
            func = [(xfin[n] - xin[n]) for n in range(len(xin))]
            e = []
            e_div = max(xin)

            for ind in range(len(func)):
                error_cum += (abs(func[ind]) / e_div) ** 2
                e.append(abs(func[ind]) / e_div)
            e_rel = max(e)
            error_local[psr_index] = e_rel

            if e_rel > e_max:
                e_max = e_rel
            tstart_local[psr_index] += dt_loc
            traverse += 1

            for num in Res_dict[r_ind]:
                g = r.thermo
                num[0].insert(g)

        except RuntimeError:
            # print marker
            """div_list.append(r_ind)
            res = Res_dict[r_ind][0]
            g = Reserv[res][0].thermo
            Reserv[res][1].insert(g)"""
            continue

    error = e_max

    # Residual analysis
    if traverse != len(PSR):
        divg = len(PSR) - traverse
        print "Diverging reactors=", divg
        error = 1.0
    print "Cummu error=", error_cum

    return error, traverse, PSRInd_updated


def objective_func(Y, PSR, species, coeffmat, massflowreactor, fvect, vectsize, wall_map):
    """
    Calculates the value of the function, net production rates of species in this case,
    for current state of all reactors. This is used as the RHS in the Newton's method.
    :param PSR:
    :param species:
    :param coeffmat:
    :param massflowreactor:
    :param fvect:
    :param wall_map:
    :return:
    """
    PSR, Yfinal = PSR_insert_gas(Y, PSR, vectsize, species)
    RatesofProduction = [0.0] * (vectsize * len(PSR))
    f_vect = [0.0] * (vectsize * len(PSR))
    Cw = [0.0] * (vectsize * len(PSR))
    J_diffu = [0.0] * (vectsize * len(PSR))

    # Obtaining rates of production
    for r_ind in range(len(PSR)):
        try:
            r = PSR[r_ind]
            rates = list(r.thermo.net_production_rates)  # kmol/m^3/s
            mw = list(r.thermo.molecular_weights)  # kg/kmol
            u_k = list(r.thermo.partial_molar_int_energies)  # J/kmol
            energy_production = 0

            for ind in range(vectsize):
                f_vect[r_ind * vectsize + ind] = fvect[r_ind * vectsize + ind]
                if ind < len(species):
                    # TODO: include turbulence diffusion effects
                    RatesofProduction[r_ind * vectsize + ind] = rates[ind] * mw[ind] * r.volume
                    energy_production -= rates[ind] * r.volume * u_k[ind]
                    for pt in range(len(PSR)):
                        # advective mass flows
                        Cw[r_ind * vectsize + ind] += -1 * coeffmat[r_ind][pt] * massflowreactor[pt] * \
                                                      PSR[pt].thermo.Y[ind]
                        # mass diffusion
                        try:
                            S_pq = wall_map[r_ind][pt].area
                            if S_pq<0.0:
                                print "Negative area!!"

                            dx = (PSR[r_ind].volume ** 0.333 + PSR[pt].volume ** 0.333) / 2.0
                            D_spec = PSR[r_ind].thermo.mix_diff_coeffs_mass[ind]
                            rho = PSR[r_ind].thermo.density_mass
                            delta_Y = PSR[r_ind].thermo.Y[ind] - PSR[pt].thermo.Y[ind]
                            J_diffu[r_ind * vectsize + ind] += -1.0 * rho * D_spec * S_pq * (delta_Y / dx)
                        except:
                            pass
                else:
                    # TODO: include thermal diffusion effects
                    RatesofProduction[r_ind * vectsize + ind] = energy_production
                    for pt in range(len(PSR)):
                        Cw[r_ind * vectsize + ind] += -1 * coeffmat[r_ind][pt] * massflowreactor[pt] * \
                                                      PSR[pt].thermo.cp_mass * PSR[pt].T
        except Warning:
            continue


    rHS_Jac = add(array(Cw), array(f_vect))
    rHS_Jac = add(rHS_Jac, array(J_diffu))
    RHS_Jac = add(array(rHS_Jac), array(RatesofProduction))

    return RHS_Jac


def source_analytic(gas, vol, mass, chem_mech, energy):
    """
    Analytic Jacobian of chemical source term
    :param gas:
    :param vol:
    :param mass:
    :param chem_mech:
    :return:
    """

    y = list(gas.concentrations)
    mw = list(gas.molecular_weights)  # kg/kmol
    mwu_k = list(gas.partial_molar_int_energies)  # J/kmol
    R = ct.Reaction.listFromFile(chem_mech)  # list of reactions in specified mechanism
    Jw = zeros((len(y), len(y)))
    if energy == 'on':
        Jw = zeros((len(y) + 1, len(y) + 1))
    forward = gas.forward_rates_of_progress  # kmol/m^3/s
    backward = gas.reverse_rates_of_progress  # kmol/m^3/s
    production_rate = gas.net_production_rates  # kmol/m^3/s

    test_rate_grad = 0
    test2_rate_grad = 0
    test3_rate_grad = 0

    for r in range(len(R)):
        # dictionary of reactant stoichiometric coeff
        react = R[r].reactants
        prod = R[r].products
        rate_grad = {}
        rate_grad_energy = 0
        try:
            E_a = R[r].rate.activation_energy  # J/kmol
            beta = R[r].rate.temperature_exponent
            Gas_constant = 8314.46261815324  # J/K/kmol
        except:
            continue

        # Calculating Reaction rate gradients
        for s in react.keys():
            ind = gas.species_index(s)
            try:
                rate_grad[ind] += react[s] * forward[r] / ((y[ind]) + 1e-15)
            except:
                rate_grad[ind] = react[s] * forward[r] / ((y[ind]) + 1e-15)
            if energy == 'on':
                rate_grad_energy += -1 * react[s] * (beta / gas.T + E_a / (Gas_constant * gas.T ** 2)) * forward[
                    r] * vol * mwu_k[ind]
                test2_rate_grad += forward[r] * vol * mwu_k[ind] / gas.T
                test3_rate_grad += -1 * react[s] * forward[r] * vol * mwu_k[ind] / gas.T

        for s in prod.keys():
            ind = gas.species_index(s)
            try:
                rate_grad[ind] += -1 * prod[s] * backward[r] / ((y[ind]) + 1e-15)
            except:
                rate_grad[ind] = -1 * prod[s] * backward[r] / ((y[ind]) + 1e-15)
            if energy == 'on':
                # TODO: include wall thermal diffusion
                rate_grad_energy -= prod[s] * (beta / gas.T + E_a / (Gas_constant * gas.T ** 2)) * backward[r] * vol * \
                                    mwu_k[ind]
                test2_rate_grad -= backward[r] * vol * mwu_k[ind] / gas.T
                test3_rate_grad -= prod[s] * backward[r] * vol * mwu_k[ind] / gas.T

        # Filling in the Jacobian
        # reactant stoichiometric coeff
        for s in react.keys():
            ind = gas.species_index(s)
            react_coeff = -1 * react[s]
            for spec in rate_grad:
                gamma = mass / (vol * mw[spec])
                Jw[ind][spec] += react_coeff * rate_grad[spec] * gamma

        # product stoichiometric coeff
        for s in prod.keys():
            ind = gas.species_index(s)
            prod_coeff = prod[s]
            for spec in rate_grad:
                gamma = mass / (vol * mw[spec])
                Jw[ind][spec] += prod_coeff * rate_grad[spec] * gamma

        if energy == 'on':
            Jw[len(y)][len(y)] += rate_grad_energy

    if energy == 'on':
        for spec in range(len(gas.species_names)):
            Jw[len(y)][len(y)] += vol * production_rate[spec] * mwu_k[spec] / gas.T
            test_rate_grad += vol * production_rate[spec] * mwu_k[spec] / gas.T

    return Jw


def Armijo(fdk, PSR, species, coeffmat, massflowreactor, fvect, fk, SpecMassFrac_update, SpecMassFrac):
    vectsize = len(species) + 1
    sigma = 1e-1
    alpha = 1
    beta = 0.1
    stop = 0
    alpha_arr = [alpha] * len(SpecMassFrac_update)
    rhs = fdk

    print SpecMassFrac
    rhs = rhs.dot(SpecMassFrac_update)

    i = 0
    for x in rhs:
        if x > 0:
            alpha_arr[i] = alpha_arr[i]
        i += 1
    while stop == 0:

        SpecMassFrac_update2 = add(SpecMassFrac, alpha * SpecMassFrac_update)

        for p in PSR:
            begin = PSR.index(p) * vectsize
            y0 = []

            for i in range(len(species)):
                y0.append(SpecMassFrac_update2[begin + i])
            g = p.thermo
            Temp = g.T

            try:
                Pr = SpecMassFrac_update2[begin + len(species)] + ct.one_atm
            except:
                Pr = g.P

            g.TPY = Temp, Pr, y0
            p.insert(g)
        fk_alpha, Spec, pr = function(0, y0, PSR, species, coeffmat, massflowreactor, fvect, vectsize)

        y_alpha = subtract(fk_alpha, fk)

        stop = 1

        rhs = rhs * sigma * alpha
        y_alpha_sum = linalg.norm(y_alpha)
        print "y_alpha=", y_alpha_sum
        rhs_sum = linalg.norm(rhs)
        print "rhs=", rhs_sum

        if y_alpha_sum > rhs_sum and alpha > 1e-7:
            alpha = alpha * beta
            print "alpha=", alpha
            stop = 0
        elif y_alpha_sum <= rhs_sum or alpha <= 1e-7:
            stop = 1

    return alpha


def Newton_sparse(g, J, x0, alpha):
    """
    Newton's method on sparse matrix
    :param g:
    :param J:
    :param x0:
    :param alpha:
    :return:
    """

    f = -1.0 * array(g)
    delta_x = spsolve(A=J, b=f, use_umfpack=True)
    exitcode = 0.0
    if exitcode != 0.0:
        print "Exit code=", exitcode
    xfinal = add(array(x0), alpha * delta_x)
    delta_x = alpha * delta_x

    return delta_x, xfinal, exitcode


def Jacobian_sparse(t, Y, PSR, coeffmat, species, vectsize, massflowreactor, wall_map, chem_mech, energy):
    """
    Sparse jacobian assembly with transport and source terms
    :param Y:
    :param PSR:
    :param coeffmat:
    :param species: PSR[0].thermo.species_names
    :param vectsize:
    :param massflowreactor:
    :param wall_map:
    :param chem_mech:
    :return:
    """

    Jac = []
    row_Jac = []
    col_Jac = []

    PSR, Y0 = PSR_insert_gas(Y, PSR, vectsize, species)
    print "Jacobian Calculation"

    for p in range(len(PSR)):
        gas = PSR[p]
        Jw = source_analytic(gas.thermo, gas.volume, gas.mass, chem_mech, energy)
        mw = list(gas.thermo.molecular_weights)

        for q in range(len(PSR)):
            coeff = coeffmat[p][q]
            if coeff != 0.0:

                for spec0 in range(vectsize):
                    for spec in range(vectsize):
                        Js = 0.0
                        Jw_add = 0.0

                        if spec == spec0:
                            if spec == len(species):
                                # energy equation
                                Js += -1.0 * coeff * massflowreactor[q] * PSR[q].thermo.cp_mass
                            else:
                                Js += -1.0 * coeff * massflowreactor[q]
                                try:
                                    S_pq = wall_map[p][q].area
                                    dx = (PSR[p].volume ** 0.333 + PSR[q].volume ** 0.333) / 2.0
                                    D_spec = PSR[p].thermo.mix_diff_coeffs_mass[spec0]
                                    rho = PSR[p].thermo.density_mass
                                    Js_diffu = -1.0 * (coeff / abs(coeff)) * rho * D_spec * S_pq / dx
                                    Js += Js_diffu
                                    # TODO: include turbulence diff
                                except:
                                    pass

                        if p == q:
                            if spec != len(species) and spec0 != len(species):
                                # Column in source term jacobian has constant denominator species
                                Jw_add += Jw[spec, spec0] * mw[spec] * PSR[p].volume

                            else:
                                # Energy equation source term
                                Jw_add -= Jw[spec, spec0]
                        J = Js + Jw_add
                        if J != 0.0:
                            Jac.append(J)
                            # Row in jacobian has constant numerator species
                            row = p * vectsize + spec
                            col = q * vectsize + spec0
                            row_Jac.append(row)
                            col_Jac.append(col)

    PSR_Jac = csc_matrix((Jac, (row_Jac, col_Jac)), shape=(len(PSR) * vectsize, len(PSR) * vectsize))
    print "Returning Jacobian"

    return PSR_Jac


def PSR_insert_gas(Y, PSR, vectsize, species):
    """
    Inserting new gas state in PSR object
    :param Y:
    :param PSR:
    :param vectsize:
    :param species:
    :return:
    """

    Yfinal = []
    for p in PSR:
        begin = PSR.index(p) * vectsize
        y0 = []
        for i in range(len(species)):
            conc = abs(Y[begin + i])
            y0.append(conc)

        g = p.thermo
        Pr = g.P
        if Pr < 1000:
            continue
        if vectsize > len(species) and 3000 > Y[begin + len(species)] > 200:
            Temp = Y[begin + len(species)]
        else:
            Temp = g.T

        g.TPY = Temp, Pr, y0
        p.insert(g)

        if vectsize > len(species):
            Yfinal += list(list(p.thermo.Y) + [p.thermo.T])
        else:
            Yfinal += list(p.thermo.Y)

    return PSR, Yfinal


def globalsolver(t0, Y0, dt, PSR, species, coeffmat, massflowreactor, fvect, vectsize, wall_map,
                 chem_mech, energy, iter, g_old):
    """
    Global solver that currently calls a global newton's method on a sparse system of equations
    :param t0:
    :param Y0:
    :param dt:
    :param PSR:
    :param species:
    :param coeffmat:
    :param massflowreactor:
    :param fvect:
    :param vectsize:
    :param wall_map: dictionary of wall objects mapped from reactor to neighbouring reactors.
    :param chem_mech:
    :return:
    """

    tol = 1e-15
    """if iter == 0:
        print "Calculating net production rates"
        RHS_Jac0 = objective_func(Y0, PSR, species, coeffmat, massflowreactor, fvect, vectsize, wall_map)
        g_old = linalg.norm(RHS_Jac0) / RHS_Jac0.size"""
    alpha = 0.01
    print "Newton Solver"
    PSR_Jac = Jacobian_sparse(t0, Y0, PSR, coeffmat, species, vectsize, massflowreactor, wall_map, chem_mech,
                              energy)
    RHS_Jac0 = objective_func(Y0, PSR, species, coeffmat, massflowreactor, fvect, vectsize, wall_map)
    Spec_Y = list(Y0)
    Yinitial = Spec_Y
    D_Y, Yfinal, exitcode = Newton_sparse(RHS_Jac0, PSR_Jac, Yinitial, alpha=alpha)
    dy = amax(absolute(D_Y))

    print "Newton dy=", dy
    PSR, Yfinal2 = PSR_insert_gas(Yfinal, PSR, vectsize, species)
    Y0 = list(Yfinal2)
    max_err = 0

    for sp in range(len(Yfinal2)):
        err = abs(Yfinal2[sp] - Yinitial[sp]) / abs(Yinitial[sp] + tol)
        if err > max_err:
            max_err = err

    print "Max error=", max_err

    psr_num = 0
    g_new = 0

    for r in range(len(RHS_Jac0)):
        if r >= (psr_num * vectsize + vectsize):
            psr_num += 1
        g_max = linalg.norm(RHS_Jac0) / RHS_Jac0.size  # abs(RHS_Jac0[r])#massflowreactor[psr_num])
        if g_max > g_new:
            g_new = g_max

    print "G_new=", g_new
    g_error = abs(g_new - g_old)
    print "G_rel=", g_error
    g_rel = g_new / g_old
    tstep = t0 + dt

    return PSR, Yfinal2, tstep, g_new, max_err


def globalsolver_odeint(t0, Y0, tmin, tmax, PSR, species, coeffmat, massflowreactor, fvect, vectsize, wall_map,
                        chem_mech,
                        energy, error_log, solver_type, error_rate_log):
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

    print "Calculating net production rates"
    Yinitial = list(Y0)

    r = BDF(lambda t, y: objective_func(y, PSR, species, coeffmat, massflowreactor, fvect, vectsize, wall_map),
            t0=t0, y0=Yinitial, t_bound=tmax, jac=lambda t, y: Jacobian_sparse(t, y, PSR, coeffmat, species, vectsize,
                                                                               massflowreactor, wall_map, chem_mech,
                                                                               energy))

    tol = 1e-15
    e_max = 0
    e_prev = -1
    e_rel = 0
    Yfinal = None

    while r.status == 'running' and r.t < tmax and abs(e_prev - e_max) > 1e-8:
        e_prev = e_rel
        xin = []
        xfin = []

        for p in PSR:
            if energy == 'on':
                xin += list(list(p.thermo.Y) + [p.thermo.T])
            else:
                xin += list(p.thermo.Y)
        r.step()
        print "Time=", r.t
        PSR, Yfinal = PSR_insert_gas(r.y, PSR, vectsize, species)
        for p in PSR:
            if energy == 'on':
                xfin += list(list(p.thermo.Y) + [p.thermo.T])
            else:
                xfin += list(p.thermo.Y)
        func = [(xfin[n] - xin[n]) for n in range(len(xin))]
        e = []
        for ind in range(len(func)):
            e.append(abs(func[ind]) / (xin[ind] + tol))

        e_rel = max(e)
        e_max = e_rel
        error_log.append(e_max)
        solver_type.append(2)
        epsilon = error_log(len(error_log) - 1)
        epsilon_omega = (e_max - epsilon) / (epsilon + tol)
        error_rate_log.append(epsilon_omega)
        """if e_rel > e_max:
            e_max = e_rel"""

        print "Error=", e_max

        if r.step_size < tmin * 1e-2:
            print "Too small time step"
            break

    print 'End global'

    if Yfinal:
        Yresult = Yfinal
    else:
        Yresult = Yinitial
    objfunc = r.fun(r.t, r.y)
    g_omega = linalg.norm(objfunc) / objfunc.size
    return PSR, Yresult, r.t, g_omega, e_max


def PSR_gas_conc(PSR, energy):
    Yfinal = []
    for p in PSR:
        if energy == 'on':
            Yfinal += list(list(p.thermo.Y) + [p.thermo.T])
        else:
            Yfinal += list(p.thermo.Y)

    return Yfinal


def DFS_VISIT(G, u, time, PSRIndex, cyclic, temp):
    time += 1
    G[u]['d'] = time
    G[u]['color'] = 'g'
    PSRIndex.append(u)
    temp.append(u)
    for v in G[u]['adjacent']:
        if G[v]['color'] == 'w':
            DFS_VISIT(G, v, time, PSRIndex, cyclic, temp)

    G[u]['color'] == 'b'
    time += 1
    G[u]['f'] = time


def DFS(G, inlet_list):
    """
    Depth First Search
    :param G:
    :param inlet_list:
    :return:
    """
    time = 0
    PSRIndex = []
    cyclic = []
    # starting DFS from inlet reactors
    for u in inlet_list:
        if G[u]['color'] == 'w':
            temp = []
            DFS_VISIT(G, u, time, PSRIndex, cyclic, temp)
    # proceeding to other reactors. As inlet reactor
    # marked 'b', the check for 'w' reactors would
    # make sure that inlet reactors are not revisited
    for u in G:
        if G[u]['color'] == 'w':
            temp = []
            DFS_VISIT(G, u, time, PSRIndex, cyclic, temp)
    return PSRIndex, cyclic


def solver(path, PSR, chem_mech, energy, mfc_gen, boundary, massimb, header, data, data_num, symm_axis):
    species = PSR[0].thermo.species_names
    # Creating DFS graph
    G = {}
    for g_ind in range(len(massimb.massflowreactor)):
        for mat_ind in range(len(massimb.massflowreactor)):
            val = abs(mfc_gen.coeffmat[mat_ind][g_ind])
            G_temp = {}
            if val != 0 and val != 1:
                G_temp[val] = mat_ind

        G_list = [0] * len(G_temp)
        d_time = []
        f_time = []
        i = 1
        g_len = len(G_temp)

        # arranging adjacent cells in descending order of mass flow exchange
        # so that highest mass flow neighbour is given priority
        for key in G_temp.keys():
            G_list[g_len - i] = G_temp[key]
            d_time.append(0)
            f_time.append(0)
            i = i + 1
        G[g_ind] = {'adjacent': list(G_list), 'color': 'w', 'd': d_time, 'f': f_time}

    PSRInd, cyclic = DFS(G, list(set(boundary.InletCondition)))
    print "Length DFS list=", len(PSRInd)
    print "Exit massflow=", massimb.m_exit  # should be equal to the inlet massflow
    print "Massflow in outlet reactors", massimb.mass_out  # should be double the total inlet massflow
    print "Number of outlet mfc=", massimb.count
    MassImb = [0] * len(PSR)
    network = []
    for i in range(len(PSR)):
        mfc_in = PSR[i].inlets
        mfc_out = PSR[i].outlets
        m_in = 0
        m_out = 0
        for j in mfc_in:
            if 'set_mass_flow_rate' in dir(j):
                m_in += j.mdot(0)
        for j in mfc_out:
            if 'set_mass_flow_rate' in dir(j):
                m_out += j.mdot(0)
        MassImb[i] += m_in - m_out
        network.append(ct.ReactorNet([PSR[i]]))

    # print "Number of Reservoirs=", res_num
    print "Saving initial state"

    sr.save_data({1: mfc_gen.cells_react}, 'ReactorCells', path)
    print "Solving network"
    print "Advancing"
    dt = min(massimb.tres) / 2  # 0.001
    print "Tres min=", dt
    t0 = time.time()
    t1 = time.clock()
    error_log = []  # logging errors
    error_rate_log = []  # logging relative rate of change of error
    solver_type = []  # 0: local solver, 1: Newton global solver, 2: Time-stepping global solver

    errorlog_co2 = []
    errorlog_no = []
    errorlog_no2 = []

    sr.save(PSR, header, data, massimb.MassImb2, boundary.OutletReact, boundary.InletReact,
            massimb.tres, path, 0, data_num, symm_axis)
    sr.save(PSR, header, data, MassImb, boundary.OutletReact, boundary.InletReact,
            massimb.tres, path, 1, data_num, symm_axis)

    tstart = 0
    changespec = []
    error = 1
    error_local = [1] * len(PSR)
    tstart_local = [0] * len(PSR)

    for i in range(len(species)):
        changespec.append(i)

    PSRInd_updated = list(PSRInd)

    dt_fac = 10  # int(tmax/tmin)
    g_old = 1.0
    epsilon = 1.0  # error
    epsilon_omega = 0.0  # rate of error
    g_omega = 1.0  # total rate of change of conserved property
    newton_count = 0
    # Solver iterations
    for iter_g in range(10):
        tol = 1e-15
        print "Iteration ", iter_g
        tmin = min(massimb.tres)
        tmax = max(massimb.tres)
        dt = tmax
        print "dt=", dt
        error_cum = 0
        dt = tmin
        if epsilon == 0:  # > 1e-02:
            # Local solver
            m_v_max, epsilon2 = localsolver_section(energy, PSR, PSRInd_updated, network, massimb.tres,
                                                    massimb.massflowreactor, error_cum,
                                                    error_local, tstart_local, mfc_gen.Res_dict, dt_fac, mfc_gen,
                                                    error_log, error, dt, solver_type)
            epsilon_omega = (epsilon2 - epsilon) / (epsilon + tol)
            error_rate_log.append(epsilon_omega)
            epsilon = epsilon2

        # elif epsilon < 1e-02:
        else:
            # Global Solver
            Y0 = PSR_gas_conc(PSR, energy)
            for loop in range(1):
                t0 = tstart
                print "Global Solver"
                # TODO: decide epsilon_omega cutoff value
                # decide Newton/Time stepping based on rate of decrease of residuals. If rate high then Newtom
                if epsilon_omega < 0 or newton_count < 4:
                    # PSR, Yfinal2, tstep, g_new, max_err
                    solver_type.append(1)
                    PSR, Yout, tnew, g_omega, epsilon2 = globalsolver(t0, Y0, dt, PSR, species, mfc_gen.coeffmat,
                                                                      massimb.massflowreactor, mfc_gen.fvect,
                                                                      mfc_gen.vectsize, mfc_gen.wall_map, chem_mech,
                                                                      energy, loop, g_old)
                    error_log.append(epsilon2)
                    epsilon_omega = (epsilon2 - epsilon) / (epsilon + tol)
                    error_rate_log.append(epsilon_omega)
                    newton_count += 1

                else:
                    # PSR, Yresult, r.t, e_max
                    print "Tmin=", tmin
                    newton_count = 0
                    PSR, Yout, tnew, g_omega, epsilon2 = globalsolver_odeint(t0, Y0, tmin, tmin, PSR, species,
                                                                             mfc_gen.coeffmat,
                                                                             massimb.massflowreactor, mfc_gen.fvect,
                                                                             mfc_gen.vectsize, mfc_gen.wall_map,
                                                                             chem_mech,
                                                                             energy, error_log, solver_type,
                                                                             error_rate_log)
                g_old = g_omega
                epsilon_omega = (epsilon2 - epsilon) / (epsilon + tol)
                epsilon = epsilon2
                tstart = tnew
                Y0 = Yout
                # Error logging
                if epsilon2 < 1e-09 and g_omega < 1e-07:
                    print "Epsilon2 = ", epsilon2
                    print "G_omega = ", g_omega
                    break

        if iter_g == 0 or iter_g == 5 or iter_g == 10 or iter_g == 15 or iter_g == 20 or iter_g == 25 or \
                iter_g == 30 or iter_g == 35 or iter_g == 40 or iter_g == 45 or iter_g == 50 or iter_g == 5000:
            sr.save(PSR, header, data, MassImb, boundary.OutletReact, boundary.InletReact,
                    massimb.tres, path, iter_g, data_num, symm_axis)

        # Residual analysis
        if epsilon < 1e-09 or g_omega < 1e-15:
            Y0 = PSR_gas_conc(PSR, energy)
            """RHS_Jac1 = objective_func(t0, Y0, PSR, species, mfc_gen.coeffmat,
                                massimb.massflowreactor, mfc_gen.fvect, mfc_gen.vectsize, mfc_gen.wall_map)
            g_new = amax(absolute(RHS_Jac1))"""
            print "Reaction rate=", g_omega
            print "tstart=", tstart
            print "error=", error
            print "max tres=", max(massimb.tres)

            break

    print "Solution completed"
    t01 = time.time() - t0
    t11 = time.clock() - t1
    print "time to Solve0=", t01
    print "time to Solve1=", t11
    sr.save(PSR, header, data, MassImb, boundary.OutletReact, boundary.InletReact,
            massimb.tres, path, 10, data_num, symm_axis)
    save_dict = {1: error_log, 2: solver_type, 3: error_rate_log}
    print "Lenght rate error", len(error_rate_log)
    print "Length error", len(error_log)
    print "Length solver type", len(solver_type)
    sr.save_data(save_dict, "ErrorLog", path)
    sr.save_data({1: errorlog_co2}, "ErrorLogCO2", path)
    sr.save_data({1: errorlog_no}, "ErrorLogNO", path)
    sr.save_data({1: errorlog_no2}, "ErrorLogNO2", path)
    with open(path + '/errorlog.pkl', 'wb') as file:
        pickle.dump(error_log, file, pickle.HIGHEST_PROTOCOL)
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
