import cantera as ct
import MassImbalance as mi
import pickle
from Generate_MFC import GenerateMFC
import os
import numpy as np
from Solver import CRN_Solver
# from CRN_Solver_OOPS import CRN_Solver
#from Solver_test import CRN_Solver
# import CRN_Solver
import CFD_Data_Handle as cfh


def gas(cell_id, data, header, key, gas_obj, data_num):
    """
    Assigns state, T+D+Y, to gas object
    :param cell_id:
    :param data:
    :param header:
    :param key:
    :param gas_obj:
    :param data_num:
    :return:
    """
    try:
        T = data[data_num][header.index('Static Temperature') + 1][cell_id - 1]
        T_corr = T
        D = data[data_num][header.index('Density') + 1][cell_id - 1]
        P = ct.one_atm
        Y = {j.upper(): data[data_num][header.index('Mass fraction of ' + j) + 1][cell_id - 1] for j in key}
        try:
            gas_obj.TDY = T_corr, D, Y
        except:
            gas_obj.TDY = 300, D, Y
    except:
        gas_obj.TPY = 1500, ct.one_atm, {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}
    return gas_obj


class ChemicalReactorNetwork:
    """
    This is a class that defines a CRN. The CRN is generated from data that
    is passed by creating 0D PSRs and interconnecting them with mass flow
    controllers and valves. Certain other operations such as mass balance
    are also performed to ensure the well-posedness of the problem.
    """

    def __init__(self, chem_mech, energy, process):
        self.PSR = []
        self.key = []
        # Defining gas object and species
        self.gas_obj = []
        self.species = []
        self.chem_mech = chem_mech
        self.chemical_mechanism(chem_mech)
        # Defining Energy equation solver(off by default)
        self.energy = ''
        self.energy_equation(energy)
        self.res_num = 0
        self.mfc_flow = 0.0
        self.interior_faces = 0
        self.volume_total = 0.0
        self.process = process

    def print_prop(self):
        print "Number of Reactors=", len(self.PSR)
        print "Total Volume=", self.volume_total

    def chemical_mechanism(self, chem_mech):
        """
        Used to update gas objects and related parameters when
        a different chemical mechanism needs to be specified
        after CRN generation is completed.
        :param chem_mech: .cti format mechanism
        :return:
        """
        self.gas_obj = ct.Solution(chem_mech)
        self.gas_obj.transport_model = 'Multi'
        self.species = self.gas_obj.species_names

    def energy_equation(self, energy):
        """
        Define energy equation solver(on/off)
        :param energy: 'on'/'off'
        :return:
        """
        self.energy = energy

    def generate_key(self, data, data_num, header):
        """
        Generates a list of species present in the input file that are also
        present in the specified reaction mechanism.

        :param data: input data dictionary
        :param data_num:
        :param header:
        :return:
        """
        for i in data[data_num].keys():
            s1 = header[i - 1].find('Mole fraction of')
            s2 = header[i - 1].find('<')
            s3 = header[i - 1].find('>')
            name = header[i - 1][len('Mass fraction of '):]
            if s1 != -1 and s2 == -1 and s3 == -1 and (name.upper() in self.species):
                self.key.append(name)
            elif s1 != -1 and s2 != -1 and s3 != -1 and (name[0].upper() in self.species):
                self.key.append(name[0])

    def reactor_create(self, graph, data, data_num, zone, header):
        """
        Creates cantera Ideal Gas Reactor objects and stores them in a list.

        :param graph: CRN connectivity graph
        :param data:
        :param data_num:
        :param zone:
        :param header:
        :return:
        """
        # creating ideal reactors
        print "Creating Reactors"
        volume_total = 0
        explored = {}
        print data[data_num].keys()
        print header
        self.PSR = np.empty(len(graph.keys()), dtype=object)
        if self.process == "PIVCRN":
            self.gas_obj.TPY = 2200, ct.one_atm, {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}
            #self.gas_obj.TPY = 200, ct.one_atm, {'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}
            psr = ct.IdealGasReactor(self.gas_obj, energy=self.energy)
            react_net = ct.ReactorNet([psr])
            react_net.advance_to_steady_state()
            """print dir(self.gas_obj)
            print psr.thermo.partial_molar_cp * 2200
            print psr.thermo.partial_molar_enthalpies-psr.thermo.partial_molar_cp * 2200"""
        for cell in graph.keys():
            explored[cell] = []
            # Modifying existing gas object state
            if self.process == "PIVCRN":
                psr = ct.IdealGasReactor(self.gas_obj, energy=self.energy)
                #vol = 1e-8*2e-3
            else:
                g = gas(cell, data, header, self.key, self.gas_obj, data_num)
                # Inserting gas object in Ideal Gas Reactor object
                psr = ct.IdealGasReactor(g, energy=self.energy)
            vol = data[data_num][header.index('Cell Volume') + 1][cell - zone[data_num][1]]
            # density = data[data_num][header.index('Density') + 1][cell - zone[data_num][1]]
            # define reactor volume
            # if vol< 0:
            # print vol
            # print cell
            psr.volume = vol
            """if density > 1:
                print density
                print cell"""
            # Todo: remove : self.PSR.append(psr)
            self.PSR[cell - 1] = psr
            # print "Pressure = ", psr.thermo.P
            self.volume_total += vol

    def crn_gen(self, path, graph, facearea, data, data_num, header, data_std, periodic, zone, bc_dict, symm_axis,
                Peclet):
        """

        :param path:
        :param graph:
        :param facearea:
        :param data:
        :param data_num:
        :param header:
        :param data_std:
        :param periodic:
        :param zone:
        :param bc_dict:
        :param symm_axis:
        :param Peclet:
        :return:
        """
        self.generate_key(data, data_num, header)
        self.reactor_create(graph, data, data_num, zone, header)  # Creating individual reactors of the network
        # Creating Mass flow controllers between reactors in network
        mfc_gen = GenerateMFC(self.PSR, self.energy, periodic, bc_dict,self.process)
        print zone
        # print data_std[20]
        print header
        # print facearea[10]

        boundary = mfc_gen.generate_mfc(graph, facearea, data, data_num, header, data_std, self.PSR, zone, Peclet)
        """for v in mfc_gen.ValveRecord:
            print "Valve flow=", v.mdot(0)"""
        # Balancing Mass Flows

        massimb = mi.MassImbalance(self.process)
        if self.process == "PIVCRN":
            massimb.PIVCRN_velocity(mfc_gen, self.PSR, path, boundary)
        else:
            massimb.PIVCRN_velocity(mfc_gen, self.PSR, path, boundary)
            #massimb.mass_imbalance(mfc_gen, boundary, self.PSR, path)
        # print "Outflow reator flowrate = ",massimb.massflowreactor[112]
        # print mfc_gen.coeffmat[112]
        """for v in mfc_gen.ValveRecord:
            print "Valve flow=", v.mdot(0)"""
        # print mfc_gen.Res_dict[0]
        # Printing properties
        # self.print_prop()
        mfc_gen.print_prop()
        boundary.print_prop()
        massimb.print_prop()
        # Solver
        """CRN_Solver.solver(path, self.PSR, self.chem_mech, self.energy, mfc_gen, boundary,
              massimb, header, data, data_num, symm_axis)"""
        crn_soln = CRN_Solver(path, self.PSR, self.chem_mech, self.energy, mfc_gen, boundary,
                              massimb, self.process)
        crn_soln.solver(header, data, data_num, symm_axis)

    def main(self, path, data_num, id_dict, symm_axis, Peclet):
        print "start"
        # CRN data
        with open(path + "/graph2plot.pkl", 'rb') as f:
            graph = pickle.load(f)
        with open(path + "/graphdata.pkl", 'rb') as f:
            data = pickle.load(f)
        with open("header.pkl", 'rb') as f:
            header = pickle.load(f)
        with open("zone.pkl", 'rb') as f:
            zone = pickle.load(f)
        print "Reading Face-Cell map"
        with open("facearea.pkl", 'rb') as f:
            facearea = pickle.load(f)
        print "Reading CFD Data"
        with open("data_std.pkl", 'rb') as f:
            data_std = pickle.load(f)
        print "Reading Periodic Faces"
        with open("pfaces.pkl", 'rb') as f:
            periodic = pickle.load(f)
        self.crn_gen(path, graph, facearea, data, data_num, header,
                     data_std, periodic, zone, id_dict, symm_axis, Peclet)

    def main_PIV(self, path, data_num, id_dict, symm_axis, Peclet):
        print "start"
        # CRN data
        with open(path + "/graph2plot_PIV.pkl", 'rb') as f:
            graph = pickle.load(f)
        with open(path + "/graphdata_PIV.pkl", 'rb') as f:
            data = pickle.load(f)
        with open("header_PIV.pkl", 'rb') as f:
            header = pickle.load(f)
        with open("zone_PIV.pkl", 'rb') as f:
            zone = pickle.load(f)
        print "Reading Face-Cell map"
        with open("facearea_PIV.pkl", 'rb') as f:
            facearea = pickle.load(f)
        print "Reading CFD Data"
        with open("data_std_PIV.pkl", 'rb') as f:
            data_std = pickle.load(f)
        periodic = []
        self.crn_gen(path, graph, facearea, data, data_num, header,
                     data_std, periodic, zone, id_dict, symm_axis, Peclet)


if __name__ == "__main__":
    process0 = "CFDCRN"  # PIVCRN//CFDCRN
    path = os.getcwd()
    #path = path + '/500_492PIVCRN_Y Velocity_X-Coordinate__'
    path = path + '/1000_968_Static Temperature_Mass fraction of o2__'
    #path = path + '/110000_100035PIVCRN_Y Velocity__'
    #path = path + '/1000000_120951PIVCRN_Y Velocity__'
    path = path.replace(os.sep, '/')
    symm_axis = 'X-Coordinate'
    chem_mech = 'gri30.cti' #'drm19.cti'/'gri30.cti'
    crn = ChemicalReactorNetwork(chem_mech=chem_mech, energy='off', process=process0)
    phi = 0.91
    xch4 = phi/((2*4.76)+phi)
    gas_obj = ct.Solution(chem_mech)
    gas_obj.TPY = 1500, ct.one_atm, {'CH4': xch4, 'O2': (1 - xch4) * 0.23, 'N2': (1 - xch4) * 0.77}
    # self.gas_obj.TPY = 200, ct.one_atm, {'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}
    psr = ct.IdealGasReactor(gas_obj, energy='on')
    react_net = ct.ReactorNet([psr])
    react_net.advance_to_steady_state()
    if process0 == "PIVCRN":
        data_num = 1
        id_dict = {3: ['air', {'CH4': xch4, 'O2': (1 - xch4) * 0.23, 'N2': (1 - xch4) * 0.77}, 273.15+20.0],
                   6: ['secondary_inlet', psr.thermo.mass_fraction_dict(), psr.thermo.T],
                   17: ['fuel', {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}, 294.0],
                   20: ['pilot', {'O2': 0.059, 'N2': 0.735, 'H2O': 0.091, 'CO2': 0.111, 'OH': 0.001, 'NO': 0.003},
                        1880.0],
                   23: ['top'],
                   24: ['side']}
        crn.main_PIV(path, data_num, id_dict, symm_axis, Peclet=False)
    else:
        data_num = 15
        id_dict = {21: ['air', {'O2': 0.23, 'N2': 0.77}, 291.0],
                   17: ['fuel', {'CH4': 0.16, 'O2': (1 - 0.16) * 0.23, 'N2': (1 - 0.16) * 0.77}, 294.0],
                   20: ['pilot', {'O2': 0.059, 'N2': 0.735, 'H2O': 0.091, 'CO2': 0.111, 'OH': 0.001, 'NO': 0.003},
                        1880.0],
                   23: ['top'],
                   24: ['side']}
        crn.main(path, data_num, id_dict, symm_axis, Peclet=True)
