import cantera as ct
import math


class BoundaryConditions:
    def __init__(self, bc_def_dict, rhsvect, gas_obj):
        # Boundary definitions
        self.bc_def = bc_def_dict
        # Mass Balance vectors
        self.rhsvect = rhsvect
        # Gas object
        self.gas_obj = gas_obj
        # inlets
        self.Inletmfc = {}
        self.InletCondition = []
        self.InletReact = []
        self.AirInletmfc = {}
        self.SecondaryInletmfc = {}
        self.PilotInletmfc = {}
        self.FuelInletmfc = {}
        self.valveAirInlet = {}
        self.valveAirInlet_rev = {}
        self.Exhaust = []
        self.mflow_tot = 0.0
        self.mflow_tot_fuel = 0.0
        self.mflow_tot_air = 0.0
        self.mflow_tot_sec_inlet= 0.0
        self.mflow_tot_pilot = 0.0
        self.inlet_faces = 0
        # outlets
        self.Outletmfc = {}
        self.OutletReact = []
        self.outlet_faces = 0
        self.mflow_per = 0.0
        self.per_face = 0.0
        self.mflow_shadow = 0.0
        self.mflow_out = 0.0
        f_area_tot = 0.0

        self.walls = 0.0
        self.heat = 0.0

        self.periodic_id = {}
        self.shadow_id = {}
        self.m_shad_direct = 0.0
        mflow_inter = 0.0

    def print_prop(self):
        print "Total inlet flow=", self.mflow_tot
        print "Total inlet air flow=", self.mflow_tot_air
        print "Total inlet fuel flow=", self.mflow_tot_fuel
        print "Total outlet flow=", self.mflow_out
        print "Total periodic flow=", self.mflow_per
        print "Total shadow flow=", self.mflow_shadow
        print "Total shadow flow direct=", self.m_shad_direct
        # print "Total interface mass flow=", self.mflow_inter
        print "Inlet faces=", self.inlet_faces
        print "Outlet faces=", self.outlet_faces
        print "Periodic faces=", self.per_face
        print "Walls=", self.walls
        # print "Inlet area total=", self.f_area_tot
        # print "Total Interior faces=", self.interior_faces
        print "Heat loss=", self.heat

    def inlet(self, mfc_gen, mflux, reactor, reactor_obj, current_key, zone_id):
        # 20: massflow inlet, 10: velocity inlet
        self.inlet_faces += 1
        """f_area = facearea[face][0]
        f_area_tot += abs(f_area)"""
        mflux = abs(mflux)
        if mflux == 0 or math.isnan(mflux):
            print mflux
            print "Inlet issues"
        self.mflow_tot += mflux

        self.InletCondition.append(reactor - 1)

        if current_key == 'air':
            self.mflow_tot_air += mflux
        elif current_key == 'secondary_inlet':
            self.mflow_tot_sec_inlet += mflux
        elif current_key == 'pilot':
            self.mflow_tot_pilot += mflux
        elif current_key == 'fuel' or 'fuel' in current_key:
            self.mflow_tot_fuel += mflux
        else:  # TODO: solve this in a better way
            raise ValueError(str(current_key) + ' is not specified in the solver as an inlet '
                                                'boundary condition. Please add it.')
        self.gas_obj.TPY = self.bc_def[zone_id][2], ct.one_atm, self.bc_def[zone_id][1]
        density = self.gas_obj.density
        res = ct.Reservoir(self.gas_obj)
        if mfc_gen.process == "PIVCRN":
            for molecule in self.bc_def[zone_id][1]:
                mfc_gen.fvect[(reactor - 1) * mfc_gen.vectsize + self.gas_obj.species_index(molecule)] \
                    += density * mflux * self.bc_def[zone_id][1][molecule]
            if mfc_gen.energy == 'on':
                mfc_gen.fvect[(reactor - 1) * mfc_gen.vectsize + len(self.gas_obj.species_names)] \
                    += density * mflux * self.gas_obj.cp_mass * self.gas_obj.T
            mfc_gen.rhsvect[reactor - 1] += density * abs(mflux)
        else:
            mfc_gen.rhsvect[reactor - 1] += abs(mflux)
            for molecule in self.bc_def[zone_id][1]:
                mfc_gen.fvect[(reactor - 1) * mfc_gen.vectsize + self.gas_obj.species_index(molecule)] \
                    += mflux * self.bc_def[zone_id][1][molecule]
            if mfc_gen.energy == 'on':
                mfc_gen.fvect[(reactor - 1) * mfc_gen.vectsize + len(self.gas_obj.species_names)] \
                    += mflux * self.gas_obj.cp_mass * self.gas_obj.T

        # Defining Massflow Controller.
        # try: add mass flow to existing mfc
        # except: define new mfc
        try:
            if current_key == 'air':
                mfc = self.AirInletmfc[reactor - 1]
            elif current_key == 'secondary_inlet':
                mfc = self.SecondaryInletmfc[reactor - 1]
            elif current_key == 'pilot':
                mfc = self.PilotInletmfc[reactor - 1]
            elif current_key == 'fuel' or 'fuel' in current_key:
                mfc = self.FuelInletmfc[reactor - 1]
            else:  # TODO: solve this in a better way
                raise ValueError(
                    str(current_key) + ' is not specified in the solver as an inlet '
                                       'boundary condition. Please add it.')
            m = mfc.mdot(0)
            m += mflux
            mfc.set_mass_flow_rate(m)

        except:
            mfc = ct.MassFlowController(res, reactor_obj, mdot=mflux)
            if current_key == 'air':
                self.AirInletmfc[reactor - 1] = mfc
            elif current_key == 'secondary_inlet':
                self.SecondaryInletmfc[reactor - 1] = mfc
            elif current_key == 'pilot':
                self.PilotInletmfc[reactor - 1] = mfc
            elif current_key == 'fuel' or 'fuel' in current_key:
                self.FuelInletmfc[reactor - 1] = mfc
            else:  # TODO: solve this in a better way
                raise ValueError(
                    str(current_key) + ' is not specified in the solver as an inlet '
                                       'boundary condition. Please add it.')
            kv = mfc_gen.valve_mult * mflux

            # Feedforward valve for inlet MFC
            if mfc_gen.ffin == 'on':
                v0 = ct.Valve(res, reactor_obj)
                self.valveAirInlet[reactor - 1] = v0
                v0.set_valve_coeff(kv)
                mfc_gen.ValveRecord[v0] = mfc

            # Feedback valve
            if mfc_gen.fbin == 'on':
                v02 = ct.Valve(reactor_obj, res)
                self.valveAirInlet_rev[reactor - 1] = v02
                v02.set_valve_coeff(kv)
                mfc_gen.ValveRecord[v02] = mfc

        self.InletReact.append(reactor - 1)


    def outlet(self, mfc_gen, mflux, reactor, reactor_obj):
        mflow = abs(mflux)
        # total outflow from reactor
        if mfc_gen.mattype == "sparse":
            mfc_gen.mass_bal(reactor - 1, reactor - 1, mflow)
        else:
            mfc_gen.coeffmat[reactor - 1][reactor - 1] += mflow
        g = reactor_obj.thermo
        exhaust = ct.Reservoir(g)
        try:
            mfc_gen.Res_dict[reactor - 1].append([exhaust, reactor_obj])
        except:
            mfc_gen.Res_dict[reactor - 1] = [[exhaust, reactor_obj]]

        try:
            mfc = self.Outletmfc[reactor - 1]
            m = mfc.mdot(0)
            m += mflow
            mfc.set_mass_flow_rate(m)
        except:
            mfc = ct.MassFlowController(reactor_obj, exhaust, mdot=mflow)

            # Feedforward valve for outlet MFC
            kv = mfc_gen.valve_mult * mflow
            if mfc_gen.ffout == 'on':
                v0 = ct.Valve(reactor_obj, exhaust)
                v0.set_valve_coeff(kv)
                mfc_gen.ValveRecord[v0] = mfc

            # Feedback valve for outlet MFC
            if mfc_gen.fbout == 'on':
                v02 = ct.Valve(exhaust, reactor_obj)
                v02.set_valve_coeff(kv)
                mfc_gen.ValveRecord[v02] = mfc
            self.Exhaust.append(exhaust)
            mfc_gen.mfc_rec[mfc] = [reactor - 1]
            self.Outletmfc[reactor - 1] = mfc
            self.OutletReact.append(reactor - 1)
        self.mflow_out += mflow
        self.outlet_faces += 1
        """reactor_out[reactor - 1] += mflow        
        """

    def periodic(self, mfc_gen, data_std, zone, PSR):
        # periodic boundary condition
        for face in self.periodic_id:
            cell = self.periodic_id[face][0]  # periodic reactor
            z = self.periodic_id[face][1]
            pair = mfc_gen.periodic[face]  # second face in pair, i.e flow towards it
            cell_check = self.shadow_id[pair][0]  # shadow reactor
            mflow = data_std[z][18][face - zone[z][1]]
            self.mflow_per += mflow
            z_shad = self.shadow_id[pair][1]
            self.mflow_shadow += data_std[z_shad][18][pair - zone[z_shad][1]]
            if cell == cell_check:
                self.Periodic_exclud.append(cell - 1)

            if mflow > 0 and cell != cell_check:
                try:
                    mfc_predef = mfc_gen.mfc_rec_inv[cell - 1][cell_check - 1]
                    m = mfc_predef[0].mdot(0) + mflow
                    mfc_predef[0].set_mass_flow_rate(m)
                    mfc_predef[1].set_mass_flow_rate(m)
                    mfc_gen.coeffmat[cell - 1][cell - 1] += mflow  # total outflow from cell
                    mfc_gen.coeffmat[cell_check - 1][cell - 1] += -1.0 * mflow  # inflow from cell to cell_check
                except:
                    mfc = mfc_gen.mfc_def(self, PSR[cell - 1], Res_temp, mflow, mfc_gen.ff, mfc_gen.fb)
                    g = PSR[cell - 1].thermo
                    Res_temp = ct.Reservoir(g)
                    try:
                        mfc_gen.Res_dict[cell - 1].append([Res_temp, PSR[cell - 1]])
                    except:
                        mfc_gen.Res_dict[cell - 1] = [[Res_temp, PSR[cell - 1]]]
                    mfc2 = mfc_gen.mfc_def(self, Res_temp, PSR[cell_check - 1], mflow, mfc_gen.ff, mfc_gen.fb)
                    mfc_gen.flow_record_per(self, cell - 1, cell_check - 1, mflow, mfc, mfc2)

            elif mflow < 0 and cell != cell_check:
                mflow = abs(mflow)

                try:
                    mfc_predef = mfc_gen.mfc_rec_inv[cell_check - 1][cell - 1]
                    m = mfc_predef[0].mdot(0) + mflow
                    mfc_predef[0].set_mass_flow_rate(m)
                    mfc_predef[1].set_mass_flow_rate(m)
                    mfc_gen.coeffmat[cell_check - 1][cell_check - 1] += mflow  # total outflow from cell_check
                    mfc_gen.coeffmat[cell - 1][cell_check - 1] += -1.0 * mflow  # inflow from cell_check to cell
                except:
                    mfc = mfc_gen.mfc_def(self, PSR[cell - 1], Res_temp, mflow, mfc_gen.ff, mfc_gen.fb)
                    g = PSR[cell - 1].thermo
                    Res_temp = ct.Reservoir(g)
                    try:
                        mfc_gen.Res_dict[cell - 1].append([Res_temp, PSR[cell - 1]])
                    except:
                        mfc_gen.Res_dict[cell - 1] = [[Res_temp, PSR[cell - 1]]]
                    mfc2 = mfc_gen.mfc_def(self, Res_temp, PSR[cell_check - 1], mflow, mfc_gen.ff, mfc_gen.fb)
                    mfc_gen.flow_record_per(self, cell_check - 1, cell - 1, mflow, mfc2, mfc)

    def resolve_bcs(self, mfc_gen, face, zone_type, zone_id, reactor, reactor_obj, mflux, heatflux):
        # Boundary conditions
        # mass flow inlet
        # requires mass flow controller and reservoir
        # print zone[z]
        if zone_type == 20 or zone_type == 1020 or zone_type == 10 or zone_type == 1010:
            current_key = self.bc_def[zone_id][0]
            if mflux != 0:
                self.inlet(mfc_gen, mflux, reactor, reactor_obj, current_key, zone_id)

        # periodic boundary
        if zone_type == 12 or zone_type == 1012:
            self.periodic_id[face] = [reactor, zone_id]
            self.per_face += 1  # counting periodic faces

        # shadow periodic boundary
        if zone_type == 8 or zone_type == 1008:
            self.shadow_id[face] = [reactor, zone_id]
            self.m_shad_direct += mflux

        # interface

        # wall
        if zone_type == 3 or zone_type == 1003:
            # print "Wall found"
            # Res_env
            self.walls += 1

            if mfc_gen.energy == 'on':
                Res_env = ct.Reservoir(reactor_obj.thermo)
                # Res_env.thermo.T = 300
                k = reactor_obj.thermo.thermal_conductivity
                w = ct.Wall(reactor_obj, Res_env, U=0, Q=abs(heatflux))
            else:
                heatflux = 0

            self.heat += heatflux

        # outlet
        if zone_type == 36 or zone_type == 1036:
            if mflux != 0:
                self.outlet(mfc_gen, mflux, reactor, reactor_obj)

        """print "number of outlets:", len(self.Outletmfc)
        print "number of inlets:", len(self.Inletmfc)
        for i in self.Inletmfc:
            print "Inlet flow=", self.Inletmfc[i].mdot(0)
        print "Relating periodic faces"
        Periodic_exclud = []
        mfc_per = {}
        mflow_shadow = 0"""
