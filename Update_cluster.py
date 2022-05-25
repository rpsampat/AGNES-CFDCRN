import math
import cantera as ct


class UpdateCluster:

    def __init__(self,cluster_obj):
        self.obj = cluster_obj


    def update_network(self, cluster, count, zid, u, i, graph):
        # Adding faces: set() allows only one occurence of common faces
        s1 = set(cluster['graph'][count]['faces'])
        s2 = set(graph[u]['faces'])

        # symmetric difference of sets, i.e elements in either s1 or s2, not in both
        # motive to remove interior faces and keep only external faces
        cluster['graph'][count]['faces'] = list(s1 ^ s2)
        # Adding cells
        c1 = cluster['graph'][count]['cells']
        c2 = graph[u]['cells']
        cluster['graph'][count]['cells'] = list(c1 + c2)
        mu = []
        x_var = []

        # updating average values of quantities to be compared
        for cri in self.obj.criteria:
            mu.append(cluster['info'][zid][self.obj.header.index(cri) + 1][count - 1])
            index = self.obj.header.index(cri) + 1
            x_var.append(self.obj.graph_dat[i][index][u - self.obj.zone[i][1]])
        for i in range(len(x_var)):
            cluster['graph'][count]['x_var'][i].append(x_var[i])
        cluster['graph'][count]['mu'] = list(mu)

    def update_values(self,cluster, zid, count, i, u, vel_list, gas):
        #gas = ct.Solution(self.obj.chem_mech)
        # Cell Volume
        vol = self.obj.header.index('Cell Volume') + 1
        v1 = cluster['info'][zid][vol][count - self.obj.zone[i][1]]
        v2 = self.obj.graph_dat[i][vol][u - self.obj.zone[i][1]]
        v = v1 + v2
        cluster['info'][zid][vol][count - self.obj.zone[i][1]] = v

        # Mass
        density = self.obj.header.index('Density') + 1
        d1 = cluster['info'][zid][density][count - self.obj.zone[i][1]]
        d2 = self.obj.graph_dat[i][density][u - self.obj.zone[i][1]]
        m1 = d1 * v1
        m2 = d2 * v2
        m = m1 + m2

        # Density
        index = self.obj.header.index('Density') + 1
        cluster['info'][zid][index][count - self.obj.zone[i][1]] = m / v

        # Progress variable
        index = self.obj.header.index('Progress Variable') + 1
        pv1 = cluster['info'][zid][index][count - 1]
        pv2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        pv = (m1 * pv1 + m2 * pv2) / m
        cluster['info'][zid][index][count - 1] = pv

        # RMS Temperature  TODO: remove?
        """rtemp = header.index('RMS Temperature') + 1
        RT1 = cluster['info'][zid][rtemp][count - zone[i][1]]
        RT2 = graph_dat[i][rtemp][u - zone[i][1]]
        RT = (m1*RT1 + m2*RT2)/m
        cluster['info'][zid][rtemp][count - zone[i][1]] = RT"""

        # Total Enthalpy
        index = self.obj.header.index('Specific Heat (Cp)') + 1
        cp1 = cluster['info'][zid][index][count - self.obj.zone[i][1]]
        cp2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]

        index = self.obj.header.index('Total Temperature') + 1
        T1 = cluster['info'][zid][index][count - self.obj.zone[i][1]]
        T2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        H1 = m1 * cp1 * T1
        H2 = m2 * cp2 * T2
        H = H1 + H2

        index = self.obj.header.index('Enthalpy') + 1
        cluster['info'][zid][index][count - 1] = H
        if H < 0:
            print H
            print count
            print H1
            print H2

        # Moles
        index = self.obj.header.index('Static Temperature') + 1
        Ts1 = cluster['info'][zid][index][count - 1]
        Ts2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        R = 8.314

        index = self.obj.header.index('Static Pressure') + 1
        ps1 = cluster['info'][zid][index][count - 1] + 101325
        ps2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]] + 101325
        n1 = ps1 * v1 / (R * Ts1)
        n2 = ps2 * v2 / (R * Ts2)
        n = n1 + n2

        index = self.obj.header.index('Mean Molecular Weight') + 1
        mw1 = cluster['info'][zid][index][count - 1]
        mw2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        mw = (n1 * mw1 + n2 * mw2) / n

        # Molecular weight
        index = self.obj.header.index('Mean Molecular Weight') + 1
        cluster['info'][zid][index][count - 1] = mw

        # Species
        for h in self.obj.header:
            """if h.find('Mole Fraction') != -1:
                index = header.index(h) + 1
                X1 = cluster['info'][zid][index][count - 1]
                X2 = graph_dat[i][index][u - zone[i][1]]
                cluster['info'][zid][index][count - 1] = (n1 * X1 + n2 * X2) / n"""

            if h.find('Mass Fraction ') != -1:
                index = self.obj.header.index(h) + 1
                Y1 = cluster['info'][zid][index][count - 1]
                Y2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
                cluster['info'][zid][index][count - 1] = (m1 * Y1 + m2 * Y2) / m

        # Cp
        index = self.obj.header.index('Specific Heat (Cp)') + 1
        cp1 = cluster['info'][zid][index][count - 1]
        cp2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        Cp0 = (m1 * cp1 + m2 * cp2) / m
        error = 1

        while error > 1e-6:
            Tt = H / (m * Cp0)

            Y = {j.upper(): cluster['info'][zid][self.obj.header.index('Mass fraction of ' + j) + 1][count - 1]
                 for j in self.obj.key}

            gas.TDY = Tt, m / v, Y
            Cp = gas.cp
            error = abs(Cp0 - Cp) / Cp0
            Cp0 = Cp

        index = self.obj.header.index('Specific Heat (Cp)') + 1
        cluster['info'][zid][index][count - 1] = Cp

        # Total Temperature
        index = self.obj.header.index('Total Temperature') + 1
        cluster['info'][zid][index][count - 1] = Tt

        """# Total Pressure
        index = header.index('Total Pressure') + 1
        P=n*R*Tt/v
        cluster['info'][zid][index][count - 1] = P - 101325"""

        # Vorticity and Velocity
        vel_u, vel_v, vel_w = 0, 0, 0
        x_vel = ['X Velocity', 'Axial Velocity']
        y_vel = ['Y Velocity', 'Radial Velocity']
        z_vel = ['Z Velocity', 'Z Velocity']
        vel_u1 = vel_list[0]
        vel_u2 = vel_list[1]
        vel_v1 = vel_list[2]
        vel_v2 = vel_list[3]
        vel_w1 = vel_list[4]
        vel_w2 = vel_list[5]

        for vel in range(len(x_vel)):
            if x_vel[vel] in self.obj.header:
                index = self.obj.header.index(x_vel[vel]) + 1
                vel_u = (m1 * vel_u1 + m2 * vel_u2) / m
                cluster['info'][zid][index][count - 1] = vel_u

            if y_vel[vel] in self.obj.header:
                index = self.obj.header.index(y_vel[vel]) + 1
                vel_v = (m1 * vel_v1 + m2 * vel_v2) / m
                cluster['info'][zid][index][count - 1] = vel_v

            if z_vel[vel] in self.obj.header:
                index = self.obj.header.index(z_vel[vel]) + 1
                vel_w = (m1 * vel_w1 + m2 * vel_w2) / m
                cluster['info'][zid][index][count - 1] = vel_w

        index = self.obj.header.index('Velocity Magnitude') + 1
        vel = math.sqrt(vel_u ** 2 + vel_v ** 2 + vel_w ** 2)
        cluster['info'][zid][index][count - 1] = vel

        index = self.obj.header.index('Radial Velocity') + 1
        vel_vr1 = cluster['info'][zid][index][count - 1]
        vel_vr2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        vel_vr = (m1 * vel_vr1 + m2 * vel_vr2) / m
        cluster['info'][zid][index][count - 1] = vel_vr

        for h in self.obj.header:
            if h.find('Vorticity') != -1:
                index = self.obj.header.index(h) + 1
                X1 = cluster['info'][zid][index][count - 1]
                X2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
                cluster['info'][zid][index][count - 1] = (m1 * X1 + m2 * X2) / m

        # Static Temperature
        index = self.obj.header.index('Static Temperature') + 1
        Ts = (H / m - (vel ** 2) / 2) / Cp
        cluster['info'][zid][index][count - 1] = Ts
        #if math.isnan(Ts):
          #  continue

        # Static pressure
        index = self.obj.header.index('Static Pressure') + 1
        ps = gas.P
        cluster['info'][zid][index][count - 1] = ps - 101325

        # Turbulence
        for h in self.obj.header:
            if h.find('Turbulent') != -1:
                index = self.obj.header.index(h) + 1
                X1 = cluster['info'][zid][index][count - 1]
                X2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
                cluster['info'][zid][index][count - 1] = (m1 * X1 + m2 * X2) / m

        # X-Coordinate
        index = self.obj.header.index('X-Coordinate') + 1
        x1 = cluster['info'][zid][index][count - 1]
        x2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        x = (m1 * x1 + m2 * x2) / m
        cluster['info'][zid][index][count - 1] = x

        # Y-Coordinate
        index = self.obj.header.index('Y-Coordinate') + 1
        y1 = cluster['info'][zid][index][count - 1]
        y2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        y = (m1 * y1 + m2 * y2) / m
        cluster['info'][zid][index][count - 1] = y

        # Z-Coordinate
        index = self.obj.header.index('Z-Coordinate') + 1
        z1 = cluster['info'][zid][index][count - 1]
        z2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        z = (m1 * z1 + m2 * z2) / m
        cluster['info'][zid][index][count - 1] = z

        # Turbulent Viscosity
        index = self.obj.header.index('Turbulent Viscosity') + 1
        mut1 = cluster['info'][zid][index][count - 1]
        mut2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        mut = (m1 * mut1 + m2 * mut2) / m
        cluster['info'][zid][index][count - 1] = mut

        # Mass Imbalance
        index = self.obj.header.index('Mass Imbalance') + 1
        imb1 = cluster['info'][zid][index][count - 1]
        imb2 = self.obj.graph_dat[i][index][u - self.obj.zone[i][1]]
        imb = imb1 + imb2
        cluster['info'][zid][index][count - 1] = imb
