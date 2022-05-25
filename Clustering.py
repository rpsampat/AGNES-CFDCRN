import cantera as ct
import math
import Update_cluster as UC
import random
import numpy as np

class Clustering:

    def __init__(self, graph_dat, zone, header, criteria, extra_num, start_ratio, criteria_derived, tol_deriv,
        chem_mech, id_dict, data_num, symm_axis, tol):

        self.zone = zone
        self.header = header
        self.criteria = criteria
        self.criteria_derived = criteria_derived
        self.tol = tol
        self.tol_deriv = tol_deriv
        self.start_loc = start_ratio
        self.id_dict = id_dict
        self.graph_dat = graph_dat
        self.chem_mech = chem_mech
        self.data_num = data_num
        self.symm_axis = symm_axis
        self.extra_num = extra_num
        self.key = []

    def startpoint(self, graph):
        if self.start_loc in self.id_dict:
            start_face = self.zone[self.id_dict[self.start_loc][0]][1]
            for cells in graph:
                if start_face in graph[cells]['faces']:
                    break
            startpt = cells

        else:
            # startpt at reactor containing cell number specified
            startpt = int(self.start_loc)
            for cells in graph:
                if startpt in graph[cells]['cells']:
                    break
            startpt = cells

        try:
            trial = 1 / startpt
        except:
            print "Failed start point=", startpt
            startpt = 1
        print "start point=", startpt

        return startpt

    def enqueue(self,Q, elem):
        Q.append(elem)

    def dequeue(self,Q):
        elem = Q[0]
        Q.remove(Q[0])
        return elem

    def cluster_initialise(self, graph, sigma):
        cluster = {}
        startpt = random.randint(1, len(graph.keys()))
        #startpt = self.startpoint(graph)
        graph[startpt]['color'] = 'gray'
        graph[startpt]['distance'] = 0
        graph[startpt]['predecessor'] = 1
        Q = []

        self.enqueue(Q, startpt)
        # searching for cell zone in which cell with index 1 is present
        strt_face = list(graph[startpt]['faces'])
        startpt_cell = graph[startpt]['cells'][0]
        z_fid_list = []
        for i in self.zone:
            if self.zone[i][0] == 12 and startpt_cell >= self.zone[i][1] and startpt_cell <= self.zone[i][2]:
                # searching for cell zone
                zid = i
            elif self.zone[i][0] == 13:
                for f_ref in strt_face:
                    if f_ref >= self.zone[i][1] and f_ref <= self.zone[i][2]:
                        # searching for face zone
                        z_fid_list.append(i)
            else:
                pass
        z_fid_list = np.unique(sorted(z_fid_list))
        if len(z_fid_list) > 2:
            z_fid = -1
        else:
            z_fid = z_fid_list[-1]

        # searching for cell zone in which cell with index 1 is present
        """for i in self.zone:
            if self.zone[i][0] == 12:
                if self.zone[i][2] >= startpt >= self.zone[i][1]:
                    break"""
        # initialising cluster data dictionary: same format as graph_dat
        #zid = i
        cluster['info'] = {zid: {}}
        for j in self.graph_dat[zid]:
            cluster['info'][zid][j] = [self.graph_dat[zid][j][startpt - self.zone[zid][1]]]

        # print data_in
        mu = []
        x_var = []
        for cri in self.criteria:
            val = cluster['info'][zid][self.header.index(cri) + 1][0]
            mu.append(val)
            x_var.append([val])

        cluster["graph"] = {1: {"adjacent": [], 'color': 'white', 'distance': 'inf', 'predecessor': '',
                                'faces': list(graph[startpt]['faces']), 'nodes': list(graph[startpt]['nodes']),
                                'cells': list(graph[startpt]['cells']),
                                'sigma': sigma,  # sigma: allowed deviation in cluster
                                'mu': list(mu),  # mu: average value of cluster
                                'n': len(graph[startpt]['cells']), 'x_var': list(x_var)}}

        gas = ct.Solution(self.chem_mech)

        # species taken from chemical mechanism
        species = gas.species_names
        for i in self.graph_dat[self.data_num].keys():
            s1 = self.header[i - 1].find('Mole fraction of')
            s2 = self.header[i - 1].find('<')
            s3 = self.header[i - 1].find('>')
            name = self.header[i - 1][len('Mole fraction of '):]
            if s1 != -1 and s2 == -1 and s3 == -1 and (name.upper() in species):
                self.key.append(name)

        return cluster, zid,z_fid, Q, startpt, gas

    def velocity_deriv_criteria(self, zid, count, u, i, cluster):
         # difference of nodes for derived criteria
         index = self.header.index('Progress Variable') + 1
         prog_var = self.graph_dat[i][index][u - self.zone[i][1]]
         pv_react = cluster['info'][zid][index][count - 1]
         pv_diff = abs(prog_var - pv_react)
         # only includes velocity if present in header, otherwise zero.
         x_vel = ['X Velocity', 'Axial Velocity']
         y_vel = ['Y Velocity', 'Radial Velocity']
         z_vel = ['Z Velocity', 'Z Velocity']

         vel_u1, vel_u2, vel_v1, vel_v2, vel_w1, vel_w2 = 0, 0, 0, 0, 0, 0

         for vel in range(len(x_vel)):
             if x_vel[vel] in self.header:
                 index = self.header.index(x_vel[vel]) + 1
                 vel_u1 = cluster['info'][zid][index][count - 1]
                 vel_u2 = self.graph_dat[i][index][u - self.zone[i][1]]

             if y_vel[vel] in self.header:
                 index = self.header.index(y_vel[vel]) + 1
                 vel_v1 = cluster['info'][zid][index][count - 1]
                 vel_v2 = self.graph_dat[i][index][u - self.zone[i][1]]

             if z_vel[vel] in self.header:
                 index = self.header.index(z_vel[vel]) + 1
                 vel_w1 = cluster['info'][zid][index][count - 1]
                 vel_w2 = self.graph_dat[i][index][u - self.zone[i][1]]

         index = self.header.index('Velocity Magnitude') + 1
         vm1 = cluster['info'][zid][index][count - 1]
         vm2 = self.graph_dat[i][index][u - self.zone[i][1]]

         try:
             sin_ang = math.sqrt(
                 (vel_v1 * vel_w2 - vel_v2 * vel_w1) ** 2 + (vel_u1 * vel_w2 - vel_u2 * vel_w1) ** 2 +
                 (vel_u1 * vel_v2 - vel_u2 * vel_v1) ** 2) / (vm1 * vm2)
         except ZeroDivisionError:
             sin_ang = 0

         # Vel dir:  TODO: remove?
         # vel_val1 = vel_u1 * vel_u2
         # vel_val2 = vel_v1 * vel_v2
         # vel_val3 = vel_w1 * vel_w2

         if self.symm_axis == 'X-Coordinate' and vel_u1 * vel_u2 >= 0:
             vel_dir = self.tol_deriv
         elif self.symm_axis == 'Y-Coordinate' and vel_v1 * vel_v2 >= 0:
             vel_dir = self.tol_deriv
         elif self.symm_axis == 'Z-Coordinate' and vel_w1 * vel_w2 >= 0:
             vel_dir = self.tol_deriv
         else:
             vel_dir = self.tol_deriv + 1

         if self.criteria_derived == 'Progress Variable':
             derived = pv_diff
         elif self.criteria_derived == 'Velocity Angle':
             derived = sin_ang
         elif self.criteria_derived == 'Velocity Dir':
             derived = vel_dir
         else:
             derived = self.tol_deriv

         vel_list = [vel_u1, vel_u2, vel_v1, vel_v2, vel_w1, vel_w2,]
         return derived, vel_list

    def vorticity_func(vor_x, vor_y):
        return math.atan(vor_y / vor_x)

    def fitness(self, quantity, dict_dat, graph, extra):
        """
        :param quantity:
        :param dict_dat:
        :param graph:
        :param extra:
        :return:
        """
        mu = list(dict_dat['mu'])
        sigma = list(dict_dat['sigma'])
        comp = 0.0
        f1 = set(graph['faces'])
        f2 = set(dict_dat['faces'])
        # list of common faces between cell and reactor to which it should be added
        f_com = list(f1 & f2)
        # common nodes
        n1 = set(graph['nodes'])
        n2 = set(dict_dat['nodes'])
        n_com = list(n1 & n2)

        if extra == 1:
            vor_theta = self.vorticity_func(quantity[0], quantity[1])
            mu_vor = self.vorticity_func(mu[0], mu[1])
            sig_vor = self.vorticity_func(sigma[0], sigma[1])
            comp = abs(vor_theta - mu_vor) / sig_vor
        else:
            for quant in range(len(quantity)):
                dev = (abs(quantity[quant] - mu[quant])) / abs(mu[quant])
                if dev <= sigma[quant]:
                    continue
                else:
                    comp += 1.0
        # checking if difference of quantity is within limits
        # and whether cells share a common face
        if comp == 0.0 and (f_com != [] or len(n_com) > 1.0):
            return True
        else:
            return False

    def bfs(self, graph, sigma):
        cluster, zid, z_fid, Q, startpt, gas = self.cluster_initialise(graph, sigma)
        neighbours = {1: []}
        within = {}
        neighbours[1] = list(graph[startpt]['adjacent'])
        within[startpt] = 1
        cells_explored = 0
        count = 1 # number of clusters created

        uc = UC.UpdateCluster(self)
        # BFS starts
        print "BFS started"
        while Q:

            u = self.dequeue(Q)

            for element in graph[u]['adjacent']:
                if element != 0:
                    if graph[element]['color'] == 'white':
                        graph[element]['color'] = 'gray'
                        graph[element]['distance'] = graph[u]['distance'] + 1
                        graph[element]['predecessor'] = u
                        self.enqueue(Q, element)
            # search for zone in which cell exists
            u_face = graph[u]['faces']
            u_cell = graph[u]['cells'][0]
            zid_curr = -1
            z_fid_curr = -1
            z_fid_list = []
            for i in self.zone:
                if self.zone[i][0] == 12 and u_cell >= self.zone[i][1] and u_cell <= self.zone[i][2]:
                    # searching for cell zone
                    zid_curr = i
                elif self.zone[i][0] == 13:
                    for f_ref in u_face:
                        if f_ref >= self.zone[i][1] and f_ref <= self.zone[i][2]:
                            # searching for face zone
                            z_fid_list.append(i)
                else:
                    pass
            # cell and face can exist in different types of zones. If so, the last zone by ascending order is chosen.
            # refer to piv_readdata.py for zoneid numbers.
            z_fid_list = np.unique(sorted(z_fid_list))
            if len(z_fid_list) > 2:
                z_fid_curr = -1
            else:
                z_fid_curr = z_fid_list[-1]
            # search for zone in which cell exists
            """for i in self.zone:
                if self.zone[i][0] == 12:
                    if self.zone[i][2] >= u >= self.zone[i][1]:
                        break"""

            q1 = []
            q_ind = []
            # extracting value of criteria for current cell
            for cri in self.criteria:
                spot = self.header.index(cri)
                q_ind.append(spot)
                q1.append(self.graph_dat[zid_curr][spot + 1][u - self.zone[zid_curr][1]])
            extra = self.extra_num

            # difference of nodes for derived criteria
            index = self.header.index('Progress Variable') + 1
            prog_var = self.graph_dat[zid_curr][index][u - self.zone[zid_curr][1]]
            pv_react = cluster['info'][zid][index][count - 1]
            pv_diff = abs(prog_var - pv_react)

            derived, vel_list = self.velocity_deriv_criteria(zid, count, u, zid_curr, cluster)

            # fitness check
            if self.fitness(q1, cluster['graph'][count], graph[u],extra) and u != startpt and zid == zid_curr and z_fid==z_fid_curr:#and derived <= self.tol_deriv
                # and check_zone(q1, avg_dict) == cluster['graph'][count]['zone']:
                #print "fit count =",count
                n = cluster['graph'][count]['n']
                cluster['graph'][count]['n'] = n + len(graph[u]['cells'])

                uc.update_values(cluster, zid, count, zid_curr, u, vel_list, gas)
                #uc.update_network(cluster, count, zid, u, i, graph)
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
                for cri in self.criteria:
                    mu.append(cluster['info'][zid][self.header.index(cri) + 1][count - 1])
                    index = self.header.index(cri) + 1
                    x_var.append(self.graph_dat[zid_curr][index][u - self.zone[zid_curr][1]])
                for i in range(len(x_var)):
                    cluster['graph'][count]['x_var'][i].append(x_var[i])
                cluster['graph'][count]['mu'] = list(mu)

            elif u != startpt:
                # creating new cluster
                count += 1
                for j in self.graph_dat[zid_curr]:
                    try:
                        cluster['info'][zid_curr][j].append(self.graph_dat[zid_curr][j][u - self.zone[zid_curr][1]])
                    except KeyError:
                        cluster['info'][zid_curr] = {j: [self.graph_dat[zid_curr][j][u - self.zone[zid_curr][1]]]}
                mu = []
                x_var = []
                for cri in self.criteria:
                    val = cluster['info'][zid_curr][self.header.index(cri) + 1][count - 1]
                    mu.append(val)
                    x_var.append([val])


                cluster['graph'][count] = {"adjacent": [], 'color': 'white', 'distance': 'inf', 'predecessor': '',
                                           'faces': list(graph[u]['faces']), 'nodes': list(graph[u]['nodes']),
                                           'cells': list(graph[u]['cells']), 'sigma': sigma,
                                           'mu': list(mu), 'n': len(graph[u]['cells']), 'x_var': list(x_var)
                                          }
                zid = zid_curr
                z_fid = z_fid_curr
                """if len(graph)<800:
                    print zid
                    print zid_curr
                    print self.fitness(q1, cluster['graph'][count-1], graph[u],extra)
                    print "Check"""""

            try:
                neighbours[count] = list(set(graph[u]['adjacent']) | set(neighbours[count]))
            except KeyError:
                neighbours[count] = list(graph[u]['adjacent'])

            if u in neighbours[count]:
                neighbours[count].remove(u)
            if 0 in neighbours[count]:
                neighbours[count].remove(0)

            within[u] = count
            graph[u]['color'] = 'black'
            cells_explored += 1

        return cluster, neighbours, within

    def cluster(self, graph, sigma):
        cluster, neighbours, within = self.bfs(graph, sigma)
        for g in neighbours:
            # g:reactor number
            for c in neighbours[g]:
                # c:cell number
                if g != within[c]:
                    # neighbouring reactors of 'g'
                    cluster["graph"][g]["adjacent"] = list(set(list(cluster["graph"][g]["adjacent"]) + [within[c]]))

        return cluster