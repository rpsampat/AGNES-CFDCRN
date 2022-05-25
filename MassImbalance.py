import pickle
import matplotlib.pyplot as plt
import numpy as np
from numpy import dot, linalg, allclose, subtract, bincount, take, sqrt, array, ones
from numpy import zeros
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix


class MassImbalance:
    def __init__(self, process):
        self.tres = []
        self.massflowreactor = []
        self.m_exit = 0  # exit mass flow from system
        self.mass_out = 0  # outlet mass flow accumulated over internal reactors
        self.count = 0
        self.MassImb2 = []
        self.process = process

    def print_prop(self):
        print "Mass Imbalance properties:"
        print "Exit mass flow = ", self.m_exit
        print "Outlet mass flow = ", self.mass_out
        print "Number of outlet mass flow controllers = ", self.count

    def mass_imbalance_sparse(self, mfc_gen, boundary, PSR, path):

        # Correcting for inherent mass imbalance from CFD cells
        num_mfc = 0
        null_react = []
        rows, cols = mfc_gen.coeffmat.nonzero()
        for pos in range(len(PSR)):
            bal = 0
            const = mfc_gen.coeffmat[pos][pos]

            # normalising coefficient matrix by dividing outflow mfcs from each reactor
            # by total outflow from reactor( values along diagonal)
            for vect in range(len(PSR)):
                val = mfc_gen.coeffmat[vect][pos]
                bal += mfc_gen.coeffmat[vect][pos]
                mfc_gen.coeffmat[vect][pos] = val / const  # reactor_out[pos]
                if mfc_gen.coeffmat[vect][pos] != 0:
                    num_mfc += 1
                if mfc_gen.coeffmat[vect][pos] > 1:
                    print "greater: ", mfc_gen.coeffmat[vect][pos]
                if const == 1 and val != 0:
                    print "const0:", val

        num_mfc = num_mfc - len(PSR)
        print "No mass flow reactor=", null_react
        print "Number of mfc=", num_mfc
        print "NUmber of mfc in record=", len(mfc_gen.mfc_rec)

        print "Solving for mfcs"
        x = linalg.solve(mfc_gen.coeffmat, mfc_gen.rhsvect)
        # Corrected mass flows
        self.massflowreactor = np.ones(np.size(self.PSR))  # list(x)

        with open(path + "/coeffmat.pkl", 'rb') as f:
            coeffmatsamp = pickle.load(f)
        """with open(path + "/massflowreactor.pkl", 'rb') as f:
            mfr_samp = pickle.load(f)"""

        for pos in range(len(PSR)):

            for vect in range(len(PSR)):

                comp = mfc_gen.coeffmat[vect][pos] - coeffmatsamp[vect][pos]
                if comp != 0:
                    print "Current coeff = ", coeffmatsamp[vect][pos]
                    print "Record coeff = ", mfc_gen.coeffmat[vect][pos]
                    print "row coeff = ", vect
                    print "pos coeff = ", pos

        """for l1 in range(len(mfr_samp)):
            if mfr_samp[l1]-self.massflowreactor[l1] !=0:
                print "Current mass flow = ", self.massflowreactor[l1]
                print "Record mass flow = ",mfr_samp[l1]
                print "Reactor no. = ",l1"""

        with open(path + '/massflowreactor.pkl', 'wb') as file:
            pickle.dump(self.massflowreactor, file, pickle.HIGHEST_PROTOCOL)
        with open(path + '/coeffmat.pkl', 'wb') as file:
            pickle.dump(mfc_gen.coeffmat, file, pickle.HIGHEST_PROTOCOL)
        with open(path + '/rhsvect.pkl', 'wb') as file:
            pickle.dump(mfc_gen.rhsvect, file, pickle.HIGHEST_PROTOCOL)

    def print_sparsematrix(self, matrix):
        """
        Creates a new matrix based on a sparse input matrix such that all nonzero values=1.
        Converts it to a dense form and prints it as a grey scale image
        :param matrix: csc matrix
        :return:
        """
        binary_list = ones(len(matrix.data))
        binary_shape = matrix.shape
        rows, cols = matrix.nonzero()
        binary_matrix = csc_matrix((binary_list, (rows, cols)), shape=binary_shape)
        binary_matrix = binary_matrix.todense()
        plt.figure()
        plt.imshow(binary_matrix, cmap="Greys")
        plt.colorbar()
        plt.show()

    def write_graph2file(self, coeffmat, mass_vect, rhs_vect, PSR, path):
        """
        writes 3 files; file_node, file_edge and file rhs to represent the nodes, edges (directed with massflow)
        and inlet boudnary conditions of the CRN.
        :param coeffmat: column normalised coefficient matrix in csc sparse format
        :param mass_vect: vector of total outlet massflow from each reactor
        :param rhs_vect: vector of inlet boundary conditions
        :param PSR:
        :param path:
        :return:
        """
        colptr = coeffmat.indptr
        rowind = coeffmat.indices
        nodes = []
        edges={}
        edge_count = 0
        file_node = open(path + '/nodes.txt', 'wb')
        file_node.write("Node\n")
        file_edge = open(path + '/edges.txt', 'wb')
        file_edge.write("Edge   From node   To node    Massflow\n")
        file_rhs = open(path + '/boundary.txt', 'wb')
        file_rhs.write("Node    Inlet flow\n")
        for p in range(np.size(PSR)):
            file_node.write(str(p)+'\n')
            if rhs_vect[p] >0.0:
                file_rhs.write(str(p)+"    "+str(rhs_vect[p])+'\n')
            rowrange = range(colptr[p], colptr[p + 1])
            for quo in rowrange:
                # elem = array([..list of elements..], dtype=...)
                coeff = coeffmat.data[quo]
                q = rowind[quo]
                # flow from reactor p to reactor q
                if p!=q:
                    edges[edge_count]={1:[p,q],2:coeff*mass_vect[p]}
                    file_edge.write(str(edge_count)+"   "+str(p)+"  "+str(q)+"  "+str(abs(coeff*mass_vect[p]))+'\n')
                    edge_count +=1

        file_node.close()
        file_edge.close()
        file_rhs.close()
    def PIVCRN_velocity(self, mfc_gen, PSR, path, boundary):
        """
        Generates record of total massflow out of a reactor and coefficient matrix of velocities between reactors for
        the PIV-CRN process.
        :param mfc_gen:
        :param PSR:
        :param path:
        :return:
        """
        # CSC can't use _cs_matrix's .nonzero method because it
        # returns the indices sorted for self transposed.
        # cols, rows = mfc_gen.coeffmat.nonzero()
        self.massflowreactor = np.zeros(np.size(mfc_gen.rhsvect))
        self.tres = 0.0001 * np.ones(np.size(mfc_gen.rhsvect))
        mass_sum =[]
        data_temp = mfc_gen.coeffmat.data
        print np.size(PSR)
        colptr = mfc_gen.coeffmat.indptr
        rowind = mfc_gen.coeffmat.indices
        for i in range(np.size(PSR)):
            rowrange = range(colptr[i], colptr[i + 1])
            mass_sum.append(0.0)
            """if i in boundary.OutletReact:
                # accounts for outflow from domain
                m_out = boundary.Outletmfc[i].mdot(0)
                mass_sum[i] = mass_sum[i]# + abs(m_out)"""
            for pt in rowrange:
                coeff = mfc_gen.coeffmat.data[pt]
                q = rowind[pt]
                if q==i:
                    mass_sum[i] = mass_sum[i] + abs(mfc_gen.coeffmat.data[pt])
                    #mfc_gen.coeffmat.data[pt] = 1.0
                else:
                    if coeff >0:
                        print coeff

            for pt2 in rowrange:
                q2 = rowind[pt2]
                mfc_gen.coeffmat.data[pt2] = mfc_gen.coeffmat.data[pt2] / mass_sum[i]
        """plt.plot(mfc_gen.rhsvect, linestyle='none', marker='.')
        plt.yscale('log')"""
        #plt.show()
        print mfc_gen.coeffmat
        print np.sum(mfc_gen.rhsvect)
        # coeffmat is normalised coefficient by total velocity flux from reactor
        # mass_sum is sum of total velocity flux from each reactor
        self.massflowreactor = spsolve(mfc_gen.coeffmat,mfc_gen.rhsvect)
        """fig = plt.figure()
        ax = fig.add_subplot(111)
        img = csc_matrix.todense(mfc_gen.coeffmat)
        im = ax.imshow(img)#, vmin=-1e-7, vmax=1e-7, interpolation='none',
                       #cmap=plt.cm.gray, aspect='auto')
        fig.colorbar(im)"""
        # imgplot = plt.imshow(img)
        # plt.clim(-1, 1)
        #plt.show()
        mfc_gen.rho_calc = np.divide(self.massflowreactor,mass_sum)
        mfc_gen.rho_calc = np.where(mfc_gen.rho_calc>1.5,1.5,mfc_gen.rho_calc)
        # Calculating total outflow mass
        for i in range(np.size(PSR)):
            if i in boundary.OutletReact:
                m_out = boundary.Outletmfc[i].mdot(0)
                self.m_exit = self.m_exit + abs(m_out)*mfc_gen.rho_calc[i]
        """plt.plot(self.massflowreactor, linestyle='none',marker='.')
        plt.yscale('log')
        plt.show()"""
        #plt.plot(mfc_gen.rho_calc, linestyle='none', marker='.')
        # plt.yscale('log')
        #plt.show()
        print np.sum(mfc_gen.rhsvect)
        print "Density=",mfc_gen.rho_calc
        with open(path + '/massflowreactor.pkl', 'wb') as file:
            pickle.dump(self.massflowreactor, file, pickle.HIGHEST_PROTOCOL)
        with open(path + '/coeffmat.pkl', 'wb') as file:
            pickle.dump(mfc_gen.coeffmat, file, pickle.HIGHEST_PROTOCOL)
        with open(path + '/rhsvect.pkl', 'wb') as file:
            pickle.dump(mfc_gen.rhsvect, file, pickle.HIGHEST_PROTOCOL)




    def mass_imbalance(self, mfc_gen, boundary, PSR, path):
        """
        Generates balanced massflow out of each reactor and coefficient matrix for CFD-CRN process where massflow
        between reactors is known.
        :param mfc_gen:
        :param boundary:
        :param PSR:
        :param path:
        :return:
        """
        # Correcting for inherent mass imbalance from CFD cells
        num_mfc = 0
        null_react = []
        if mfc_gen.mattype == "sparse":
            rows, cols = mfc_gen.coeffmat.nonzero()
            colptr = mfc_gen.coeffmat.indptr
            rowind = mfc_gen.coeffmat.indices
            # bin_wts = np.abs(mfc_gen.coeffmat.data)
            # col_norm = bincount(cols, weights=bin_wts)
            col_norm = zeros(len(PSR))
            # searching for diagonal elements
            for p in range(len(PSR)):
                rowrange = range(colptr[p], colptr[p + 1])
                for quo in rowrange:
                    # elem = array([..list of elements..], dtype=...)
                    coeff = mfc_gen.coeffmat.data[quo]
                    q = rowind[quo]
                    if p == q:
                        col_norm[p] = coeff
            # normalising every column
            for p in range(len(PSR)):
                rowrange = range(colptr[p], colptr[p + 1])
                for quo in rowrange:
                    q = rowind[quo]
                    mfc_gen.coeffmat.data[quo]/=col_norm[p]
            print "Solving for mfcs"
            x = spsolve(A=mfc_gen.coeffmat, b=mfc_gen.rhsvect, use_umfpack=True)
            print np.allclose(mfc_gen.coeffmat.dot(x), mfc_gen.rhsvect)
            #print x
            self.write_graph2file(mfc_gen.coeffmat, x, mfc_gen.rhsvect, PSR, path)
        elif mfc_gen.mattype == "dense":
            for pos in range(len(PSR)):
                bal = 0
                const = mfc_gen.coeffmat[pos][pos]

                # normalising coefficient matrix by dividing outflow mfcs from each reactor
                # by total outflow from reactor( values along diagonal)
                for vect in range(len(PSR)):
                    val = mfc_gen.coeffmat[vect][pos]
                    bal += mfc_gen.coeffmat[vect][pos]
                    mfc_gen.coeffmat[vect][pos] = val / const  # reactor_out[pos]
                    if mfc_gen.coeffmat[vect][pos] != 0:
                        num_mfc += 1
                    if mfc_gen.coeffmat[vect][pos] > 1:
                        print "greater: ", mfc_gen.coeffmat[vect][pos]
                    if const == 1 and val != 0:
                        print "const0:", val
            print "Solving for mfcs"
            x = linalg.solve(mfc_gen.coeffmat, mfc_gen.rhsvect)
        num_mfc = num_mfc - len(PSR)
        print "No mass flow reactor=", null_react
        print "Number of mfc=", num_mfc
        print "NUmber of mfc in record=", len(mfc_gen.mfc_rec)
        # Corrected mass flows
        self.massflowreactor = np.array(x)
        with open(path + '/massflowreactor.pkl', 'wb') as file:
            pickle.dump(self.massflowreactor, file, pickle.HIGHEST_PROTOCOL)
        with open(path + '/coeffmat.pkl', 'wb') as file:
            pickle.dump(mfc_gen.coeffmat, file, pickle.HIGHEST_PROTOCOL)
        with open(path + '/rhsvect.pkl', 'wb') as file:
            pickle.dump(mfc_gen.rhsvect, file, pickle.HIGHEST_PROTOCOL)

        print "Mass correction files written"
        for m_veri in self.massflowreactor:
            if m_veri == 0.0:
                null_react.append(m_veri)
            if m_veri <= 0:
                print "Negative mass flow=", m_veri
                print "Reactor=", self.massflowreactor.index(m_veri)

        # dotprod = dot(mfc_gen.coeffmat, x)
        dotprod = mfc_gen.coeffmat.dot(x)
        # print subtract(dotprod, mfc_gen.rhsvect)
        # answer = allclose(dot(mfc_gen.coeffmat, x), mfc_gen.rhsvect)
        answer = allclose(mfc_gen.coeffmat.dot(x), mfc_gen.rhsvect)
        print "Answer is closed:", answer
        mass_null = []

        for m in null_react:
            mass_null.append(x[m])

        print "Null flow reactor=", null_react
        self.MassImb2 = [0] * len(PSR)

        # Calculating Initial mass imbalance in CRN
        for i in range(len(PSR)):
            mfc_in = PSR[i].inlets
            mfc_out = PSR[i].outlets
            m_in = 0
            m_out = 0
            for j in mfc_in:
                m_in += j.mdot(0)
            for j in mfc_out:
                m_out += j.mdot(0)
            self.MassImb2[i] += m_in - m_out

        for r in mfc_gen.rhsvect:
            self.count += r
        print "Total inflow =", self.count

        self.count = 0

        # Reassigning mass flows to mass flow controllers after mass balance correction
        exception_flows = 0
        for p in PSR:
            # calculating residence time of reactor
            elem = np.where(PSR == p)
            self.tres.append(p.mass / self.massflowreactor[elem[0][0]])
            mfc_out = p.outlets  # might include valves also
            mfc_in = p.inlets  # might include valves also
            mfc_list = mfc_out + mfc_in
            elem = np.where(PSR == p)
            for mfc in mfc_list:
                try:
                    c1 = mfc_gen.mfc_rec[mfc][0]  # 'from' reactor
                except:  # for bc inlet mfcs not on list
                    try:
                        if mfc != boundary.FuelInletmfc[p] and mfc != boundary.AirInletmfc[p]:
                            print "Something's wrong with Periodic"
                            print mfc_gen.mfc_rec[mfc]
                    except:
                        # print mfc
                        # print "Mass flow=", mfc.mdot(0)
                        exception_flows += mfc.mdot(0)
                        try:
                            print mfc_gen.mfc_rec[mfc]
                        except:
                            continue
                    continue

                # outlet reactors
                if len(mfc_gen.mfc_rec[mfc]) == 1:
                    if c1 != elem[0][0]:
                        print "Error index =", c1
                        print "outlet reactor num =", PSR.index(p)
                    self.count += 1
                    m = 0
                    mass = 0
                    M_k = self.massflowreactor[c1]

                    for cell in range(len(PSR)):
                        # Todo: "dense": mass += self.massflowreactor[cell] * mfc_gen.coeffmat[c1][cell]
                        mass += self.massflowreactor[cell] * mfc_gen.coeffmat[c1, cell]
                        # m is sum of total outflow from outlet reactor and inlet flows from outlet reactor
                        # to other reactors(inlet flows denoted with minus sign
                        # Todo: "dense": m += mfc_gen.coeffmat[cell][c1] * M_k
                        m += mfc_gen.coeffmat[cell, c1] * M_k

                    # print "Massflow in outer cell =", mass
                    if m <= 0:
                        print "Cell num negative:", c1
                    self.m_exit += m
                    mfc.set_mass_flow_rate(m)
                    # print M_k-m

                # other internal reactors
                else:
                    M_k = self.massflowreactor[c1]  # outflow from kth reactor
                    c2 = mfc_gen.mfc_rec[mfc][1]  # 'to' reactor
                    if c1 != elem[0][0] and c2 != elem[0][0]:
                        print "Error index=", c1
                        print "reactor num=", elem[0][0]
                    # Todo: "dense"
                    if (c2 in boundary.OutletReact) and (c1 not in boundary.OutletReact):
                        self.mass_out += abs(mfc_gen.coeffmat[c2, c1] * M_k)
                    if (c1 in boundary.OutletReact) and (c2 not in boundary.OutletReact):
                        self.mass_out += -abs(mfc_gen.coeffmat[c2, c1] * M_k)

                    # reassigning mass flow rate of MFC, assuming only one MFC between two reactors in one direction!!
                    # Todo: "dense"
                    m_set = abs(mfc_gen.coeffmat[c2, c1] * M_k)
                    mfc.set_mass_flow_rate(m_set)

        print "Total exception flow=", exception_flows
        # Correcting Valve Coefficients based on related MFCs
        for v in mfc_gen.ValveRecord:
            mfc_ref = mfc_gen.ValveRecord[v]
            kv_valve = mfc_ref.mdot(0) * mfc_gen.valve_mult
            # TODO: check valve coefficient
            v.set_valve_coeff(kv_valve)  # lambda dP: mfc_ref.mdot(0)*valve_mult*(dP))
