import time
import Clustering
import SaveReactors as sr
from BFS_pchanged_monaghan import bfs
from Para_plot import *


def generate(num, criteria, criteria_derived, extra, addon, startpt, chem_mech, tol_lim, id_dict, data_num, symm_axis,
             avg_diction):
    """
    :param num: threshold number of reactors
    :param criteria: list of criteria in string format
    :param criteria_derived: user defined criteria, not present in Fluent
    :param extra: 1/0, to activate extra parameter function check for clustering
    :param addon: addition to Folder name
    :param startpt: starting point of  clustering, options mentioned in startlist
    :param chem_mech: chemical mechanism used
    :param tol_lim:
    :param id_dict:
    :param data_num:
    :param symm_axis:
    :return:
    """

    nCRN = num
    t0 = time.time()
    t1 = time.clock()
    print "Reading Mesh and Data required for Clustering"
    with open("mesh.pkl", 'rb') as f:
        graph = pickle.load(f)
    with open("zone.pkl", 'rb') as f:
        zone = pickle.load(f)
    print "Completed reading mesh"
    with open("data.pkl", 'rb') as f:
        graph_dat = pickle.load(f)
    with open("header.pkl", 'rb') as f:
        header = pickle.load(f)
    with open("facearea.pkl", 'rb') as f:
        facelist = pickle.load(f)
    print "Completed reading data"

    # tolerance for clustering
    tol = [0.1]*len(criteria)
    # increment of tolerance
    increment = [0.01]*len(criteria)
    # tolerance for derived criteria
    tol_deriv = 0.258  # sin(15)=0.2588
    # increment of tolerance for derived criteria
    incr_deriv = 0.01
    print 'criteria:', criteria
    avg_base = []
    stdev_base = []

    # the following section of code was an attempt to generate gradients of selected criteria in the entire
    # cfd domain before clustering. The intentino was to use these gradients as a clustering criteria.
    # It turned out to be too inefficient in its current form, hence was abandoned, but this section of the code
    # is kept to facilitate future modifications.
    data_univ = {}
    gradient = []
    # for i in zone:
    #     if i==data_num:  #zone[i][0] == 12:
    #         data_univ[header.index('X-Coordinate') + 1] = graph_dat[i][header.index('X-Coordinate') + 1]
    #         data_univ[header.index('Y-Coordinate') + 1] = graph_dat[i][header.index('Y-Coordinate') + 1]
    #         # data_univ[header.index('Z-Coordinate') + 1] = graph_dat[i][header.index('Z-Coordinate') + 1]
    #         data_univ[header.index('Static Pressure') + 1] = graph_dat[i][header.index('Static Pressure') + 1]
    #         gradient = graph_dat[i][header.index('dX-Velocity/dx') + 1]
    #         numer = graph_dat[i][header.index('X Velocity') + 1]
    #         grad_max=max(gradient)/(max(numer)-min(numer))
    #         break
    # cell_volume=graph_dat[i][header.index('Cell Surface Area') + 1]
    # grad_max=1/((min(cell_volume))**(1/3))
    # end of abandoned code
    for cri in criteria:
        q_max = []
        q_min = []
        for z in graph_dat:
            if (header.index(cri) + 1 in graph_dat[z]) and zone[z][0] == 12:
                quantities = graph_dat[z][header.index(cri) + 1]
                q_max.append(max(quantities))
                q_min.append(min(quantities))
                data_univ[header.index(cri) + 1] = quantities
                stdev = np.std(quantities)
                avg = np.average(quantities)
        # avg_base.append((max(q_max) - min(q_min)))
        avg_base.append(avg)
        stdev_base.append(stdev)
    #import matplotlib.pyplot as plt
    #xlist = np.linspace(1,len(quantities),num=len(quantities))
    #plt.scatter(xlist,quantities)
    #plt.hist(quantities,bins ='auto')#int(len(quantities)/2.0))
    #plt.show()
    print "Average=", avg_base
    print "Standard dev=", stdev_base
    print "Tol range=", stdev_base[0] / avg_base[0]

    num_reactor = [len(graph)]
    tol_conv = []

    # list of available start points
    startlist = list(id_dict.keys() + [1])
    red_count = 0
    ind_start = startlist.index(startpt)
    data_std = graph_dat
    avg_diction = []
    for cri in criteria:
        avg_diction.append(tol[criteria.index(cri)])  # * avg_base[criteria.index(cri)])
    print avg_diction

    # Starting Clustering iterations
    for loop in range(300000):
        """if avg_diction == {}:
            avg_diction[tuple(['min', 'max'] * len(criteria))] = []
            for cri in criteria:
                avg_diction[tuple(['min', 'max'] * len(criteria))].append(tol[criteria.index(cri)] *
                                                                          avg_base[criteria.index(cri)])"""
        for cri in range(len(criteria)):
            avg_diction[cri] = (tol[cri])# * avg_base[cri])
        length = len(graph)
        tol_conv.append(tol[0])
        print "Number of Clusters=", length
        if length < nCRN or tol[0] > tol_lim:
            break
        # Rotating start point: the starting point of clustering is
        # rotated from amongst the elements of startlist,
        # the point for the first iteration being specified by the variable 'startpt'
        startpt = startlist[ind_start]
        """cluster = bfs(graph, graph_dat, zone, header, criteria, extra, startpt, avg_diction, criteria_derived,
                      tol_deriv,
                      chem_mech, id_dict, data_num, symm_axis, tol)"""
        cs = Clustering.Clustering(graph_dat, zone, header, criteria, extra, startpt, criteria_derived,
                                   tol_deriv, chem_mech, id_dict, data_num, symm_axis, tol)
        cluster = cs.cluster(graph, avg_diction)
        ind_start += 1
        if ind_start >= len(startlist):
            ind_start = 0
        graph = cluster["graph"]
        graph_dat = cluster["info"]
        num_reactor.append(len(graph))
        if loop > 0:
            red2 = length-len(graph)
            if red2 < 0.01*length:
                red_count += 1
                if red_count>0.01*length:
                    tol = [t + increment[tol.index(t)] for t in tol]
                    tol_deriv += incr_deriv
                    print "tolerance =", tol
                    red_count = 0
            else:
                continue

    t01 = time.time()-t0
    t11 = time.clock() - t1
    print "time to cluster0 =", t01
    print "time to cluster1 =", t11
    print "tolerance =", tol

    # Creating folder
    loc = os.getcwd()
    path = loc + '/%i_%i' % (nCRN, len(graph))
    for title in criteria:
        path = path + '_' + title
    path = path + '_' + criteria_derived + '_' + addon
    try:
        os.stat(path)
    except:
        os.mkdir(path)

    with open(path+'/readme.txt', 'w') as file:
        file.write("Number of reactors specified=%i\n" % nCRN)
        file.write("Number of reactors created=%i\n" % len(graph))
        if isinstance(startpt, basestring):
            file.write("Start point="+startpt+'\n')
        else:
            file.write("Start point=%f\n" % startpt)
        file.write("Criteria of clustering=")
        for c in criteria:
            file.write(c)
        file.write("\nExtra=%i\n"%extra)
        file.write("Add on="+addon+"\n")
        file.write("Time to cluster(s)=%f\n" % t01)
        file.write("Time solving started=%s\n" % time.asctime())

    with open(path+'/graph2plot.pkl', 'wb') as file:
        pickle.dump(graph, file, pickle.HIGHEST_PROTOCOL)

    with open(path+'/graphdata.pkl', 'wb') as file:
        pickle.dump(graph_dat, file, pickle.HIGHEST_PROTOCOL)
    sr.save_data({1:num_reactor}, "NumberOfReactors", path)
    sr.save_data({1:tol_conv}, "Tolerance", path)
    para3D(graph_dat, header, path)
    t02 = time.time() - t01
    t12 = time.clock() - t11
    print "time to para0=", t02
    print "time to para1=", t12
   
    return path


if __name__ == "__main__":
    generate(1, ['Static Temperature'], '', 0, 'individual_run', 'fuel', 'gri30.cti', 3, {})
