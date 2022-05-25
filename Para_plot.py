from pyevtk.hl import*
import numpy
import pickle
import os


def graph2para(dict, header):
    """
    This function converts a graph=={reactor0:{'cells':{..},
    'info':{'temperature':...,'pressure':..},reactor1:..}
    to a paraview compatible dictionary=={'temperature':numpyarray[R1_temp,R2_temp..]
    'pressure':numpyarray[R1p,R2p..]}
    :param dict: graph
    :param header:
    :return: dictionary of nummpy arrays
    """
    para = {}
    for i in range(len(header)):
        print header[i]
        if header[i].find("c<s>") == -1:
            temp = []
            for j in dict.keys():
                temp.append(dict[j][i])
            para[header[i]] = numpy.array(temp, dtype=numpy.float64)

    return para


def para3D(dict, header, path):
    """
    converts point data to VTK format
    :param dict:
    :param header:
    :param path:
    :return: 
    """
    print dict.keys()
    for zone in dict:
        print dict[zone].keys()
        x = header.index('X-Coordinate')+1
        y = header.index('Y-Coordinate')+1
        z = header.index('Z-Coordinate')+1
        # checking if zone is cell zone, indicated by presence of coordinate data
        if x in dict[zone]:
            X = numpy.array(dict[zone][x], dtype=numpy.float64)
            Y = numpy.array(dict[zone][y], dtype=numpy.float64)
            Z = numpy.array(dict[zone][z], dtype=numpy.float64)
            dat = {}
            for index in dict[zone]:
                if header[index-1].find('<') == -1 and header[index-1].find('>') == -1:
                    dat[header[index-1]] = numpy.array(dict[zone][index], dtype=numpy.float64)
            pointsToVTK(path+"/para_orig"+str(zone), X, Y, Z, data=dat)
        else:
            continue


if __name__ == "__main__":
    path = os.getcwd()
    with open("data.pkl", 'rb') as f:
        data = pickle.load(f)
    with open("header.pkl", 'rb') as f:
        header = pickle.load(f)
    para3D(data, header, path)
