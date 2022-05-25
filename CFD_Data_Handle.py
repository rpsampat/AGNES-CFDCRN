import pickle

class CfdDataHandle:
    def __init__(self, data_std, zone, facearea):
        self.data_std = data_std
        self.zone = zone
        self.facearea = facearea
        self.a = 10
        self.mflux_std = {}
        self.heatflux_std = {}


    def mflux(self):
        """
        Converts mass flux data extracted from CFD into a format independent of zones.
        :return:
        """

        face = self.facearea.keys()
        for f in face:
            for z in self.zone:
                if self.zone[z][0] == 13 and self.zone[z][2] >= f >= self.zone[z][1]:
                    mflux = self.data_std[z][18][f - self.zone[z][1]]
                    self.mflux_std[f] = mflux
    def heatflux(self):
        face = self.facearea.keys()
        for f in face:
            for z in self.zone:
                if self.zone[z][0] == 13 and self.zone[z][2] >= f >= self.zone[z][1]:
                    heatflux = self.data_std[z][20][f - self.zone[z][1]]# 20 index in xfile.h for boundary heat flux
                    self.heatflux_std[f] = heatflux

    def main(self):
        self.mflux()
        self.heatflux()
