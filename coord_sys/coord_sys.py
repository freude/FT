#!/usr/bin/python

import numpy as np
import myconst as MyConst
from sys import stdout


class CoordSys(object):

    def __init__(self, **kw):

        self.num_cells = kw.get('num_cells', 10)
        self.array_sizes = kw.get('array_sizes', 8)
        units = kw.get('units', 'au')

        if units == 'au':
            self.lattice_const = MyConst.a_Si/MyConst.ab
        else:
            self.lattice_const = MyConst.a_Si

        self.coord_sizes = self.num_cells*self.array_sizes+1

        self.length = self.num_cells*self.lattice_const

        self.coord_stps=self.length/(self.coord_sizes-1)

        self.origin_cells = self.num_cells/2
        self.origin_inds=self.array_sizes*(self.origin_cells-1)+1;
        self.coord_limits = [-self.origin_cells*self.lattice_const,\
                             (self.num_cells-self.origin_cells)*self.lattice_const]

    #----------------------------------------------------------------------

    def set_origin_cells(self, origin_cells):
        self.origin_cells = self.num_cells/2
        self.origin_inds=self.array_sizes*(self.origin_cells-1)+1;
        self.coord_limits = [-self.origin_cells*self.lattice_const,\
                             (self.num_cells-self.origin_cells+1)*self.lattice_const]

    def x(self):
        return np.linspace(self.coord_limits[0],self.coord_limits[1],self.coord_sizes, endpoint = True)
