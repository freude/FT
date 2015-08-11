#!/usr/bin/python

import numpy as np
import myconst as MyConst


class CoordSys(object):

    def __init__(self, num_cell, array_sizes, units='au'):

        if units == 'au':
            self.lattice_const=0.5431*1e-9/MyConst.ab;
        else
            self.lattice_const=0.5431*1e-9;

        self.num_cells = num_cells;
        self.arrays_sizes = arrays_sizes;
        self.coord_sizes = obj.num_cells*obj.arrays_sizes;
        dist=self.num_cells*obj.lattice_const;
        self.coord_stps=dist/self.coord_sizes;
        self.origin_inds=self.arrays_sizes*(self.origin_cells-1)+1;
        self.coord_limits(1)=-self.origin_inds(1)*self.coord_stps+self.coord_stps;
        self.coord_limits(2)=-self.origin_inds(1)*self.coord_stps+self.num_cells*self.lattice_const;

    #----------------------------------------------------------------------

