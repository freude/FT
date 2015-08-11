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

    def set_origin_cells(self, origin_cells):
            self.origin_cells = origin_cells;
            self.origin_inds=self.arrays_sizes*(self.origin_cells-1)+1;
            self.coord_limits(1)=-self.origin_inds(1)*self.coord_stps+self.coord_stps;
            self.coord_limits(2)=-self.origin_inds(1)*self.coord_stps+self.num_cells*self.lattice_const;
        end

        [kx ,ky ,kz] = global_basis_rec(obj,j1,j2,j3)
        [x ,y ,z] = global_basis(obj,j1,j2,j3)

        function [x ,y ,z] = global_basis_array(obj,X,Y,Z)
            [x ,y ,z] = arrayfun(@obj.global_basis,X,Y,Z);
        end;

        function [x ,y ,z] = global_basis_rec_array(obj,X,Y,Z)
            [x ,y ,z] = arrayfun(@obj.global_basis_rec,X,Y,Z);
        end;

        function [x ,y ,z] = global_coords_gen(obj)
            x=1:obj.coord_sizes;
            [X,Y,Z]=meshgrid(x,x,x);
            [x ,y ,z] = arrayfun(@obj.global_basis,X,Y,Z);
        end;

        function [x ,y ,z] = global_coords_rec_gen(obj)
            x=1:obj.coord_sizes;
            [X,Y,Z]=meshgrid(x,x,x);
            [x ,y ,z] = arrayfun(@obj.global_basis_rec,X,Y,Z);
        end;

        function x = x(obj)
            x = obj.coord_limits(1):obj.coord_stps:obj.coord_limits(2);
        end;
