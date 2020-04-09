# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 16:10:15 2020

@author: ghaza
"""

from math import pi
import pygmsh
import meshio

def test():
    geom = pygmsh.opencascade.Geometry()
    geom.add_cylinder([0.0,0.0,0.0], [0.0,0.0,1.0], 20, 0.25*pi,char_length=1)
  #  ref=0.097625512963
    mesh = pygmsh.generate_mesh(geom)
    #meshio.write("final.vtk", mesh)
    return mesh
#meshio.write("opencascade_cylinder.vtk",test())
    


