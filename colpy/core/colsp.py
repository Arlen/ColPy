#!/usr/bin/env python3
################################################################################
#
#  Copyright (C) 2011 Arlen Avakian
#
#  This file is part of ColPy.
#
#  ColPy is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ColPy is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ColPy.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


import abc
import math


_all_colour_spaces = ['XYZ', 'xyY', 'Lab', 'LCHab', 'Luv', 'LCHuv']

_CIE_EPSILON = 216 / 24389
_CIE_KAPPA = 24389 / 27


def cie_epsilon():
    return _CIE_EPSILON


def cie_kappa():
    return _CIE_KAPPA


def cie_ke():
    return cie_kappa() * cie_epsilon()


def _compute_uo_vo(rw):
    d = rw.X + 15.0 * rw.Y + 3.0 * rw.Z
    return ((4.0 * rw.X) / d, (9.0 * rw.Y) / d)


class ColourSpace(metaclass=abc.ABCMeta):
    
    @abc.abstractmethod
    def __init__(self, tri):
        self.__tri = tuple(float(n) for n in tri[0:3])


    def __repr__(self):
        return "{0}{1}".format(self.__class__.__name__, self.tri)

    
    def __str__(self):
        return "{0}{1}".format(self.__class__.__name__.strip("Colour_"),
                               self.tri)    
    
    
    @property
    def tri(self):
        """
        A 3-tuple representing a position in the colour space.
        """
        return self.__tri

    
    @tri.setter
    def tri(self, tri):
        assert len(tri) is 3, "tri must contain 3 elements"
        self.__tri = tuple(float(n) for n in tri[0:3])


    @staticmethod
    def supported_colour_spaces():
        for cs in _all_colour_spaces:
            yield cs

            
    
class Colour_XYZ(ColourSpace):
    
    def __init__(self, X=1.0, Y=1.0, Z=1.0):
        super().__init__((X, Y, Z))

        
    @property
    def X(self):
        return self.tri[0]


    @property
    def Y(self):
        return self.tri[1]


    @property
    def Z(self):
        return self.tri[2]

        
    def from_xyY(self, col):
        X = col.x / col.y * col.Y
        Y = col.Y
        Z = (1.0 - col.x - col.y) / col.y * col.Y
        self.tri = (X, Y, Z)

        
    def from_Lab(self, col, rw):
        fy = (col.L + 16.0) / 116.0
        fx = (col.a / 500.0) + fy
        fz = fy - (col.b / 200.0)
        fxCube = fx * fx * fx
        fzCube = fz * fz * fz

        if fxCube > cie_epsilon():
            xr = fxCube
        else:
            xr = ((116.0 * fx) - 16.0) / cie_kappa()

        if col.L > cie_ke():
            yr = math.pow( (col.L + 16.0) / 116.0, 3.0 )
        else:
            yr = col.L / cie_kappa()

        if fzCube > cie_epsilon():
            zr = fzCube
        else:
            zr = ((116.0 * fz) - 16.0) / cie_kappa()

        X = xr * rw.X
        Y = yr * rw.Y
        Z = zr * rw.Z
        self.tri = (X, Y, Z)


    def from_LCHab(self, col, rw):
        self.from_Lab(col.to_Lab(), rw)


    def from_Luv(self, col, rw):
        (uo, vo) = _compute_uo_vo(rw)
        c = -1/3
        a = (((52 * col.L) / (col.u + 13 * col.L * uo)) - 1) / 3
        if col.L > cie_ke():
            Y = math.pow( (col.tri[0] + 16.0) / 116.0, 3.0 )
        else:
            Y = col.L / cie_kappa()

        b = -5.0 * Y
        d = (((39 * col.L) / (col.v + 13 * col.L * vo)) - 5) * Y
        X = (d - b) / (a - c)
        Z = X * a + b
        self.tri = (X, Y, Z)


    def from_LCHuv(self, col, rw):
        self.from_Luv(col.to_Luv(), rw)

    
    for cs in ColourSpace.supported_colour_spaces():
        if cs is not 'XYZ':
            if cs in ('xyY'):
                exec("def to_{0}(self):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self)\n"
                     "  return rt".format(cs, 'XYZ'))
            else:
                exec("def to_{0}(self, rw):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self, rw)\n"
                     "  return rt".format(cs, 'XYZ'))
    del cs
                 


class Colour_xyY(ColourSpace):

    def __init__(self, x=1/3, y=1/3, Y=1.0):
        super().__init__((x, y, Y))


    @property
    def x(self):
        return self.tri[0]


    @property
    def y(self):
        return self.tri[1]


    @property
    def Y(self):
        return self.tri[2]

    
    def from_XYZ(self, col):
        s = sum(col.tri)
        x = col.X / s
        y = col.Y / s
        Y = col.Y
        self.tri = (x, y, Y)

        
    def from_Lab(self, col, rw):
        self.from_XYZ(col.to_XYZ(rw))


    def from_LCHab(self, col, rw):
        self.from_XYZ(col.to_XYZ(rw))


    def from_Luv(self, col, rw):
        self.from_XYZ(col.to_XYZ(rw))


    def from_LCHuv(self, col, rw):
        self.from_XYZ(col.to_XYZ(rw))
        

    for cs in ColourSpace.supported_colour_spaces():
        if cs is not 'xyY':
            if cs in ('XYZ'):
                exec("def to_{0}(self):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self)\n"
                     "  return rt".format(cs, 'xyY'))
            else:
                exec("def to_{0}(self, rw):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self, rw)\n"
                     "  return rt".format(cs, 'xyY'))
    del cs


class Colour_Lab(ColourSpace):

    def __init__(self, L=100.0, a=0.0, b=0.0):
        super().__init__((L, a, b))


    @property
    def L(self):
        return self.tri[0]

    
    @property
    def a(self):
        return self.tri[1]


    @property
    def b(self):
        return self.tri[2]


    def from_XYZ(self, col, rw):
        xr = col.X / rw.X
        yr = col.Y / rw.Y
        zr = col.Z / rw.Z

        if xr > cie_epsilon():
            fx = math.pow(xr, 1/3)
        else:
            fx = ((cie_kappa() * xr) + 16.0) / 116.0

        if yr > cie_epsilon():
            fy = math.pow(yr, 1/3)
        else:
            fy = ((cie_kappa() * yr) + 16.0) / 116.0

        if zr > cie_epsilon():
            fz = math.pow(zr, 1/3)
        else:
            fz = ((cie_kappa() * zr) + 16.0) / 116.0

        L = 116.0 * fy - 16.0
        a = 500.0 * (fx - fy)
        b = 200.0 * (fy - fz)
        self.tri = (L, a, b)


    def from_xyY(self, col, rw):
        self.from_XYZ(col.to_XYZ(), rw)


    def from_LCHab(self, col):
        L = col.L
        a = col.C * math.cos(math.radians(col.H))
        b = col.C * math.sin(math.radians(col.H))
        self.tri = (L, a, b)


    def from_Luv(self, col, rw):
        self.from_XYZ(col.to_XYZ(rw), rw)


    def from_LCHuv(self, col, rw):
        self.from_XYZ(col.to_XYZ(rw), rw)
        
        
    for cs in ColourSpace.supported_colour_spaces():
        if cs is not 'Lab':
            if cs in ('LCHab'):
                exec("def to_{0}(self):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self)\n"
                     "  return rt".format(cs, 'Lab'))
            else:
                exec("def to_{0}(self, rw):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self, rw)\n"
                     "  return rt".format(cs, 'Lab'))
    del cs
        

class Colour_LCHab(ColourSpace):

    def __init__(self, L=100.0, C=0.0, H=0.0):
        super().__init__((L, C, H))


    @property
    def L(self):
        return self.tri[0]

    
    @property
    def C(self):
        return self.tri[1]


    @property
    def H(self):
        return self.tri[2]

    
    def from_XYZ(self, col, rw):
        self.from_Lab(col.to_Lab(rw))


    def from_xyY(self, col, rw):
        self.from_Lab(col.to_Lab(rw))


    def from_Lab(self, col):
        L = col.L
        C = math.sqrt(col.a * col.a + col.b * col.b)
        h = math.angles(math.atan2(col.tri[2], col.tri[1]))
        if h < 0.0:
            H = h + 360
        elif h > 360.0:
            H = h - 360
        else:
            H = h
        self.tri = (L, C, H)


    def from_Luv(self, col, rw):
        self.from_Lab(col.to_Lab(rw))


    def from_LCHuv(self, col, rw):
        self.from_Lab(col.to_Lab(rw))
        

    for cs in ColourSpace.supported_colour_spaces():
        if cs is not 'LCHab':
            if cs in ('Lab'):
                exec("def to_{0}(self):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self)\n"
                     "  return rt".format(cs, 'LCHab'))
            else:
                exec("def to_{0}(self, rw):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self, rw)\n"
                     "  return rt".format(cs, 'LCHab'))
    del cs
        

class Colour_Luv(ColourSpace):

    def __init__(self, L=100.0, u=0.0, v=0.0):
        super().__init__((L, u, v))

        
    @property
    def L(self):
        return self.tri[0]

    
    @property
    def u(self):
        return self.tri[1]


    @property
    def v(self):
        return self.tri[2]

    
    def from_XYZ(self, col, rw):
        yr = col.Y / rw.Y
        (up, vp) = _compute_uo_vo(col)
        (urp, vrp) = _compute_uo_vo(rw)
        if yr > cie_epsilon():
            L = 116 * math.pow(yr, 1/3) - 16
        else:
            L = cie_kappa() * yr
            
        u = 13 * L * (up - urp)
        v = 13 * L * (vp - vrp)
        self.tri = (L, u, v)


    def from_xyY(self, col, rw):
        self.from_XYZ(col.to_XYZ(), rw)


    def from_Lab(self, col, rw):
        self.from_XYZ(col.to_XYZ(rw), rw)
        
        
    def from_LCHab(self, col, rw):
        self.from_XYZ(col.to_XYZ(rw), rw)


    def from_LCHuv(self, col):
        L = col.L
        u = col.C * math.cos(math.radians(col.H))
        v = col.C * math.sin(math.radians(col.H))
        self.tri = (L, u, v)
        
        
    for cs in ColourSpace.supported_colour_spaces():
        if cs is not 'Luv':
            if cs in ('LCHuv'):
                exec("def to_{0}(self):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self)\n"
                     "  return rt".format(cs, 'Luv'))
            else:
                exec("def to_{0}(self, rw):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self, rw)\n"
                     "  return rt".format(cs, 'Luv'))
    del cs




class Colour_LCHuv(ColourSpace):

    def __init__(self, L=100.0, C=0.0, H=0.0):
        super().__init__((L, C, H))

        
    @property
    def L(self):
        return self.tri[0]

    
    @property
    def C(self):
        return self.tri[1]


    @property
    def H(self):
        return self.tri[2]

    
    def from_XYZ(self, col, rw):
        self.from_Luv(col.to_Luv(rw))


    def from_xyY(self, col, rw):
        self.from_Luv(col.to_Luv(rw))


    def from_Lab(self, col, rw):
        self.from_Luv(col.to_Luv(rw))

        
    def from_LCHab(self, col, rw):
        self.from_Luv(col.to_Luv(rw), rw)


    def from_Luv(self, col):
        L = col.L
        C = math.sqrt(col.u * col.u + col.v * col.v)
        h = math.angles(math.ata2(col.v, col.u))
        if h < 0.0:
            H = h + 360
        elif h >= 360:
            H = h - 360
        else:
            H = h
        self.tri = (L, C, H)
        

    for cs in ColourSpace.supported_colour_spaces():
        if cs is not 'LCHuv':
            if cs in ('Luv'):
                exec("def to_{0}(self):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self)\n"
                     "  return rt".format(cs, 'LCHuv'))
            else:
                exec("def to_{0}(self, rw):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self, rw)\n"
                     "  return rt".format(cs, 'LCHuv'))
    del cs

                

def colour(t1, t2, t3, cs):
    if cs not in ColourSpace.supported_colour_spaces():
        raise ValueError("unsupported colour space")
    return eval("Colour_{0}(t1, t2, t3)".format(cs))

