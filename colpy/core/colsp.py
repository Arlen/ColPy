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
    """
    The Abstract Base Class which all the other colour classes inherit.
    """
    
    @abc.abstractmethod
    def __init__(self, *tri):
        """
        Initializes the tristimulus values of a colour.

        :param tri: tristimulus
        :type tri: floats
        """
        self.t1, self.t2, self.t3 = tri


    def __repr__(self):
        return "{0}{1}".format(self.__class__.__name__, self.tri)

    
    def __str__(self):
        return "{0}{1}".format(self.__class__.__name__.strip("Colour_"),
                               self.tri)    


    @property
    def t1(self):
        return self._t1

    
    @t1.setter
    def t1(self, t1):
        self._t1 = t1

        
    @property
    def t2(self):
        return self._t2

    
    @t2.setter
    def t2(self, t2):
        self._t2 = t2

        
    @property
    def t3(self):
        return self._t3

    
    @t3.setter
    def t3(self, t3):
        self._t3 = t3

        
    @property
    def tri(self):
        """
        Returns 3-tuple of floats representing a position in the colour space.
        """
        return (self.t1, self.t2, self.t3)

    
    @tri.setter
    def tri(self, tri):
        """
        Sets the tristimulus values.

        :param tri: tristimulus
        :type tri: floats
        """
        self.t1, self.t2, self.t3 = tri


    @staticmethod
    def supported_colour_spaces():
        """
        A generator that returns all the supported colour spaces in ColPy.

        :return: name of each supported colour space.
        :rtype: string
        """
        for cs in _all_colour_spaces:
            yield cs

            
    
class Colour_XYZ(ColourSpace):
    
    def __init__(self, X=1.0, Y=1.0, Z=1.0):
        """
        Initializes the tristimulus values of a colour in XYZ colour space.

        :param X: The X coordinate
        :type X: float
        :param Y: The Y coordinate
        :type Y: float
        :param Z: The Z coordinate
        :type Z: float
        """
        super().__init__(X, Y, Z)


    @property
    def XYZ(self):
        """
        Returns 3-tuple of floats representing a position in XYZ colour space.
        """
        return self.tri

    
    @XYZ.setter
    def XYZ(self, XYZ):
        """
        Sets the tristimulus values of a colour in XYZ colour space.

        :param XYZ: tristimulus values
        :type XYZ: floats
        """       
        self.tri = XYZ
        
    
    @property
    def X(self):
        """
        Returns the `X` coordinate.

        :rtype: float
        """
        return self.t1


    @X.setter
    def X(self, X):
        """
        Sets the `X` coordinate.

        :param X: The X coordinate
        :type X: float
        """
        self.t1 = X


    @property
    def Y(self):
        """
        Returns the `Y` coordinate.

        :rtype: float
        """
        return self.t2


    @Y.setter
    def Y(self, Y):
        """
        Sets the `Y` coordinate.

        :param Y: The Y coordinate
        :type Y: float
        """
        self.t2 = Y


    @property
    def Z(self):
        """
        Returns the `Z` coordinate.

        :rtype: float
        """
        return self.t3


    @Z.setter
    def Z(self, Z):
        """
        Sets the `Z` coordinate.

        :param Z: The Z coordinate
        :type Z: float
        """
        self.t3 = Z

        
    def from_xyY(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `xyY` colour space to `XYZ` colour space.

        :param col: colour
        :type col: :class:`Colour_xyY`
        """
        self.X = col.x / col.y * col.Y
        self.Y = col.Y
        self.Z = (1.0 - col.x - col.y) / col.y * col.Y

        
    def from_Lab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Lab` colour space to `XYZ` colour space.

        :param col: colour
        :type col: :class:`Colour_Lab`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """
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

        self.X = xr * rw.X
        self.Y = yr * rw.Y
        self.Z = zr * rw.Z


    def from_LCHab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHab` colour space to `XYZ` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHab`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """
        self.from_Lab(col.to_Lab(), rw)


    def from_Luv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `XYZ` colour space.

        :param col: colour
        :type col: :class:`Colour_Luv`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """
        (uo, vo) = _compute_uo_vo(rw)
        c = -1/3
        a = (((52 * col.L) / (col.u + 13 * col.L * uo)) - 1) / 3
        if col.L > cie_ke():
            self.Y = math.pow( (col.tri[0] + 16.0) / 116.0, 3.0 )
        else:
            self.Y = col.L / cie_kappa()

        b = -5.0 * Y
        d = (((39 * col.L) / (col.v + 13 * col.L * vo)) - 5) * Y
        self.X = (d - b) / (a - c)
        self.Z = X * a + b


    def from_LCHuv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHuv` colour space to `XYZ` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHuv`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """
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
        """
        Initializes the tristimulus values of a colour in xyY colour space.

        :param x: The x coordinate
        :type x: float
        :param y: The y coordinate
        :type y: float
        :param Y: The Y coordinate
        :type Y: float
        """
        super().__init__(x, y, Y)


    @property
    def xyY(self):
        """
        Returns 3-tuple of floats representing a position in xyY colour space.
        """        
        return self.tri

    
    @xyY.setter
    def xyY(self, xyY):
        """
        Sets the tristimulus values of a colour in xyY colour space.

        :param tri: tristimulus
        :type tri: floats
        """               
        self.tri = xyY
        
        
    @property
    def x(self):
        """
        Returns the `x` coordinate.

        :rtype: float
        """        
        return self.t1

    
    @x.setter
    def x(self, x):
        """
        Sets the `x` coordinate.

        :param x: The x coordinate
        :type x: float
        """        
        self.t1 = x


    @property
    def y(self):
        """
        Returns the `y` coordinate.

        :rtype: float
        """        
        return self.t2

    
    @y.setter
    def y(self, y):
        """
        Sets the `y` coordinate.

        :param y: The y coordinate
        :type y: float
        """         
        self.t2 = y


    @property
    def Y(self):
        """
        Returns the `Y` coordinate.

        :rtype: float
        """        
        return self.t3

    
    @Y.setter
    def Y(self, Y):
        """
        Sets the `Y` coordinate.

        :param Y: The Y coordinate
        :type Y: float
        """        
        self.t3 = Y

    
    def from_XYZ(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        """        
        s = sum(col.XYZ)
        self.x = col.X / s
        self.y = col.Y / s
        self.Y = col.Y

        
    def from_Lab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Lab` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_Lab`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_XYZ(col.to_XYZ(rw))


    def from_LCHab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHab` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHab`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """       
        self.from_XYZ(col.to_XYZ(rw))


    def from_Luv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_Luv`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_XYZ(col.to_XYZ(rw))


    def from_LCHuv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHuv` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHuv`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
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
        """
        Initializes the tristimulus values of a colour in Lab colour space.

        :param L: The L coordinate
        :type L: float
        :param a: The a coordinate
        :type a: float
        :param b: The b coordinate
        :type b: float        
        """
        super().__init__(L, a, b)


    @property
    def Lab(self):
        """
        Returns 3-tuple of floats representing a position in Lab colour space.
        """        
        return self.tri

    
    @Lab.setter
    def Lab(self, Lab):
        """
        Sets the tristimulus values of a colour in Lab colour space.

        :param Lab: tristimulus values
        :type Lab: floats
        """       
        self.tri = Lab
        
        
    @property
    def L(self):
        """
        Returns the `L` coordinate.

        :rtype: float
        """
        return self.t1

    
    @L.setter
    def L(self, L):
        """
        Sets the `L` coordinate.

        :param L: The L coordinate
        :type L: float
        """
        self.t1 = L

    
    @property
    def a(self):
        """
        Returns the `a` coordinate.

        :rtype: float
        """       
        return self.t2

    
    @a.setter
    def a(self, a):
        """
        Sets the `a` coordinate.

        :param a: The a coordinate
        :type a: float
        """
        self.t2 = a
        

    @property
    def b(self):
        """
        Returns the `b` coordinate.

        :rtype: float
        """        
        return self.t3

    
    @b.setter
    def b(self, b):
        """
        Sets the `b` coordinate.

        :param b: The b coordinate
        :type b: float
        """
        self.t3 = b


    def from_XYZ(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `Lab` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """         
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

        self.L = 116.0 * fy - 16.0
        self.a = 500.0 * (fx - fy)
        self.b = 200.0 * (fy - fz)


    def from_xyY(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `xyY` colour space to `Lab` colour space.

        :param col: colour
        :type col: :class:`Colour_xyY`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """          
        self.from_XYZ(col.to_XYZ(), rw)


    def from_LCHab(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHab` colour space to `Lab` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHab`
        """        
        self.L = col.L
        self.a = col.C * math.cos(math.radians(col.H))
        self.b = col.C * math.sin(math.radians(col.H))


    def from_Luv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `Lab` colour space.

        :param col: colour
        :type col: :class:`Colour_Luv`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_XYZ(col.to_XYZ(rw), rw)


    def from_LCHuv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHuv` colour space to `Lab` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHuv`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
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
        """
        Initializes the tristimulus values of a colour in LCHab colour space.

        :param L: The L coordinate
        :type L: float
        :param C: The C coordinate
        :type C: float
        :param H: The H coordinate
        :type H: float        
        """
        super().__init__(L, C, H)

        
    @property
    def LCH(self):
        """
        Returns 3-tuple of floats representing a position in LCHab colour space.
        """        
        return self.tri

    
    @LCH.setter
    def LCH(self, LCH):
        """
        Sets the tristimulus values of a colour in LCHab colour space.

        :param LCH: tristimulus values
        :type LCH: floats
        """       
        self.tri = LCH
        
        
    @property
    def L(self):
        """
        Returns the `L` coordinate.

        :rtype: float
        """        
        return self.t1

    
    @L.setter
    def L(self, L):
        """
        Sets the `L` coordinate.

        :param L: The L coordinate
        :type L: float
        """
        self.t1 = L
        
    
    @property
    def C(self):
        """
        Returns the `C` coordinate.

        :rtype: float
        """        
        return self.t2

    
    @C.setter
    def C(self, C):
        """
        Sets the `C` coordinate.

        :param C: The C coordinate
        :type C: float
        """
        self.t2 = C


    @property
    def H(self):
        """
        Returns the `H` coordinate.

        :rtype: float
        """       
        return self.t3

    
    @H.setter
    def H(self, H):
        """
        Sets the `H` coordinate.

        :param H: The H coordinate
        :type H: float
        """
        self.t3 = H

    
    def from_XYZ(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `LCHab` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_Lab(col.to_Lab(rw))


    def from_xyY(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `xyY` colour space to `LCHab` colour space.

        :param col: colour
        :type col: :class:`Colour_xyY`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_Lab(col.to_Lab(rw))


    def from_Lab(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Lab` colour space to `LCHab` colour space.

        :param col: colour
        :type col: :class:`Colour_Lab`
        """        
        self.L = col.L
        self.C = math.sqrt(col.a * col.a + col.b * col.b)
        h = math.angles(math.atan2(col.b, col.a))
        if h < 0.0:
            self.H = h + 360
        elif h > 360.0:
            self.H = h - 360
        else:
            self.H = h


    def from_Luv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `LCHab` colour space.

        :param col: colour
        :type col: :class:`Colour_Luv`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_Lab(col.to_Lab(rw))


    def from_LCHuv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHuv` colour space to `LCHab` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHuv`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
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
        """
        Initializes the tristimulus values of a colour in Luv colour space.

        :param L: The L coordinate
        :type L: float
        :param u: The u coordinate
        :type u: float
        :param v: The v coordinate
        :type v: float        
        """
        super().__init__(L, u, v)

        
    @property
    def Luv(self):
        """
        Returns 3-tuple of floats representing a position in Luv colour space.
        """        
        return self.tri

    
    @Luv.setter
    def Luv(self, Luv):
        """
        Sets the tristimulus values of a colour in Luv colour space.

        :param Luv: tristimulus values
        :type Luv: floats
        """       
        self.tri = Luv
        
        
    @property
    def L(self):
        """
        Returns the `L` coordinate.

        :rtype: float
        """
        return self.t1

    
    @L.setter
    def L(self, L):
        """
        Sets the `L` coordinate.

        :param L: The L coordinate
        :type L: float
        """
        self.t1 = L

    
    @property
    def u(self):
        """
        Returns the `u` coordinate.

        :rtype: float
        """       
        return self.t2

    
    @u.setter
    def u(self, u):
        """
        Sets the `u` coordinate.

        :param u: The u coordinate
        :type u: float
        """
        self.t2 = u
        

    @property
    def v(self):
        """
        Returns the `v` coordinate.

        :rtype: float
        """        
        return self.t3

    
    @v.setter
    def v(self, v):
        """
        Sets the `v` coordinate.

        :param v: The v coordinate
        :type v: float
        """
        self.t3 = v

    
    def from_XYZ(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `Luv` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        yr = col.Y / rw.Y
        (up, vp) = _compute_uo_vo(col)
        (urp, vrp) = _compute_uo_vo(rw)
        if yr > cie_epsilon():
            self.L = 116 * math.pow(yr, 1/3) - 16
        else:
            self.L = cie_kappa() * yr
            
        self.u = 13 * L * (up - urp)
        self.v = 13 * L * (vp - vrp)


    def from_xyY(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `xyY` colour space to `Luv` colour space.

        :param col: colour
        :type col: :class:`Colour_xyY`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_XYZ(col.to_XYZ(), rw)


    def from_Lab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Lab` colour space to `Luv` colour space.

        :param col: colour
        :type col: :class:`Colour_Lab`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_XYZ(col.to_XYZ(rw), rw)
        
        
    def from_LCHab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHab` colour space to `Luv` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHab`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_XYZ(col.to_XYZ(rw), rw)


    def from_LCHuv(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHuv` colour space to `Luv` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHuv`
        """         
        self.L = col.L
        self.u = col.C * math.cos(math.radians(col.H))
        self.v = col.C * math.sin(math.radians(col.H))
        
        
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
        """
        Initializes the tristimulus values of a colour in LCHuv colour space.

        :param L: The L coordinate
        :type L: float
        :param C: The C coordinate
        :type C: float
        :param H: The H coordinate
        :type H: float        
        """
        super().__init__(L, C, H)

        
    @property
    def LCH(self):
        """
        Returns 3-tuple of floats representing a position in LCHuv colour space.
        """        
        return self.tri

    
    @LCH.setter
    def LCH(self, LCH):
        """
        Sets the tristimulus values of a colour in LCHuv colour space.

        :param LCH: tristimulus values
        :type LCH: floats
        """       
        self.tri = LCH
        
        
    @property
    def L(self):
        """
        Returns the `L` coordinate.

        :rtype: float
        """        
        return self.t1

    
    @L.setter
    def L(self, L):
        """
        Sets the `L` coordinate.

        :param L: The L coordinate
        :type L: float
        """
        self.t1 = L
        
    
    @property
    def C(self):
        """
        Returns the `C` coordinate.

        :rtype: float
        """        
        return self.t2

    
    @C.setter
    def C(self, C):
        """
        Sets the `C` coordinate.

        :param C: The C coordinate
        :type C: float
        """
        self.t2 = C


    @property
    def H(self):
        """
        Returns the `H` coordinate.

        :rtype: float
        """       
        return self.t3

    
    @H.setter
    def H(self, H):
        """
        Sets the `H` coordinate.

        :param H: The H coordinate
        :type H: float
        """
        self.t3 = H
        
    
    def from_XYZ(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_Luv(col.to_Luv(rw))


    def from_xyY(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `xyY` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_xyY`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_Luv(col.to_Luv(rw))


    def from_Lab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Lab` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_Lab`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`        
        """         
        self.from_Luv(col.to_Luv(rw))

        
    def from_LCHab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHab`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.from_Luv(col.to_Luv(rw), rw)


    def from_Luv(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_Luv`
        :param rw: reference white
        :type rw: :class:`Colour_XYZ`
        """        
        self.L = col.L
        self.C = math.sqrt(col.u * col.u + col.v * col.v)
        h = math.angles(math.ata2(col.v, col.u))
        if h < 0.0:
            self.H = h + 360
        elif h >= 360:
            self.H = h - 360
        else:
            self.H = h
        

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

