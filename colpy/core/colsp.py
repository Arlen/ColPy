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
import numpy as np


_RGB = ['AdobeRGB', 'AppleRGB', 'BestRGB', 'BetaRGB','BruceRGB', 'CIERGB',
        'ColorMatchRGB', 'DonRGB4', 'ECIRGB', 'EktaSpacePS5', 'NTSCRGB',
        'PAL_SECAMRGB', 'ProPhotoRGB', 'SMPTE_CRGB', 'WideGamutRGB', 'sRGB']

_non_RGB = ['XYZ', 'xyY', 'Lab', 'LCHab', 'Luv', 'LCHuv']

_all_colour_spaces = _RGB + _non_RGB


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
        self._t1, self._t2, self._t3 = tri


    def __repr__(self):
        return "{0}{1}".format(self.__class__.__name__, self.tri)


    def __str__(self):
        return "{0}{1}".format(self.__class__.__name__.strip("Colour_"),
                               self.tri)


    @property
    def _t1(self):
        return self.__t1


    @_t1.setter
    def _t1(self, t1):
        self.__t1 = t1


    @property
    def _t2(self):
        return self.__t2


    @_t2.setter
    def _t2(self, t2):
        self.__t2 = t2


    @property
    def _t3(self):
        return self.__t3


    @_t3.setter
    def _t3(self, t3):
        self.__t3 = t3


    @property
    def tri(self):
        """
        Returns 3-tuple of floats representing the tristimulus values of
        a colour in the colour space.
        """
        return (self._t1, self._t2, self._t3)


    @tri.setter
    def tri(self, tri):
        """
        Sets the tristimulus values.

        :param tri: tristimulus
        :type tri: floats
        """
        self._t1, self._t2, self._t3 = tri


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
    def X(self):
        """
        Returns the `X` coordinate.

        :rtype: float
        """
        return self._t1


    @X.setter
    def X(self, X):
        """
        Sets the `X` coordinate.

        :param X: The X coordinate
        :type X: float
        """
        self._t1 = X


    @property
    def Y(self):
        """
        Returns the `Y` coordinate.

        :rtype: float
        """
        return self._t2


    @Y.setter
    def Y(self, Y):
        """
        Sets the `Y` coordinate.

        :param Y: The Y coordinate
        :type Y: float
        """
        self._t2 = Y


    @property
    def Z(self):
        """
        Returns the `Z` coordinate.

        :rtype: float
        """
        return self._t3


    @Z.setter
    def Z(self, Z):
        """
        Sets the `Z` coordinate.

        :param Z: The Z coordinate
        :type Z: float
        """
        self._t3 = Z


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
        :type rw: :class:`BaseIlluminant`
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
        :type rw: :class:`BaseIlluminant`
        """
        self.from_Lab(col.to_Lab(), rw)


    def from_Luv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `XYZ` colour space.

        :param col: colour
        :type col: :class:`Colour_Luv`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        (uo, vo) = _compute_uo_vo(rw)
        c = -1/3
        a = (((52 * col.L) / (col.u + 13 * col.L * uo)) - 1) / 3
        if col.L > cie_ke():
            self.Y = math.pow( (col.L + 16.0) / 116.0, 3.0 )
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
        :type rw: :class:`BaseIlluminant`
        """
        self.from_Luv(col.to_Luv(), rw)


    def from_RGB(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in any `RGB` colour space to `XYZ` colour space.

        :param col: colour
        :type col: :class:`Colour_RGB`
        """
        t = col.inverse_companding(col.gamma, col.r, col.g, col.b)
        self.tri = (col.m_adapted * t).flat      


    for cs in ColourSpace.supported_colour_spaces():
        reference_white_not_required = ['xyY'] + _RGB
        if cs is not 'XYZ':
            if cs in reference_white_not_required:
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
    del reference_white_not_required


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
    def x(self):
        """
        Returns the `x` coordinate.

        :rtype: float
        """
        return self._t1


    @x.setter
    def x(self, x):
        """
        Sets the `x` coordinate.

        :param x: The x coordinate
        :type x: float
        """
        self._t1 = x


    @property
    def y(self):
        """
        Returns the `y` coordinate.

        :rtype: float
        """
        return self._t2


    @y.setter
    def y(self, y):
        """
        Sets the `y` coordinate.

        :param y: The y coordinate
        :type y: float
        """
        self._t2 = y


    @property
    def Y(self):
        """
        Returns the `Y` coordinate.

        :rtype: float
        """
        return self._t3


    @Y.setter
    def Y(self, Y):
        """
        Sets the `Y` coordinate.

        :param Y: The Y coordinate
        :type Y: float
        """
        self._t3 = Y


    def from_XYZ(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        """
        s = sum(col.tri)
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
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw))


    def from_LCHab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHab` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHab`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw))


    def from_Luv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_Luv`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw))


    def from_LCHuv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHuv` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHuv`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw))

        
    def from_RGB(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in any `RGB` colour space to `xyY` colour space.

        :param col: colour
        :type col: :class:`Colour_RGB`
        """
        self.from_XYZ(col.to_XYZ())


    for cs in ColourSpace.supported_colour_spaces():
        reference_white_not_required = ['XYZ'] + _RGB
        if cs is not 'xyY':
            if cs in reference_white_not_required:
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
    del reference_white_not_required


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
    def L(self):
        """
        Returns the `L` coordinate.

        :rtype: float
        """
        return self._t1


    @L.setter
    def L(self, L):
        """
        Sets the `L` coordinate.

        :param L: The L coordinate
        :type L: float
        """
        self._t1 = L


    @property
    def a(self):
        """
        Returns the `a` coordinate.

        :rtype: float
        """
        return self._t2


    @a.setter
    def a(self, a):
        """
        Sets the `a` coordinate.

        :param a: The a coordinate
        :type a: float
        """
        self._t2 = a


    @property
    def b(self):
        """
        Returns the `b` coordinate.

        :rtype: float
        """
        return self._t3


    @b.setter
    def b(self, b):
        """
        Sets the `b` coordinate.

        :param b: The b coordinate
        :type b: float
        """
        self._t3 = b


    def from_XYZ(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `Lab` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
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
        :type rw: :class:`BaseIlluminant`
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
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw), rw)


    def from_LCHuv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHuv` colour space to `Lab` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHuv`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw), rw)


    def from_RGB(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in any `RGB` colour space to `Lab` colour space.

        :param col: colour
        :type col: :class:`Colour_RGB`
        """
        self.from_XYZ(col.to_XYZ(), rw)


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
    def L(self):
        """
        Returns the `L` coordinate.

        :rtype: float
        """
        return self._t1


    @L.setter
    def L(self, L):
        """
        Sets the `L` coordinate.

        :param L: The L coordinate
        :type L: float
        """
        self._t1 = L


    @property
    def C(self):
        """
        Returns the `C` coordinate.

        :rtype: float
        """
        return self._t2


    @C.setter
    def C(self, C):
        """
        Sets the `C` coordinate.

        :param C: The C coordinate
        :type C: float
        """
        self._t2 = C


    @property
    def H(self):
        """
        Returns the `H` coordinate.

        :rtype: float
        """
        return self._t3


    @H.setter
    def H(self, H):
        """
        Sets the `H` coordinate.

        :param H: The H coordinate
        :type H: float
        """
        self._t3 = H


    def from_XYZ(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `LCHab` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_Lab(col.to_Lab(rw))


    def from_xyY(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `xyY` colour space to `LCHab` colour space.

        :param col: colour
        :type col: :class:`Colour_xyY`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
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
        :type rw: :class:`BaseIlluminant`
        """
        self.from_Lab(col.to_Lab(rw))


    def from_LCHuv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHuv` colour space to `LCHab` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHuv`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_Lab(col.to_Lab(rw))


    def from_RGB(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in any `RGB` colour space to `LCHab` colour
        space.

        :param col: colour
        :type col: :class:`Colour_RGB`
        """
        self.from_XYZ(col.to_XYZ(), rw)


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
    def L(self):
        """
        Returns the `L` coordinate.

        :rtype: float
        """
        return self._t1


    @L.setter
    def L(self, L):
        """
        Sets the `L` coordinate.

        :param L: The L coordinate
        :type L: float
        """
        self._t1 = L


    @property
    def u(self):
        """
        Returns the `u` coordinate.

        :rtype: float
        """
        return self._t2


    @u.setter
    def u(self, u):
        """
        Sets the `u` coordinate.

        :param u: The u coordinate
        :type u: float
        """
        self._t2 = u


    @property
    def v(self):
        """
        Returns the `v` coordinate.

        :rtype: float
        """
        return self._t3


    @v.setter
    def v(self, v):
        """
        Sets the `v` coordinate.

        :param v: The v coordinate
        :type v: float
        """
        self._t3 = v


    def from_XYZ(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `Luv` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
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
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(), rw)


    def from_Lab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Lab` colour space to `Luv` colour space.

        :param col: colour
        :type col: :class:`Colour_Lab`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw), rw)


    def from_LCHab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `LCHab` colour space to `Luv` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHab`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
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


    def from_RGB(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in any `RGB` colour space to `Luv` colour space.

        :param col: colour
        :type col: :class:`Colour_RGB`
        """
        self.from_XYZ(col.to_XYZ(), rw)


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
    def L(self):
        """
        Returns the `L` coordinate.

        :rtype: float
        """
        return self._t1


    @L.setter
    def L(self, L):
        """
        Sets the `L` coordinate.

        :param L: The L coordinate
        :type L: float
        """
        self._t1 = L


    @property
    def C(self):
        """
        Returns the `C` coordinate.

        :rtype: float
        """
        return self._t2


    @C.setter
    def C(self, C):
        """
        Sets the `C` coordinate.

        :param C: The C coordinate
        :type C: float
        """
        self._t2 = C


    @property
    def H(self):
        """
        Returns the `H` coordinate.

        :rtype: float
        """
        return self._t3


    @H.setter
    def H(self, H):
        """
        Sets the `H` coordinate.

        :param H: The H coordinate
        :type H: float
        """
        self._t3 = H


    def from_XYZ(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `XYZ` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_Luv(col.to_Luv(rw))


    def from_xyY(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `xyY` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_xyY`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_Luv(col.to_Luv(rw))


    def from_Lab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Lab` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_Lab`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_Luv(col.to_Luv(rw))


    def from_LCHab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_LCHab`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_Luv(col.to_Luv(rw), rw)


    def from_Luv(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in `Luv` colour space to `LCHuv` colour space.

        :param col: colour
        :type col: :class:`Colour_Luv`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
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


    def from_RGB(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in any `RGB` colour space to `LCHuv` colour
        space.

        :param col: colour
        :type col: :class:`Colour_RGB`
        """
        self.from_XYZ(col.to_XYZ(), rw)


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


_xyz_scaling = np.matrix(np.identity(3))

_von_kries = np.matrix(' 0.40024 -0.2263  0.0;     \
                         0.7076   1.16532 0.0;     \
                        -0.08081  0.0457  0.91822')

_bradford = np.matrix(' 0.8951 -0.7502  0.0389; \
                        0.2664  1.7135 -0.0685; \
                       -0.1614  0.0367  1.0296' )

_adaptation_methods = {'xyz_scaling': _xyz_scaling,
                       'von_kries': _von_kries,
                       'bradford': _bradford}

class BaseRGB(ColourSpace):

    def __init__(self, gamma, rw, r_primary, g_primary, b_primary, r, g, b):
        """
        Computes the XYZ-to-RGB and RGB-to-XYZ matrices, and initializes the
        tristimulus values of a colour based on the RGB colour model.

        :param gamma: The gamma value defined in the specification
        :type gamma: float
        :param rw: The reference white chosen by the specification
        :type rw: :class:`BaseIlluminant`
        :param r_primary: The red primary 
        :type r_primary: :class:`Colour_xyY`
        :param g_primary: The green primary
        :type g_primary: :class:`Colour_xyY`
        :param b_primary: The blue primary
        :type b_primary: :class:`Colour_xyY`
        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        self.__gamma = gamma
        self.__rw = rw
        self.__m = __class__.__computeConversionMatrix(r_primary,
                                                       g_primary,
                                                       b_primary,
                                                       rw)
        self.__m_1 = self.__m.I
        self.adapt(self.__rw, 'bradford')
        super().__init__(r, g, b)


    @property
    def m(self):
        """
        Returns the RGB to XYZ matrix.

        :rtype: numpy.matrix
        """
        return self.__m


    @property
    def m_adapted(self):
        """
        Returns the chromatic adaptation matrix for RGB to XYZ conversion.

        :rtype: numpy.matrix
        """
        return self.__m_adapted


    @property
    def m_1(self):
        """
        Returns the XYZ to RGB matrix.

        :rtype: numpy.matrix
        """
        return self.__m_1


    @property
    def m_1_adapted(self):
        """
        Returns the chromatic adaptation matrix for XYZ to RGB conversion.

        :rtype: numpy.matrix
        """
        return self.__m_1_adapted


    @property
    def gamma(self):
        """
        Returns the `gamma`.

        :rtype: float
        """
        return self.__gamma


    @property
    def r(self):
        """
        Returns the `r` coordinate.
        """
        return self._t1


    @r.setter
    def r(self, r):
        """
        Sets the `r` coordinate.

        :param r: The r coordinate
        :rtype r: float
        """
        self._t1 = r


    @property
    def g(self):
        """
        Returns the `g` coordinate.
        """
        return self._t2


    @g.setter
    def g(self, g):
        """
        Sets the `g` coordinate.

        :param g: The g coordinate
        :rtype g: float
        """
        self._t2 = g


    @property
    def b(self):
        """
        Returns the `b` coordinate.
        """
        return self._t3


    @b.setter
    def b(self, b):
        """
        Sets the `b` coordinate.

        :param b: The b coordinate
        :rtype b: float
        """
        self._t3 = b


    def adapt(self, target_rw, am='bradford'):
        """
        If the user wants XYZ relative to a different reference white, then
        chromatic adaptation transform must be applied to the XYZ colour to
        convert it from the reference white of the RGB system to the desired
        reference white.  Given a target reference white, `target_rw`, this
        method computes the chromatic adaptation matrix that is needed to
        convert RGB to XYZ and XYZ to RGB.

        :param target_rw: The desired reference white
        :type target_rw: :class:`Colour_xyY`
        :param am: The name of adaptation method to be used; 'xyz_scaling',
        'von_kries', and 'bradford' methods are supported.  'bradford' is used
        by default.
        :type am: string
        """
        if am not in _adaptation_methods:
            raise ValueError("unsupported adaptation method")
        method = _adaptation_methods.get(am)
        self.__m_adapted = __class__. \
            __computeChromaticAdaptationMatrix(self.__rw,
                                               target_rw,
                                               method) * self.m
        self.__m_1_adapted = self.__m_adapted.I


    @staticmethod
    def __computeConversionMatrix(red, green, blue, rw):
        """
        Computes the conversion matrix used in linear-RGB to XYZ conversion.

        :param red: the primary red
        :type red: :class:`Colour_xyY`
        :param green: the primary green
        :type green: :class:`Colour_xyY`
        :param blue: the primary blue
        :type blue: :class:`Colour_xyY`
        :return: linear-RGB to XYZ conversion matrix
        :rtype: numpy.matrix
        """
        xyzs = np.matrix(np.zeros((3, 3)))
        xyzs[:,0] = np.asarray(red.to_XYZ().tri).reshape(3, 1)
        xyzs[:,1] = np.asarray(green.to_XYZ().tri).reshape(3, 1)
        xyzs[:,2] = np.asarray(blue.to_XYZ().tri).reshape(3, 1)
        S = xyzs.I * np.array([[rw.X], [rw.Y], [rw.Z]])
        M = np.matrix(np.zeros((3, 3)))
        M[:,0] = S[0,0] * xyzs[:,0]
        M[:,1] = S[1,0] * xyzs[:,1]
        M[:,2] = S[2,0] * xyzs[:,2]
        return M


    @staticmethod
    def __computeChromaticAdaptationMatrix(source, target, method):
        """
        Computes the chromatic adaptation matrix

        :param source: source reference white
        :type source: :class:`Colour_BaseIlluminant`
        :param target: target reference white
        :type target: :class:`Colour_BaseIlluminant`
        :param method: adaptation method
        :type method: numpy.matrix
        :return: chromatic adaptation matrix
        :rtype: numpy.matrix
        """
        S = method * np.array([[source.X], [source.Y], [source.Z]]);
        D = method * np.array([[target.X], [target.Y], [target.Z]]);
        tmp = np.matrix(np.zeros((3, 3)))
        tmp[0, 0] = D[0] / S[0]
        tmp[1, 1] = D[1] / S[1]
        tmp[2, 2] = D[2] / S[2]
        return method.I * tmp * method


    # used when converting XYZ to RGB
    @staticmethod
    def _gamma_companding(gamma, t1, t2, t3):
        p = 1.0 / gamma
        if t1 < 0.0:
            t1 = math.pow(-t1, p) * -1.0
        else:
            t1 = math.pow(t1, p)

        if t2 < 0.0:
            t2 = math.pow(-t2, p) * -1.0
        else:
            t2 = math.pow(t2, p)

        if t3 < 0.0:
            t3 = math.pow(-t3, p) * -1.0
        else:
            t3 = math.pow(t3, p)
        return (t1, t2, t3)


    # used when converting RGB to XYZ
    @staticmethod
    def _inverse_gamma_companding(gamma, r, g, b):
        import numpy as np
        t = np.zeros((3, 1))
        if r < 0.0:
            t[0, 0] = math.pow(-r, gamma) * -1.0
        else:
            t[0, 0] = math.pow(r, gamma)

        if g < 0.0:
            t[1, 0] = math.pow(-g, gamma) * -1.0
        else:
            t[1, 0] = math.pow(g, gamma)

        if b < 0.0:
            t[2, 0] = math.pow(-b, gamma) * -1.0
        else:
            t[2, 0] = math.pow(b, gamma)
        return t


    # used when converting XYZ to sRGB
    @staticmethod
    def _sRGB_companding(gamma, t1, t2, t3):
        p = 1.0 / gamma
        if t1 > 0.0031308:
            t1 = 1.055 * math.pow(t1, p) - 0.055
        else:
            t1 = 12.92 * t1

        if t2 > 0.0031308:
            t2 = 1.055 * math.pow(t2, p) - 0.055
        else:
            t2 = 12.92 * t2

        if t3 > 0.0031308:
            t3 = 1.055 * math.pow(t3, p) - 0.055
        else:
            t3 = 12.92 * t3
        return (t1, t2, t3)


    # used when converting sRGB to XYZ
    @staticmethod
    def _inverse_sRGB_companding(gamma, r, g, b):
        import numpy as np
        t = np.zeros((3, 1))
        if r > 0.04045:
            t[0, 0] = math.pow((r + 0.055) / 1.055, gamma)
        else:
            t[0, 0] = r / 12.92

        if g > 0.04045:
            t[1, 0] = math.pow((g + 0.055) / 1.055, gamma)
        else:
            t[1, 0] = g / 12.92

        if b > 0.04045:
            t[2, 0] = math.pow((b + 0.055) / 1.055, gamma)
        else:
            t[2, 0] = b / 12.92
        return t        


class RGB(BaseRGB):

    def __init__(self, gamma, rw, rp, gp, bp, r, g, b):
        """
        see BaseRGB __init__().
        """
        if(self.__class__.__name__.strip('Colour_') == 'sRGB'):
            self.companding = BaseRGB._sRGB_companding
            self.inverse_companding = BaseRGB._inverse_sRGB_companding
        else:
            self.companding = BaseRGB._gamma_companding
            self.inverse_companding = BaseRGB._inverse_gamma_companding
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


    def from_XYZ(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in 'XYZ' colour space to a 'RGB' colour space.

        :param col: colour
        :type col: :class:`Colour_XYZ`
        """
        self.tri = (self.m_1_adapted * np.array([[col._t1],
                                                 [col._t2],
                                                 [col._t3]])).flat
        self.tri = self.companding(self.gamma, self._t1, self._t2, self._t3)


    def from_xyY(self, col):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in 'xyY' colour space to a 'RGB' colour space.

        :param col: colour
        :type col: :class:`Colour_xyY`
        """        
        self.from_XYZ(col.to_XYZ())


    def from_Lab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in 'Lab' colour space to a 'RGB' colour space.

        :param col: colour
        :type col: :class:`Colour_Lab`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw))


    def from_LCHab(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in 'LCHab' colour space to a 'RGB' colour space.

        :param col: colour
        :type col: :class:`Colour_LCHab`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw))


    def from_Luv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in 'Luv' colour space to a 'RGB' colour space.

        :param col: colour
        :type col: :class:`Colour_Luv`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw))


    def from_LCHuv(self, col, rw):
        """
        Sets the tristimulus values to new values computed by transforming
        tristimulus values in 'LCHuv' colour space to a 'RGB' colour space.

        :param col: colour
        :type col: :class:`Colour_LCHuv`
        :param rw: reference white
        :type rw: :class:`BaseIlluminant`
        """
        self.from_XYZ(col.to_XYZ(rw))


    for cs in ColourSpace.supported_colour_spaces():
        reference_white_not_required = ['XYZ', 'xyY']
        if cs is not 'RGB':
            if cs in reference_white_not_required:
                exec("def to_{0}(self):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self)\n"
                     "  return rt".format(cs, 'RGB'))
            else:
                exec("def to_{0}(self, rw):\n"
                     "  rt = Colour_{0}()\n"
                     "  rt.from_{1}(self, rw)\n"
                     "  return rt".format(cs, 'RGB'))
    del cs
    del reference_white_not_required


class Colour_AdobeRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in Adobe RGB colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.0 + 51.0 / 256.0
        rw = illum.IlluminantD65('1931_2')
        rp = Colour_xyY(0.64, 0.33)
        gp = Colour_xyY(0.21, 0.71)
        bp = Colour_xyY(0.15, 0.06)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_AppleRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in Apple RGB colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 1.8
        rw = illum.IlluminantD65('1931_2')
        rp = Colour_xyY(0.625, 0.340)
        gp = Colour_xyY(0.280, 0.595)
        bp = Colour_xyY(0.155, 0.070)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_BestRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in Best RGB colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantD50('1931_2')
        rp = Colour_xyY(0.7347, 0.2653)
        gp = Colour_xyY(0.2150, 0.7750)
        bp = Colour_xyY(0.1300, 0.0350)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_BetaRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in Beta RGB colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantD50('1931_2')
        rp = Colour_xyY(0.6888, 0.3112)
        gp = Colour_xyY(0.1986, 0.7551)
        bp = Colour_xyY(0.1265, 0.0352)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_BruceRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in Bruce RGB colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantD65('1931_2')
        rp = Colour_xyY(0.64, 0.33)
        gp = Colour_xyY(0.28, 0.65)
        bp = Colour_xyY(0.15, 0.06)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_CIERGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in CIE RGB colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantE('1931_2')
        rp = Colour_xyY(0.735, 0.265)
        gp = Colour_xyY(0.274, 0.717)
        bp = Colour_xyY(0.167, 0.009)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_ColorMatchRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in ColorMatch RGB colour
        space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 1.8
        rw = illum.IlluminantD50('1931_2')
        rp = Colour_xyY(0.630, 0.340)
        gp = Colour_xyY(0.295, 0.605)
        bp = Colour_xyY(0.150, 0.075)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_DonRGB4(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in Don RGB 4 colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantD50('1931_2')
        rp = Colour_xyY(0.696, 0.300)
        gp = Colour_xyY(0.215, 0.765)
        bp = Colour_xyY(0.130, 0.035)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_ECIRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in ECI RGB colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 1.8
        rw = illum.IlluminantD50('1931_2')
        rp = Colour_xyY(0.67, 0.33)
        gp = Colour_xyY(0.21, 0.71)
        bp = Colour_xyY(0.14, 0.08)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_EktaSpacePS5(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in Ekta Space PS5 colour
        space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantD50('1931_2')
        rp = Colour_xyY(0.695, 0.305)
        gp = Colour_xyY(0.260, 0.700)
        bp = Colour_xyY(0.110, 0.005)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_NTSCRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in NTSC RGB colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantC('1931_2')
        rp = Colour_xyY(0.67, 0.33)
        gp = Colour_xyY(0.21, 0.71)
        bp = Colour_xyY(0.14, 0.08)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_PAL_SECAMRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in PAL/SECAM RGB colour
        space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantD65('1931_2')
        rp = Colour_xyY(0.64, 0.33)
        gp = Colour_xyY(0.29, 0.60)
        bp = Colour_xyY(0.15, 0.06)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_ProPhotoRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in ProPhoto RGB colour
        space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 1.8
        rw = illum.IlluminantD50('1931_2')
        rp = Colour_xyY(0.7347, 0.2653)
        gp = Colour_xyY(0.1596, 0.8404)
        bp = Colour_xyY(0.0366, 0.0001)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_SMPTE_CRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in SMPTE-C RGB colour
        space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantD65('1931_2')
        rp = Colour_xyY(0.630, 0.340)
        gp = Colour_xyY(0.310, 0.595)
        bp = Colour_xyY(0.155, 0.070)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_WideGamutRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in Wide Gamut RGB colour
        space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.2
        rw = illum.IlluminantD50('1931_2')
        rp = Colour_xyY(0.735, 0.265)
        gp = Colour_xyY(0.115, 0.826)
        bp = Colour_xyY(0.157, 0.018)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


class Colour_sRGB(RGB):

    def __init__(self, r=1.0, g=1.0, b=1.0):
        """
        Initializes the tristimulus values of a colour in sRGB colour space.

        :param r: The r coordinate
        :type r: float
        :param g: The g coordinate
        :type g: float
        :param b: The b coordinate
        :type b: float
        """
        import colpy.core.illum as illum
        gamma = 2.4
        rw = illum.IlluminantD65('1931_2')
        rp = Colour_xyY(0.64, 0.33)
        gp = Colour_xyY(0.30, 0.60)
        bp = Colour_xyY(0.15, 0.06)
        super().__init__(gamma, rw, rp, gp, bp, r, g, b)


def colour(t1, t2, t3, cs):
    if cs not in ColourSpace.supported_colour_spaces():
        raise ValueError("unsupported colour space")
    return eval("Colour_{0}(t1, t2, t3)".format(cs))


class converter():

    def __init__(self, source, target):
        self.__src = eval("Colour_{0}()".format(source))
        tar = eval("Colour_{0}()".format(target))
        if source in _RGB:
            self.__func = eval("tar.from_RGB")
        else:
            self.__func = eval("tar.from_{0}".format(source))

    def __call__(self, t1, t2, t3, rw=None):
        self.__src.tri = (t1, t2, t3)
        if rw is None:
            self.__func(self.__src)
        return self.__func.__self__.tri


#if __name__ == "__main__":
#    import doctest
#    doctest.testmod()
