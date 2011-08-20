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
import colpy.core.colsp as colsp


class BaseIlluminant(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def __init__(self, x, y):
        self.__XYZ = colsp.Colour_xyY(x, y, 1.0).to_XYZ()
        
        # set all the 'from_' methods to None to prevent overwrite
        for cs in colsp.ColourSpace.supported_colour_spaces():
            attrname = "from_" + cs
            if hasattr(self.XYZ, attrname):
                setattr(self.XYZ, attrname, None)
                

    def __repr__(self):
        return "{0}{1}".format(self.__class__.__name__, self.XYZ.tri)

    
    def __str__(self):
        return "{0}{1}".format(self.__class__.__name__.strip("Illuminant"),
                               self.XYZ.tri)
    

    @property
    def XYZ(self):
        return self.__XYZ


    @property
    def X(self):
        return self.XYZ.X


    @property
    def Y(self):
        return self.XYZ.Y


    @property
    def Z(self):
        return self.XYZ.Z

    
_A_definition   = {'1964_10': (0.45117, 0.40594), '1931_2': (0.44757, 0.40745)}
_B_definition   = {'1964_10': (0.3498, 0.3527),   '1931_2': (0.34842, 0.35161)}
_C_definition   = {'1964_10': (0.31039, 0.31905), '1931_2': (0.31006, 0.31616)}
_D50_definition = {'1964_10': (0.34773, 0.35952), '1931_2': (0.34567, 0.3585)}
_D55_definition = {'1964_10': (0.33411, 0.34877), '1931_2': (0.33242, 0.34743)}
_D65_definition = {'1964_10': (0.31382, 0.331),   '1931_2': (0.31271, 0.32902)}
_D75_definition = {'1964_10': (0.29968, 0.3174),  '1931_2': (0.29902, 0.31485)}
_E_definition   = {'1964_10': (1/3, 1/3),         '1931_2': (1/3, 1/3)}
_F1_definition  = {'1964_10': (0.31811, 0.33559), '1931_2': (0.3131, 0.33727)}
_F2_definition  = {'1964_10': (0.37925, 0.36733), '1931_2': (0.37208, 0.37529)}
_F3_definition  = {'1964_10': (0.41761, 0.38324), '1931_2': (0.4091, 0.3943)}
_F4_definition  = {'1964_10': (0.4492, 0.39074),  '1931_2': (0.44018, 0.40329)}
_F5_definition  = {'1964_10': (0.31975, 0.34246), '1931_2': (0.31379, 0.34531)}
_F6_definition  = {'1964_10': (0.3866, 0.37847),  '1931_2': (0.3779, 0.38835)}
_F7_definition  = {'1964_10': (0.31569, 0.3296),  '1931_2': (0.31292, 0.32933)}
_F8_definition  = {'1964_10': (0.34902, 0.35939), '1931_2': (0.34588, 0.35875)}
_F9_definition  = {'1964_10': (0.37829, 0.37045), '1931_2': (0.37417, 0.37281)}
_F10_definition = {'1964_10': (0.3509, 0.35444),  '1931_2': (0.34609, 0.35986)}
_F11_definition = {'1964_10': (0.38541, 0.37123), '1931_2': (0.38052, 0.37713)}
_F12_definition = {'1964_10': (0.44256, 0.39717), '1931_2': (0.43695, 0.40441)}


class IlluminantA(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _A_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_A_definition.get(std_obs))
        

class IlluminantB(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _B_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_B_definition.get(std_obs))


class IlluminantC(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _C_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_C_definition.get(std_obs))

      
class IlluminantD50(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _D50_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_D50_definition.get(std_obs))

        
class IlluminantD55(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _D55_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_D55_definition.get(std_obs))

        
class IlluminantD65(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _D65_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_D65_definition.get(std_obs))

        
class IlluminantD75(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _D75_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_D75_definition.get(std_obs))

        
class IlluminantE(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _E_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_E_definition.get(std_obs))

        
class IlluminantF1(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F1_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F1_definition.get(std_obs))

        
class IlluminantF2(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F2_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F2_definition.get(std_obs))

        
class IlluminantF3(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F3_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F3_definition.get(std_obs))

       
class IlluminantF4(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F4_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F4_definition.get(std_obs))

        
class IlluminantF5(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F5_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F5_definition.get(std_obs))

        
class IlluminantF6(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F6_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F6_definition.get(std_obs))

        
class IlluminantF7(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F7_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F7_definition.get(std_obs))

        
class IlluminantF8(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F8_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F8_definition.get(std_obs))

        
class IlluminantF9(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F9_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F9_definition.get(std_obs))

        
class IlluminantF10(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F10_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F10_definition.get(std_obs))

        
class IlluminantF11(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F11_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F11_definition.get(std_obs))

        
class IlluminantF12(BaseIlluminant):

    def __init__(self, std_obs='1964_10'):
        if std_obs not in _F12_definition.keys():
            raise ValueError("unsupported CIE standard observer")
        super().__init__(*_F12_definition.get(std_obs))


class IlluminantCustom(BaseIlluminant):

    def __init__(self, x, y):
        super().__init__(x, y)
