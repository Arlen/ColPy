ColPy Tutorial
==============

supported colour spaces
-----------------------
* non-RGB colour spaces:  ``XYZ``, ``xyY``, ``Lab``, ``LCHab``, ``Luv``, ``LCHuv``
* RGB colour spaces: ``AdobeRGB``, ``AppleRGB``, ``BestRGB``, ``BetaRGB``, ``BruceRGB``, ``CIERGB``, ``ColorMatchRGB``, ``DonRGB4``, ``ECIRGB``, ``EktaSpacePS5``, ``NTSCRGB``, ``PAL_SECAMRGB``, ``ProPhotoRGB``, ``SMPTE_CRGB``, ``WideGamutRGB``, ``sRGB``.

creating colours
----------------
Each colour spaces has its own specific Python class, and all the colour space classes are prefixed with ``Colour_``.  For example, to create a colour in ``XYZ``:

>>> # First make sure you import ColPy
>>> import colpy as cp
>>> # Then create a new colour
>>> col = cp.Colour_XYZ()
>>> col
Colour_XYZ(1.0, 1.0, 1.0)


ColPy will initialize the tristimulus values of a colour to a default value if none given.

>>> # Initialize the tristimulus values of a colour
>>> col = cp.Colour_AppleRGB(r=.6, g=.8, b=.4)
>>> col
Colour_AppleRGB(0.6, 0.8, 0.4)
>>> # or for short:
>>> col = cp.Colour_AppleRGB(.6, .8, .4)
>>> col
Colour_AppleRGB(0.6, 0.8, 0.4)


.. note::
    For non-RGB colour spaces the named parameters for initialization are the firs three letters of the colour space.  For RGB colour spaces the named parameters are ``r``, ``g``, and ``b``.


The tristimulus values can be accessed individually:

>>> col = cp.Colour_xyY(.35, .9, .39)
>>> col.x
.35
>>> col.Y = .6
>>> col
Colour_xyY(0.35, 0.9, 0.6)


If you need to set or get all the tristimulus values of a colour, use :meth:`tri`:

>>> col = cp.Colour_xyY(.35, .9, .39)
>>> col.tri
(.35, .9, .39)
>>> col.tri = .1, .2, .6
>>> col
Colour_xyY(0.1, 0.2, 0.6)


In certain situations it might be convenient to use a convenience function to create colour objects.  :func:`colour` does just that:

>>> col = cp.colour(.2, .5, .8, 'Lab')
>>> col
Colour_Lab(.2, .5, .8)



supported illuminants
---------------------
The following standardized illuminants are supported:

* ``A``, ``B``, ``C``, ``D50``, ``D55``, ``D65``, ``D75``, ``E``, ``F1``, ``F2``, ``F3``, ``F4``, ``F5``, ``F6``, ``F7``, ``F8``, ``F9``, ``F10``, ``F11``, ``F12``.

.. note::
    Both 2 degree field of view (1931) and the 10 degree field of view (1964) are available.

creating illuminants
--------------------
Each illuminant has its own specific Python class, and all illuminant classes are prefixed with ``Illuminant_``.

>>> rw = cp.IlluminantD65()
>>> rw
lluminantD65(0.94809667673716, 1.0, 1.0730513595166162)
>>> # 1964_10 is the default; if you need 1931_2, do as such:
>>> rw = cp.IlluminantD65(std_obs='1931_2')
>>> rw
IlluminantD65(0.9504285453771807, 1.0, 1.0889003707981277)

Internally, a :class:`Colour_XYZ` object is used to store the data in an illuminant class.  It can not be overwritten, but you may access:

>>> rw = cp.IlluminantD50()
>>> rw.XYZ
Colour_XYZ(0.9672062750333777, 1.0, 0.8142801513128616)
>>> # or any of XYZ's coordinates
>>> rw.Y
1.0
>>> rw.X
0.9672062750333777

colour conversion
-----------------
Any colour value can be converted from a colour space to another colour space, and there are two basic ways to do this.  The first is to use the *from_** methods:

>>> # to convert xyY to XYZ, create the target object:
>>> tar = cp.Colour_XYZ()
>>> tar
Colour_XYZ(1.0, 1.0, 1.0)
>>> # create the source
>>> src = cp.Colour_xyY(.33, .4, .9)
>>> src
Colour_xyY(0.33, 0.4, 0.9)
>>> # then call the appropriate from_* method
>>> tar.from_xyY(src)
>>> tar
Colour_XYZ(0.7424999999999999, 0.9, 0.6074999999999997)
>>> # the from_* methods do not return anything.  The object holds the new converted values.

The second method is to use the *to_** methods.  Similar to *from_**, except they return new objects:

>>> in this case only the source is needed.
>>> src = cp.Colour_xyY(.33, .4, .9)
>>> src
Colour_xyY(0.33, 0.4, 0.9)
>>> tar = src.to_XYZ()
>>> tar
Colour_XYZ(0.7424999999999999, 0.9, 0.6074999999999997)

Converting from any RGB is the same as above, but when converting to RGB there is no specific *to_** method:

>>> # call the specific to_ method to convert to any RGB:
>>> src.to_AdobeRGB()
Colour_AdobeRGB(0.9024592551310153, 0.9972363370485827, 0.7429555634969075)
>>> # call the from_RGB() method to convert from any RGB
>>> tar = cp.Colour_xyY()
>>> col = cp.Colour_AdobeRGB(.5, .5, .8)
>>> col
Colour_AdobeRGB(0.5, 0.5, 0.8)
>>> tar.from_RGB(col)
>>> tar
Colour_xyY(0.24309894987593703, 0.21392710648162114, 0.24744698754582328)

Depending on the source or the target, conversion may require a reference white.

>>> # converting xyY to Lab or Luv will require a reference white
>>> src = cp.Colour_xyY(.3, .4, .9)
>>> rw = cp.IlluminantD65()
>>> src.to_Lab(rw)
Colour_Lab(95.99676861425304, -36.28006732658423, 21.731972691334754)


In certain situations it might be convenient to use a convenience class to convert colour from one colour space to another.  :class:`converter` does just that:

>>> c = cp.converter(source='ProPhotoRGB', target='xyY')
>>> # then simply call c() with colour values in source colour space.
>>> # the functor will return a 3-tuple of floats in the target colour space.
>>> c(.8, .8, .3)
(0.45121136299250814, 0.48088659364064756, 0.66916179607839854)
