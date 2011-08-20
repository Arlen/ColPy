#!/usr/bin/env python3
"""ColPy: conversion and computation of colour.

"""

DOCLINES = __doc__.split("\n")

from distutils.core import setup

CLASSIFIERS = [
    'Development Status :: 2 - Pre-Alpha',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License (GPL)',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Operating System :: MacOS',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Topic :: Software Development',
    'Topic :: Scientific/Engineering',
    ]

NAME                = 'colpy'
PACKAGES            = ['colpy']
MAINTAINER          = "Arlen Aghakians"
MAINTAINER_EMAIL    = "arlen.prime@gmail.com"
DESCRIPTION         = DOCLINES[0]
#LONG_DESCRIPTION   =
URL                 = "https://github.com/Arlen/ColPy"
#DOWNLOAD_URL        = ""
LICENSE             = 'GPLv3'
AUTHOR              = "Arlen Aghakians"
AUTHOR_EMAIL        = "arlen.prime@gmail.com"
#PLATFORMS           = ["Windows", "Linux", "Mac OS-X"]
MAJOR               = 0
MINOR               = 1
MICRO               = 0
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


setup(
    name=NAME,
    packages=PACKAGES,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    description=DESCRIPTION,
    #long_description=LONG_DESCRIPTION
    url=URL,
    #download_url=DOWNLOAD_URL,
    license=LICENSE,
    classifiers=CLASSIFIERS,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    #platforms=PLATFORMS,
    version=VERSION
    )
