from distutils.core import setup
from Cython.Build import cythonize
from os import system

system('cp modPolyCy.spyx modPolyCy.pyx')

setup(
    name='ModPoly',
    ext_modules=cythonize('modPolyCy.pyx'),
)
