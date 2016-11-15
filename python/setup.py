#from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages

description_multiline = '''
This package creates a coordinate system to compare gene trees based on a set of reference species trees, allowing for
Machine Learning and Multivariate Analysis of these resulting signals.
It includes low-level C functions to calculate distances between gene family trees and sets of species trees.
'''

module_c = Extension('_treesignalc', 
                     define_macros = [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
                     include_dirs = ['/home/leo/local/include'],
#                     libraries = ['tcl83'],
                     library_dirs = ['/home/leo/local/lib'],
                     runtime_library_dirs = ['/home/leo/local/lib'],
                     libraries = ['genefamdist'],
                    sources = ['src/treesignalmodule.c'])

print(find_packages())

requirements = ["dendropy", "numpy"]

setup (name = 'treesignal', 
       version = '1.0', 
       description = 'Calculates the signal of gene family trees',
       author = 'Leonardo de Oliveira Martins',
       author_email = 'leomrtns@gmail.com',
       url='https://github.com/leomrtns/genefam-dist',
       long_description=description_multiline,
#       py_modules=['treesignal'],
       packages=find_packages(),
       install_requires=requirements,
       ext_modules = [module_c])

