from distutils.core import setup, Extension

module1 = Extension('treesignal', 
                     define_macros = [('MAJOR_VERSION', '1'), ('MINOR_VERSION', '0')],
#                     include_dirs = ['/home/leo/local/include'],
#                     libraries = ['tcl83'],
                     library_dirs = ['/home/leo/local/lib'],
                     runtime_library_dirs = ['/home/leo/local/lib'],
                     libraries = ['genefamdist'],
                    sources = ['treesignalmodule.c'])

setup (name = 'treesignal', 
       version = '1.0', 
       description = 'Calculates the signal of gene family trees',
       author = 'Leonardo de Oliveira Martins',
       author_email = 'leomrtns@gmail.com',
       url='https://github.com/leomrtns/genefam-dist',
       long_description='''
This package provide C functions to calculate distances between gene family trees and sets of species trees.
''',
       ext_modules = [module1])

