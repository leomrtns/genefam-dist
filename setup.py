"""Build the native CPython extension.

The C library sources are compiled into the extension so a separate system-wide
libgenefamdist installation is not required.
"""

import sys

from setuptools import Extension, setup


library_sources = [
    "hashtable.c",
    "lowlevel.c",
    "random_number_gen.c",
    "random_number.c",
    "nexus_common.c",
    "topology_common.c",
    "topology_space.c",
    "topology_mrca.c",
    "bipartition.c",
    "topology_build.c",
    "topology_splitset.c",
    "prob_distribution.c",
    "empirical_frequency.c",
    "genefam_dist.c",
]

define_macros = [("_GNU_SOURCE", "1")]
extra_compile_args = []
libraries = []

if sys.platform != "win32":
    extra_compile_args.append("-std=c11")
    libraries.append("m")

extension = Extension(
    "treesignal._treesignalc",
    sources=[
        "python/treesignal/_treesignalc.c",
        *(f"lib/{source}" for source in library_sources),
    ],
    include_dirs=["lib"],
    define_macros=define_macros,
    extra_compile_args=extra_compile_args,
    libraries=libraries,
)

setup(ext_modules=[extension])
