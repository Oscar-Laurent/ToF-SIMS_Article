from setuptools import setup, Extension
import pybind11
import numpy as np
from pybind11.setup_helpers import Pybind11Extension, build_ext

# Get the include paths
pybind11_include = pybind11.get_include()
numpy_include = np.get_include()
eigen_include = '/opt/homebrew/opt/eigen/include/eigen3/'  # Adjust this path if Eigen is installed elsewhere

ext_modules = [
    Pybind11Extension(
        'nnls_admm',
        sources=['nnls_admm.cpp'],
        include_dirs=[
            pybind11_include,  # Path to pybind11 headers
            numpy_include,  # Path to NumPy headers
            eigen_include  # Path to Eigen headers
        ],
        language='c++'
    ),
]

setup(
    name='nnls_admm',
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
)