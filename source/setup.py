from setuptools import setup

from Cython.Build import cythonize

from distutils.core import setup, Extension

# setup(ext_modules=cythonize("rect.pyx"))

import glob

cpp_files = glob.glob('**/*.cpp', recursive=True)

'''
ext = [Extension("PythonMain",
    sources=["rect.pyx"] + cpp_files,
    language="c++",
    include_dirs=['..'])]
setup(ext_modules=cythonize(ext))
'''
#if "rect.cpp" in cpp_files:
    #cpp_files.remove("rect.cpp")

if "main.cpp" in cpp_files:
    cpp_files.remove("main.cpp")

if "Mesh.cpp" in cpp_files:
    cpp_files.remove("Mesh.cpp")

if "Vertex3D.cpp" in cpp_files:
    cpp_files.remove("Vertex3D.cpp")

if "Triangle.cpp" in cpp_files:
    cpp_files.remove("Triangle.cpp")

if "Reader.cpp" in cpp_files:
    cpp_files.remove("Reader.cpp")

if "FormanGradientVector.cpp" in cpp_files:
    cpp_files.remove("FormanGradientVector.cpp")

# if "Lib/Rectangle.cpp" in cpp_files:
#     cpp_files.remove("Lib/Rectangle.cpp")

print("cpp_files:",cpp_files)

# sources_print = ["Rectangle_pyx.pyx", "Square_pyx.pyx"] + cpp_files
# sources_print = ["Rectangle_pyx.pyx"] + cpp_files
sources_print = ["Mesh_Forman.pyx"] + cpp_files

include_dir_print = ['..','Dag','Eigen','Ig','Lib_Tri']

print("sources:", sources_print)
print("include_dirs:", include_dir_print)

ext = [Extension("PythonMain",
    sources=sources_print,
    language="c++",
    include_dirs=include_dir_print)]
setup(ext_modules=cythonize(ext))