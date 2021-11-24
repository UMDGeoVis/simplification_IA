from Vertex3D cimport Vertex3D
from Triangle cimport Triangle
from Mesh cimport Mesh
from libcpp.string cimport string
from libcpp cimport bool

# Declare the class Mesh with cdef
cdef extern from "Lib_Tri/Reader.h":
    cdef cppclass Reader:
        @staticmethod
        bool readMeshFile(Mesh[Vertex3D,Triangle]& mesh, string path);        