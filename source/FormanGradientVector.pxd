#cdef extern from "Lib_Tri/formangradientvector.cpp":
#    pass

#include <utility>

from libcpp.pair cimport *
from Vertex3D cimport Vertex3D
from Triangle cimport Triangle
from Mesh cimport Mesh

# Declare the class with cdef
cdef extern from "Lib_Tri/formangradientvector.h":
    cdef cppclass FormanGradientVector:
        #FormanGradientVector() except +
        FormanGradientVector(Mesh[Vertex3D,Triangle]*, double) except +
        inline void initial_filtering()
        void build()
        void change_vtstar_mesh()
        void simplify_geometry(bool, double)
        void descending_2cells_extraction(bool)
        void descending_1cells_extraction(bool)
        void ascending_2cells_extraction(bool)
        void ascending_1cells_extraction(bool)
        pair[double, double] compute_incidence_graph()
        void writeVTK_IG(char*)