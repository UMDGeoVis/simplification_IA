#cdef extern from "Lib_Tri/Vertex3D.cpp":
#    pass

# template<class V>
# Declare the class Mesh with cdef
cdef extern from "Lib_Tri/Vertex3D.h":
    cdef cppclass Vertex3D:
        Vertex3D() except +