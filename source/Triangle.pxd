#cdef extern from "Lib_Tri/Triangle.cpp":
#    pass

# template<class V>
# Declare the class Mesh with cdef
cdef extern from "Lib_Tri/Triangle.h":
    cdef cppclass Triangle:
        Triangle() except +