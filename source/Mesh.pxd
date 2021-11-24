cdef extern from "Lib_Tri/Mesh.h":
    pass

# template<class V>
# Declare the class Mesh with cdef
cdef extern from "Lib_Tri/Mesh.h":
    cdef cppclass Mesh[V,T]:
        Mesh() except +
        # Mesh(const Mesh) except +
        #V& getVertex(int)
        #double getX()
        #double getY()
        #double getZ()
        #int getNumVertex()
        #int getTopSimplexesNum()
        void build()

'''
cdef extern from "Lib_Tri/Vertex2D.cpp":
    pass
cdef extern from "Lib_Tri/Vertex2D.h":
    cdef cppclass Vertex2D:
        Vertex2D() except + 
        double getX()
        double getY()

cdef extern from "Lib_Tri/Vertex3D.cpp":
    pass
cdef extern from "Lib_Tri/Vertex3D.h":
    cdef cppclass Vertex3D:
        Vertex3D() except +
        double getZ()
'''