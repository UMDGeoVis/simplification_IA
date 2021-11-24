# line 3 is to indicate to Cython that this .pyx file 
# has to be compiled to C++
# distutils: language = c++

from Vertex3D cimport Vertex3D
cdef class PyVertex3D:
    cdef Vertex3D c_Vertex3D
    def __cinit__(self):
        self.c_Vertex3D = Vertex3D()

from Triangle cimport Triangle
cdef class PyTriangle:
    cdef Triangle c_Triangle
    def __cinit__(self):
        self.c_Triangle = Triangle()


from Mesh cimport Mesh
#cimport Mesh

# Implement Python extension type
cdef class PyMesh:
    cdef Mesh[Vertex3D,Triangle] *c_Mesh  # Hold a pointer to the C++ instance which we're wrapping

    def __cinit__(self):
        self.c_Mesh = new Mesh[Vertex3D,Triangle]()

    def Build(self):
        self.c_Mesh.build()
        
    #def __cinit__(self, const Mesh& orig):
    #    self.c_Mesh = new Mesh(orig)

    def __dealloc__(self):
        del self.c_Mesh

from Reader cimport Reader
from cython.operator cimport dereference as deref

cdef class PyReader:
    cdef Reader *c_Reader

    #cdef read_MeshFile(self, Mesh[Vertex3D,Triangle]& mesh, str path):
    #    path_new = bytes(path, encoding='utf8')
    #    
    #    return Reader.readMeshFile(mesh, path_new)

    def PyreadMeshFile(self, PyMesh mesh, str path):
        path_new = bytes(path, encoding='utf8')
        
        return self.c_Reader.readMeshFile(deref(mesh.c_Mesh), path_new)

    #‘virtual Reader::~Reader()’ is private within this context
    #def __dealloc__(self):
    #    del self.c_Reader

from FormanGradientVector cimport FormanGradientVector
from libcpp cimport bool
# Implement Python extension type

cdef class PyFormanGradientVector:
    cdef FormanGradientVector *c_FormGraVec
    
    def __cinit__(self, PyMesh mesh, double epsilon):
        self.c_FormGraVec = new FormanGradientVector(mesh.c_Mesh, epsilon)
    
    def initial_filtering(self):
        return self.c_FormGraVec.initial_filtering()

    def build(self):
        return self.c_FormGraVec.build()

    def change_vtstar_mesh(self):
        return self.c_FormGraVec.change_vtstar_mesh()

    def simplify_geometry(self, bool tf, double dou):
        return self.c_FormGraVec.simplify_geometry(tf,dou)

    def descending_2cells_extraction(self, bool tf):
        return self.c_FormGraVec.descending_2cells_extraction(tf)

    def descending_1cells_extraction(self, bool tf):
        return self.c_FormGraVec.descending_1cells_extraction(tf)

    def ascending_2cells_extraction(self, bool tf):
        return self.c_FormGraVec.ascending_2cells_extraction(tf)

    def ascending_1cells_extraction(self, bool tf):
        return self.c_FormGraVec.ascending_1cells_extraction(tf)
    
    def compute_incidence_graph(self):
        return self.c_FormGraVec.compute_incidence_graph()

    def writeVTK_IG(self, str path):
        path_new = bytes(path, encoding='utf8')
        return self.c_FormGraVec.writeVTK_IG(path_new)
    
    def __dealloc__(self):
        del self.c_FormGraVec