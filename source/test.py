import PythonMain
import time

# create a mesh object
py_mesh = PythonMain.PyMesh()

# input data file
# py_path = "../Data/points_yh_0112.off"
py_path = "../Data/sin_sum.tri"

# read a mesh
py_reader = PythonMain.PyReader()
py_reader.PyreadMeshFile(py_mesh,py_path)
py_mesh.Build()

# create a gradient object
t1_GraCom = time.time()
py_gradient = PythonMain.PyFormanGradientVector(py_mesh, 0.0)
py_gradient.initial_filtering()
py_gradient.build()
py_gradient.change_vtstar_mesh()
t2_GraCom = time.time()
print("- time gradient computation:", t2_GraCom - t1_GraCom)
py_gradient.compute_incidence_graph()
py_gradient.write_mesh_VTK("original_mesh")
py_gradient.writeVTK_IG("original_IG.vtk")

# simplify geometry
t1_SimGeo = time.time()
# test with different threholds
py_gradient.simplify_geometry(True, -1)
t2_SimGeo = time.time()
print("- time simplifications:", t2_SimGeo - t1_SimGeo)
py_gradient.initial_filtering()
py_gradient.descending_2cells_extraction(True)
py_gradient.descending_1cells_extraction(True)
py_gradient.ascending_2cells_extraction(True)
py_gradient.ascending_1cells_extraction(True)
py_gradient.compute_incidence_graph()
# output VTK file
py_gradient.write_mesh_VTK("simplified_mesh")
py_gradient.writeVTK_IG("simplified_IG.vtk")