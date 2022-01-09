\import numpy as np
import os
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk import (vtkUnstructuredGridReader, vtkDataSetMapper, vtkActor,vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor)
from mayavi import mlab

def vtk_plot(file_name):
    '''
    
    '''
    # Read the source file.
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()
    output_port = reader.GetOutputPort()
    scalar_range = output.GetScalarRange()

    
    # obtain vtk_array and then convert them into numpy_array
    nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
    # print(nodes_vtk_array)
    
    # read txt file
    #file_name_txt = "C:/PhD_study/Programming/data/morse/pts-dt_kv_20_desc2cells.txt"  
    file_name_txt = os.path.dirname(file_name) + '/' + os.path.basename(file_name).split('.')[0] + '.txt'
    with open(file_name_txt, 'r') as f:
        lines = f.readlines()
        last_line = lines[-1]
    #convert str to number
    last_line = str.split(last_line)
    color_real = []
    for i in range(len(last_line)):
        color_real.append(float(last_line[i]))
    
    lookUpTable = vtk.vtkLookupTable()
    #lookUpTable.SetTableRange(z_min,z_max)
    if min(color_real) > 0:
        lookUpTable.SetTableRange(int(min(color_real)/(max(color_real)) * 255),255)
    if min(color_real) < 0:
        lookUpTable.SetTableRange(0,255)
    lookUpTable.Build()

    # Create a unsigned char array to color
    uCharArray = vtk.vtkUnsignedCharArray()
    uCharArray.SetNumberOfComponents(3)
    uCharArray.SetName("colors")

    # Assign color by extracting each color
    for i in range(output.GetNumberOfPoints()):
        color = [0]*3
        color[0] = 150
        color[2] = 255
        #color[1] = int(color_real[i]/(max(color_real)) * 255)
        color[1] = float((255-0)/(max(color_real) - min(color_real)) * (color_real[i]-min(color_real)) + 0)
        color[1] = int(color[1])
        #print("i:",i)
        #print("color:",color)
        uCharArray.InsertTypedTuple(i,color)
        
    # # Set Scalars
    output.GetPointData().SetScalars(uCharArray)

    # Create the mapper that corresponds the objects of the vtk file
    # into graphics elements
    mapper = vtkDataSetMapper()
    mapper.SetInputConnection(output_port)
    mapper.SetColorModeToDefault()
    mapper.SetScalarRange(scalar_range)

    # Create the Actor
    actor = vtkActor()
    actor.SetMapper(mapper)

    colors = vtk.vtkNamedColors()
    actor.GetProperty().SetColor(colors.GetColor3d("Bisque"))

    # Create the Renderer
    renderer = vtkRenderer()
    renderer.AddActor(actor)
    # renderer.AddActor(colorbar)
    renderer.SetBackground(colors.GetColor3d('Silver')) # Set background to white
    renderer.ResetCamera()

    # Create the RendererWindow
    renderer_window = vtkRenderWindow()
    renderer_window.AddRenderer(renderer)

    # Create the RendererWindowInteractor and display the vtk_file
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)
    interactor.Initialize()
    interactor.Start()
