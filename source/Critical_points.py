import sys
import os
import numpy as np
import pandas as pd
from mayavi import mlab
from mayavi.modules.surface import Surface

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk import (vtkUnstructuredGridReader, vtkDataSetMapper, vtkActor,vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor)

def VTK_CriticalPoints(file_name):
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
    numpy_nodes = vtk_to_numpy(nodes_vtk_array)
    
    # read txt file
    file_name_txt = os.path.dirname(file_name) + '/' + os.path.basename(file_name).split('.')[0] + '.txt'
    with open(file_name_txt, 'r') as f:
        lines = f.readlines()
        #if last_line is empty, this command should be changed to last_line = lines[-2]
        last_line = lines[-1] 
        
    #convert str to numberz
    last_line = str.split(last_line)

    # read and store critical points
    points = np.empty(shape=[0, 4])
    
    for i in range(len(numpy_nodes)):
        x = float(numpy_nodes[i,0])
        y = float(numpy_nodes[i,1])
        z = float(numpy_nodes[i,2])
        point_type = int(last_line[i])
        points = np.append(points, [[x,y,z,point_type]], axis=0)

    # get point type of each critical point
    all_point_type = points[:,3]
    all_point_type = pd.Series(all_point_type)
    
    # get how many kinds of color are included
    all_point_type_set = set(all_point_type)
    
    # show the number of points of different type
    all_point_type.value_counts()
    
    # store the color types in a list
    all_point_type_unique = list(all_point_type_set)
    
    # set different color for different points
    all_point_color = np.zeros(points.shape[0])
    
    index = 0
    for i in range(points.shape[0]):
        for j in range(len(all_point_type_unique)):
            if points[i,3] == all_point_type_unique[j]:
                #all_point_color[i] = float(j/(len(all_point_type_unique)))
                all_point_color[i] = float(all_point_type_unique[j]/max(all_point_type_unique))
    
    crit_point_x = points[:,0]
    crit_point_y = points[:,1]
    crit_point_z = points[:,2]
    crit_point_color = all_point_color
    
    # for function mlab.points3d, scale_factor should be set according to different scales of actual points
    # scale_factor should be set as a larger number if you don't see any points when showing it
    p3d = mlab.points3d(crit_point_x,crit_point_y,crit_point_z,crit_point_color,colormap="copper",scale_factor=6,scale_mode='none')
    mlab.draw()
    mlab.show()