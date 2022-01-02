import sys
import os
import time
import numpy as np
from scipy.spatial import Delaunay

def TxtPoints2TIN(points_input_file,tin_output_folder):
    '''
    Triangulation from input points.
    Return the triangulation results with .ff file
    
    points_input_file: point cloud data in .txt format
    tin_output_folder: generated TIN in .off format
    '''
    points = np.empty(shape=[0, 2])
    z_old=list()
    
    # count the lines of the whole txt file
    lines_num = 0
    for lines_num,line in enumerate(open(points_input_file,'rU').readlines()):
        lines_num += 1
    
    with open(points_input_file) as infile:
        for l in range(lines_num):
            line = (infile.readline()).split()
            v1 = float(line[0])
            v2 = float(line[1])
            existed_p = False
            for p in points:
                if p[0]==v1 and p[1]==v2:
                    existed_p = True
                    break
            if existed_p == False:
                points = np.append(points, [[v1,v2]], axis=0)
                z_old.append(float(line[2]))
    
    triangles = Delaunay(points) 
    
    used_pts = set()
    for tri in triangles.simplices:
        used_pts.add(tri[0])
        used_pts.add(tri[1])
        used_pts.add(tri[2])
    old2new = dict()
    new2old = list()
    z_new = list()
    points1= np.empty(shape=[0, 2])
    for pid in used_pts:
        points1 = np.append(points1,[[points[pid][0], points[pid][1]]],axis=0)
        old2new[pid] = len(points1)-1
        new2old.append(pid)
        z_new.append(z_old[pid])
        
    tris = np.empty(shape=[0, 3])
    for tri in triangles.simplices:
        v1 = tri[0]
        v2 = tri[1]
        v3 = tri[2]
        nv1 = old2new[v1]
        nv2 = old2new[v2]
        nv3 = old2new[v3]
        tris = np.append(tris,[[nv1, nv2, nv3]],axis=0)
   
   # output TIN file, the finename is the same as input file but with different suffix
    TIN = tin_output_folder + '/' + os.path.basename(points_input_file).split('.')[0]  + '.off'
    # TIN = "pts-dt.off" 
    with open(TIN,'w') as ofs:
        ofs.write("OFF\n")
        ofs.write("{} {} 0\n".format(len(points1), len(tris)))
        for pid in range(len(points1)):
            ofs.write("{} {} {}\n".format(points[pid][0], points[pid][1],z_new[pid]))
        for tri in tris:
            ofs.write("3 {} {} {}\n".format(int(tri[0]), int(tri[1]), int(tri[2])))
    
    return TIN