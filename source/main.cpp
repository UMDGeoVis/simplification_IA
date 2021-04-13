#include <iostream>
#include <string.h>
#include "source/Lib_Tri/formangradientvector.h"
#include "source/Lib_Tri/Mesh.h"
#include "source/Lib_Tri/Timer.h"

#include <stdio.h>
//#include <vlCore/VisualizationLibrary.hpp>
//#include <vlGLUT/GLUTWindow.hpp>

#define BOLD  "\033[1m\033[33m" //for dark background shell
//#define BOLD "\033[1m\033[31m"  //for white background shell
#define RESET   "\033[0m"


using namespace std;

void print_help();

int main(int argc, char* argv[])
{

    bool output_required = true;

    Mesh<Vertex3D,Triangle> mesh = Mesh<Vertex3D,Triangle>();
    FormanGradientVector* gradient;

    if(argc < 2){
        print_help();
        return 0;
    }

    assert(strcmp("-i",argv[1])==0);

    Reader::readMeshFile(mesh,argv[2]);
    mesh.build();

    gradient = new FormanGradientVector(&mesh, 0.0);
    gradient->initial_filtering();
    gradient->build();
    gradient->change_vtstar_mesh();

    //-----------------------------------------------------------------------------------------------------------------------------------------------
    double minx,miny,minz,maxx,maxy,maxz;
    minx=mesh.getVertex(0).getX();
    miny=mesh.getVertex(0).getY();
    minz=mesh.getVertex(0).getZ();
    maxx=mesh.getVertex(0).getX();
    maxy=mesh.getVertex(0).getY();
    maxz=mesh.getVertex(0).getZ();

    for(int i=1; i<mesh.getNumVertex(); i++){
        if(mesh.getVertex(i).getX() > maxx) maxx=mesh.getVertex(i).getX();
        if(mesh.getVertex(i).getY() > maxy) maxy=mesh.getVertex(i).getY();
        if(mesh.getVertex(i).getZ() > maxz) maxz=mesh.getVertex(i).getZ();

        if(mesh.getVertex(i).getX() < minx) minx=mesh.getVertex(i).getX();
        if(mesh.getVertex(i).getY() < miny) miny=mesh.getVertex(i).getY();
        if(mesh.getVertex(i).getZ() < minz) minz=mesh.getVertex(i).getZ();
    }

    cout << "Vertices:  " << mesh.getNumVertex() << endl;
    cout << "Triangles: " << mesh.getTopSimplexesNum() << endl;
    cout << "(" << minx << "," << miny << "," << minz << ")   (" << maxx << "," << maxy << "," << maxz << ")" << endl;
    //-----------------------------------------------------------------------------------------------------------------------------------------------


    if(argc > 3){

        if(strcmp("-s",argv[3])==0){
            gradient->refined_topo = 0;
            //simplification version
            Timer time_refine_topo;

            gradient->compute_incidence_graph();
            gradient->writeVTK_IG("ig_orig.vtk");
            time_refine_topo.start();
            cout << "Simplification Version:" << endl;
            cout << "   - simplification step" << endl;
            gradient->simplify_persistence(atof(argv[4]));
            cout << "   - simplification ended" << endl;
            time_refine_topo.stop();
            cout << "       - performed simplifications: " << gradient->refined_topo << endl;
            cout << "       - time simplifications:      " << time_refine_topo.getElapsedTime() << endl;

        }
        else if(strcmp("-m",argv[3])==0){
            //multiresolution version
            cout << "Multiresolution Version:" << endl;
            cout << "   - simplification step" << endl;
            gradient->simplify(true, argv[4]);
            cout << "   - simplification ended" << endl;
{
//            Timer time_refine_geom, time_refine_topo;
//            if(strcmp(argv[4],"uniform")==0){

//                cout << "   - uniform resolution refinement" << endl;
//                time_refine_topo.start();
//                gradient->refine_topology(atof(argv[5]));
//                time_refine_topo.stop();
//                cout << "       - topology refined" << endl;

//                time_refine_geom.start();
//                gradient->refine_geometry(atof(argv[5]));
//                time_refine_geom.stop();
//                cout << "       - geometry refined" << endl;
//            }
//            else if(strcmp(argv[2],"variable")==0){

//                cout << "   - variable resolution refinement" << endl;
//                vector<double> box(6);

//                if(atof(argv[3]) < 0 || atof(argv[3]) > 1){
//                    cout << "percentuale di raffinamento non valida: deve essere compresa fra 0 e 1" << endl;
//                    return 1;
//                }

//                box[0]= ((minx+maxx) - (maxx-minx)*atof(argv[3]))/2.0;
//                box[1]= ((miny+maxy) - (maxy-miny)*atof(argv[3]))/2.0;
//                box[2]= minz;

//                box[3]= ((minx+maxx) + (maxx-minx)*atof(argv[3]))/2.0;
//                box[4]= ((miny+maxy) + (maxy-miny)*atof(argv[3]))/2.0;
//                box[5]= maxz;

//                time_refine_topo.start();
//                gradient->refine_topology_box(0,box);
//                time_refine_topo.stop();
//                cout << "       - topology refined" << endl;

//                time_refine_geom.start();
//                gradient->refine_geometry_box(0,box);
//                time_refine_geom.stop();
//                cout << "       - geometry refined" << endl;

//                if(output_required) gradient->writeBoxVTK(box);
//            }

//            cout << "Performed refinement: " << gradient->refined_topo << " " << gradient->refined_geometry << endl;
//            cout << "Tempo per i raffinamenti " << endl;
//            cout << "-topologiche: " << time_refine_topo.getElapsedTime() << endl;
//            cout << "-geometriche: " << time_refine_geom.getElapsedTime() << endl;
       }
        }
        else if(strcmp("-c",argv[3])==0){
            output_required=false;
             cout << "   Geometry simplification" << endl;
       // gradient->write_mesh_VTK("orig_mesh");
        if(strcmp("-q",argv[4])==0){
            Timer simplify_timer;
            output_required=true;
            gradient->writeVTK_gradient("orig_gradiente_doporefine.vtk");
            simplify_timer.start();
            gradient->simplify_geometry(true,atof(argv[5]));
            simplify_timer.stop();
             cout << "       - time simplifications:      " << simplify_timer.getElapsedTime() << endl;

        }
        else if(strcmp("-l",argv[4])==0){
        gradient->simplify_geometry(false, atof(argv[5]));  
        }
        gradient->write_mesh_VTK("simplified_mesh");

        }
    }

    if(output_required){
        //multiresolution version
        cout << "Printing output:" << endl;
       
        gradient->writeVTK_gradient("prove_gradiente_doporefine.vtk");

        gradient->descending_2cells_extraction(true);
        gradient->descending_1cells_extraction(true);

        gradient->ascending_2cells_extraction(true);
        gradient->ascending_1cells_extraction(true);

        gradient->compute_incidence_graph();
        gradient->writeVTK_IG("ig_simplified.vtk");

        cout << "   - done" << endl;
    }



//    if(true){

//        /* init GLUT */
//          int pargc = 1;
//          glutInit( &pargc, argv );

//          /* init Visualization Library */
//          VisualizationLibrary::init();

//          /* install Visualization Library shutdown function */
//          atexit( vlGLUT::atexit_visualization_library_shutdown );

//          /* setup the OpenGL context format */
//          OpenGLContextFormat format;
//          format.setDoubleBuffer(true);
//          format.setRGBABits( 8,8,8,8 );
//          format.setDepthBufferBits(24);
//          format.setStencilBufferBits(8);
//          format.setFullscreen(false);
////          format.setMultisampleSamples(8);
////          format.setMultisample(true);

//          ref<Effect> fx = new Effect;
//          fx->shader()->enable(EN_BLEND);
//          fx->shader()->enable(EN_POINT_SMOOTH);
//          fx->shader()->enable(EN_LINE_SMOOTH);
//          fx->shader()->enable(EN_POLYGON_SMOOTH);

//          /* create the applet to be run */
//          ref<Applet> applet = new Vis_TriangleMesh;
//          applet->initialize();
//          /* create a native GLUT window */
//          ref<vlGLUT::GLUTWindow> glut_window = new vlGLUT::GLUTWindow;
//          /* bind the applet so it receives all the GUI events related to the OpenGLContext */
//          glut_window->addEventListener(applet.get());
//          /* target the window so we can render on it */
//          applet->rendering()->as<Rendering>()->renderer()->setFramebuffer( glut_window->framebuffer() );
//          /* black background */
//          applet->rendering()->as<Rendering>()->camera()->viewport()->setClearColor( white );
//          /* define the camera position and orientation */
//          vec3 eye    = vec3(0,10,35); // camera position
//          vec3 center = vec3(0,0,0);   // point the camera is looking at
//          vec3 up     = vec3(0,1,0);   // up direction
//          mat4 view_mat = mat4::getLookAt(eye, center, up);
//          applet->rendering()->as<Rendering>()->camera()->setViewMatrix( view_mat );
//          /* Initialize the OpenGL context and window properties */
//          int x = 0;
//          int y = 0;
//          int width = 512;
//          int height= 512;
//          glut_window->initGLUTWindow( "Visualization Library on GLUT - Rotating Cube", format, x, y, width, height );

//          /* ... you can open more than one GLUT window! */

//          /* enter the GLUT main loop */
//          glutMainLoop();

//          /* this point is never reached since glutMainLoop() never returns! */

//          return 0;

//    }



    return 0;
}

void print_paragraph(string stringa, int cols){
    if(stringa.size() < cols-20){
        printf("          %s\n\n",stringa.c_str());
    }
    else{
        float dim = (float)(stringa.size()/((float)cols-20));
        int dim_int = dim;
        if(dim > dim_int) dim_int++;
        for(int i=0; i<dim_int; i++){
            if(stringa.at((cols-20)*i) == ' ')
                printf("         %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
            else printf("          %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
        }
        printf("\n");
    }
}

void print_help(){

  //annoying stuff to get the dimension of the output shell (!!! not sure it works on Mac,
  //everything based on the command tput. If it doesn't work the dimension si setted to 80 by default)
  FILE* fp;
  char path[1035];

  int cols;
  fp = popen("tput cols", "r");
  if(fp != NULL){
    fgets(path, sizeof(path)-1, fp);
    cols = atoi(path);
  }
  else{
    cols = 80;
  }

  //printf start
    printf(BOLD "\n  NAME:\n\n" RESET);
    printf("\tMulti-resolution Forman Graidnet \n\n" RESET);

    printf(BOLD "  USAGE: \n\n" RESET);
    printf(BOLD "    For the simplification algorithm" RESET);
    printf(BOLD "    ./MM_gradient -i [input file] -s threshold \n\n" RESET);
    print_paragraph("the [input file] in .tri format and the persistence threshold", cols);


    printf(BOLD "  IMPLEMENTATION:\n\n" RESET);
    printf("          Author: Ciccio\n");
    printf("          Group: G3 Geometry and Graphics Group\n");
    printf("          Last Update: 13/06/2014\n\n");

    printf(BOLD "  DESCRIPTION: \n\n" RESET);
    print_paragraph("Implementation of Forman Gradient simplification and multi-resolution model [based on edge collapses] for triangle meshes.", cols);
}


