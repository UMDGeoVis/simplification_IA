#include "Reader.h"
#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>

Reader::Reader() {}

Reader::Reader(const Reader& ) {}

Reader::~Reader() {}

bool Reader::readMeshFile(Mesh<Vertex2D, Triangle> &mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int num_vertices;
    input >> num_vertices;

    if (num_vertices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_Vertices(num_vertices);

    //legge i vertici aggiustando il dominio..
    for (int i = 0; i < num_vertices; i++) {
        double x, y;

        input >> x;
        input >> y;
        if (input.eof())
            break;

        Vertex2D v = Vertex2D(x, y);
        mesh.addVertex(v);
    }

    int num_topSimplexes;
    input >> num_topSimplexes;

    if(num_topSimplexes == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_TopSimplexes(num_topSimplexes);

    //legge i top simplessi
    for (int i = 0; i < num_topSimplexes; i++) {
        int v[3];
        for (int j = 0; j < 3; j++)
            input >> v[j];
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.addTopSimplex(t);
    }

    return true;
}

bool Reader::readMeshFile(Mesh<Vertex3D, Triangle> &mesh, string path)
{
    string extension = string_management::get_file_extension(path);
    if(extension == "tri")
        return Reader::read_mesh_tri(mesh,path);
    else if(extension == "off")
        return Reader::read_mesh_off(mesh,path);
    else
    {
        cerr << "[ERROR] unsopported file format. " << endl;
        return false;
    }
}


bool Reader::read_mesh_tri(Mesh<Vertex3D, Triangle> &mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    int num_vertices;
    input >> num_vertices;

    if (num_vertices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_Vertices(num_vertices);

    //legge i vertici aggiustando il dominio..
    for (int i = 0; i < num_vertices; i++) {
        double x, y, z;

        input >> x;
        input >> y;
        input >> z;
        if (input.eof())
            break;

        Vertex3D v = Vertex3D(x, y, z);
        mesh.addVertex(v);
    }

    int num_topSimplexes;
    input >> num_topSimplexes;

    if(num_topSimplexes == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_TopSimplexes(num_topSimplexes);

    //legge i top simplessi
    for (int i = 0; i < num_topSimplexes; i++) {
        int v[3];
        for (int j = 0; j < 3; j++)
            input >> v[j];
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.addTopSimplex(t);
    }

    return true;
}


bool Reader::read_mesh_off(Mesh<Vertex3D, Triangle> &mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string l;
    getline(input,l); // trow away the first line
    getline(input,l);
    vector<string> lt;
    string_management::tokenize(l,lt," ");
    int num_vertices = atol(lt.at(0).c_str());
    int num_triangles = atol(lt.at(1).c_str());

    if (num_vertices == 0 || num_triangles == 0)
    {
        cerr << "This is not a valid .off file: " << path << endl;
        return false;
    }

    mesh.reserveVectorSpace_Vertices(num_vertices);
    mesh.reserveVectorSpace_TopSimplexes(num_triangles);
        //legge i vertici aggiustando il dominio..
    for (int i = 0; i < num_vertices; i++) {
        double x, y, z;

        input >> x;
        input >> y;
        input >> z;
        if (input.eof())
            break;

        Vertex3D v = Vertex3D(x, y, z);
        mesh.addVertex(v);
    }
    
    int v[3], index;
    int num_v;
    for (int i = 0; i < num_triangles; i++)
    {
        input >> num_v;

        if(num_v != 3)
        {
            cerr << "[ERROR] the input mesh must be a pure triangle mesh. read a simplex with "<< num_v << "vertices." << endl;
            return false;
        }

        for (int j = 0; j < num_v; j++)
        {
            input >> index;
            v[j] = index;
        }
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.addTopSimplex(t);
    }

    return true;
}