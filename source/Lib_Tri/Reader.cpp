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
