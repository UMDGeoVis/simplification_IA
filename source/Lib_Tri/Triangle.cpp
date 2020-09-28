#include "Triangle.h"

Triangle::Triangle()
{
}

Triangle::Triangle(int v1, int v2, int v3)
{
    this->vertices[0] = v1;
    this->vertices[1] = v2;
    this->vertices[2] = v3;
}

int Triangle::TV(int pos)
{
    return this->vertices[pos];
}

Edge* Triangle::TE(int pos)
{
   int v1=(vertices[(pos+1)%3]<vertices[(pos+2)%3])?vertices[(pos+1)%3]:vertices[(pos+2)%3];
    int v2=(vertices[(pos+1)%3]<vertices[(pos+2)%3])?vertices[(pos+2)%3]:vertices[(pos+1)%3];
    return new Edge(v1,v2);
}

int Triangle::TT(int pos)
{
    return this->adj[pos];
}

void Triangle::setTT(int pos, int adjId)
{
    this->adj[pos]=adjId;
}

int Triangle::getVerticesNum()
{
    return 3;
}
