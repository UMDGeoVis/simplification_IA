#include "Mesh.h"
#include "Vertex3D.h"
#include "../Eigen/Eigen"
#include "../Eigen/SVD"
#include "../Eigen/Dense"
#include "../Eigen/LU"
#include "../Dag/Dag.h"


template<class V, class T>
DAG_GeomNode* Mesh<V,T>::half_edge_collapse(int v1, int v2, int t1, int t2, vector<double> new_v){

    DAG_GeomNode* nodegeom=NULL;

    removed_vertex[v2]=true;
    removed_triangle[t1]=true;
    removed_triangle[t2]=true;

    vector<int> vt = VT(v2);
    vector<int> vv = VV(v2);

    int v_sinistro, v_destro, v_sinistro_sopra, v_destro_sopra;
    v_sinistro = v_destro = v_sinistro_sopra = v_destro_sopra = -1;

    //modifiche geometriche per t1
    int t3_sin = getTopSimplex(t1).TT(getTopSimplex(t1).vertex_index(v1));
    int t3_adj_sin = getTopSimplex(t1).TT(getTopSimplex(t1).vertex_index(v2));


    for(int i=0; i<3; i++){
        if(getTopSimplex(t1).TV(i) != v1 && getTopSimplex(t1).TV(i) != v2){
            v_sinistro = getTopSimplex(t1).TV(i);
            break;
        }
    }

    int v5_sin=-1;
    for(int i=0; i<3; i++){
        if(getTopSimplex(t3_sin).TV(i) != v2 && getTopSimplex(t3_sin).TV(i) != v_sinistro){
            v_sinistro_sopra = getTopSimplex(t3_sin).TV(i);
            break;
        }
    }

    for(int i=0; i<3; i++){
        if(getTopSimplex(t3_adj_sin).TV(i) != v1 && getTopSimplex(t3_adj_sin).TV(i) != v_sinistro){
            v5_sin = getTopSimplex(t3_adj_sin).TV(i);
            break;
        }
    }

    int t3_des = getTopSimplex(t2).TT(getTopSimplex(t2).vertex_index(v1));
    int t3_adj_des = getTopSimplex(t2).TT(getTopSimplex(t2).vertex_index(v2));

    for(int i=0; i<3; i++){
        if(getTopSimplex(t2).TV(i) != v1 && getTopSimplex(t2).TV(i) != v2){
            v_destro = getTopSimplex(t2).TV(i);
            break;
        }
    }

    for(int i=0; i<3; i++){
        if(getTopSimplex(t3_des).TV(i) != v2 && getTopSimplex(t3_des).TV(i) != v_destro){
            v_destro_sopra = getTopSimplex(t3_des).TV(i);
            break;
        }
    }

    int v5_des=-1;
    for(int i=0; i<3; i++){
        if(getTopSimplex(t3_adj_des).TV(i) != v1 && getTopSimplex(t3_adj_des).TV(i) != v_destro){
            v5_des = getTopSimplex(t3_adj_des).TV(i);
            break;
        }
    }

    //assert(t3_sin != t3_des);
    //assert(t3_adj_sin != t3_adj_des);

    /*t3 adiacente a t1 e t2 si becca gli adiacenti di t1 e t2 dal corrispettivo lato*/
    getTopSimplex(t3_sin).setTT(getTopSimplex(t3_sin).vertex_index(v_sinistro_sopra), t3_adj_sin);
    getTopSimplex(t3_adj_sin).setTT(getTopSimplex(t3_adj_sin).vertex_index(v5_sin), t3_sin);

    getTopSimplex(t3_des).setTT(getTopSimplex(t3_des).vertex_index(v_destro_sopra), t3_adj_des);
    getTopSimplex(t3_adj_des).setTT(getTopSimplex(t3_adj_des).vertex_index(v5_des), t3_des);


    for(int i=0; i<vt.size(); i++){
        if(vt[i] != t1 && vt[i] != t2){
            //assert(getTopSimplex(vt[i]).contains(v2) && !getTopSimplex(vt[i]).contains(v1));
            getTopSimplex(vt[i]).setTV(getTopSimplex(vt[i]).vertex_index(v2),v1);
        }
    }

    if(getVertex(v1).VTstar() == t1){
        getVertex(v1).VTstar(t3_sin);
    }

    if(getVertex(v1).VTstar() == t2){
        getVertex(v1).VTstar(t3_des);
    }


    if(getVertex(v_sinistro).VTstar() == t1){
        getVertex(v_sinistro).VTstar(t3_sin);
    }

    if(getVertex(v_destro).VTstar() == t2){
        getVertex(v_destro).VTstar(t3_des);
    }



    vector<double> coordsv2(3);
    coordsv2[0] = (getVertex(v2).getX());
    coordsv2[1] = (getVertex(v2).getY());
    coordsv2[2] = (getVertex(v2).getZ());

    vector<double> coordsv1(3);
    coordsv1[0] = (getVertex(v1).getX());
    coordsv1[1] = (getVertex(v1).getY());
    coordsv1[2] = (getVertex(v1).getZ());


    getVertex(v1).setX((new_v)[0]);
    getVertex(v1).setY((new_v)[1]);
    getVertex(v1).setZ((new_v)[2]);

    nodegeom = new DAG_GeomNode(v_sinistro, v_sinistro_sopra, v_destro, v_destro_sopra, v1, coordsv1, coordsv2);


    nodegeom->add_vv(vv);

    return nodegeom;
}


template<class V, class T>
bool Mesh<V,T>::convex_neighborhood(int v1, int v2, int t1, int t2){

    //first check for t1
    int v3,v4; v3=v4=-1;
    for(int i=0; i<3; i++){
        if(getTopSimplex(t1).TV(i) != v1 && getTopSimplex(t1).TV(i) != v2){
            v3 = getTopSimplex(t1).TV(i);
            break;
        }
    }
    //assert(v3 != -1);

    int t3 = getTopSimplex(t1).TT(getTopSimplex(t1).vertex_index(v1));

    if(t3 == -1) return false;

    for(int i=0; i<3; i++){
        if(getTopSimplex(t3).TV(i) != v2 && getTopSimplex(t3).TV(i) != v3){
            v4 = getTopSimplex(t1).TV(i);
            break;
        }
    }
    //assert(v4 != -1);

    //svolta di v1v3v4 == svolta v1v3v2
    Eigen::Matrix3f v1v2v3;
    Eigen::Matrix3f v1v3v4;

    Vertex3D V1 = getVertex(v1);
    Vertex3D V2 = getVertex(v2);
    Vertex3D V3 = getVertex(v3);
    Vertex3D V4 = getVertex(v4);

    v1v2v3(0,0) = V1.getZ() - V3.getZ();
    v1v2v3(1,0) = V1.getY() - V3.getY();
    v1v2v3(2,0) = V1.getX() - V3.getX();

    v1v2v3(0,1) = V2.getZ() - V3.getZ();
    v1v2v3(1,1) = V2.getY() - V3.getY();
    v1v2v3(2,1) = V2.getX() - V3.getX();

    v1v2v3(0,2) = 1;
    v1v2v3(1,2) = 1;
    v1v2v3(2,2) = 1;

    v1v3v4(0,0) = V1.getZ() - V3.getZ();
    v1v3v4(1,0) = V1.getY() - V3.getY();
    v1v3v4(2,0) = V1.getX() - V3.getX();

    v1v3v4(0,1) = V4.getZ() - V3.getZ();
    v1v3v4(1,1) = V4.getY() - V3.getY();
    v1v3v4(2,1) = V4.getX() - V3.getX();

    v1v3v4(0,2) = 1;
    v1v3v4(1,2) = 1;
    v1v3v4(2,2) = 1;

    if(((v1v2v3.determinant()) >0) != ((v1v3v4.determinant())>0) ) return false;


    //stessa cosa con altro tetraedro
    //------------------------------------------
    for(int i=0; i<3; i++){
        if(getTopSimplex(t2).TV(i) != v1 && getTopSimplex(t2).TV(i) != v2){
            v3 = getTopSimplex(t2).TV(i);
            break;
        }
    }

    t3 = getTopSimplex(t2).TT(getTopSimplex(t2).vertex_index(v1));

    if(t3 == -1) return false;

    for(int i=0; i<3; i++){
        if(getTopSimplex(t3).TV(i) != v2 && getTopSimplex(t3).TV(i) != v3){
            v4 = getTopSimplex(t1).TV(i);
            break;
        }
    }

    //svolta di v1v3v4 == svolta v1v3v2
    V3 = getVertex(v3);
    V4 = getVertex(v4);

    v1v2v3(0,0) = V1.getZ() - V3.getZ();
    v1v2v3(1,0) = V1.getY() - V3.getY();
    v1v2v3(2,0) = V1.getX() - V3.getX();

    v1v2v3(0,1) = V2.getZ() - V3.getZ();
    v1v2v3(1,1) = V2.getY() - V3.getY();
    v1v2v3(2,1) = V2.getX() - V3.getX();

    v1v2v3(0,2) = 1;
    v1v2v3(1,2) = 1;
    v1v2v3(2,2) = 1;

    v1v3v4(0,0) = V1.getZ() - V3.getZ();
    v1v3v4(1,0) = V1.getY() - V3.getY();
    v1v3v4(2,0) = V1.getX() - V3.getX();

    v1v3v4(0,1) = V4.getZ() - V3.getZ();
    v1v3v4(1,1) = V4.getY() - V3.getY();
    v1v3v4(2,1) = V4.getX() - V3.getX();

    v1v3v4(0,2) = 1;
    v1v3v4(1,2) = 1;
    v1v3v4(2,2) = 1;

    if(((v1v2v3.determinant()) >0) != ((v1v3v4.determinant())>0) ) return false;

    return true;
}

