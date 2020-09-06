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

        //Set the TV of VT(v2) to replace v2 with v1
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

 template<class V, class T> void Mesh<V,T>::half_edge_collapse_simple(int v1, int v2, int t1, int t2, vector<double> new_v,priority_queue<Geom_Sempl*, vector<Geom_Sempl*>, sort_arcs_geom>& queue,double limit){
    cout<<"[EDGE CONTRACTION] v1 and v2:"<<v1<<", "<<v2<<endl;
    removed_vertex[v2]=true;
    removed_triangle[t1]=true;
    if(t2!=-1)
    removed_triangle[t2]=true;
   cout<<"=====Extract relations====="<<endl;
    vector<int> vt = VT(v2);
    vector<int> vv = VV(v2);
 

    int v_sinistro, v_destro, v_sinistro_sopra, v_destro_sopra;
    v_sinistro = v_destro = v_sinistro_sopra = v_destro_sopra = -1;
  cout<<"=====Extract neighboring triangles====="<<endl;
    //modifiche geometriche per t1
    int t3_sin = getTopSimplex(t1).TT(getTopSimplex(t1).vertex_index(v1));
    int t3_adj_sin = getTopSimplex(t1).TT(getTopSimplex(t1).vertex_index(v2));
   cout<<"t3 sin and t3 adj sin: "<<t3_sin<<"; "<<t3_adj_sin<<endl;

    for(int i=0; i<3; i++){
        if(getTopSimplex(t1).TV(i) != v1 && getTopSimplex(t1).TV(i) != v2){
            v_sinistro = getTopSimplex(t1).TV(i);
            break;
        }
    }
   cout<<"v_sinistro: "<<v_sinistro<<endl;
    int v5_sin=-1;
    for(int i=0; i<3; i++){
        if(t3_sin!=-1 && getTopSimplex(t3_sin).TV(i) != v2 && getTopSimplex(t3_sin).TV(i) != v_sinistro){
            v_sinistro_sopra = getTopSimplex(t3_sin).TV(i);
            break;
        }
    }
  cout<<"v_sinistro_sopra: "<<v_sinistro_sopra<<endl;
    for(int i=0; i<3; i++){
        if(t3_adj_sin!=-1&&getTopSimplex(t3_adj_sin).TV(i) != v1 && getTopSimplex(t3_adj_sin).TV(i) != v_sinistro){
            v5_sin = getTopSimplex(t3_adj_sin).TV(i);
            break;
        }
    }
   cout<<"v5_sin: "<<v5_sin<<endl;
   int t3_des=-1,t3_adj_des=-1;
    int v5_des=-1;
    if(t2!=-1){
    t3_des = getTopSimplex(t2).TT(getTopSimplex(t2).vertex_index(v1));
    t3_adj_des = getTopSimplex(t2).TT(getTopSimplex(t2).vertex_index(v2));
   cout<<"t3 des and t3 adj des: "<<t3_des<<"; "<<t3_adj_des<<endl;
    for(int i=0; i<3; i++){
        if( getTopSimplex(t2).TV(i) != v1 && getTopSimplex(t2).TV(i) != v2){
            v_destro = getTopSimplex(t2).TV(i);
            break;
        }
    }
 cout<<"v_destro: "<<v_destro<<endl;
    for(int i=0; i<3; i++){
        if(t3_des != -1 && getTopSimplex(t3_des).TV(i) != v2 && getTopSimplex(t3_des).TV(i) != v_destro){
            v_destro_sopra = getTopSimplex(t3_des).TV(i);
            break;
        }
    }
   cout<<"v_destro_sopra: "<<v_destro_sopra<<endl;
   
    for(int i=0; i<3; i++){
        if(t3_adj_des != -1 &&getTopSimplex(t3_adj_des).TV(i) != v1 && getTopSimplex(t3_adj_des).TV(i) != v_destro){
            v5_des = getTopSimplex(t3_adj_des).TV(i);
            break;
        }
    }
   cout<<"v5_des: "<<v5_des <<endl;
    //assert(t3_sin != t3_des);
    //assert(t3_adj_sin != t3_adj_des);
//cout<<"Set TT relations"<<endl;
    /*t3 adiacente a t1 e t2 si becca gli adiacenti di t1 e t2 dal corrispettivo lato*/
    }
    if(t3_sin!=-1)
        getTopSimplex(t3_sin).setTT(getTopSimplex(t3_sin).vertex_index(v_sinistro_sopra), t3_adj_sin);
    if(t3_adj_sin!=-1)
        getTopSimplex(t3_adj_sin).setTT(getTopSimplex(t3_adj_sin).vertex_index(v5_sin), t3_sin);

    if(t3_des!=-1)
    getTopSimplex(t3_des).setTT(getTopSimplex(t3_des).vertex_index(v_destro_sopra), t3_adj_des);

    if(t3_adj_des!=-1)
    getTopSimplex(t3_adj_des).setTT(getTopSimplex(t3_adj_des).vertex_index(v5_des), t3_des);
    
        //Set the TV of VT(v2) to replace v2 with v1
    for(int i=0; i<vt.size(); i++){
        if(vt[i] != t1 && vt[i] != t2){
            //assert(getTopSimplex(vt[i]).contains(v2) && !getTopSimplex(vt[i]).contains(v1));
            int v2_id=getTopSimplex(vt[i]).vertex_index(v2);

            getTopSimplex(vt[i]).setTV(v2_id,v1);    }
    }
   //  cout<<"Update New edges.  VV size:"<<vv.size()<<endl;
   
    for(int i=0; i<vv.size(); i++){
       // cout<<vv[i]<<endl;
        if(vv[i] != v1 && vv[i] != v_sinistro && vv[i] != v_destro){
            vector<double> new_vertex(3,0);
            
            Vertex3D va=getVertex(v1);
            Vertex3D vb=getVertex(vv[i]);
            
            //assert(getTopSimplex(vt[i]).contains(v2) && !getTopSimplex(vt[i]).contains(v1));
            Edge * e=new Edge(v1,vv[i]);
           // if(is_QEM)
            //if(va.getZ()>vb.getZ())
            { new_vertex[0] = va.getX(); new_vertex[1] = va.getY(), new_vertex[2] = va.getZ(); }
          //  else
            //{ new_vertex[0] = vb.getX(); new_vertex[1] = vb.getY(), new_vertex[2] = vb.getZ(); }

            vector<double> dif = {va.getX()-vb.getX(),va.getY()-vb.getY(),va.getZ()-vb.getZ()};
            double length = sqrt(dif[0]*dif[0]+dif[1]*dif[1]+dif[2]*dif[2]);
            if(length<limit){
            queue.push(new Geom_Sempl(e,length,new_vertex));
       //     cout<<"New edge updated"<<endl;
            }
              }
    }

    if(getVertex(v1).VTstar() == t1||getVertex(v1).VTstar() == t2){
        if(t3_sin!=-1)
        getVertex(v1).VTstar(t3_sin);
        else if(t3_adj_sin!=-1)
        getVertex(v1).VTstar(t3_adj_sin);
        else if(t3_des!=-1)
        getVertex(v1).VTstar(t3_des);
        else if(t3_adj_des!=-1)
        getVertex(v1).VTstar(t3_adj_des);
    }

    // if(getVertex(v1).VTstar() == t2){
    //     if(t3_des!=-1)
    //     getVertex(v1).VTstar(t3_des);
    //     else if(t3_adj_des!=-1)
    //     getVertex(v1).VTstar(t3_adj_des);
    // }

    if(v_sinistro!=-1)
    if(getVertex(v_sinistro).VTstar() == t1){
        if(t3_sin!=-1)
        getVertex(v_sinistro).VTstar(t3_sin);
        else if(t3_adj_sin!=-1)
        getVertex(v_sinistro).VTstar(t3_adj_sin);
    }
    if(v_destro!=-1)
    if(getVertex(v_destro).VTstar() == t2){
        if(t3_des!=-1)
        getVertex(v_destro).VTstar(t3_des);
        else if(t3_adj_des!=-1)
        getVertex(v_destro).VTstar(t3_adj_des);
    }


  cout<<"FINISHED EDGE CONTRACTION"<<endl;

    //return 0;
    }

 template<class V, class T> void Mesh<V,T>::half_edge_collapse_QEM(int v1, int v2, int t1, int t2, vector<double> new_v,priority_queue<Geom_Sempl*, vector<Geom_Sempl*>, sort_arcs_geom>& queue,double limit,vector<Matrix>* vQEM){
    cout<<"[EDGE CONTRACTION] v1 and v2:"<<v1<<", "<<v2<<endl;
    removed_vertex[v2]=true;
    removed_triangle[t1]=true;
    if(t2!=-1)
    removed_triangle[t2]=true;
   cout<<"=====Extract relations====="<<endl;
    vector<int> vt = VT(v2);
    vector<int> vv = VV(v2);

    int v_sinistro, v_destro, v_sinistro_sopra, v_destro_sopra;
    v_sinistro = v_destro = v_sinistro_sopra = v_destro_sopra = -1;
  cout<<"=====Extract neighboring triangles====="<<endl;
    //modifiche geometriche per t1
    int t3_sin = getTopSimplex(t1).TT(getTopSimplex(t1).vertex_index(v1));
    int t3_adj_sin = getTopSimplex(t1).TT(getTopSimplex(t1).vertex_index(v2));
   cout<<"t3 sin and t3 adj sin: "<<t3_sin<<"; "<<t3_adj_sin<<endl;

    for(int i=0; i<3; i++){
        if(getTopSimplex(t1).TV(i) != v1 && getTopSimplex(t1).TV(i) != v2){
            v_sinistro = getTopSimplex(t1).TV(i);
            break;
        }
    }
   cout<<"v_sinistro: "<<v_sinistro<<endl;
    int v5_sin=-1;
    for(int i=0; i<3; i++){
        if(t3_sin!=-1 && getTopSimplex(t3_sin).TV(i) != v2 && getTopSimplex(t3_sin).TV(i) != v_sinistro){
            v_sinistro_sopra = getTopSimplex(t3_sin).TV(i);
            break;
        }
    }
  cout<<"v_sinistro_sopra: "<<v_sinistro_sopra<<endl;
    for(int i=0; i<3; i++){
        if(t3_adj_sin!=-1&&getTopSimplex(t3_adj_sin).TV(i) != v1 && getTopSimplex(t3_adj_sin).TV(i) != v_sinistro){
            v5_sin = getTopSimplex(t3_adj_sin).TV(i);
            break;
        }
    }
   cout<<"v5_sin: "<<v5_sin<<endl;
   int t3_des=-1,t3_adj_des=-1;
    int v5_des=-1;
    if(t2!=-1){
    t3_des = getTopSimplex(t2).TT(getTopSimplex(t2).vertex_index(v1));
    t3_adj_des = getTopSimplex(t2).TT(getTopSimplex(t2).vertex_index(v2));
   cout<<"t3 des and t3 adj des: "<<t3_des<<"; "<<t3_adj_des<<endl;
    for(int i=0; i<3; i++){
        if( getTopSimplex(t2).TV(i) != v1 && getTopSimplex(t2).TV(i) != v2){
            v_destro = getTopSimplex(t2).TV(i);
            break;
        }
    }
 cout<<"v_destro: "<<v_destro<<endl;
    for(int i=0; i<3; i++){
        if(t3_des != -1 && getTopSimplex(t3_des).TV(i) != v2 && getTopSimplex(t3_des).TV(i) != v_destro){
            v_destro_sopra = getTopSimplex(t3_des).TV(i);
            break;
        }
    }
   cout<<"v_destro_sopra: "<<v_destro_sopra<<endl;
   
    for(int i=0; i<3; i++){
        if(t3_adj_des != -1 &&getTopSimplex(t3_adj_des).TV(i) != v1 && getTopSimplex(t3_adj_des).TV(i) != v_destro){
            v5_des = getTopSimplex(t3_adj_des).TV(i);
            break;
        }
    }
   cout<<"v5_des: "<<v5_des <<endl;
    //assert(t3_sin != t3_des);
    //assert(t3_adj_sin != t3_adj_des);
//cout<<"Set TT relations"<<endl;
    /*t3 adiacente a t1 e t2 si becca gli adiacenti di t1 e t2 dal corrispettivo lato*/
    }
    if(t3_sin!=-1)
        getTopSimplex(t3_sin).setTT(getTopSimplex(t3_sin).vertex_index(v_sinistro_sopra), t3_adj_sin);
    if(t3_adj_sin!=-1)
        getTopSimplex(t3_adj_sin).setTT(getTopSimplex(t3_adj_sin).vertex_index(v5_sin), t3_sin);

    if(t3_des!=-1)
    getTopSimplex(t3_des).setTT(getTopSimplex(t3_des).vertex_index(v_destro_sopra), t3_adj_des);

    if(t3_adj_des!=-1)
    getTopSimplex(t3_adj_des).setTT(getTopSimplex(t3_adj_des).vertex_index(v5_des), t3_des);
    
        //Set the TV of VT(v2) to replace v2 with v1
    for(int i=0; i<vt.size(); i++){
        if(vt[i] != t1 && vt[i] != t2){
            //assert(getTopSimplex(vt[i]).contains(v2) && !getTopSimplex(vt[i]).contains(v1));
            int v2_id=getTopSimplex(vt[i]).vertex_index(v2);

            getTopSimplex(vt[i]).setTV(v2_id,v1);    }
    }
   //  cout<<"Update New edges.  VV size:"<<vv.size()<<endl;
   ///TODO： Update the QEM before computing the error 


    for(int i=0; i<vv.size(); i++){
       // cout<<vv[i]<<endl;
        if(vv[i] != v1 && vv[i] != v_sinistro && vv[i] != v_destro){
            vector<double> new_vertex(3,0);
            int new_vertex_pos=-1;
            Edge * e=new Edge(v1,vv[i]);
            double error = compute_error(v1,vv[i],vQEM,new_vertex_pos);
            assert(new_vertex_pos!=-1);
            if(new_vertex_pos==1)
                {
                *e= Edge(vv[i],v1); 
                }

            if(error<limit){
            cout<<"["<<e->EV(0)<<","<<e->EV(1)<<"]  Error will be introduced: "<<error<<endl; 
            queue.push(new Geom_Sempl(e,error,new_vertex));
       //     cout<<"New edge updated"<<endl;
            }
              }
    }

    if(getVertex(v1).VTstar() == t1||getVertex(v1).VTstar() == t2){
        if(t3_sin!=-1)
        getVertex(v1).VTstar(t3_sin);
        else if(t3_adj_sin!=-1)
        getVertex(v1).VTstar(t3_adj_sin);
        else if(t3_des!=-1)
        getVertex(v1).VTstar(t3_des);
        else if(t3_adj_des!=-1)
        getVertex(v1).VTstar(t3_adj_des);
    }


    if(v_sinistro!=-1)
    if(getVertex(v_sinistro).VTstar() == t1){
        if(t3_sin!=-1)
        getVertex(v_sinistro).VTstar(t3_sin);
        else if(t3_adj_sin!=-1)
        getVertex(v_sinistro).VTstar(t3_adj_sin);
    }
    if(v_destro!=-1)
    if(getVertex(v_destro).VTstar() == t2){
        if(t3_des!=-1)
        getVertex(v_destro).VTstar(t3_des);
        else if(t3_adj_des!=-1)
        getVertex(v_destro).VTstar(t3_adj_des);
    }
    getVertex(v1).setX((new_v)[0]);
    getVertex(v1).setY((new_v)[1]);
    getVertex(v1).setZ((new_v)[2]);

  cout<<"FINISHED EDGE CONTRACTION"<<endl;

    //return 0;
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

