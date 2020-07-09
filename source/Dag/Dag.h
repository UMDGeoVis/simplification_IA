#ifndef DAG_H
#define DAG_H

#include <set>
#include <vector>
#include <list>

#include "../Ig/ig.h"

using namespace std;

class DAG_Node{

protected:


    bool refined, visited;



    bool is_topo; //per discriminare fra operazioni topologiche e geometriche (per comodita', non indispensabile)

};


class DAG_GeomNode : public DAG_Node{


    //coppia di triangoli di sinistra
    int v_sinistro, v_sinistro_sopra;
    int v_destro, v_destro_sopra;

    //vertice centrale (quello che ricevera' la freccia del gradiente)
    int v1;
    bool to_switch_sin,to_switch_des;

    vector<int> vv;

    double edge_lenght;
    //coordinate del vecchio vertice;
    vector<double> coordsv2;
    vector<double> coordsv1;

    list<DAG_GeomNode*> dep;

public:

    inline ~DAG_GeomNode(){dep.clear();}
    inline DAG_GeomNode(int v_sin, int v_sin_up, int v_des, int v_des_up, int v1, vector<double> coordsv1,vector<double> coordsv2){
        this->v_sinistro = v_sin;
        this->v_destro = v_des;
        this->v_sinistro_sopra = v_sin_up;
        this->v_destro_sopra = v_des_up;

        this->v1 = v1;
        this->coordsv1 = coordsv1;
        this->coordsv2 = coordsv2;

        refined = visited = is_topo = false;
    }

    inline pair<bool,bool> get_to_switch(){return pair<bool,bool>(to_switch_sin, to_switch_des);}
    inline bool set_to_switch(bool to_switch_sin,bool to_switch_des){this->to_switch_sin=to_switch_sin; this->to_switch_des=to_switch_des;}
    inline void set_edge_lenght(double edg){edge_lenght=edg;}

    inline void add_vv(vector<int> vv){this->vv=vv;}
    inline vector<int> getVV(){return vv;}

    inline void add_dep(DAG_GeomNode* node){ dep.push_back(node);}
    inline vector<DAG_GeomNode*> get_dep(){return vector<DAG_GeomNode*>(dep.begin(), dep.end());}

    inline void setV1(int v1){this->v1=v1;}//una volta che l'operazione e' stata eseguita v1 memorizza il vertice inserito. Prima memorizzava il vertice sopravvissuto
    inline int inserted_vertex(){return v1;}

    inline bool isRefined(){return refined;}
    inline void setRefined(){refined = true;}

    inline bool isVisited(){return visited;}
    inline void setVisited(){visited=true;}
    inline void setNotVisited(){visited=false;}

    inline double getEdgeLenght(){return edge_lenght;}

    inline vector<int> getVertices(){
        vector<int> vert(5);
        vert[0]=(v_sinistro);
        vert[1]=(v_sinistro_sopra);
        vert[2]=(v_destro);
        vert[3]=(v_destro_sopra);
        vert[4]=(v1);

        return vert;
    }

    inline vector<double> getCoordsv1(){return coordsv1;}
    inline vector<double> getCoordsv2(){return coordsv2;}

    inline vector<int> getUpdatedVertices(vector<DAG_GeomNode*>* dag_per_vertex){
        vector<int> vert(5);
        if((*dag_per_vertex)[v_sinistro] != NULL){
            //assert((*dag_per_vertex)[v_sinistro]->isRefined());
            vert[0]=(*dag_per_vertex)[v_sinistro]->inserted_vertex();
        }
        else{
            vert[0] = -v_sinistro;
        }

        if((*dag_per_vertex)[v_sinistro_sopra] != NULL){
            //assert((*dag_per_vertex)[v_sinistro_sopra]->isRefined());
            vert[1]=((*dag_per_vertex)[v_sinistro_sopra]->inserted_vertex());
        }
        else{
            vert[1]=(-v_sinistro_sopra);
        }

        if((*dag_per_vertex)[v_destro] != NULL){
            //assert((*dag_per_vertex)[v_destro]->isRefined());
            vert[2]=((*dag_per_vertex)[v_destro]->inserted_vertex());
        }
        else{
            vert[2]=(-v_destro);
        }

        if((*dag_per_vertex)[v_destro_sopra] != NULL){
            //assert((*dag_per_vertex)[v_destro_sopra]->isRefined());
            vert[3]=((*dag_per_vertex)[v_destro_sopra]->inserted_vertex());
        }
        else{
            vert[3]=(-v_destro_sopra);
        }

        if((*dag_per_vertex)[v1] != NULL){
            //assert((*dag_per_vertex)[v1]->isRefined());
            vert[4]=((*dag_per_vertex)[v1]->inserted_vertex());
        }
        else{
            vert[4]=(-v1);
        }

        return vert;
    }

    inline void update_vertex_index(int j, int v_index){
        switch(j){
            case 0: v_sinistro = v_index; break;
            case 1: v_sinistro_sopra = v_index; break;
            case 2: v_destro = v_index; break;
            case 3: v_destro_sopra = v_index; break;
            case 4: v1 = v_index; break;
            default: cout << "Error" << endl; break;

        }
    }

    inline bool isRefinable(vector<DAG_GeomNode*>* dag_per_vertex){

        vector<int> vert = this->getVV();
        for(int i=0; i<vert.size(); i++){
            if(((*dag_per_vertex)[vert[i]] != NULL && !(*dag_per_vertex)[vert[i]]->isRefined()))
                return false;
        }
        return true;
    }



};

class DAG_TopoNode : public DAG_Node{

    bool is_con;
    int label_pp;

    double height;

    Node* q;
    Node* p;
    Node* p1;

    list<int> to_p;
    list<int> to_q;

    list<iNode*> set_p; //nodi che vanno riconnessi con p (di indice > di p)
    list<nNode*> set_q; //nodi che vanno riconnessi con q (di indice < di q)

    //parte geometrica
    vector<int> simplex_vert;
    pair<int,int> edge; //edge da raggiungere nel cammino

    int start_vertex;

    set<DAG_TopoNode*> parents;
    set<DAG_TopoNode*> children;
    set<DAG_GeomNode*> parents_geom;


public:

    inline ~DAG_TopoNode(){children.clear(); parents.clear();}
    inline DAG_TopoNode(){
        refined = false;
        visited =false;
        height = 0;
    }

    inline DAG_TopoNode(nNode* p, iNode* q, nNode* p1,
                        list<iNode*> set_p, list<nNode*> set_q,
                        list<int> to_p, list<int> to_q,
                        int label_p1,
                        bool is_con,
                        vector<int> simplex,
                        pair<int,int> edge,
                        int start_vertex){

        this->p = p;
        this->p1 = p1;
        this->q = q;

        this->set_p = set_p;
        this->set_q = set_q;

        this->to_p = to_p;
        this->to_q = to_q;

        this->is_con = is_con;

        label_pp = label_p1;
        refined =false;
        visited =false;

        children = set<DAG_TopoNode*>();
        parents = set<DAG_TopoNode*>();

        height = 0;

        simplex_vert=simplex;
        this->edge=edge;

        is_topo = true;
        this->start_vertex=start_vertex;
    }

    inline vector<int> getVertices(){
        vector<int> ret;

        for(int i=0; i<simplex_vert.size(); i++)
            ret.push_back(simplex_vert[i]);

        ret.push_back(edge.first);
        ret.push_back(edge.second);
        ret.push_back(start_vertex);
        return ret;
    }

    inline int get_dim_q(){return 1;}

    inline void add_child(DAG_TopoNode* node){children.insert(node);}
    inline void add_father(DAG_TopoNode* node){parents.insert(node);}
    inline void add_geom(DAG_GeomNode* node){parents_geom.insert(node);}

    inline set<DAG_TopoNode*>* getChild(){return &children;}
    inline set<DAG_TopoNode*>* getFathers(){return &parents;}
    inline set<DAG_GeomNode*>* getFathersGeom(){return &parents_geom;}

    inline Node * get_q(){return q;}
    inline Node * get_p(){return p;}
    inline Node * get_p1(){return p1;}

    inline list<iNode*>* get_to_p(){return &set_p;}
    inline list<nNode*>* get_to_q(){return &set_q;}

    inline int get_label_pp(){return label_pp;}

    inline list<int>* get_label_to_p(){return &to_p;}
    inline list<int>* get_label_to_q(){return &to_q;}

    inline void setHeight(double h){height = h;}
    inline double getHeight(){return height;}

    inline bool isRefined(){return refined;}
    inline void setRefined(){refined = true;}

    inline bool isVisited(){return visited;}
    inline void setVisited(){visited=true;}
    inline void setNotVisited(){visited=false;}

    inline list<nNode*> getMaximaList(){
        list<nNode*> lista;
        if(is_con){
            lista.insert(lista.end(), set_q.begin(), set_q.end());
        }
        else{
            lista.push_back((nNode*)p1);
        }
        return lista;
    }

    inline list<iNode*> getSaddleList(){
        return list<iNode*>(set_p.begin(), set_p.end());
    }


    inline list<nNode*> getMinimaList(){
        list<nNode*> lista;
        if(is_con){
            lista.push_back((nNode*)p1);
        }
        else{
            lista.insert(lista.end(), set_q.begin(), set_q.end());
        }
        return lista;
    }


    inline bool isExpansion(){
        return this->is_con;
    }

    inline bool isRefinable(){
        for(set<DAG_TopoNode*>::iterator it = parents.begin(); it != parents.end(); it++){
            if(!(*it)->isRefined())
                return false;
        }

        for(set<DAG_GeomNode*>::iterator it = parents_geom.begin(); it != parents_geom.end(); it++){
            if(!(*it)->isRefined())
                return false;
        }

        return true;
    }

    inline vector<DAG_TopoNode*> get_dep(){
        return vector<DAG_TopoNode*>(children.begin(), children.end());
    }

    inline vector<int> getUpdatedVertices(vector<DAG_GeomNode*>* dag_per_vertex){
        vector<int> vert= vector<int>();

        for(int i=0; i<simplex_vert.size(); i++){
            if((*dag_per_vertex)[simplex_vert[i]] != NULL){
                //assert((*dag_per_vertex)[simplex_vert[i]]->isRefined());
                vert.push_back((*dag_per_vertex)[simplex_vert[i]]->inserted_vertex());
            }
            else{
                vert.push_back(-simplex_vert[i]);
            }
        }

        if((*dag_per_vertex)[edge.first] != NULL){
            //assert((*dag_per_vertex)[edge.first]->isRefined());
            vert.push_back((*dag_per_vertex)[edge.first]->inserted_vertex());
        }
        else{
            vert.push_back(-edge.first);
        }

        if((*dag_per_vertex)[edge.second] != NULL){
            //assert((*dag_per_vertex)[edge.second]->isRefined());
            vert.push_back((*dag_per_vertex)[edge.second]->inserted_vertex());
        }
        else{
            vert.push_back(-edge.second);
        }

        if((*dag_per_vertex)[start_vertex] != NULL){
            //assert((*dag_per_vertex)[start_vertex]->isRefined());
            vert.push_back((*dag_per_vertex)[start_vertex]->inserted_vertex());
        }
        else{
            vert.push_back(-start_vertex);
        }

        return vert;
    }

};




#endif // DAG_H
