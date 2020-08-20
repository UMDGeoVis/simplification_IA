#ifndef FORMANGRADIENTVECTOR_H
#define FORMANGRADIENTVECTOR_H

#include "./Reader.h"
#include "forman_arrow.h"
#include "aux_forman_gradient_vector.h"
#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>
#include "assert.h"

#include "Mesh.h"
#include "Vertex3D.h"
#include "../Ig/ig.h"
#include "../Dag/Dag.h"

#define NON_VALID_CONF 1000

struct Simplex_Graph{
    vector<int> simplex;
    vector<Simplex_Graph*> coboundary;
};

struct Geom_Sempl{

public:
    Edge* edge;
    double val;
    vector<double> new_v;

    inline Geom_Sempl(Edge* edge, double val, vector<double> new_v){ this->edge = edge; this->val=val; this->new_v = new_v;}
};

struct sort_arcs_geom{
    bool operator()(Geom_Sempl* a1, Geom_Sempl* a2){

        return a1->val > a2->val;

    }
};


struct Topo_Sempl{

public:
    Arc* arc;
    double val;
    int lvl;
    int filt_s0,filt_s1;
    vector<int> filt_ex;

    inline Topo_Sempl(Arc* arc, double val, int lvl,int filt_s0,int filt_s1,vector<int> filt_ex){ this->arc = arc; this->val=val; this->lvl=lvl;this->filt_s0=filt_s0;this->filt_s1=filt_s1;this->filt_ex=filt_ex;}
};

struct sort_arcs_topo{
    bool operator()(Topo_Sempl* s1, Topo_Sempl* s2){
        if(s1->val != s2->val)
            return s1->val > s2->val;
         else if(s1->lvl!=s2->lvl)
             return s1->lvl >s2->lvl;
       /* else if(s1->arc->getNode_i()->getCriticalIndex()!=s2->arc->getNode_i()->getCriticalIndex())
            return s1->arc->getNode_i()->getCriticalIndex()>s2->arc->getNode_i()->getCriticalIndex();
        else
            return s1->arc->getNode_j()->getCriticalIndex()>s2->arc->getNode_j()->getCriticalIndex();
        */
         else if(s1->filt_s0!=s2->filt_s0)
            return s1->filt_s0>s2->filt_s0;
        else if(s1->filt_s1!=s2->filt_s1)
            {//cout<<"All the same"<<"s1 index"<<s1->arc->getNode_j()->getCriticalIndex()<<" s2 index"<<s2->arc->getNode_j()->getCriticalIndex()<<endl;
                   return s1->filt_s1>s2->filt_s1;
            }
        else
        {
            for (int i=0;i<s1->filt_ex.size();i++)
            {
                   if( s1->filt_ex[i]!=s2->filt_ex[i])
                    return s1->filt_ex[i] > s2->filt_ex[i];
            }
            return s1->filt_ex[s1->filt_ex.size()-1] > s2->filt_ex[s1->filt_ex.size()-1];
        }
        
    }
};


typedef std::vector<Arrows>				CompressedToExpandedLUT;
typedef std::map<Arrows, unsigned int>	ExpandedToCompressedLUT;


class FormanGradientVector
{

    Mesh<Vertex3D,Triangle>* mesh;

    vector<unsigned int> forman_gradient;
    IG forman_ig;

    CompressedToExpandedLUT compressed_to_expanded;
    ExpandedToCompressedLUT expanded_to_compressed;

    int min, selle1, max;
    double epsilon;

    vector<double> field;


    //multiresolution here
    list<DAG_GeomNode*> geom_root;
    list<DAG_TopoNode*> topo_root;

    vector<DAG_GeomNode*>* dag_per_vertex;

    map<Node*, DAG_TopoNode*> min_map;
    map<Node*, DAG_TopoNode*> saddle_map;
    map<Node*, DAG_TopoNode*> max_map;

    //For filtering. NEW
    vector<unsigned> filtration;

    bool QEM_based;

public:
    int refined_topo;
    int refined_geometry;


    int tot_topo;
    int tot_geom;

public:
    FormanGradientVector(Mesh<Vertex3D,Triangle>* mesh, double);
    void build();

    void change_vtstar_mesh();
    TriGradient convert_compressed_to_expand(unsigned int ga);
    unsigned int convert_expand_to_compressed(Arrows ga);
    void setVE(int v, int v2);
    void setEF(int v, int f);
    Edge* getVE(int v);
    int getEF(Edge*);
    void freeVE(int v1, int v2);
    void freeEF(int v, int f);
    bool isValidTetCase(Arrows arrow);
    bool is_vertex_critical(int v);
    bool is_edge_critical(int v, int v1);
    bool is_face_critical(int f);

    void descending_2cells_extraction(bool);
    void descending_1cells_extraction(bool);

    void ascending_2cells_extraction(bool);
    void ascending_1cells_extraction(bool);
    pair<double, double> compute_incidence_graph();
    void compute_critical_simplexes(map<int, nNode*>* min, map<pair<int,int>, iNode*>* sad, map<int, nNode*>* max);
    void compute_critical_simplexes_new(map<int, nNode *> *min, map<pair<int, int>, iNode *> *sad, map<int, nNode *> *max);
    
    pair<double, double> compute_incidence_graph(map<int, nNode*>* min, map<pair<int,int>, iNode*>* sad, map<int, nNode*>* max);

    void writeVTK_1cells(char*, set<int>, set<pair<int, bool> >, set<pair<int,int> >);
    void writeVTK_2cells(char* nome_file_output, set<pair<int, bool> > critici, vector<int> triangles);
    void writeVTK_IG(char* file_name);
    void writeBoxVTK(vector<double> box);
    void writeVTK_gradient(char* nomeFile);
    void writeVTK_2cells_on_vert(char* nome_file_output, set<pair<int, bool> > critici, vector<int> triangles);
    void writeVTK_1cells_on_tri(char*,  set<pair<int, bool> >, map<int,int>);

    void output_mm(list<DAG_TopoNode*>* , list<DAG_GeomNode*>* , char* filename);

    //simplification process
    void simplify_persistence(float pers);
    void contraction(nNode*, iNode*, priority_queue<Topo_Sempl*, vector<Topo_Sempl*>, sort_arcs_topo>*);
    void removal(nNode*, iNode*, priority_queue<Topo_Sempl*, vector<Topo_Sempl*>, sort_arcs_topo>*);

    //simplify for multiresolution
    void simplify(bool, char*);
    list<DAG_GeomNode*>* simplify_geometry(vector<DAG_GeomNode*>*);
    list<DAG_TopoNode*>* simplify_persistence(map<Node*, DAG_TopoNode*>* ,map<Node*, DAG_TopoNode*>* ,map<Node*, DAG_TopoNode*>*);
    void build_persistence_queue(priority_queue<Topo_Sempl*, vector<Topo_Sempl*>, sort_arcs_topo>*);
    bool valid_gradient_configuration(int v1,int v2,int t1,int t2, bool*,bool*);
    bool valid_gradient();

    DAG_TopoNode* contraction_mr(nNode*, iNode*, priority_queue<Topo_Sempl*, vector<Topo_Sempl*>, sort_arcs_topo>*);
    DAG_TopoNode* removal_mr(nNode*, iNode*, priority_queue<Topo_Sempl*, vector<Topo_Sempl*>, sort_arcs_topo>*);

    map<vector<int>, Simplex_Graph*> compute_lower_star(int v, vector<int>* vt);
    gradient* compute_gradient(vector<int>& v, map<vector<int>, Simplex_Graph*>* lower_link_structure, map<vector<int>, int>* crtical_points);

    void push_in_coboundary(Simplex_Graph* sface, Simplex_Graph* co_sface);
    int num_unpared_faces(vector<int>& co_face, gradient* gradient_v, map<vector<int>, int>* critical_points);
    vector<int> unique_pairable_face(vector<int>& s_face, gradient* gradient_v, map<vector<int>, int>* critical_points);

    //refinement
    void build_topo_hierarchy(list<DAG_TopoNode*>*,map<Node*,DAG_TopoNode*>* ,map<Node*, DAG_TopoNode*>* ,map<Node*, DAG_TopoNode*>*);
    int refine_mesh(int adj_t1_1, int adj_t1_2, int adj_t2_1, int adj_t2_2, int v1, vector<double> coordsv2, vector<double> coordsv1, pair<bool,bool> to_switch);

    void getRefinable(DAG_GeomNode* node, vector<DAG_GeomNode*>* dag_per_vertex);
    void getRefinable(DAG_TopoNode* node);

    int insert(DAG_TopoNode*);
    int expand(DAG_TopoNode*);

    //refinement query
    void refine_geometry(double edge_lenght);
    int refine_topology(double pers);

    void refine_geometry_box(double edge_lenght,vector<double> box);
    int refine_topology_box(double pers,vector<double> box);

    void refine_geometry_sphere(double edge_lenght,vector<double> center, double radius);
    int refine_topology_sphere(double pers,vector<double> center, double radius);

    int real_index(int v_index);

    //parte su saliency
    void computeFinalSaliencies();
    bool isLocalMaximum(int v);
    void init(double eps);
    double computeCCurvature(int v);

    inline double euclidean_distance(int v1, int v2){

        return sqrt(pow(mesh->getVertex(v1).getX() - mesh->getVertex(v2).getX(),2) +
        pow(mesh->getVertex(v1).getY() - mesh->getVertex(v2).getY(),2) +
        pow(mesh->getVertex(v1).getZ() - mesh->getVertex(v2).getZ(),2));
//        pow(field[v1] - field[v2],2));
    }

    inline void set_alive_simplexes(){

        for(int i=0; i<mesh->getTopSimplexesNum(); i++)
            mesh->set_alive(i);

        for(int i=0; i<mesh->getNumVertex(); i++)
            mesh->set_v_alive(i);
    }

    inline bool is_inside_box(DAG_GeomNode* node, vector<double> box){

        vector<int> vert = node->getVertices();
        for(int i=0; i<vert.size(); i++){
            int v = vert[i];

            if((*dag_per_vertex)[v]==NULL){
                int index = mesh->get_new_index(v);
                if(!(mesh->getVertex(index).getX() > box[0] &&
                    mesh->getVertex(index).getY() > box[1] &&
                    mesh->getVertex(index).getZ() > box[2] &&
                    mesh->getVertex(index).getX() < box[3] &&
                    mesh->getVertex(index).getY() < box[4] &&
                    mesh->getVertex(index).getZ() < box[5]))
                    return false;

            }
            else if((*dag_per_vertex)[v]->isRefined()){
                int index = (*dag_per_vertex)[v]->inserted_vertex();
                if(!(mesh->getVertex(index).getX() > box[0] &&
                    mesh->getVertex(index).getY() > box[1] &&
                    mesh->getVertex(index).getZ() > box[2] &&
                    mesh->getVertex(index).getX() < box[3] &&
                    mesh->getVertex(index).getY() < box[4] &&
                    mesh->getVertex(index).getZ() < box[5]))
                    return false;
            }
            else{
                vector<double> coords = (*dag_per_vertex)[v]->getCoordsv2();
                if(!(coords[0] > box[0] &&
                    coords[1] > box[1] &&
                    coords[2] > box[2] &&
                    coords[0] < box[3] &&
                    coords[1] < box[4] &&
                    coords[2] < box[5]))
                    return false;
            }
        }

        return true;
    }

    inline bool is_inside_box(DAG_TopoNode* node, vector<double> box){

        vector<int> vert = node->getVertices();
        for(int i=0; i<vert.size(); i++){
            int v = vert[i];

            if((*dag_per_vertex)[v]==NULL){
                int index = mesh->get_new_index(v);
                if(!(mesh->getVertex(index).getX() > box[0] &&
                    mesh->getVertex(index).getY() > box[1] &&
                    mesh->getVertex(index).getZ() > box[2] &&
                    mesh->getVertex(index).getX() < box[3] &&
                    mesh->getVertex(index).getY() < box[4] &&
                    mesh->getVertex(index).getZ() < box[5]))
                    return false;

            }
            else if((*dag_per_vertex)[v]->isRefined()){
                int index = (*dag_per_vertex)[v]->inserted_vertex();
                if(!(mesh->getVertex(index).getX() > box[0] &&
                    mesh->getVertex(index).getY() > box[1] &&
                    mesh->getVertex(index).getZ() > box[2] &&
                    mesh->getVertex(index).getX() < box[3] &&
                    mesh->getVertex(index).getY() < box[4] &&
                    mesh->getVertex(index).getZ() < box[5]))
                    return false;
            }
            else{
                vector<double> coords = (*dag_per_vertex)[v]->getCoordsv2();
                if(!(coords[0] > box[0] &&
                    coords[1] > box[1] &&
                    coords[2] > box[2] &&
                    coords[0] < box[3] &&
                    coords[1] < box[4] &&
                    coords[2] < box[5]))
                    return false;
            }
        }

        return true;
    }

    inline bool is_inside_sphere(DAG_GeomNode* node, vector<double> center, double radius){

        vector<int> vert = node->getVertices();
        for(int i=0; i<vert.size(); i++){
            int v = vert[i];

            if((*dag_per_vertex)[v]==NULL){
                int index = mesh->get_new_index(v);
                float distance = pow(mesh->getVertex(index).getX()-center[0],2) +
                                 pow(mesh->getVertex(index).getY()-center[1],2) +
                                 pow(mesh->getVertex(index).getZ()-center[2],2);

                if(distance <= pow(radius,2))
                    return false;

            }
            else if((*dag_per_vertex)[v]->isRefined()){
                int index = (*dag_per_vertex)[v]->inserted_vertex();
                float distance = pow(mesh->getVertex(index).getX()-center[0],2) +
                                 pow(mesh->getVertex(index).getY()-center[1],2) +
                                 pow(mesh->getVertex(index).getZ()-center[2],2);

                if(distance <= pow(radius,2))
                    return false;
            }
            else{
                vector<double> coords = (*dag_per_vertex)[v]->getCoordsv2();

                float distance = pow(coords[0]-center[0], 2) +
                                 pow(coords[1]-center[1], 2) +
                                 pow(coords[2]-center[2], 2);

                if(distance <= pow(radius,2))
                    return false;
            }
        }

        return true;
    }

    inline bool is_inside_sphere(DAG_TopoNode* node, vector<double> center, double radius){

        vector<int> vert = node->getVertices();
        for(int i=0; i<vert.size(); i++){
            int v = vert[i];

            if((*dag_per_vertex)[v]==NULL){
                int index = mesh->get_new_index(v);
                float distance = pow(mesh->getVertex(index).getX()-center[0],2) +
                                 pow(mesh->getVertex(index).getY()-center[1],2) +
                                 pow(mesh->getVertex(index).getZ()-center[2],2);

                if(distance <= pow(radius,2))
                    return false;

            }
            else if((*dag_per_vertex)[v]->isRefined()){
                int index = (*dag_per_vertex)[v]->inserted_vertex();
                float distance = pow(mesh->getVertex(index).getX()-center[0],2) +
                                 pow(mesh->getVertex(index).getY()-center[1],2) +
                                 pow(mesh->getVertex(index).getZ()-center[2],2);

                if(distance <= pow(radius,2))
                    return false;
            }
            else{
                vector<double> coords = (*dag_per_vertex)[v]->getCoordsv2();

                float distance = pow(coords[0]-center[0], 2) +
                                 pow(coords[1]-center[1], 2) +
                                 pow(coords[2]-center[2], 2);

                if(distance <= pow(radius,2))
                    return false;
            }
        }

        return true;
    }


    inline pair<double, double> update_persistence(int v1, int v2, pair<double, double> pers){
        double pers2 = abs(mesh->getVertex(v1).getZ()-mesh->getVertex(v2).getZ());
        if(pers.first > pers2)
            pers.first = pers2;

        if(pers.second < pers2)
            pers.second = pers2;

        return pers;
    }


    inline void reorder_forman_gradient(){

        vector<unsigned int> new_forman_gradient;

        for(int i=0; i<mesh->getTopSimplexesNum(); i++)
            if(mesh->is_alive(i))
                new_forman_gradient.push_back(forman_gradient[i]);

        forman_gradient = new_forman_gradient;

    }

    inline int vector_to_index(vector<int> simplex){

        vector<int> vt = mesh->VT(simplex.front());

        for(unsigned int i=0; i<vt.size(); i++){

            int t = vt.at(i);
            int founded=0;

            for(int j=0; j<3; j++){

                int find = simplex.at(j);
                for(int k=0; k<3; k++){
                    if(find == mesh->getTopSimplex(t).TV(k)){
                        founded++;
                        break;
                    }
                }

                if(founded != j+1) break;
            }

            if(founded == 3) return t;
        }

        return -1;
    }


    inline vector<int> sort_simplex(vector<int>* simplex){
        list<int> ordered;
        for(int i=0; i<simplex->size(); i++){
            if(i==0) ordered.push_back(simplex->at(0));
            else{
                list<int>::iterator it = ordered.begin();
                for(; it!=ordered.end(); it++){
                    // if(mesh->getVertex(simplex->at(i)).getZ() > mesh->getVertex(*it).getZ())
                    //     break;
                    if (filtration[simplex->at(i)]>filtration[*it])
                        break;

                }

                if(it == ordered.begin()) ordered.push_front(simplex->at(i));
                else if(it == ordered.end()) ordered.push_back(simplex->at(i));
                else ordered.insert(it, simplex->at(i));
            }
        }

        return vector<int>(ordered.begin(), ordered.end());
    }

    inline bool is_lower(const vector<int>& simplex1,const vector<int>& simplex2, bool* same){

        if(simplex1.size() < simplex2.size()) return true;
        if(simplex1.size() > simplex2.size()) return false;

        int similar_pair=-1;
        for(int i=1; i<simplex1.size(); i++){
           // if(mesh->getVertex(simplex1.at(i)).getZ() == mesh->getVertex(simplex2.at(i)).getZ()){
               if(filtration[simplex1.at(i)]==filtration[simplex2.at(i)]){
                    if(simplex1[i]!=simplex2[i])similar_pair=i;
                    continue;
            }
            //else if(mesh->getVertex(simplex1.at(i)).getZ() > mesh->getVertex(simplex2.at(i)).getZ()) return false;
            else if(filtration[simplex1.at(i)]>filtration[simplex2.at(i)]) return false;
            else return true;
        }

        if(simplex1.back() == simplex2.back() && similar_pair == -1)
            *same=true;
        else
            *same=false;

        return (simplex1[similar_pair]<simplex2[similar_pair]);
    }

    inline void push_ordered(list<vector<int> >* coda, vector<int>& simplesso){

        bool same=false;

        list<vector<int> >::iterator it = coda->begin();
        for(; it != coda->end(); it++){
            if(!is_lower(*it, simplesso, &same)) break;
        }

        if(!same){
            if(it == coda->begin()) coda->push_front(simplesso);
            else if(it == coda->end()) coda->push_back(simplesso);
            else coda->insert(it, simplesso);
        }
    }

//    inline void check_triangles(){

//        for(int i=0; i<mesh->getTopSimplexesNum(); i++){

//            vector<int> simplex;
//            for(int j=0; j<3; i++){
//                simplex.push_back(mesh->getTopSimplex(i).TV(j));
//            }

//            simplex = sort_simplex(&simplex);
//            for(int j=0; j<3; i++){
//                cout << simplex.at(j) << " ";
//            }
//            cout << endl;


//            bool trovato = false;
//            for(int j=0; j<3; i++){
//                vector<int> edge = simplex;
//                edge.erase(edge.begin()+j);

//                simplices_map::iterator it_grad = gradient_vector.find(edge);
//                if(it_grad != gradient_vector.end(edge) /*&& it_grad->second == simplex*/){
//                    vector<int> tri = edge;
//                    tri.push_back(it_grad->second);
//                    tri = sort_simplex(&tri);
//                    if(tri == simplex)
//                        trovato = true;
//                }
//            }

//            //assert(critical_simplexes.find(simplex) != critical_simplexes.end() || trovato);

//        }
//    }

    inline int extract_missing_id(vector<int> &smaller, vector<int> &bigger)
    {

        int ret=-1;

        for(vector<int>::iterator it=bigger.begin(); it!=bigger.end(); ++it)
        {
            bool found=false;
            for(vector<int>::iterator it2=smaller.begin(); it2!=smaller.end(); ++it2)
            {
                if(*it==*it2)
                {
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                ret = *it;
                break;
            }
        }
        return ret;
    }

    inline bool in_the_box(int v, vector<double> bl, vector<double> tr){

        if(mesh->getVertex(v).getX() >= bl[0] &&
           mesh->getVertex(v).getY() >= bl[1] &&

           mesh->getVertex(v).getX() < tr[0] &&
           mesh->getVertex(v).getY() < tr[1] ) return true;

        else return false;

    }

    inline void initial_filtering(){
        unsigned int num_v = mesh->getNumVertex();
        map<double, vector<int>> vert;
cout<<"Vertices number:"<<num_v<<endl;
filtration =vector<unsigned int>(num_v,0) ;

      for(int i=0; i< num_v; i++){
            vert[mesh->getVertex(i).getZ()].push_back(i);
        }

        int count=0;
        for(auto m : vert){
            for(auto v : m.second){
                filtration[v]=count++;
            }
        }


    cout<<filtration[100]<<"; "<<filtration[59]<<endl;
    }

inline int getTriangleHighestVertex(int t){
        double z = filtration[mesh->getTopSimplex(t).TV(0)];
        int l =mesh->getTopSimplex(t).TV(0);
        for(int i=1; i<3; i++){
            if(z < filtration[mesh->getTopSimplex(t).TV(i)]){
                l = mesh->getTopSimplex(t).TV(i);
                z = filtration[mesh->getTopSimplex(t).TV(i)];
            }
        }

        return l;
}


};

#endif // FORMANGRADIENTVECTOR_H
