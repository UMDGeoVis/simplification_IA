#ifndef AUX_FORMAN_GRADIENT_VECTOR_H
#define AUX_FORMAN_GRADIENT_VECTOR_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>

using namespace std;

typedef map<vector<int>, int> simplices_map;
typedef set<pair<vector<int>,int> > simplices_label_set;
typedef map<int,simplices_label_set > leaf_cells_map;


//questa struct permette di memorizzare in una sola lista tutti i simplessi
//con un overhead (posizioni dell'array inutilizzate) che varia a seconda del tipo di simplesso memorizzato
//il tipo di simplesso si ottiene gratis dal numero di posizioni utilizzate
struct simplex
{
    vector<int> simpl_id;
    vector<simplex*> coboundary;

    simplex(vector<int> &ids)
    {
        coboundary = vector<simplex*>();
        simpl_id.assign(ids.begin(),ids.end());
    }

    ~simplex()
    {
        simpl_id.clear();
        coboundary.clear();
    }

    void init(vector<int> &ids)
    {
        simpl_id.assign(ids.begin(),ids.end());

    }

};

struct gradient
{
    simplices_map vertices_vector;
    simplices_map edges_vector;

    gradient()
    {
        vertices_vector = map<vector<int>, int>();
        edges_vector = map<vector<int>, int>();
    }

    pair<simplices_map::iterator,bool> insert(const vector<int> &key, int v)
    {
        if(key.size() == 1)
            return vertices_vector.insert(make_pair(key,v));//key is  vertex id of tail, v is the head 
        else if(key.size() == 2)
            return edges_vector.insert(make_pair(key,v)); //key is edge, v is the id of the third vertex in the paired triangle

        cout<<"[insert] non dovrei mai arrivarci"<<endl;
        int a; cin>>a;
        return edges_vector.insert(make_pair(key,v));
    }


    void erase(const vector<int> &key, map<vector<int>,int>::iterator it)
    {
        if(key.size() == 1)
            vertices_vector.erase(it);
        else if(key.size() == 2)
            edges_vector.erase(it);
    }

    simplices_map::iterator find(const vector<int> &key)
    {
        if(key.size() == 1)
            return vertices_vector.find(key);
        else if(key.size() == 2)
            return edges_vector.find(key);

        cout<<"[find] non dovrei mai arrivarci"<<endl;
        int a; cin>>a;
        return edges_vector.find(key);
    }

    simplices_map::iterator begin(const vector<int> &key) //mi serve il vector per ritornare l'iteratore corretto
    {
        if(key.size() == 1)
            return vertices_vector.begin();
        else if(key.size() == 2)
            return edges_vector.begin();

        cout<<"[end] non dovrei mai arrivarci"<<endl;
        int a; cin>>a;
        return edges_vector.end();
    }

    simplices_map::iterator end(const vector<int> &key) //mi serve il vector per ritornare l'iteratore corretto
    {
        if(key.size() == 1)
            return vertices_vector.end();
        else if(key.size() == 2)
            return edges_vector.end();

        cout<<"[end] non dovrei mai arrivarci"<<endl;
        int a; cin>>a;
        return edges_vector.end();
    }


    void add_local_gradients(gradient &vecs)
    {
        vertices_vector.insert(vecs.vertices_vector.begin(),vecs.vertices_vector.end());
        edges_vector.insert(vecs.edges_vector.begin(),vecs.edges_vector.end());
    }



    void clear()
    {
        vertices_vector.clear();
        edges_vector.clear();
    }

    bool is_empty()
    {
        return ((vertices_vector.size() == 0) && (edges_vector.size() == 0) );
    }
};

#endif // AUX_FORMAN_GRADIENT_VECTOR_H
