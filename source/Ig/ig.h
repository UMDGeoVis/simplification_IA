#ifndef IG_H
#define IG_H

#include <set>
#include <vector>

using namespace std;

class Node{
protected:
    int index; //DEBUG
    int critical_index;

public:
    inline int getCriticalIndex(){return critical_index;}
};

class Arc{

    Node* node_i;
    Node* node_j; //j=i+1;

    int simplex_to_node_i;
    int simplex_to_node_j;

    int label;

public:
    inline Arc();
    inline Arc(Node* node1, int simplex_i, Node* node2, int simplex_j){node_i=node1; node_j=node2;label=1;simplex_to_node_i=simplex_i; simplex_to_node_j=simplex_j;}
    inline Node* getNode_i(){return node_i;}
    inline Node* getNode_j(){return node_j;}
    inline int getLabel(){return label;}
    inline void setLabel(int lab){label=lab;}
    inline int getSimplexi(){return simplex_to_node_i;}
    inline int getSimplexj(){return simplex_to_node_j;}
};

class nNode : public Node{

  set<Arc*> arcs;

public:

  inline nNode(int ci){critical_index=ci;}
  inline bool addArc(Arc*  arc){return arcs.insert(arc).second;}
  inline vector<Arc*> getArcs(){return vector<Arc*>(arcs.begin(), arcs.end());}
  inline void removeArc(Arc* arc_to_remove){arcs.erase(arc_to_remove);}
  inline void clear_arcs(){arcs=set<Arc*>();}

};

class iNode : public Node{


    set<Arc*> arcs_up;
    set<Arc*> arcs_down;

    pair<int,int> edge_id;

public:
    inline iNode(int edge){critical_index = edge;}
    inline void add_edge_id(int t1, int t2){edge_id=pair<int,int>(t1,t2);}
    inline pair<int,int> get_edge_id(){return edge_id;}
    inline bool addArc(bool up, Arc* arc){
        if(up) return arcs_up.insert(arc).second;
        else return arcs_down.insert(arc).second;
    }
    inline vector<Arc*> getArcs(bool up){
        if(up) return vector<Arc*>(arcs_up.begin(), arcs_up.end());
        else return vector<Arc*>(arcs_down.begin(), arcs_down.end());
    }
    inline void removeArc(bool up, Arc* arc_to_remove){
        if(up) arcs_up.erase(arc_to_remove);
        else arcs_down.erase(arc_to_remove);
    }
    inline void clear_arcs(){arcs_up=set<Arc*>(); arcs_down=set<Arc*>();}

};


class IG
{

private:
    set<nNode*> minima;
    set<iNode*> saddle;
    set<nNode*> maxima;

    vector<set<Arc*> > arcs;

public:
    inline IG(){arcs =vector<set<Arc*> >(2,set<Arc*>());}
    inline vector<nNode*> getMinima(){return vector<nNode*>(minima.begin(), minima.end());}
    inline vector<iNode*> getSaddle(){return vector<iNode*>(saddle.begin(), saddle.end());}
    inline vector<nNode*> getMaxima(){return vector<nNode*>(maxima.begin(), maxima.end());}

    inline set<Arc*> getLevelArcs(int lvl){return arcs[lvl];}

    inline bool addNode(Node* node, int lvl){
        if(lvl == 0) return minima.insert((nNode*)node).second;
        else if(lvl == 1) return saddle.insert((iNode*)node).second;
        else return maxima.insert((nNode*)node).second;
    }

    inline void removeArc(int lvl, Arc* arc_to_remove){
        arcs[lvl].erase(arc_to_remove);
    }

    inline void removeNode(Node* node, int lvl){
        if(lvl==0) minima.erase((nNode*)node);
        else if(lvl==1) saddle.erase((iNode*)node);
        else maxima.erase((nNode*)node);

        delete node;
    }

    inline Arc* addArc(Node* n1, int s_n1, Node* n2, int s_n2, int lvl){

        Arc* new_arc = new Arc(n1,s_n1,n2,s_n2);
        arcs[lvl].insert(new_arc);

        if(lvl == 0){
            ((nNode*)n1)->addArc(new_arc);
            ((iNode*)n2)->addArc(true, new_arc);
        }
        else{
            //assert(lvl==1);
            ((iNode*)n1)->addArc(false, new_arc);
            ((nNode*)n2)->addArc(new_arc);
        }


        return new_arc;
    }

    inline Arc* already_connected(nNode* extrema, iNode* saddle){
        vector<Arc*> arcs = extrema->getArcs();
        for(int i=0; i<arcs.size(); i++)
            if(((arcs[i]->getNode_i() == extrema && arcs[i]->getNode_j() == saddle) ||
               (arcs[i]->getNode_j() == extrema && arcs[i]->getNode_i() == saddle))
                    && arcs[i]->getLabel() > 0 )
                return arcs[i];

        return NULL;
    }
};

#endif // IG_H
