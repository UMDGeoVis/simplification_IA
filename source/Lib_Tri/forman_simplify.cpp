#include "formangradientvector.h"
#include "Timer.h"

void FormanGradientVector::simplify_persistence(float pers_limit){

    map<int, nNode*>* minima_nodes = new map<int, nNode*>();
    map<pair<int,int>, iNode*>* saddle_nodes = new map<pair<int,int>, iNode*>();
    map<int, nNode*>* maxima_nodes = new map<int, nNode*>();

    compute_critical_simplexes(minima_nodes, saddle_nodes, maxima_nodes);
    pair<double, double> pers = compute_incidence_graph(minima_nodes, saddle_nodes, maxima_nodes);//pers: persistence range


    priority_queue<Topo_Sempl*, vector<Topo_Sempl*>, sort_arcs_topo>* queue = new priority_queue<Topo_Sempl*, vector<Topo_Sempl*>, sort_arcs_topo>();
    build_persistence_queue(queue);

    Topo_Sempl* sempl;
    iNode* saddle=NULL;
    nNode* extrema=NULL;
    int count=0;
    int lvl0=0,lvl1=0,lvl2=0;
    while(!queue->empty() && queue->top()->val <= pers_limit){

        count++;
        sempl = queue->top();
        queue->pop();

        if(sempl->arc->getLabel() != 1 ){
            if(sempl->arc->getLabel() == -1)
                delete sempl->arc;

            lvl2++;
            delete sempl; // if label==2? don't simplify
            continue;
        }

        if(sempl->lvl == 0){
            lvl0++;
            saddle  = sempl->arc->getNode_j();
            extrema = sempl->arc->getNode_i();

            if(saddle->getArcs(true).size() != 2){ //In which case? Sometimes edge on the border. And some other cases
                delete sempl;
                continue;
            }
            contraction(extrema, saddle, queue);
            refined_topo++;

        }
        else{
            lvl1++;
            saddle  = sempl->arc->getNode_i();
            extrema = sempl->arc->getNode_j();

            if(saddle->getArcs(false).size() != 2){
                delete sempl;
                continue;
            }
            removal(extrema,saddle, queue);
            refined_topo++;
        }
        delete sempl;
//        delete extrema;
//        delete saddle;
    }

    // cout<<"original queue size:"<<count<<endl;
    cout<<"Level 0:"<<lvl0<<", Level 1:"<<lvl1<<", Level 2:"<<lvl2<<endl;
    delete queue;

    delete minima_nodes;
    delete saddle_nodes;
    delete maxima_nodes;
}


void FormanGradientVector::contraction(nNode *extrema, iNode *saddle, priority_queue<Topo_Sempl *, vector<Topo_Sempl *>, sort_arcs_topo> *queue){

    vector<Arc*> arcs = saddle->getArcs(true);
    //assert(arcs.size() == 2);
    nNode* other_extrema;
    int vertex,next_vertex;
    int ending_path_simplex;


    int ex_minimum = extrema->getCriticalIndex();
    // Find critical edge from two adjacent triangles.
    pair<int,int> critical_edge_tetra=saddle->get_edge_id();
    Edge* critical_edge=NULL;
    if(critical_edge_tetra.second<0){ //Only one side
        critical_edge = mesh->getTopSimplex(critical_edge_tetra.first).TE(-critical_edge_tetra.second-1);
    }
    else{
        for(int i=0; i<3; i++){
            if(!(mesh->getTopSimplex(critical_edge_tetra.second).contains(mesh->getTopSimplex(critical_edge_tetra.first).TV(i)))){
                critical_edge = mesh->getTopSimplex(critical_edge_tetra.first).TE(i);
                break;
            }
        }
    }


    //SETUP per dopo
    if(arcs[0]->getNode_i() == extrema){  // if arcs[0] is the minimum to be deleted, find the other minimum
        other_extrema = (nNode*)arcs[1]->getNode_i();
        next_vertex = mesh->getTopSimplex(saddle->get_edge_id().first).TV(arcs[0]->getSimplexj());
        vertex = mesh->getTopSimplex(saddle->get_edge_id().first).TV(arcs[1]->getSimplexj());
        ending_path_simplex = arcs[1]->getSimplexi();
    }
    else{
        other_extrema = (nNode*)arcs[0]->getNode_i();
        next_vertex = mesh->getTopSimplex(saddle->get_edge_id().first).TV(arcs[1]->getSimplexj());
        vertex = mesh->getTopSimplex(saddle->get_edge_id().first).TV(arcs[0]->getSimplexj());
        ending_path_simplex = arcs[0]->getSimplexi();
    }

    Edge* old_edge = critical_edge;

    //modifiche sul gradiente
    while(next_vertex != ex_minimum){

        //trovo il nuovo edge;
        Edge* edge = getVE(next_vertex);
        vertex=next_vertex;
        edge->EV(0) == vertex ? next_vertex = edge->EV(1) : next_vertex = edge->EV(0);

        //azzero la sua adiacenza;
        freeVE(vertex,next_vertex);

        //accoppio il vecchio edge;
        if(old_edge->EV(0) == vertex){
            setVE(old_edge->EV(0), old_edge->EV(1));
        }
        else{
            setVE(old_edge->EV(1), old_edge->EV(0));
        }

        delete old_edge;
        old_edge = edge;
    }
    setVE(next_vertex, vertex);


    //modifiche sul IG 
    //// Only delete arc, no need to change directions? Because the edge is no longer a saddle
    arcs[0]->setLabel(-1);
    arcs[1]->setLabel(-1);
    vector<Arc*> maxima_arcs = saddle->getArcs(false);
    nNode* maxima;
    for(int i=0;i<maxima_arcs.size(); i++){

        maxima_arcs[i]->setLabel(-1);
        maxima = ((nNode*)maxima_arcs[i]->getNode_j());

        forman_ig.removeArc(1, maxima_arcs[i]);
        maxima->removeArc(maxima_arcs[i]);
        saddle->removeArc(false, maxima_arcs[i]);

    }


    vector<Arc*> saddle_arcs = extrema->getArcs();
    iNode* node_saddle1=NULL;
    for(int i=0; i<saddle_arcs.size(); i++){

        saddle_arcs[i]->setLabel(-1);
        int starting_index = saddle_arcs[i]->getSimplexj();

        node_saddle1 = ((iNode*)saddle_arcs[i]->getNode_j());

        forman_ig.removeArc(0,saddle_arcs[i]);
        node_saddle1->removeArc(true, saddle_arcs[i]);
        extrema->removeArc(saddle_arcs[i]);

        if(node_saddle1 != saddle){
            if(forman_ig.already_connected(other_extrema,node_saddle1)==NULL){
                Arc* arco = forman_ig.addArc(other_extrema, ending_path_simplex, node_saddle1, starting_index, 0);

                if(arco->getLabel() == 1){
                    double val = abs(field[arco->getNode_i()->getCriticalIndex()] - field[arco->getNode_j()->getCriticalIndex()]);
                }//Q:Should arco be added to persistence queue? val here is not used. 
            }   //A: If next step is to simplify the geometry, then here no need to update the queue. A new queue will be generated when simplify the gradient again. 
            else{
                forman_ig.already_connected(other_extrema,node_saddle1)->setLabel(2);
                Arc* arco = forman_ig.addArc(other_extrema, ending_path_simplex, node_saddle1, starting_index, 0);
                arco->setLabel(2);  
            }
        }
    }

    vector<Arc*> minima_arcs = saddle->getArcs(true);
    nNode* minima;
    for(int i=0;i<minima_arcs.size(); i++){
        minima_arcs[i]->setLabel(-1);
        minima = ((nNode*)minima_arcs[i]->getNode_i());

        forman_ig.removeArc(0, minima_arcs[i]);
        minima->removeArc(minima_arcs[i]);
        saddle->removeArc(true, minima_arcs[i]);
    }

    forman_ig.removeNode(saddle,1);
    forman_ig.removeNode(extrema,0);

}

void FormanGradientVector::removal(nNode *extrema, iNode *saddle, priority_queue<Topo_Sempl *, vector<Topo_Sempl *>, sort_arcs_topo> *queue){

    vector<Arc*> arcs = saddle->getArcs(false);

    nNode* other_extrema;
    int triangle;
    int vertex=-1;
    int ending_path_simplex;


    if(arcs[0]->getNode_j() == extrema){
        other_extrema = (nNode*)arcs[1]->getNode_j();
        triangle = arcs[0]->getSimplexi();
        ending_path_simplex = arcs[1]->getSimplexj();
    }
    else{
        other_extrema = (nNode*)arcs[0]->getNode_j();
        triangle = arcs[1]->getSimplexi();
        ending_path_simplex = arcs[0]->getSimplexj();
    }

    Edge* old_edge=NULL;
    Edge* edge=NULL;
    pair<int,int> twotriangles= saddle->get_edge_id();
    if(twotriangles.second >= 0){
        for(int j=0; j<3;j++){
            if(twotriangles.second == mesh->getTopSimplex(twotriangles.first).TT(j)){
                old_edge = mesh->getTopSimplex(twotriangles.first).TE(j);
                break;
            }
        }
    }
    else{
        old_edge = mesh->getTopSimplex(twotriangles.first).TE(-twotriangles.second -1);
    }


    //modifiche sul gradiente
    int next_triangle=triangle;

    while(!is_face_critical(next_triangle)){

        //trovo il nuovo edge;
        int i=0;
        for(; i<3; i++){
            edge = mesh->getTopSimplex(next_triangle).TE(i);
            int tri = getEF(edge);
            if(tri == next_triangle) break;
            delete edge;
        }
        // cout<<old_edge->EV(0)<<", "<<old_edge->EV(1)<<endl;
        // cout<<edge->EV(0)<<", "<<edge->EV(1)<<endl;
        triangle=next_triangle;
        next_triangle=mesh->getTopSimplex(triangle).TT(i); // i is determined in the for loop above
        // cout<<"ET1:"<<triangle<<"  ET2:"<<next_triangle<<endl;
        //cout<<mesh->getTopSimplex(triangle).TV(0)<<endl;
        vertex = mesh->getTopSimplex(triangle).TV(i);
        //cout<<"Vertex "<<vertex<<": "<<mesh->getVertex(vertex).getX()<<", "<<mesh->getVertex(vertex).getY()<<endl;
        // cout<<"Vertex index:"<<vertex<<endl;
        //azzero la sua adiacenza;
        freeEF(vertex,triangle);  //Guess: E is the edge opposite to vertex

        //accoppio il vecchio edge;
        for(int i=0;i<3; i++){
            if(old_edge->EV(0) != mesh->getTopSimplex(triangle).TV(i) &&
               old_edge->EV(1) != mesh->getTopSimplex(triangle).TV(i)   ){
                vertex = mesh->getTopSimplex(triangle).TV(i);
                break;}
        }
        setEF(vertex,triangle); //revert the direction of gradient

        delete old_edge;
        old_edge = edge;
    }


    for(int i=0;i<3; i++){
        if(old_edge->EV(0) != mesh->getTopSimplex(next_triangle).TV(i) &&
           old_edge->EV(1) != mesh->getTopSimplex(next_triangle).TV(i)   ){
            vertex = mesh->getTopSimplex(next_triangle).TV(i);
            break;
        }
    }
    setEF(vertex,next_triangle);

    //modifiche sul IG (aggiustare situazioni qui!!!!)
    arcs[0]->setLabel(-1);
    arcs[1]->setLabel(-1);


    vector<Arc*> minima_arcs = saddle->getArcs(true);
    nNode* minima;
    for(int i=0;i<minima_arcs.size(); i++){
        minima_arcs[i]->setLabel(-1);
        minima = ((nNode*)minima_arcs[i]->getNode_i());

        forman_ig.removeArc(0, minima_arcs[i]);
        minima->removeArc(minima_arcs[i]);
        saddle->removeArc(true, minima_arcs[i]);
    }


    vector<Arc*> saddle_arcs = extrema->getArcs();
    iNode* node_saddle1=NULL;
    for(int i=0; i<saddle_arcs.size(); i++){

        node_saddle1 = ((iNode*)saddle_arcs[i]->getNode_i());

        if(node_saddle1 != saddle){

            saddle_arcs[i]->setLabel(-1);
            int starting_path_simplex = saddle_arcs[i]->getSimplexi();
            forman_ig.removeArc(1,saddle_arcs[i]);
            node_saddle1->removeArc(false, saddle_arcs[i]);
            extrema->removeArc(saddle_arcs[i]);

            if(forman_ig.already_connected(other_extrema,node_saddle1)==NULL){
                Arc* arco = forman_ig.addArc(node_saddle1, starting_path_simplex, other_extrema, ending_path_simplex, 1);
                if(arco->getLabel() == 1){  //when will the label be 1 and when is 0
                    double val = abs(field[arco->getNode_i()->getCriticalIndex()] - field[mesh->getTopSimplexHighestVertex(arco->getNode_j()->getCriticalIndex())]);
                }
            }
            else{
                forman_ig.already_connected(other_extrema,node_saddle1)->setLabel(2);
            }
        }
    }

    vector<Arc*> maxima_arcs = saddle->getArcs(false);
    nNode* maxima;
    for(int i=0;i<maxima_arcs.size(); i++){
        maxima_arcs[i]->setLabel(-1);
        maxima = ((nNode*)maxima_arcs[i]->getNode_j());

        forman_ig.removeArc(1, maxima_arcs[i]);
        maxima->removeArc(maxima_arcs[i]);
        saddle->removeArc(false, maxima_arcs[i]);

    }

    forman_ig.removeNode(saddle,1);
    forman_ig.removeNode(extrema,2);


}
