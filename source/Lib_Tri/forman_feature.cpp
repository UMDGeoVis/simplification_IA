#include "formangradientvector.h"
#include "Timer.h"

void FormanGradientVector::descending_2cells_extraction(bool with_geometry){


    Timer time;
    time.start();

    set<pair<int, bool> > critici;
    vector<int> labeling = vector<int>(mesh->getTopSimplexesNum(),-1);
    for(int i=0; i<mesh->getTopSimplexesNum(); i++){

        if(is_face_critical(i)){

            queue<Edge* > coda;
            labeling[i]=i;

            vector<int> simplexi;
            int v=-1;
            for(int j=0; j<3; j++){
                    simplexi.push_back(mesh->getTopSimplex(i).TV(j));
                    if(with_geometry){
                        if(v==-1 || mesh->getVertex(v).getZ() < mesh->getVertex(mesh->getTopSimplex(i).TV(j)).getZ())
                        v=mesh->getTopSimplex(i).TV(j);
                    }
            }

            if(with_geometry)
                critici.insert(pair<int,bool>(v,true));

            for(int j=0; j<3; j++){
                vector<int> simplex = simplexi;
                simplex.erase(simplex.begin()+j);
                coda.push(new Edge(simplex[0], simplex[1]));
            }
            simplexi.clear();

            while(!coda.empty()){
                Edge* edge = coda.front();
                coda.pop();

                int tri = getEF(edge);
                if(tri != -1){

                    labeling[tri]=i;

                    for(int j=0; j<3; j++){
                        Edge* new_e = mesh->getTopSimplex(tri).TE(j);
                        if(*new_e != *edge){
                            coda.push(new_e);
                        }
                        else delete new_e;
                    }

                    delete edge;
                }
                else if(is_edge_critical(edge->EV(0), edge->EV(1))){
                    //here for critical edges;
                    if(with_geometry)
                        mesh->getVertex(edge->EV(0)).getZ() > mesh->getVertex(edge->EV(1)).getZ() ? critici.insert(pair<int,bool>(edge->EV(0),true)) : critici.insert(pair<int,bool>(edge->EV(1),true));
                }
            }
        }
    }


    time.stop();
    cout << "2-cell extraction terminated in " << time.getElapsedTimeInSec() << endl;

    if(with_geometry) writeVTK_2cells("prove_2cell.vtk", critici, labeling);
}


void FormanGradientVector::descending_1cells_extraction(bool with_geometry){

    Timer time;
    time.start();
    queue<int> coda;

    set<int> vertici;
    set<pair<int, bool> > critici;
    set<pair<int,int> > edges;

    vector<bool> visited = vector<bool>(mesh->getTopSimplexesNum(),false);
    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        visited[i]=true;
        for(int j=0; j<3; j++){
            if(mesh->getTopSimplex(i).TT(j) == -1 || !visited[mesh->getTopSimplex(i).TT(j)]){

                Edge* edge = mesh->getTopSimplex(i).TE(j);
                if(is_edge_critical(edge->EV(0), edge->EV(1))){


                    if(with_geometry){
                        mesh->getVertex(edge->EV(0)).getZ() > mesh->getVertex(edge->EV(1)).getZ() ? edges.insert(pair<int,int>(edge->EV(0),edge->EV(1))) : edges.insert(pair<int,int>(edge->EV(1),edge->EV(0)));
                        mesh->getVertex(edge->EV(0)).getZ() > mesh->getVertex(edge->EV(1)).getZ() ? critici.insert(pair<int,bool>(edge->EV(0),true)) : critici.insert(pair<int,bool>(edge->EV(1),true));
                        vertici.insert(edge->EV(0));
                        vertici.insert(edge->EV(1));
                    }

                    coda.push(edge->EV(0));
                    coda.push(edge->EV(1));
                    delete edge;

                    while(!coda.empty()){

                        int vert = coda.front();
                        coda.pop();

                        edge = getVE(vert);
                        if(edge != NULL){
                            if(with_geometry){
                                mesh->getVertex(edge->EV(0)).getZ() > mesh->getVertex(edge->EV(1)).getZ() ? edges.insert(pair<int,int>(edge->EV(0),edge->EV(1))) : edges.insert(pair<int,int>(edge->EV(1),edge->EV(0)));
                                vertici.insert(edge->EV(0));
                                vertici.insert(edge->EV(1));
                            }
                            int v2 = edge->EV(0) == vert ? edge->EV(1) : edge->EV(0);
                            delete edge;
                            coda.push(v2);
                        }
                        else if(is_vertex_critical(vert)){
                            //here critical vertex
                            if(with_geometry)
                                critici.insert(pair<int,bool>(vert, false));
                        }
                    }
                }
            }
        }
    }

    time.stop();
    cout << "1-cell extraction terminated in " << time.getElapsedTimeInSec() << endl;

    if(with_geometry)writeVTK_1cells("prove_descending.vtk", vertici, critici, edges);
}

void FormanGradientVector::ascending_2cells_extraction(bool with_geometry){

    Timer time;
    time.start();

    queue<int> coda;
    vector<int> label = vector<int>(mesh->getNumVertex(), -1);
    set<pair<int,bool> > critici;

    int critical_edge =0;
    int minima=0;

    for(int i=0; i<mesh->getNumVertex(); i++){
        if(mesh->is_v_alive(i) && is_vertex_critical(i)){
            minima++;
            label[i]=i;

            if(with_geometry)
                critici.insert(pair<int,bool>(i,true));

            vector<Edge*> ve = mesh->VE(i);
            for(unsigned int j=0; j<ve.size(); j++){
                int v2 = ve[j]->EV(0) == i ? ve[j]->EV(1) : ve[j]->EV(0);
                Edge* edge = getVE(v2);
                if(edge != NULL && *edge == *(ve[j])){
                    coda.push(v2);
                }
                else if(edge != NULL && is_edge_critical(edge->EV(0), edge->EV(1))){
                    critical_edge++;
                    if(with_geometry)
                        mesh->getVertex(edge->EV(0)).getZ() > mesh->getVertex(edge->EV(1)).getZ() ? critici.insert(pair<int,bool>(edge->EV(0),true)) : critici.insert(pair<int,bool>(edge->EV(1),true));
                }
                if(edge != NULL) delete edge;
                delete ve[j];
            }


            while(!coda.empty()){

                int v = coda.front();
                coda.pop();
                label[v]=i;

                ve = mesh->VE(v);
                for(unsigned int j=0; j<ve.size(); j++){
                    int v2 = ve[j]->EV(0) == v ? ve[j]->EV(1) : ve[j]->EV(0);
                    Edge* edge = getVE(v2);
                    if(edge != NULL && *edge == *(ve[j])){
                        coda.push(v2);
                    }
                    else if(edge != NULL && is_edge_critical(edge->EV(0), edge->EV(1))){
                        critical_edge++;
                        if(with_geometry)
                            mesh->getVertex(edge->EV(0)).getZ() > mesh->getVertex(edge->EV(1)).getZ() ? critici.insert(pair<int,bool>(edge->EV(0),true)) : critici.insert(pair<int,bool>(edge->EV(1),true));
                    }
                    if(edge != NULL) delete edge;
                    delete ve[j];
                }
            }

        }
    }

    time.stop();
    cout << "2-cell ascending extraction terminated in " << time.getElapsedTimeInSec() << endl;

    if(with_geometry) writeVTK_2cells_on_vert("prove_2cell_vert.vtk", critici, label);
}

void FormanGradientVector::ascending_1cells_extraction(bool with_geometry){

    Timer time;
    time.start();

    vector<bool> visited = vector<bool>(mesh->getTopSimplexesNum(),false);
    set<pair<int,bool> > critici;
    map<int,int> visited_triangle;

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        visited[i]=true;
        for(int j=0; j<3; j++){
            if(mesh->getTopSimplex(i).TT(j) == -1 || !visited[mesh->getTopSimplex(i).TT(j)]){

                queue<int> coda;
                Edge* edge = mesh->getTopSimplex(i).TE(j);
                if(is_edge_critical(edge->EV(0), edge->EV(1))){
                    if(with_geometry)
                        mesh->getVertex(edge->EV(0)).getZ() > mesh->getVertex(edge->EV(1)).getZ() ? critici.insert(pair<int,bool>(edge->EV(0),true)) : critici.insert(pair<int,bool>(edge->EV(1),true));


                    vector<int> et = mesh->ET(*edge);
                    for(int j=0; j<et.size(); j++){
                        coda.push(et[j]);
                    }

                    while(!coda.empty()){

                        int t = coda.front();
                        coda.pop();

                        visited_triangle[t]=i;

                        if(is_face_critical(t)){
                            //here critical face
                            if(with_geometry){
                                int v=-1;
                                for(int j=0; j<3; j++){
                                    if(v==-1 || filtration[v] < filtration[mesh->getTopSimplex(i).TV(j)])
                                        v=mesh->getTopSimplex(i).TV(j);
                                }
                            }
                        }
                        else{
                            for(int j=0; j<3; j++){
                                int t1 = getEF(mesh->getTopSimplex(t).TE(j));
                                int t_adj = mesh->getTopSimplex(t).TT(j);
                                if(t == t1 && t_adj != -1){
                                    coda.push(t_adj);
                                }
                            }
                        }
                    }

                }
            }
        }
    }

    time.stop();
    cout << "2-cell ascending extraction terminated in " << time.getElapsedTimeInSec() << endl;

    if(with_geometry)
        writeVTK_1cells_on_tri("prove_1cell_asc.vtk", critici, visited_triangle);
}

pair<double, double> FormanGradientVector::compute_incidence_graph(){

    forman_ig=IG();
    int arc_saddle_minimum=0, arc_saddle_maximum=0,level2=0;
    pair<double, double> pers = pair<double, double>(100,0);

    vector<bool> visited = vector<bool>(mesh->getTopSimplexesNum(),false);

    map<int,nNode*> minima_nodes;
    map<int,nNode*> maxima_nodes;

    int critical_edge=0;

    vector<bool> vertices_visited = vector<bool>(mesh->getNumVertex(),false);

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        visited[i]=true;
        if(!mesh->is_alive(i))continue; //ATTENZIONE, QUALCHE TRIANGOLO RIMANE NON VIVO E INVALIDA LA RICERCA
        for(int j=0; j<3; j++){
            if(mesh->getTopSimplex(i).TT(j) == -1 || !visited[mesh->getTopSimplex(i).TT(j)]){

                Edge* edge = mesh->getTopSimplex(i).TE(j);
                int vert_saddle;
                if(is_edge_critical(edge->EV(0), edge->EV(1))){
                    critical_edge++;
                    iNode* saddle_node;
                    if(filtration[edge->EV(0)] > filtration[edge->EV(1)]){
                        vert_saddle = edge->EV(0);
                    }
                    else{
                        vert_saddle = edge->EV(1);
                    }


//                    //assert(!vertices_visited[vert_saddle]);
//                    vertices_visited[vert_saddle]=true;

                    saddle_node = new iNode(vert_saddle);
                    forman_ig.addNode(saddle_node,1);

                    vector<int> et = mesh->ET(*edge);
                    for(int k=0; k<et.size(); k++){
                        list<int> stack;
                        int last_t=et[k];
                        stack.push_front(et[k]);

                        while(!stack.empty()){

                            int t = stack.front();
                            stack.pop_front();


                            if(is_face_critical(t)){
                                if(maxima_nodes.find(t) == maxima_nodes.end()){
                                    maxima_nodes[t]=new nNode(mesh->getTopSimplexHighestVertex(t));
                                  //  cout<<"[FOR DEBUG] tid:"<<t<<"; highest vertex"<<mesh->getTopSimplexHighestVertex(t)<<endl;
                                    forman_ig.addNode(maxima_nodes[t],2);
                                }
                                Arc* arc = forman_ig.already_connected(maxima_nodes[t],saddle_node);
                                if(arc == NULL){
                                    forman_ig.addArc(saddle_node, et[k],maxima_nodes[t], last_t, 1);
                                    pers = update_persistence(saddle_node->getCriticalIndex(), maxima_nodes[t]->getCriticalIndex() ,pers);
                                arc_saddle_maximum++;
                                }
                                else
                                   { arc->setLabel(2);
                                   level2++;
                                   }
                            }
                            else{
                                for(int k=0; k<3; k++){
                                    int t1 = getEF(mesh->getTopSimplex(t).TE(k));
                                    int t_adj = mesh->getTopSimplex(t).TT(k);
                                    if(t == t1 && t_adj != -1){
                                        stack.push_front(t_adj);
                                        last_t = t;
                                    }
                                }

                            }
                        }
                    }

                    for(int k=0; k<2; k++){
                        list<int> stack;
                        int last_v=edge->EV(k);
                        stack.push_front(edge->EV(k));

                        while(!stack.empty()){

                            int vert = stack.front();
                            stack.pop_front();
                            ////assert(mesh->is_v_alive(vert));

                            Edge* edge2 = getVE(vert);
                            if(edge2 != NULL){
                                int v2 = edge2->EV(0) == vert ? edge2->EV(1) : edge2->EV(0);
                                stack.push_front(v2);
                                last_v = vert;
                            }
                            else if(is_vertex_critical(vert)){
                                if(minima_nodes.find(vert) == minima_nodes.end()){
                                    minima_nodes[vert]=new nNode(vert);
                                    forman_ig.addNode(minima_nodes[vert],0);
                                }
                                Arc* arc = forman_ig.already_connected(minima_nodes[vert],saddle_node);
                                if(arc == NULL){
                                    forman_ig.addArc(minima_nodes[vert], last_v,saddle_node, edge->EV(k), 0);
                                    pers = update_persistence(minima_nodes[vert]->getCriticalIndex(), saddle_node->getCriticalIndex() ,pers);
                                arc_saddle_minimum++;
                                }
                                else
                                   { arc->setLabel(2);
                                   level2++;
                                   }
                            }
                            delete edge2;
                        }
                    }


                    delete edge;
                }
            }
        }
    }
    cout<<"Number of critical Points:"<<endl;
    cout << minima_nodes.size() << " " << critical_edge << " " << maxima_nodes.size() << endl;

    
    cout<<"Number of arcs:"<<endl;
    cout<<"Saddle - Maximum: "<<arc_saddle_maximum<<endl;
    cout<<"Saddle - Minimum: "<< arc_saddle_minimum<<endl;
    cout<<"Level 2 arcs:"<<level2<<endl;

    return pers;
}


pair<double, double> FormanGradientVector::compute_incidence_graph(map<int, nNode*>* minima_nodes, map<pair<int,int>, iNode*>* saddle_nodes, map<int, nNode*>* maxima_nodes){

//    for(map<int, nNode*>::iterator it=minima_nodes->begin(); it != minima_nodes->end(); it++)
//        (*it).second->clear_arcs();

//    for(map<int, nNode*>::iterator it=maxima_nodes->begin(); it != maxima_nodes->end(); it++)
//        (*it).second->clear_arcs();

//    for(map<pair<int,int>, iNode*>::iterator it=saddle_nodes->begin(); it != saddle_nodes->end(); it++)
//        (*it).second->clear_arcs();

    forman_ig=IG();
    int arc_saddle_minimum=0, arc_saddle_maximum=0, level2=0;

    //assert(forman_ig.getLevelArcs(0).size() == 0 && forman_ig.getLevelArcs(1).size() ==0 );

    pair<double, double> pers = pair<double, double>(100,0);

    vector<bool> visited = vector<bool>(mesh->getTopSimplexesNum(),false);

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        visited[i]=true;
        if(!mesh->is_alive(i))continue;
        for(int j=0; j<3; j++){
            if(mesh->getTopSimplex(i).TT(j) == -1 || !visited[mesh->getTopSimplex(i).TT(j)]){

                Edge* edge = mesh->getTopSimplex(i).TE(j);
                if(is_edge_critical(edge->EV(0), edge->EV(1))){

                    iNode* saddle_node=NULL;
                    //qui dobbiamo recuperare la sella
                    if(mesh->getTopSimplex(i).TT(j) != -1){
                        //assert(saddle_nodes->find(pair<int,int>(i,mesh->getTopSimplex(i).TT(j))) != saddle_nodes->end() || saddle_nodes->find(pair<int,int>(mesh->getTopSimplex(i).TT(j),i)) != saddle_nodes->end() );
                        saddle_node = mesh->getTopSimplex(i).TT(j) < i ? (*saddle_nodes)[pair<int,int>(i,mesh->getTopSimplex(i).TT(j))] : (*saddle_nodes)[pair<int,int>(mesh->getTopSimplex(i).TT(j),i)];
                    }
                    else{
                        //assert(saddle_nodes->find(pair<int,int>(i,-j-1)) != saddle_nodes->end());
                        saddle_node = (*saddle_nodes)[pair<int,int>(i,-j-1)];
                    }
                    //assert(saddle_node!=NULL);
                    //assert(saddle_node->getCriticalIndex() == edge->EV(0) || saddle_node->getCriticalIndex() == edge->EV(1));

                    //
                    //
                    //////explore_asc1cell_mig() in Terrain trees
                    vector<int> et = mesh->ET(*edge);
                    for(int k=0; k<et.size(); k++){
                        list<int> stack;
                        int last_t=et[k];
                        stack.push_front(et[k]);

                        while(!stack.empty()){

                            int t = stack.front();
                            stack.pop_front();


                            if(is_face_critical(t)){ //TODO: add an arc count
                                nNode* maximum = (*maxima_nodes)[t];
                                //assert(maximum != NULL);
                                //  cout<<"[FOR DEBUG] tid:"<<t<<"; highest vertex"<<mesh->getTopSimplexHighestVertex(t)<<endl;
                                //  cout<<"Critical index"<<maximum->getCriticalIndex()<<endl;
                                //  cout<<"Number of vertices:"<<mesh->getNumVertex()<<endl;
                                Arc* arc = forman_ig.already_connected(maximum,saddle_node);
                                if(arc == NULL){
                                    forman_ig.addArc(saddle_node, et[k],maximum, last_t, 1);
                                    pers = update_persistence(saddle_node->getCriticalIndex(), getTriangleHighestVertex(maximum->getCriticalIndex()) ,pers);
                                arc_saddle_maximum++;
                                }
                                else
                                   { arc->setLabel(2);
                                    level2++;
                                   }
                                
                            }
                            else{
                                for(int k=0; k<3; k++){
                                    int t1 = getEF(mesh->getTopSimplex(t).TE(k));
                                    int t_adj = mesh->getTopSimplex(t).TT(k);
                                    if(t == t1 && t_adj != -1){
                                        stack.push_front(t_adj);
                                        last_t = t;
                                    }
                                }

                            }
                        }
                    }
                    //////
                    //
                    //
                    //
                    //////explore_desc1cell_mig() in Terrain trees

                    for(int k=0; k<2; k++){
                        list<int> stack;
                        int last_v=edge->EV(k);
                        stack.push_front(edge->EV(k));

                        while(!stack.empty()){

                            int vert = stack.front();
                            stack.pop_front();
                            ////assert(mesh->is_v_alive(vert));

                            Edge* edge2 = getVE(vert);
                            if(edge2 != NULL){
                                int v2 = edge2->EV(0) == vert ? edge2->EV(1) : edge2->EV(0);
                                stack.push_front(v2);
                                last_v = vert;
                            }
                            else if(is_vertex_critical(vert)){ //TODO: add an arc count
                                nNode* minimum = (*minima_nodes)[vert];
                                //assert(minimum != NULL);
                                Arc* arc = forman_ig.already_connected(minimum,saddle_node);
                                if(arc == NULL){
                                    int v_index = mesh->getTopSimplex(saddle_node->get_edge_id().first).vertex_index(edge->EV(k));
                                    //assert(v_index != -1);

                                    forman_ig.addArc(minimum, last_v,saddle_node, v_index, 0);
                                    pers = update_persistence(minimum->getCriticalIndex(), saddle_node->getCriticalIndex() ,pers);
                                 arc_saddle_minimum++;
                                }
                                else
                                    {arc->setLabel(2);
                                    level2++;
                                    }
                               
                            }
                            delete edge2;
                        }
                    }


                    delete edge;
                }
            }
        }
    }

    cout<<"Number of arcs:"<<endl;
    cout<<"Saddle - Maximum: "<<arc_saddle_maximum<<endl;
    cout<<"Saddle - Minimum: "<< arc_saddle_minimum<<endl;
        cout<<"Level 2 arcs:"<<level2<<endl;
    return pers;
}

void FormanGradientVector::compute_critical_simplexes(map<int, nNode *> *min, map<pair<int, int>, iNode *> *sad, map<int, nNode *> *max){

    vector<bool> visited = vector<bool>(mesh->getTopSimplexesNum(),false);

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        visited[i]=true;
        if(!mesh->is_alive(i))continue;
        for(int j=0; j<3; j++){
            if(mesh->getTopSimplex(i).TT(j) == -1 || !visited[mesh->getTopSimplex(i).TT(j)]){

                Edge* edge = mesh->getTopSimplex(i).TE(j);
                int vert_saddle;
                //assert(mesh->is_v_alive(edge->EV(0)) && mesh->is_v_alive(edge->EV(1)));
                if(is_edge_critical(edge->EV(0), edge->EV(1))){
                    if(filtration[edge->EV(0)] > filtration[edge->EV(1)]){
                        vert_saddle = edge->EV(0);
                    }
                    else{
                        vert_saddle = edge->EV(1);
                    }


                    iNode* node = new iNode(vert_saddle);
                    if(mesh->getTopSimplex(i).TT(j) != -1){
                        if(mesh->getTopSimplex(i).TT(j) < i){
                            (*sad)[pair<int,int>(i,mesh->getTopSimplex(i).TT(j))]=node;
                            node->add_edge_id(i,mesh->getTopSimplex(i).TT(j)); //So edge id are the two adjacent triangle?
                        }
                        else{
                            (*sad)[pair<int,int>(mesh->getTopSimplex(i).TT(j),i)]=node;
                            node->add_edge_id(mesh->getTopSimplex(i).TT(j),i);
                        }
                    }
                    else{
                        (*sad)[pair<int,int>(i,-j-1)]=node;
                        node->add_edge_id(i,-j-1);
                    }
                }
                delete edge;
            }
        }

        if(is_face_critical(i)){
            (*max)[i]=new nNode(i);
        }
    }

    for(int i=0;i<mesh->getNumVertex();i++){
        if(mesh->is_v_alive(i) && is_vertex_critical(i))
            (*min)[i]=new nNode(i);
    }

}
