#include "formangradientvector.h"


int FormanGradientVector::real_index(int v_index){
     if(v_index < 0)
        return mesh->getUpdateVertexIndex(-v_index);
    else return v_index;
}

void FormanGradientVector::getRefinable(DAG_GeomNode* node, vector<DAG_GeomNode*>* dag_per_vertex){
    vector<int> vert = node->getVV();

    for(int i=0; i<vert.size(); i++){
        if((*dag_per_vertex)[vert[i]] != NULL && !(*dag_per_vertex)[vert[i]]->isRefined()){
            if(!(*dag_per_vertex)[vert[i]]->isRefinable(dag_per_vertex)){
                getRefinable((*dag_per_vertex)[vert[i]], dag_per_vertex);
            }

            vector<int> vertices = (*dag_per_vertex)[vert[i]]->getUpdatedVertices(dag_per_vertex);

            int adj_t1_1,adj_t1_2,adj_t2_1,adj_t2_2;
            adj_t1_1 = adj_t1_2 = adj_t2_1 = adj_t2_2 = -1;

            vector<int> vt = mesh->VT(real_index(vertices[4]));

            for(int i=0; i<vt.size(); i++){


                if(mesh->getTopSimplex(vt[i]).contains(real_index(vertices[0])) && mesh->getTopSimplex(vt[i]).contains(real_index(vertices[1]))){
                    adj_t1_1 = vt[i];
                }
                else if(mesh->getTopSimplex(vt[i]).contains(real_index(vertices[0]))){
                    adj_t1_2 = vt[i];
                }

                if(mesh->getTopSimplex(vt[i]).contains(real_index(vertices[2])) && mesh->getTopSimplex(vt[i]).contains(real_index(vertices[3]))){
                    adj_t2_1 = vt[i];
                }
                else if(mesh->getTopSimplex(vt[i]).contains(real_index(vertices[2]))){
                    adj_t2_2 = vt[i];
                }
            }

            refined_geometry++;
            (*dag_per_vertex)[vert[i]]->setV1(refine_mesh(adj_t1_1,adj_t1_2,adj_t2_1,adj_t2_2, real_index(vertices[4]), (*dag_per_vertex)[vert[i]]->getCoordsv2(), (*dag_per_vertex)[vert[i]]->getCoordsv1(), (*dag_per_vertex)[vert[i]]->get_to_switch()));
            (*dag_per_vertex)[vert[i]]->setRefined();
        }
    }
}


int FormanGradientVector::refine_mesh(int adj_t1_1, int adj_t1_2, int adj_t2_1, int adj_t2_2, int v1, vector<double> coordsv2, vector<double> coordsv1, pair<bool,bool> to_switch){


    int new_t1=-1;
    int new_t2=-1;
    int new_v=-1;

    bool up_direction=false;
    int tri_up=-1;

    int v3_sinistra;
    int v3_destra;

    int t_pointed_sinistra=0;
    int t_pointed_destra=0;

    int tail_sin=-1;
    int tail_des=-1;
    int head_sin=-1;
    int head_des=-1;

    int pointed_by_v1=-1;
    if(getVE(v1)!=NULL){
        pointed_by_v1=getVE(v1)->EV(1);
    }


    for(int i=0; i<3; i++){
        if(mesh->getTopSimplex(adj_t1_1).contains(mesh->getTopSimplex(adj_t1_2).TV(i)) && mesh->getTopSimplex(adj_t1_2).TV(i) != v1){
            v3_sinistra = mesh->getTopSimplex(adj_t1_2).TV(i);
            break;
        }
    }

    for(int i=0; i<3; i++){
        if(mesh->getTopSimplex(adj_t2_1).contains(mesh->getTopSimplex(adj_t2_2).TV(i)) && mesh->getTopSimplex(adj_t2_2).TV(i) != v1){
            v3_destra = mesh->getTopSimplex(adj_t2_2).TV(i);
            break;
        }
    }

//    //assert(mesh->getTopSimplex(adj_t1_1).contains(v3_sinistra) &&
//           mesh->getTopSimplex(adj_t1_1).contains(v1));
//    //assert(mesh->getTopSimplex(adj_t1_2).contains(v3_sinistra) &&
//           mesh->getTopSimplex(adj_t1_2).contains(v1));

//    //assert(mesh->getTopSimplex(adj_t2_1).contains(v3_destra) &&
//           mesh->getTopSimplex(adj_t2_1).contains(v1));
//    //assert(mesh->getTopSimplex(adj_t2_2).contains(v3_destra) &&
//           mesh->getTopSimplex(adj_t2_2).contains(v1));


    t_pointed_sinistra = getEF(new Edge(v1,v3_sinistra));
    t_pointed_destra = getEF(new Edge(v1,v3_destra));

//    //assert(getEF(new Edge(v1,v3_sinistra)) != -1 ||
//           (getVE(v1) != NULL && getVE(v1)->EV(1) == v3_sinistra) ||
//           (getVE(v3_sinistra) != NULL && getVE(v3_sinistra)->EV(1) == v1) ||
//            is_edge_critical(v1,v3_sinistra));

//    //assert(getEF(new Edge(v1,v3_destra)) != -1 ||
//           (getVE(v1) != NULL && getVE(v1)->EV(1) == v3_destra) ||
//           (getVE(v3_destra) != NULL && getVE(v3_destra)->EV(1) == v1) ||
//            is_edge_critical(v1,v3_destra));

//    cout << adj_t1_1 << " " << adj_t1_2 << " " << adj_t2_1 << " " << adj_t2_2 << " " << endl;


    if(t_pointed_sinistra == -1){
        if(!is_edge_critical(v1,v3_sinistra)){
            Edge* ar = getVE(v3_sinistra);
            if(ar==NULL || ar->EV(1) != v1){
                ar = getVE(v1);
                //assert(ar->EV(1) == v3_sinistra);
            }
            else{
                //assert(ar->EV(1) == v1);
            }
            tail_sin = ar->EV(0);
            head_sin = ar->EV(1);
            delete ar;
        }
    }

//    cout << getEF(new Edge(v1,v3_destra)) << " "
//         << is_edge_critical(v1, v3_destra) << " "
//         << getVE(v1) << " "
//         << getVE(v3_destra) << endl;

    if(t_pointed_destra == -1){
        if(!is_edge_critical(v1,v3_destra)){
            Edge* ar=getVE(v3_destra);
            if(ar==NULL || ar->EV(1) != v1){
                ar = getVE(v1);
                //assert(ar->EV(1) == v3_destra);
            }
            else{
                //assert(ar->EV(1) == v1);
            }
//            cout << "la freccia coinvolge " << ar->EV(0) << " " << ar->EV(1) << endl;
            tail_des = ar->EV(0);
            head_des = ar->EV(1);
            delete ar;
        }
    }

    //nuovo vertice che viene inserito
    new_v = mesh->getNumVertex();
    Vertex3D v = Vertex3D(coordsv2[0], coordsv2[1], coordsv2[2]);
    mesh->addVertex(v);


    mesh->getVertex(v1).setX(coordsv1[0]);
    mesh->getVertex(v1).setY(coordsv1[1]);
    mesh->getVertex(v1).setZ(coordsv1[2]);


//    //assert(mesh->getTopSimplex(adj_t1_1).contains(v3_sinistra));
//    //assert(mesh->getTopSimplex(adj_t1_2).contains(v3_sinistra));
//    //assert(mesh->getTopSimplex(adj_t2_1).contains(v3_destra));
//    //assert(mesh->getTopSimplex(adj_t2_2).contains(v3_destra));

//    cout << getVE(v3_sinistra) << " " << getVE(v1) << " " << getEF(new Edge(v1,v3_sinistra)) << endl;
//    cout << getVE(v3_destra) << " " << getVE(v1) << " " << getEF(new Edge(v1,v3_destra)) << endl;
//    if(getVE(v1) != NULL) cout << getVE(v1)->EV(0) << " " << getVE(v1)->EV(1) << endl;


//    cout << "fine setup" << endl;

    if(t_pointed_sinistra == -1 && !is_edge_critical(v1,v3_sinistra)){
        ////assert(tail_sin != -1 && head_sin != -1);
        freeVE(tail_sin,head_sin);
    }

    if(t_pointed_destra == -1  && !is_edge_critical(v1,v3_destra)){
        ////assert(tail_des != -1 && head_des != -1);
        freeVE(tail_des,head_des);
    }

//    cout << "fine free relations" << endl;
    //----------------------------------------------------------------------------------------------------------
    //operazioni per il lato sinistro
    //operazioni da gradiente

    //modifiche geometriche
    //nuovo triangolo che viene inserito a sinistra
    new_t1=mesh->getTopSimplexesNum();
    new_t2=new_t1+1;
    Triangle t = Triangle(v1,v3_sinistra,new_v);
    mesh->addTopSimplex(t);

    mesh->getTopSimplex(new_t1).setTT(0,adj_t1_1);
    mesh->getTopSimplex(new_t1).setTT(1,new_t2);
    mesh->getTopSimplex(new_t1).setTT(2,adj_t1_2);

    //modifiche su adj_t1_1;
    mesh->getTopSimplex(adj_t1_1).setTV(mesh->getTopSimplex(adj_t1_1).vertex_index(v1),new_v);
    int pos = 3-mesh->getTopSimplex(adj_t1_1).vertex_index(new_v)-mesh->getTopSimplex(adj_t1_1).vertex_index(v3_sinistra);
    mesh->getTopSimplex(adj_t1_1).setTT(pos,new_t1);

    //modifiche su adj_t1_2;
    pos = 3-mesh->getTopSimplex(adj_t1_2).vertex_index(v1)-mesh->getTopSimplex(adj_t1_2).vertex_index(v3_sinistra);
    mesh->getTopSimplex(adj_t1_2).setTT(pos,new_t1);


    //----------------------------------------------------------------------------------------------------------
    //operazioni per il lato destro
    //operazioni da gradiente

    //modifiche geometriche
    //nuovo triangolo che viene inserito a sinistra

    Triangle t1 = Triangle(v3_destra,v1,new_v);
    mesh->addTopSimplex(t1);


    mesh->getTopSimplex(new_t2).setTT(0,new_t1);
    mesh->getTopSimplex(new_t2).setTT(1,adj_t2_1);
    mesh->getTopSimplex(new_t2).setTT(2,adj_t2_2);


    //modifiche su adj_t2_1;
    mesh->getTopSimplex(adj_t2_1).setTV(mesh->getTopSimplex(adj_t2_1).vertex_index(v1),new_v);
    pos = 3-mesh->getTopSimplex(adj_t2_1).vertex_index(new_v)-mesh->getTopSimplex(adj_t2_1).vertex_index(v3_destra);
    mesh->getTopSimplex(adj_t2_1).setTT(pos,new_t2);

    //modifiche su adj_t2_2;
    pos = 3-mesh->getTopSimplex(adj_t2_2).vertex_index(v1)-mesh->getTopSimplex(adj_t2_2).vertex_index(v3_destra);
    mesh->getTopSimplex(adj_t2_2).setTT(pos,new_t2);

    //----------------------------------------------------------------------------------------------------------


    int current = adj_t1_1;
    if(mesh->getTopSimplex(current).contains(pointed_by_v1) && pointed_by_v1 != v3_destra && pointed_by_v1 != v3_sinistra){
        up_direction=true;
        tri_up=current;
    }

    int next = mesh->getTopSimplex(current).TT(mesh->getTopSimplex(current).vertex_index(v3_sinistra));
    while(next != -1 && next != adj_t2_1){

        if(mesh->getTopSimplex(next).contains(pointed_by_v1) && pointed_by_v1 != v3_destra && pointed_by_v1 != v3_sinistra){
            up_direction=true;
            tri_up=next;
        }


        pos = mesh->getTopSimplex(next).vertex_index(v1);
        mesh->getTopSimplex(next).setTV(pos, new_v);

        if(mesh->getTopSimplex(next).TT((pos+1)%3) == current){
            current = next;
            next = mesh->getTopSimplex(next).TT((pos+2)%3);
        }
        else{
            current = next;
            next = mesh->getTopSimplex(next).TT((pos+1)%3);
        }
    }

    if(next == -1){
        //assert(false);
        int current = adj_t2_1;
        int next = mesh->getTopSimplex(current).TT(mesh->getTopSimplex(current).vertex_index(v3_destra));
        while(next != -1 && next != adj_t1_1){
            pos = mesh->getTopSimplex(next).vertex_index(v1);
            mesh->getTopSimplex(next).setTV(pos, new_v);
            if(mesh->getTopSimplex(next).TT((pos+1)%3) == current){
                current = next;
                next = mesh->getTopSimplex(next).TT((pos+2)%3);
            }
            else{
                current = next;
                next = mesh->getTopSimplex(next).TT((pos+1)%3);
            }
        }
    }

    forman_gradient.push_back(0);
    forman_gradient.push_back(0);

    if(up_direction){
        //assert(mesh->getTopSimplex(tri_up).contains(new_v));
        mesh->getVertex(new_v).VTstar(tri_up);
        mesh->getVertex(v1).VTstar(new_t1);
        setVE(v1,new_v);
    }
    else{
        mesh->getVertex(new_v).VTstar(new_t1);
        setVE(new_v,v1);
    }

    if(!mesh->getTopSimplex(mesh->getVertex(v1).VTstar()).contains(v1)){
        mesh->getVertex(v1).VTstar(new_t1);
    }


    if(t_pointed_sinistra == -1){
        if(head_sin==-1 && tail_sin==-1){
            setEF(v1,new_t1);
        }
        else if(to_switch.first){
            //assert(false);
            //assert(tail_sin==v3_sinistra);
            mesh->getVertex(tail_sin).VTstar(new_t1);
            setVE(tail_sin,new_v);
            setEF(new_v,new_t1);
        }
        else{
            mesh->getVertex(tail_sin).VTstar(new_t1);
            setVE(tail_sin,head_sin);
            setEF(v1,new_t1);
        }

    }
    else if(t_pointed_sinistra == adj_t1_1){
        setEF(new_v,new_t1);
    }
    else{
        //assert(t_pointed_sinistra == adj_t1_2);
        setEF(v1,new_t1);
    }


    if(t_pointed_destra == -1){
        if(head_des==-1 && tail_des==-1){
            setEF(v1,new_t2);
        }
        else if(to_switch.second){
            //assert(false);
            //assert(tail_des==v3_destra);
            mesh->getVertex(tail_des).VTstar(new_t2);
            setVE(tail_des,new_v);
            setEF(new_v,new_t2);
        }
        else{
            mesh->getVertex(tail_des).VTstar(new_t2);
            setVE(tail_des,head_des);
            setEF(v1,new_t2);
        }
    }
    else if(t_pointed_destra == adj_t2_1){
        setEF(new_v,new_t2);
    }
    else{
        //assert(t_pointed_destra == adj_t2_2);
        setEF(v1,new_t2);
    }

//    //assert(!is_edge_critical(v1,v3_sinistra) && !is_edge_critical(new_v,v3_sinistra) &&
//           !is_edge_critical(v1,v3_destra) && !is_edge_critical(new_v,v3_destra));

    ////assert(valid_gradient());
    return new_v;

}

//---------------------------------------------------------------------------------------------------------------------------------------------
//qui le operazioni per topologia
void FormanGradientVector::getRefinable(DAG_TopoNode* node){

    set<DAG_GeomNode*>* geom_nodes = node->getFathersGeom();
    for(set<DAG_GeomNode*>::iterator it = geom_nodes->begin(); it != geom_nodes->end(); it++){
        if((*it)->isRefined()) continue;

        if(!(*it)->isRefinable(dag_per_vertex)){
            getRefinable(*it,dag_per_vertex);
        }

        vector<int> vertices = (*it)->getUpdatedVertices(dag_per_vertex);

        int adj_t1_1,adj_t1_2,adj_t2_1,adj_t2_2;
        adj_t1_1 = adj_t1_2 = adj_t2_1 = adj_t2_2 = -1;

        vector<int> vt = mesh->VT(real_index(vertices[4]));

        for(int i=0; i<vt.size(); i++){

            if(mesh->getTopSimplex(vt[i]).contains(real_index(vertices[0])) && mesh->getTopSimplex(vt[i]).contains(real_index(vertices[1]))){
                adj_t1_1 = vt[i];
            }
            else if(mesh->getTopSimplex(vt[i]).contains(real_index(vertices[0]))){
                adj_t1_2 = vt[i];
            }

            if(mesh->getTopSimplex(vt[i]).contains(real_index(vertices[2])) && mesh->getTopSimplex(vt[i]).contains(real_index(vertices[3]))){
                adj_t2_1 = vt[i];
            }
            else if(mesh->getTopSimplex(vt[i]).contains(real_index(vertices[2]))){
                adj_t2_2 = vt[i];
            }
        }

        (*it)->setV1(refine_mesh(adj_t1_1,adj_t1_2,adj_t2_1,adj_t2_2, real_index(vertices[4]), (*it)->getCoordsv2(), (*it)->getCoordsv1(), (*it)->get_to_switch()));
        (*it)->setRefined();
        refined_geometry++;
    }
    //qui sopra semplifico la geometria se ne ho bisogno

    set<DAG_TopoNode*>* topo_nodes = node->getFathers();
    for(set<DAG_TopoNode*>::iterator it = topo_nodes->begin(); it != topo_nodes->end(); it++){
        if((*it)->isRefined()) continue;

        if(!(*it)->isRefinable()){
            getRefinable(*it);
        }
        if((*it)->isExpansion()){
            expand(*it);
        }
        else{
            insert(*it);
        }

        refined_topo++;
        (*it)->setRefined();
    }
}

int FormanGradientVector::expand(DAG_TopoNode* node){

    //espansione
    //lavoriamo fra minimi e selle
    //quindi partiamo da un vertice di minimo e navighiamo fino all'edge critico;
    vector<int> all_vert= node->getUpdatedVertices(this->dag_per_vertex);

    //assert(all_vert.size()==4);

    int min_vertex= real_index(all_vert[0]);
    Edge last = Edge(real_index(all_vert[1]), real_index(all_vert[2]));

    //assert(mesh->ET(last).size() > 0);

    //assert(!is_vertex_critical(min_vertex));
    //assert(!is_edge_critical(last.EV(0), last.EV(1)));

    int current_vertex = min_vertex;
    Edge* next_edge = getVE(current_vertex);
    //assert(next_edge != NULL);
    int next_vertex = current_vertex == next_edge->EV(0) ? next_edge->EV(1) : next_edge->EV(0);
    freeVE(current_vertex, next_vertex);


    while(*next_edge != last){

        delete next_edge;
        next_edge = getVE(next_vertex);

        freeVE(next_edge->EV(0),next_edge->EV(1));
        setVE(next_vertex, current_vertex);
        current_vertex = next_vertex;
        next_vertex = current_vertex == next_edge->EV(0) ? next_edge->EV(1) : next_edge->EV(0);
    }

    //assert(is_edge_critical(last.EV(0), last.EV(1)));
    //assert(is_vertex_critical(min_vertex));

    ////assert(valid_gradient());

}

int FormanGradientVector::insert(DAG_TopoNode* node){


    vector<int> all_vert = node->getUpdatedVertices(this->dag_per_vertex);
    //assert(all_vert.size()==6);

    Edge last = Edge(real_index(all_vert[3]), real_index(all_vert[4]));

    //assert(mesh->VT(last.EV(0)).size() > 0);
    //assert(mesh->VT(last.EV(1)).size() > 0);
    //assert(mesh->ET(last).size() > 0);
    //assert(!is_edge_critical(last.EV(0), last.EV(1)));

    int actual_t = -1;
    vector<int> vt = mesh->VT(real_index(all_vert[0]));
    for(int i=0; i<vt.size(); i++){
        if(mesh->getTopSimplex(vt[i]).contains(real_index(all_vert[1])) && mesh->getTopSimplex(vt[i]).contains(real_index(all_vert[2]))){
            actual_t=vt[i];
            break;
        }
    }
    int original_max=actual_t;
    //assert(actual_t != -1);
    //assert(!is_face_critical(original_max));

    int i=0;
    Edge* old_edge;
    for(; i<3; i++){
        old_edge = mesh->getTopSimplex(actual_t).TE(i);
        int tri = getEF(old_edge);
        if(tri == actual_t) break;
    }

    //assert(i!=3);
    int next_triangle = mesh->getTopSimplex(actual_t).TT(i);
    int vert_to_navigate=mesh->getTopSimplex(actual_t).TV(i);
    int vert_to_pair;

    freeEF(vert_to_navigate,actual_t);

    Edge* find_edge;
    while(*old_edge != last){

        int i=0;
        for(; i<3; i++){
            find_edge = mesh->getTopSimplex(next_triangle).TE(i);
            int tri = getEF(find_edge);
            if(tri == next_triangle) break;
            delete find_edge;
        }
        //assert(i!=3);

        actual_t=next_triangle;
        next_triangle=mesh->getTopSimplex(actual_t).TT(i);
        vert_to_navigate = mesh->getTopSimplex(actual_t).TV(i);

        freeEF(vert_to_navigate, actual_t);

        for(int i=0;i<3; i++){
            if(old_edge->EV(0) != mesh->getTopSimplex(actual_t).TV(i) &&
               old_edge->EV(1) != mesh->getTopSimplex(actual_t).TV(i)   ){
                vert_to_pair = mesh->getTopSimplex(actual_t).TV(i);
                break;}
        }

        setEF(vert_to_pair, actual_t);

        delete old_edge;
        old_edge = find_edge;

    }

    //assert(is_edge_critical(old_edge->EV(0),old_edge->EV(1)));
    //assert(is_face_critical(original_max));
    ////assert(valid_gradient());
}




