#include "formangradientvector.h"

void FormanGradientVector::refine_geometry(double edge_lenght){

    queue<DAG_GeomNode*> ref_queue;

    //creo una coda di possibili raffinamenti iniziata con i valori di root
    for(list<DAG_GeomNode*>::iterator it=geom_root.begin(); it != geom_root.end(); it++){
            ref_queue.push(*it);
    }

    //inizio a raffinarli se mi servono
    while(!ref_queue.empty()){

        DAG_GeomNode* node = ref_queue.front();
        ref_queue.pop();

        if(node->isVisited()) continue;
        node->setVisited();



        if(!node->isRefined() && node->getEdgeLenght() >= edge_lenght){

            if(!node->isRefinable(dag_per_vertex)){
                getRefinable(node, dag_per_vertex);
            }

            vector<int> vertices = node->getUpdatedVertices(dag_per_vertex);

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

            node->setV1(refine_mesh(adj_t1_1,adj_t1_2,adj_t2_1,adj_t2_2, real_index(vertices[4]), node->getCoordsv2(), node->getCoordsv1(), node->get_to_switch()));
            node->setRefined();
        }

        vector<DAG_GeomNode*> dep_op = node->get_dep();
        for(int i=0; i<dep_op.size(); i++)
            if(!dep_op[i]->isVisited())
                ref_queue.push(dep_op[i]);


        if(refined_geometry >= tot_geom*edge_lenght)
            break;
    }


    set_alive_simplexes();
}

int FormanGradientVector::refine_topology(double pers){

    queue<DAG_TopoNode*> ref_queue;

    refined_topo=refined_geometry=0;

    //creo una coda di possibili raffinamenti iniziata con i valori di root
    cout << "Siamo qui " << pers << endl;
    for(list<DAG_TopoNode*>::iterator it=topo_root.begin(); it != topo_root.end(); it++){
            ref_queue.push(*it);
    }

    //inizio a raffinarli se mi servono
    while(!ref_queue.empty()){

        DAG_TopoNode* node = ref_queue.front();
        ref_queue.pop();

        if(node->isVisited()) continue;
        node->setVisited();

        if(!node->isRefined() /*&& node->getHeight() >= pers*/){
            if(!node->isRefinable()){
                getRefinable(node);
            }

            if(node->isExpansion()){
                expand(node);
            }
            else{
                insert(node);
            }
            refined_topo++;
            node->setRefined();
        }

        vector<DAG_TopoNode*> dep_op = node->get_dep();
        for(int i=0; i<dep_op.size(); i++)
            if(!dep_op[i]->isVisited())
                ref_queue.push(dep_op[i]);


        if(refined_topo >= tot_topo*pers)
            break;
    }


    //set il numero di triangoli e vertici vivi correttamente
    set_alive_simplexes();

    cout << refined_geometry << " " << refined_topo << endl;

}

void FormanGradientVector::refine_geometry_box(double edge_lenght, vector<double> box){

    queue<DAG_GeomNode*> ref_queue;

    //creo una coda di possibili raffinamenti iniziata con i valori di root
    for(list<DAG_GeomNode*>::iterator it=geom_root.begin(); it != geom_root.end(); it++){
        //if((*it)->getEdgeLenght() >= edge_lenght &&  is_inside_box(*it, box))
            ref_queue.push(*it);
    }

    //inizio a raffinarli se mi servono
    while(!ref_queue.empty()){

        DAG_GeomNode* node = ref_queue.front();
        ref_queue.pop();

        if(node->isVisited()) continue;
        node->setVisited();

        if(!node->isRefined()  && node->getEdgeLenght() >= edge_lenght && is_inside_box(node, box) ){
            if(!node->isRefinable(dag_per_vertex)){
                getRefinable(node, dag_per_vertex);
            }

            vector<int> vertices = node->getUpdatedVertices(dag_per_vertex);

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

            node->setV1(refine_mesh(adj_t1_1,adj_t1_2,adj_t2_1,adj_t2_2, real_index(vertices[4]), node->getCoordsv2(), node->getCoordsv1(), node->get_to_switch()));
            node->setRefined();
            refined_geometry++;
        }

        vector<DAG_GeomNode*> dep_op = node->get_dep();
        for(int i=0; i<dep_op.size(); i++)
            if(!dep_op[i]->isVisited())
                ref_queue.push(dep_op[i]);

    }

}

int FormanGradientVector::refine_topology_box(double pers, vector<double> box){


    refined_topo=0;
    refined_geometry=0;
    queue<DAG_TopoNode*> ref_queue;

    //creo una coda di possibili raffinamenti iniziata con i valori di root
    for(list<DAG_TopoNode*>::iterator it=topo_root.begin(); it != topo_root.end(); it++){
        //if((*it)->getHeight() >= pers && is_inside_box(*it, box))
            ref_queue.push(*it);
    }

    //inizio a raffinarli se mi servono
    while(!ref_queue.empty()){

        DAG_TopoNode* node = ref_queue.front();
        ref_queue.pop();

        if(node->isVisited()) continue;

        node->setVisited();

        if(!node->isRefined() && node->getHeight() >= pers && is_inside_box(node, box)){
            if(!node->isRefinable()){
                getRefinable(node);
            }
            if(node->isExpansion()){
                expand(node);
            }
            else{
                insert(node);
            }
            node->setRefined();
            refined_topo++;
        }

        vector<DAG_TopoNode*> dep_op = node->get_dep();
        for(int i=0; i<dep_op.size(); i++)
            if(!dep_op[i]->isVisited())
                ref_queue.push(dep_op[i]);
    }
}

void FormanGradientVector::refine_geometry_sphere(double edge_lenght, vector<double> center, double radius){

    queue<DAG_GeomNode*> ref_queue;

    //creo una coda di possibili raffinamenti iniziata con i valori di root
    for(list<DAG_GeomNode*>::iterator it=geom_root.begin(); it != geom_root.end(); it++){
        if((*it)->getEdgeLenght() >= edge_lenght &&  is_inside_sphere(*it, center, radius))
            ref_queue.push(*it);
    }

    //inizio a raffinarli se mi servono
    while(!ref_queue.empty()){

        DAG_GeomNode* node = ref_queue.front();
        ref_queue.pop();

        if(node->isVisited()) continue;
        node->setVisited();

        if(node->isRefined()){
            //qui mi serve avere una lista di operazioni che dipendono da questa e che vorrei poter fare
            //queste vengono aggiunte qualora necessario alla coda;
        }
        else{
            if(!node->isRefinable(dag_per_vertex)){
                getRefinable(node, dag_per_vertex);
            }

            vector<int> vertices = node->getUpdatedVertices(dag_per_vertex);

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

            node->setV1(refine_mesh(adj_t1_1,adj_t1_2,adj_t2_1,adj_t2_2, real_index(vertices[4]),  node->getCoordsv2(), node->getCoordsv1(), node->get_to_switch()));
            node->setRefined();
        }

        vector<DAG_GeomNode*> dep_op = node->get_dep();
        for(int i=0; i<dep_op.size(); i++)
            if(dep_op[i]->getEdgeLenght() >= edge_lenght && !dep_op[i]->isVisited() &&  is_inside_sphere(dep_op[i], center, radius))
                ref_queue.push(dep_op[i]);

    }

}

int FormanGradientVector::refine_topology_sphere(double pers, vector<double> center, double radius){

    queue<DAG_TopoNode*> ref_queue;

    //creo una coda di possibili raffinamenti iniziata con i valori di root
    for(list<DAG_TopoNode*>::iterator it=topo_root.begin(); it != topo_root.end(); it++){
        if((*it)->getHeight() >= pers && is_inside_sphere(*it, center, radius))
            ref_queue.push(*it);
    }

    //inizio a raffinarli se mi servono
    while(!ref_queue.empty()){

        DAG_TopoNode* node = ref_queue.front();
        ref_queue.pop();

        if(node->isVisited()) continue;

        node->setVisited();

        if(node->isRefined()){

        }
        else{
            if(!node->isRefinable()){
                getRefinable(node);
            }
            if(node->isExpansion()){
                expand(node);
            }
            else{
                insert(node);
            }
            node->setRefined();
        }

        vector<DAG_TopoNode*> dep_op = node->get_dep();
        for(int i=0; i<dep_op.size(); i++)
            if(dep_op[i]->getHeight() >= pers && !dep_op[i]->isVisited() && is_inside_sphere(dep_op[i], center, radius))
                ref_queue.push(dep_op[i]);
    }
}
