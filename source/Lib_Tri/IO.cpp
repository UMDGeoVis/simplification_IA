#include "formangradientvector.h"
#include <stdio.h>

void FormanGradientVector::writeVTK_1cells(char* file_name, set<int> vertici, set<pair<int, bool> > vertici_critici, set<pair<int,int> > edges){

    int vertex_number = vertici.size();
    int edge_number = edges.size();

    vector<int> new_vertex_index = vector<int>(mesh->getNumVertex(), -1);
    vector<int> critical_index = vector<int>(mesh->getNumVertex(), -1);

    FILE* file;
    file = fopen(file_name, "w");


    int i=0;
    for(set<int>::iterator it = vertici.begin(); it!=vertici.end(); it++){
            new_vertex_index[*it] = i++;
    }

    for(set<pair<int, bool> >::iterator it = vertici_critici.begin(); it!=vertici_critici.end(); it++){
        if(it->second)   critical_index[new_vertex_index[it->first]] = 1;
        else critical_index[new_vertex_index[it->first]] = 0;

    }

    cout << vertex_number << endl;

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vertex_number);

    for(int i=0; i<mesh->getNumVertex(); i++){
        if(new_vertex_index[i] != -1){
            fprintf(file, "%f %f %f\n", mesh->getVertex(i).getX(), mesh->getVertex(i).getY(), mesh->getVertex(i).getZ());
        }

    }
    fprintf(file, "\n\n");

    //QUA HO FINITO DI SCRIVERE I VERTICI

    fprintf(file, "CELLS %d %d\n", edge_number, edge_number*3);

    for(set<pair<int,int> >::iterator it = edges.begin(); it != edges.end(); it++){
        fprintf(file, "2 %d %d \n", new_vertex_index[it->first], new_vertex_index[it->second]);
    }
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", edge_number);

    for(int i=0; i<edge_number; i++)
        fprintf(file, "%d ", 3);
    fprintf(file, "\n\n");


    fprintf(file, "POINT_DATA %d \n", vertex_number);
    fprintf(file, "FIELD FieldData 2\n");
    fprintf(file, "originalField 1 %d float\n", vertex_number);

    int j=0;
    for(int i=0; i<mesh->getNumVertex(); i++){
        if(new_vertex_index[i] != -1){
            fprintf(file, "%f ", mesh->getVertex(i).getZ());
         j++;
        }
    }

    fprintf(file, "\n\n");
    fprintf(file, "criticalPoint 1 %d int\n", vertex_number);

    for(int i=0; i<critical_index.size(); i++)
        fprintf(file, "%d ", critical_index[i]);

    fprintf(file, "\n\n");


    fclose(file);

}

void FormanGradientVector::writeVTK_2cells(char* nome_file_output, set<pair<int, bool> > critici, vector<int> triangles){

    cout << "arrivo qui" << endl;

    FILE* file;
    file = fopen(nome_file_output, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", mesh->getNumVertex());

    for(int i=0; i<mesh->getNumVertex(); i++)
        fprintf(file, "%f %f %f\n", mesh->getVertex(i).getX(), mesh->getVertex(i).getY(), mesh->getVertex(i).getZ());
    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", mesh->getTopSimplexesNum(), mesh->getTopSimplexesNum()*4);

    for(int i=0; i<mesh->getTopSimplexesNum(); i++)
        fprintf(file, "3 %d %d %d\n", mesh->getTopSimplex(i).TV(0), mesh->getTopSimplex(i).TV(1), mesh->getTopSimplex(i).TV(2));
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", mesh->getTopSimplexesNum());

    for(int i=0; i<mesh->getTopSimplexesNum(); i++)
        fprintf(file, "%d ", 5);
    fprintf(file, "\n\n");


    fprintf(file, "POINT_DATA %d \n", mesh->getNumVertex());
    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "originalfield 1 %d float\n", mesh->getNumVertex());

    for(int i=0; i<mesh->getNumVertex(); i++)
        fprintf(file, "%f ", mesh->getVertex(i).getZ());

    fprintf(file, "\n\n");


    fprintf(file, "CELL_DATA %d \n", mesh->getTopSimplexesNum());
    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "descending3cells 1 %d float\n", mesh->getTopSimplexesNum());

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        fprintf(file, "%d ", triangles[i]);
    }
    fprintf(file, "\n");


    fclose(file);
}

void FormanGradientVector::writeVTK_2cells_on_vert(char* nome_file_output, set<pair<int, bool> > critici, vector<int> triangles){


    FILE* file;
    file = fopen(nome_file_output, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", mesh->getNumVertex());

    for(int i=0; i<mesh->getNumVertex(); i++)
        fprintf(file, "%f %f %f\n", mesh->getVertex(i).getX(), mesh->getVertex(i).getY(), mesh->getVertex(i).getZ());
    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", mesh->getTopSimplexesNum(), mesh->getTopSimplexesNum()*4);

    for(int i=0; i<mesh->getTopSimplexesNum(); i++)
        fprintf(file, "3 %d %d %d\n", mesh->getTopSimplex(i).TV(0), mesh->getTopSimplex(i).TV(1), mesh->getTopSimplex(i).TV(2));
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", mesh->getTopSimplexesNum());

    for(int i=0; i<mesh->getTopSimplexesNum(); i++)
        fprintf(file, "%d ", 5);
    fprintf(file, "\n\n");


    fprintf(file, "POINT_DATA %d \n", mesh->getNumVertex());
    fprintf(file, "FIELD FieldData 3\n");
    fprintf(file, "originalfield 1 %d float\n", mesh->getNumVertex());

    for(int i=0; i<mesh->getNumVertex(); i++)
        fprintf(file, "%f ", mesh->getVertex(i).getZ());

    fprintf(file, "\n");

    fprintf(file, "critic 1 %d float\n", mesh->getNumVertex());

    for(int i=0; i<mesh->getNumVertex(); i++){
        if(critici.find(pair<int,bool>(i,true)) != critici.end())
            fprintf(file, "0 ");
        else if(critici.find(pair<int,bool>(i,false)) != critici.end())
            fprintf(file, "1 ");
        else
            fprintf(file, "-1 ");

    }
    fprintf(file, "\n");

    fprintf(file, "ascending 1 %d float\n", mesh->getNumVertex());

    for(int i=0; i<mesh->getNumVertex(); i++)
        fprintf(file, "%d ", triangles[i]);

    fprintf(file, "\n");

    fclose(file);
}

void FormanGradientVector::writeVTK_1cells_on_tri(char* nome_file_output,  set<pair<int, bool> > critici, map<int,int> triangles){

    FILE* file;
    file = fopen(nome_file_output, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", mesh->getNumVertex());

    for(int i=0; i<mesh->getNumVertex(); i++)
        fprintf(file, "%f %f %f\n", mesh->getVertex(i).getX(), mesh->getVertex(i).getY(), mesh->getVertex(i).getZ());
    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", mesh->getTopSimplexesNum(), mesh->getTopSimplexesNum()*4);

    for(int i=0; i<mesh->getTopSimplexesNum(); i++)
        fprintf(file, "3 %d %d %d\n", mesh->getTopSimplex(i).TV(0), mesh->getTopSimplex(i).TV(1), mesh->getTopSimplex(i).TV(2));
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", mesh->getTopSimplexesNum());

    for(int i=0; i<mesh->getTopSimplexesNum(); i++)
        fprintf(file, "%d ", 5);
    fprintf(file, "\n\n");


    fprintf(file, "POINT_DATA %d \n", mesh->getNumVertex());
    fprintf(file, "FIELD FieldData 2\n");
    fprintf(file, "originalfield 1 %d float\n", mesh->getNumVertex());

    for(int i=0; i<mesh->getNumVertex(); i++)
        fprintf(file, "%f ", mesh->getVertex(i).getZ());

    fprintf(file, "\n\n");

    fprintf(file, "critic 1 %d float\n", mesh->getNumVertex());

    for(int i=0; i<mesh->getNumVertex(); i++){
        if(critici.find(pair<int,bool>(i,true)) != critici.end())
            fprintf(file, "0 ");
        else if(critici.find(pair<int,bool>(i,false)) != critici.end())
            fprintf(file, "1 ");
        else
            fprintf(file, "-1 ");

    }
    fprintf(file, "\n");


    fprintf(file, "CELL_DATA %d \n", mesh->getTopSimplexesNum());
    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "ascending1cells 1 %d float\n", mesh->getTopSimplexesNum());

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        if(triangles.find(i) == triangles.end())
            fprintf(file, "-1 ");
        else
            fprintf(file, "%d ", triangles[i]);
    }
    fprintf(file, "\n");


    fclose(file);

}

void FormanGradientVector::writeVTK_IG(char* file_name){

    vector<int> new_vertex_index = vector<int>(mesh->getNumVertex(), -1);
    vector<bool> connected = vector<bool>(mesh->getNumVertex(), false);
    set<pair<int,int> > arcs = set<pair<int,int> >();

//    compute_incidence_graph();
//    cout << forman_ig.getMaxima().size() << endl;
//    cout << forman_ig.getSaddle().size() << endl;
//    cout << forman_ig.getMinima().size() << endl;

    int vertex_number = 0;
    int edge_number =0;

    vector<int> critici = vector<int>(mesh->getNumVertex(), -1);

    vector<nNode*> minima = forman_ig.getMinima();
    for(int i=0; i<minima.size(); i++){
        nNode* node = minima[i];
        critici[node->getCriticalIndex()]=0;
    }

    vector<nNode*> maxima = forman_ig.getMaxima();
    for(int i=0; i<maxima.size(); i++){
        nNode* node = maxima[i];
        critici[node->getCriticalIndex()]=2;
    }

    vector<iNode*> saddle = forman_ig.getSaddle();
    for(int i=0; i<saddle.size(); i++){
        iNode* node = saddle[i];
        int sad_vertex = node->getCriticalIndex();
        critici[sad_vertex]=1;

        vector<Arc*> arcs_up = node->getArcs(true);
        for(int j=0; j<arcs_up.size(); j++){
            arcs.insert(pair<int,int>(((nNode*)arcs_up[j]->getNode_i())->getCriticalIndex() ,sad_vertex));
        }

        vector<Arc*> arcs_down = node->getArcs(false);
        for(int j=0; j<arcs_down.size(); j++){
            arcs.insert(pair<int,int>(sad_vertex, ((nNode*)arcs_down[j]->getNode_j())->getCriticalIndex()));
        }
    }




    for(set<pair<int,int> >::iterator it = arcs.begin(); it!=arcs.end(); it++){
        connected[it->first] = true;
        connected[it->second] = true;
    }

    for(int i=0; i<mesh->getNumVertex(); i++){
        if(critici[i] != -1 && connected[i]){
            new_vertex_index[i] = vertex_number;
            vertex_number++;
        }
    }

    edge_number = arcs.size();

    FILE* file;
    file = fopen(file_name, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vertex_number);

    for(int i=0; i<mesh->getNumVertex(); i++){
        if(critici[i] != -1 && connected[i]){
            fprintf(file, "%f %f %f\n", mesh->getVertex(i).getX(), mesh->getVertex(i).getY(), mesh->getVertex(i).getZ());
        }

    }
    fprintf(file, "\n\n");

    //QUA HO FINITO DI SCRIVERE I VERTICI

    fprintf(file, "CELLS %d %d\n", edge_number, edge_number*3);

    for(set<pair<int,int> >::iterator it = arcs.begin(); it != arcs.end(); it++){
        fprintf(file, "2 %d %d \n", new_vertex_index[it->first], new_vertex_index[it->second]);
    }
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", edge_number);

    for(int i=0; i<edge_number; i++)
        fprintf(file, "%d ", 3);
    fprintf(file, "\n\n");


    fprintf(file, "POINT_DATA %d \n", vertex_number);
    fprintf(file, "FIELD FieldData 2\n");
    fprintf(file, "originalField 1 %d float\n", vertex_number);

    int j=0;
    for(int i=0; i<mesh->getNumVertex(); i++){
        if(critici[i] != -1 && connected[i]){
            fprintf(file, "%f ", mesh->getVertex(i).getZ());
         j++;
        }
    }

    fprintf(file, "\n\n");
    fprintf(file, "criticalPoint 1 %d int\n", vertex_number);

    for(int i=0; i<mesh->getNumVertex(); i++)
        if(critici[i] != -1 && connected[i])
            fprintf(file, "%d ", critici[i]);

    fprintf(file, "\n\n");


    fclose(file);
}


void FormanGradientVector::writeVTK_gradient(char* nomeFile){

    vector<int> critical_vertexes = vector<int>(mesh->getNumVertex(),-1);

    vector<vector<float> > new_vertexes_triangles;

    map<pair<int,int>, vector<float> > edges;

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        if(edges.find(pair<int,int>(mesh->getTopSimplex(i).TV(0), mesh->getTopSimplex(i).TV(1))) == edges.end() && edges.find(pair<int,int>(mesh->getTopSimplex(i).TV(1), mesh->getTopSimplex(i).TV(0))) == edges.end())
            edges[pair<int,int>(mesh->getTopSimplex(i).TV(0), mesh->getTopSimplex(i).TV(1))] = vector<float>();
        if(edges.find(pair<int,int>(mesh->getTopSimplex(i).TV(1), mesh->getTopSimplex(i).TV(2))) == edges.end() && edges.find(pair<int,int>(mesh->getTopSimplex(i).TV(2), mesh->getTopSimplex(i).TV(1))) == edges.end())
            edges[pair<int,int>(mesh->getTopSimplex(i).TV(1), mesh->getTopSimplex(i).TV(2))] = vector<float>();
        if(edges.find(pair<int,int>(mesh->getTopSimplex(i).TV(2), mesh->getTopSimplex(i).TV(0))) == edges.end() && edges.find(pair<int,int>(mesh->getTopSimplex(i).TV(0), mesh->getTopSimplex(i).TV(2))) == edges.end())
            edges[pair<int,int>(mesh->getTopSimplex(i).TV(2), mesh->getTopSimplex(i).TV(0))] = vector<float>();
    }

    for(map<pair<int,int>, vector<float> >::iterator it = edges.begin(); it != edges.end(); it++){

        float x = (mesh->getVertex(it->first.first).getX() + mesh->getVertex(it->first.second).getX())/2.0;
        float y = (mesh->getVertex(it->first.first).getY() + mesh->getVertex(it->first.second).getY())/2.0;
        float z = (mesh->getVertex(it->first.first).getZ() + mesh->getVertex(it->first.second).getZ())/2.0;

        vector<float> bar;
        bar.push_back(x);
        bar.push_back(y);
        bar.push_back(z);

        edges[it->first] = bar;
    }


    vector<vector<float> > tri_baricenter;

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){

        float x = (mesh->getVertex(mesh->getTopSimplex(i).TV(0)).getX() + mesh->getVertex(mesh->getTopSimplex(i).TV(1)).getX()+ mesh->getVertex(mesh->getTopSimplex(i).TV(2)).getX())/3.0;
        float y = (mesh->getVertex(mesh->getTopSimplex(i).TV(0)).getY() + mesh->getVertex(mesh->getTopSimplex(i).TV(1)).getY()+ mesh->getVertex(mesh->getTopSimplex(i).TV(2)).getY())/3.0;
        float z = (mesh->getVertex(mesh->getTopSimplex(i).TV(0)).getZ() + mesh->getVertex(mesh->getTopSimplex(i).TV(1)).getZ()+ mesh->getVertex(mesh->getTopSimplex(i).TV(2)).getZ())/3.0;

        vector<float> bar;
        bar.push_back(x);
        bar.push_back(y);
        bar.push_back(z);

        tri_baricenter.push_back(bar);


        if(is_face_critical(i))
            new_vertexes_triangles.push_back(tri_baricenter[i]);
    }


    vector<vector<float> > vectors_gradient;
    for(int i=0; i<mesh->getNumVertex(); i++ ){

        if(is_vertex_critical(i)){
            critical_vertexes[i] = 0;
        }

        vector<int> vert;
        vert.push_back(i);

        vector<float> vect;
        Edge* e = getVE(i);
        if(e != NULL){

            vector<int> edge;
            edge.push_back(e->EV(0));
            edge.push_back(e->EV(1));

            if(edges.find(pair<int,int>(edge[0], edge[1])) != edges.end()){

                vector<float> bar_edge = edges.find(pair<int,int>(edge[0], edge[1]))->second;
                vector<float> bar_point;
                bar_point.push_back(mesh->getVertex(i).getX());
                bar_point.push_back(mesh->getVertex(i).getY());
                bar_point.push_back(mesh->getVertex(i).getZ());

                vect.push_back(bar_edge[0] - bar_point[0]);
                vect.push_back(bar_edge[1] - bar_point[1]);
                vect.push_back(bar_edge[2] - bar_point[2]);
            }
            else{
                //assert(edges.find(pair<int,int>(edge[1], edge[0])) != edges.end());

                vector<float> bar_edge = edges.find(pair<int,int>(edge[1], edge[0]))->second;
                vector<float> bar_point;
                bar_point.push_back(mesh->getVertex(i).getX());
                bar_point.push_back(mesh->getVertex(i).getY());
                bar_point.push_back(mesh->getVertex(i).getZ());

                vect.push_back(bar_edge[0] - bar_point[0]);
                vect.push_back(bar_edge[1] - bar_point[1]);
                vect.push_back(bar_edge[2] - bar_point[2]);
            }

        }
        else{

            vect.push_back(0.0);
            vect.push_back(0.0);
            vect.push_back(0.0);
        }

        vectors_gradient.push_back(vect);
    }


    vector<vector<float> > new_vertexes;
    vector<vector<float> > new_vectors;
    vector<int> from_edges;
    for(map<pair<int,int>, vector<float> >::iterator it = edges.begin(); it != edges.end(); it++){

        vector<int> edge;
        edge.push_back(it->first.first);
        edge.push_back(it->first.second);

        edge = sort_simplex(&edge);
        int tri = getEF(new Edge(edge[0],edge[1]));
        if(tri != -1){
            vector<float> new_vert = it->second;

            int t_ind = tri;
            vector<float> vert_vector;

            new_vertexes.push_back(new_vert);

            vert_vector.push_back(tri_baricenter[t_ind][0] -new_vert[0]);
            vert_vector.push_back(tri_baricenter[t_ind][1] -new_vert[1]);
            vert_vector.push_back(tri_baricenter[t_ind][2] -new_vert[2]);

            new_vectors.push_back(vert_vector);

            from_edges.push_back(-1);
        }
        else{
            //qua faccio il controllo nel caso in cui siamo su una 1-sella
            if(is_edge_critical(edge[0], edge[1])){

                vector<float> new_vert = it->second;
                new_vertexes.push_back(new_vert);

                from_edges.push_back(1);
                vector<float> vert_vector;

                vert_vector.push_back(0);
                vert_vector.push_back(0);
                vert_vector.push_back(0);

                new_vectors.push_back(vert_vector);
            }
        }
    }

    //assert(new_vertexes.size() == new_vectors.size() && new_vertexes.size() == from_edges.size());

    FILE* file;
    file = fopen(nomeFile, "w");

    int vertex_number = mesh->getNumVertex()+new_vertexes.size()+new_vertexes_triangles.size();
    int triangles = mesh->getTopSimplexesNum();

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vertex_number);

    for(int i=0; i<mesh->getNumVertex(); i++){
            fprintf(file, "%f %f %f\n", mesh->getVertex(i).getX(), mesh->getVertex(i).getY(), mesh->getVertex(i).getZ());
    }
    for(int i=0; i<new_vertexes.size(); i++){
        fprintf(file, "%f %f %f\n", new_vertexes[i][0], new_vertexes[i][1], new_vertexes[i][2]);
    }
    for(int i=0; i<new_vertexes_triangles.size(); i++){
        fprintf(file, "%f %f %f\n", new_vertexes_triangles[i][0], new_vertexes_triangles[i][1], new_vertexes_triangles[i][2]);
    }
    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", triangles, triangles*4);

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        fprintf(file, "3 %d %d %d \n", mesh->getTopSimplex(i).TV(0), mesh->getTopSimplex(i).TV(1), mesh->getTopSimplex(i).TV(2));
    }
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", triangles);

    for(int i=0; i<triangles; i++)
        fprintf(file, "%d ", 5);
    fprintf(file, "\n\n");


    fprintf(file, "POINT_DATA %d \n", vertex_number);
    fprintf(file, "VECTORS vector float\n");

    for(int i=0; i<mesh->getNumVertex(); i++){
        fprintf(file, "%f %f %f   ", vectors_gradient[i][0], vectors_gradient[i][1], vectors_gradient[i][2]);
    }
    for(int i=0; i<new_vectors.size(); i++){
        fprintf(file, "%f %f %f   ", new_vectors[i][0], new_vectors[i][1], new_vectors[i][2]);
    }
    for(int i=0; i<new_vertexes_triangles.size(); i++){
        fprintf(file, "%f %f %f   ", 0, 0, 0);
    }
    fprintf(file, "\n\n");

    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "critical 1 %d int\n", vertex_number);

    for(int i=0; i<mesh->getNumVertex(); i++){
        fprintf(file, "%d ", critical_vertexes[i]);
    }
    for(int i=0; i<from_edges.size(); i++){
        fprintf(file, "%d ", from_edges[i]);
    }
    for(int i=0; i<new_vertexes_triangles.size(); i++){
        fprintf(file, "%d ", 2);
    }
    fprintf(file, "\n\n");


//    int j=0;
//    for(int i=0; i<mesh->getNumVertex(); i++){
//        if(new_vertex_index[i] != -1){
//            fprintf(file, "%f ", mesh->getVertex(i).getZ());
//         j++;
//        }
//    }
//    cout << j << " " << vertex_number << endl;

//    fprintf(file, "\n\n");
//    fprintf(file, "criticalPoint 1 %d int\n", vertex_number);

//    for(int i=0; i<new_critici.size(); i++)
//        fprintf(file, "%d ", new_critici[i]);

//    fprintf(file, "\n\n");

//    fprintf(file, "CELL_DATA %d \n", triangles.size());
//    fprintf(file, "FIELD FieldData 1\n");
//    fprintf(file, "descending_2_cells 1 %d int \n", triangles.size());

//    for(set<pair<vector<int>,int > >::iterator it = triangles.begin(); it != triangles.end(); it++){
//        fprintf(file, "%d ", it->second);
//    }

    fclose(file);

}

void FormanGradientVector::output_mm(list<DAG_TopoNode*>* topodags, list<DAG_GeomNode*>* geomdags,char* filename){

    FILE* file;
    file = fopen(filename, "w");

    map<DAG_TopoNode*,int> topo_to_int;
    //qua metto a posto la parte topologica
    int topon=0;
    for(list<DAG_TopoNode*>::iterator it = topodags->begin(); it != topodags->end(); it++){
        topo_to_int[*it]=topon++;
        cout <<  (*it)->getFathers()->size() << " ";
    }
    cout << endl;


    map<DAG_GeomNode*,int> geom_to_int;
    //qua metto a posto la parte geometrica
    for(int i=0; i<dag_per_vertex->size(); i++)
        geom_to_int[(*dag_per_vertex)[i]]=i;

    //stampe dei DAG nodes
    fprintf(file, "%d %d\n", dag_per_vertex->size(), topodags->size());

    //stampe dei DAG nodes geom
    for(int i=0; i<dag_per_vertex->size(); i++){
        if((*dag_per_vertex)[i] == NULL)
            fprintf(file, "0\n");
        else
        {
            DAG_GeomNode* node = (*dag_per_vertex)[i];
            fprintf(file, "1\n");

            vector<int> verts=node->getVertices();
            for(int j=0; j<verts.size(); j++){
                if((*dag_per_vertex)[verts[j]] == NULL){
                    fprintf(file, "%d ", (-mesh->getUpdateVertexIndex(verts[j]))-1);
                }
                else
                    fprintf(file, "%d ", verts[j]);
            }
            fprintf(file, "\n");

            fprintf(file, "%lf %d %d\n", node->getEdgeLenght(), node->get_to_switch().first,node->get_to_switch().second);
            for(int i=0; i<node->getCoordsv1().size(); i++)
                fprintf(file, "%lf ", node->getCoordsv1()[i]);
            fprintf(file, "\n");
            for(int i=0; i<node->getCoordsv2().size(); i++)
                fprintf(file, "%lf ", node->getCoordsv2()[i]);
            fprintf(file, "\n");

            vector<int> fathers = node->getVV();
            fprintf(file, "%d ",fathers.size());
            for(int i=0; i<fathers.size(); i++){
                if((*dag_per_vertex)[fathers[i]] == NULL)
                    fprintf(file, "%d ", -1);
                else
                    fprintf(file, "%d ", geom_to_int[(*dag_per_vertex)[fathers[i]]]);
            }
            fprintf(file, "\n");

            vector<DAG_GeomNode*> children = node->get_dep();
            fprintf(file, "%d ",children.size());
            for(int i=0; i<children.size(); i++){
                fprintf(file, "%d ", geom_to_int[children[i]]);
            }
            fprintf(file, "\n");
        }
    }

    //stampe dei DAG nodes topo
    for(list<DAG_TopoNode*>::iterator it = topodags->begin(); it != topodags->end(); it++){
        vector<int> verts = (*it)->getVertices();
        fprintf(file, "%d ", verts.size());
        for(int i=0; i<verts.size(); i++){
            if((*dag_per_vertex)[verts[i]] == NULL)
                fprintf(file, "%d ", (-mesh->getUpdateVertexIndex(verts[i]))-1);
            else
                fprintf(file, "%d ", verts[i]);
        }
        fprintf(file, "%lf\n", (*it)->getHeight());

        set<DAG_TopoNode*>* topofathers = (*it)->getFathers();
        fprintf(file, "%d ", topofathers->size());
        for(set<DAG_TopoNode*>::iterator it = topofathers->begin(); it != topofathers->end(); it++){
            fprintf(file, "%d ", topo_to_int[*it]);
        }
        fprintf(file, "\n");

        set<DAG_TopoNode*>* topochildren = (*it)->getChild();
        fprintf(file, "%d ", topochildren->size());
        for(set<DAG_TopoNode*>::iterator it = topochildren->begin(); it != topochildren->end(); it++){
            fprintf(file, "%d ", topo_to_int[*it]);
        }
        fprintf(file, "\n");

        set<DAG_GeomNode*>* geomfathers = (*it)->getFathersGeom();
        fprintf(file, "%d ", geomfathers->size());
        for(set<DAG_GeomNode*>::iterator it = geomfathers->begin(); it != geomfathers->end(); it++){
            fprintf(file, "%d ", geom_to_int[*it]);
        }
        fprintf(file, "\n");
    }

    //stampe della geometria
    fprintf(file, "%d %d\n", mesh->getNumVertex(), mesh->getTopSimplexesNum());
    for(int i=0; i<mesh->getNumVertex(); i++){
        fprintf(file, "%lf %lf %lf\n", mesh->getVertex(i).getX(),mesh->getVertex(i).getY(),mesh->getVertex(i).getZ());
    }

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        fprintf(file, "%d %d %d %d %d %d\n", mesh->getTopSimplex(i).TV(0),mesh->getTopSimplex(i).TV(1),mesh->getTopSimplex(i).TV(2),mesh->getTopSimplex(i).TT(0),mesh->getTopSimplex(i).TT(1),mesh->getTopSimplex(i).TT(2));
    }

    //stampe del gradiente

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        fprintf(file,"%u ", forman_gradient[i]);
    }
    fprintf(file, "\n");

    fclose(file);
}



void FormanGradientVector::write_mesh_VTK(string mesh_name ){
    stringstream stream;
    stream<<mesh_name<<".vtk";
    ofstream output(stream.str().c_str());
    output.unsetf( std::ios::floatfield ); // floatfield not set
    output.precision(15);

     output<<"# vtk DataFile Version 2.0" << endl << endl
         << "ASCII" << endl << "DATASET UNSTRUCTURED_GRID " <<  endl << endl;

    output<< "POINTS " << mesh->getNumVertex() << " float" << endl;

 for(int v=0; v<mesh->getNumVertex(); v++)
    {
        Vertex3D& vert = mesh->getVertex(v);
        output<<vert.getX()<<" "<<vert.getY()<<" "<<vert.getZ()<<endl;
    }

    output<<endl << "CELLS " << mesh->getTopSimplexesNum() << " " << (mesh->getTopSimplexesNum()*4) << endl;

    for(int t=0; t<mesh->getTopSimplexesNum(); t++)
    {
        output<<"3 ";
        for(int i=0; i< 3; i++)
            output<<mesh->getTopSimplex(t).TV(i)<<" ";
        output<<endl;
    }

    output<< endl << "CELL_TYPES " <<  mesh->getTopSimplexesNum()<< endl;
    for (int i = 0; i <  mesh->getTopSimplexesNum(); ++i)
        output<< "6 ";
    output<< endl;

    output<< "POINT_DATA " << mesh->getNumVertex() << endl << endl;
    output<< "FIELD FieldData 1" << endl << endl;
    output<< "fieldvalue 1 " << mesh->getNumVertex() << " float" << endl;

    for (int i=0; i <mesh->getNumVertex(); ++i)
        output<<mesh->getVertex(i).getZ()<< " ";
    output<<endl;

    output.close();





}















