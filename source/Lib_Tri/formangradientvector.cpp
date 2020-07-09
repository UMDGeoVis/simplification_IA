#include "formangradientvector.h"
#include <stdio.h>
#include "assert.h"
#include "Timer.h"
#include <climits>

FormanGradientVector::FormanGradientVector(Mesh<Vertex3D,Triangle>* mesh, double epsilon)
{
    this->mesh = mesh;

    min = selle1 = max = 0;
    this->epsilon = epsilon;

    forman_ig = IG();

    forman_gradient = vector<unsigned int>(mesh->getTopSimplexesNum(), NON_VALID_CONF);
    unsigned int tetCaseCountOuterAllowed=0;
    unsigned int tot_case=0;

    field = vector<double>(mesh->getNumVertex());
    for(int i=0; i<mesh->getNumVertex(); i++){
        field[i]=mesh->getVertex(i).getZ();
    }


    CompressedToExpandedLUT compressed_to_expanded_local = CompressedToExpandedLUT(512);
    ExpandedToCompressedLUT expanded_to_compressed_local = ExpandedToCompressedLUT();
    compressed_to_expanded = CompressedToExpandedLUT(512);
    expanded_to_compressed = ExpandedToCompressedLUT();


//      #pragma omp parallel shared(compressed_to_expanded_local, expanded_to_compressed_local)
//      {

        //#pragma omp for
        for(unsigned int i = 0; i < 511; i++)
        {

            Arrows ga = static_cast<Arrows>( static_cast<unsigned int>( i ));
            if( isValidTetCase(ga) )
            {
//                 #pragma omp critical
//                 {
                    compressed_to_expanded_local[tetCaseCountOuterAllowed]       = ga;
                    expanded_to_compressed_local[ga]                             = tetCaseCountOuterAllowed;
                    tetCaseCountOuterAllowed++;
                 //}
             }

            tot_case++;
        }
    //}

    compressed_to_expanded = compressed_to_expanded_local;
    expanded_to_compressed = expanded_to_compressed_local;
    cout << "Case " << tetCaseCountOuterAllowed << " " << expanded_to_compressed.size() << " " << tot_case << endl;
}

bool FormanGradientVector::isValidTetCase(Arrows arrow)
{
    // This disables the outer four bits (corresponding to the FT relations for the opposite T to each F) from being valid

    // test vertices
    if( BitTwiddle::popcount( arrow & CONTAINS_V0 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_V1 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_V2 ) > 1 )		return false;

    // test edges
    if( BitTwiddle::popcount( arrow & CONTAINS_E01 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_E02 ) > 1 )		return false;
    if( BitTwiddle::popcount( arrow & CONTAINS_E12 ) > 1 )		return false;

    // test faces
    if( BitTwiddle::popcount( arrow & CONTAINS_F012 ) > 1 )		return false;

    // valid!
    return true;
}

void FormanGradientVector::change_vtstar_mesh(){

//#pragma omp parallel
//    {
        vector<int> vt;
        TriGradient grad;

//        #pragma omp for
        for(int i=0;i<mesh->getNumVertex(); i++){

            vt = mesh->VT(i);
            for(int j=0; j<vt.size(); j++){

                grad = convert_compressed_to_expand(forman_gradient[vt[j]]);
                if(!grad.is_vertex_unpaired(mesh->getTopSimplex(vt[j]).vertex_index(i))){
                    mesh->getVertex(i).VTstar(vt[j]);
                    break;
                }
            }
        }
//    }
}

TriGradient FormanGradientVector::convert_compressed_to_expand(unsigned int ga){
    if(ga == NON_VALID_CONF) return TriGradient();
    return TriGradient(compressed_to_expanded[ga]);
}

unsigned int FormanGradientVector::convert_expand_to_compressed(Arrows ga){
    return expanded_to_compressed[ga];
}

void FormanGradientVector::setVE(int v, int v2){

    vector<int> et = mesh->ET(Edge(v,v2));
    //assert(et.size() > 0);
    mesh->getVertex(v).VTstar(et[0]);
    for(unsigned int i=0; i<et.size(); i++){
        TriGradient ga = convert_compressed_to_expand(forman_gradient[et[i]]);
        ga.setVE(mesh->getTopSimplex(et[i]).vertex_index(v), mesh->getTopSimplex(et[i]).vertex_index(v2));
        forman_gradient[et[i]] = convert_expand_to_compressed(ga.getArrow());
    }
}

void FormanGradientVector::setEF(int v, int f){

    TriGradient ga = convert_compressed_to_expand(forman_gradient[f]);
    ga.setEF(mesh->getTopSimplex(f).vertex_index(v));
    forman_gradient[f] = convert_expand_to_compressed(ga.getArrow());

}

Edge* FormanGradientVector::getVE(int v){

    int t = mesh->getVertex(v).VTstar();
    int v2 = convert_compressed_to_expand(forman_gradient[t]).get_vertex_pair(mesh->getTopSimplex(t).vertex_index(v));

    if(v2 != -1){
        return new Edge(v, mesh->getTopSimplex(t).TV(v2));
    }

    return NULL;
}

int FormanGradientVector::getEF(Edge* e){

    vector<int> et = mesh->ET(*e);
    for(int i=0; i<et.size(); i++){
        int t = convert_compressed_to_expand(forman_gradient[et[i]]).get_edge_pair(mesh->getTopSimplex(et[i]).vertex_index(e->EV(0)),mesh->getTopSimplex(et[i]).vertex_index(e->EV(1)));
        if(t == 3) return et[i];
    }

    return -1;
}

void FormanGradientVector::freeVE(int v1, int v2){
    vector<int> et = mesh->ET(Edge(v1,v2));
    for(unsigned int i=0; i<et.size(); i++){
        TriGradient ga = convert_compressed_to_expand(forman_gradient[et[i]]);
        ga.clearVE(mesh->getTopSimplex(et[i]).vertex_index(v1), mesh->getTopSimplex(et[i]).vertex_index(v2));
        forman_gradient[et[i]] = convert_expand_to_compressed(ga.getArrow());
    }
}

void FormanGradientVector::freeEF(int v, int f){
    TriGradient ga = convert_compressed_to_expand(forman_gradient[f]);
    ga.clearEF(mesh->getTopSimplex(f).vertex_index(v));
    forman_gradient[f] = convert_expand_to_compressed(ga.getArrow());
}

bool FormanGradientVector::is_vertex_critical(int v){
    return convert_compressed_to_expand(forman_gradient[mesh->getVertex(v).VTstar()]).is_vertex_unpaired(mesh->getTopSimplex(mesh->getVertex(v).VTstar()).vertex_index(v));
}

bool FormanGradientVector::is_edge_critical(int v1, int v2){

    int t1 = mesh->getVertex(v1).VTstar();
    int t2 = mesh->getVertex(v2).VTstar();

    if((mesh->getTopSimplex(t1).contains(v2) && !(convert_compressed_to_expand(forman_gradient[t1]).is_edge_unpaired(mesh->getTopSimplex(t1).vertex_index(v1),mesh->getTopSimplex(t1).vertex_index(v2)))) ||
       (mesh->getTopSimplex(t2).contains(v1) && !(convert_compressed_to_expand(forman_gradient[t2]).is_edge_unpaired(mesh->getTopSimplex(t2).vertex_index(v2),mesh->getTopSimplex(t2).vertex_index(v1))))) {
        return false;
    }

    vector<int> et = mesh->ET(Edge(v1,v2));
    for(unsigned int i=0; i<et.size(); i++){
        if(!(convert_compressed_to_expand(forman_gradient[et[i]]).is_edge_unpaired(mesh->getTopSimplex(et[i]).vertex_index(v1),mesh->getTopSimplex(et[i]).vertex_index(v2))))
            return false;
    }

    return true;
}

bool FormanGradientVector::is_face_critical(int f){
    return convert_compressed_to_expand(forman_gradient[f]).is_face_unpaired();
}

void FormanGradientVector::build(){
    forman_ig=IG();
    for(int i=0; i< mesh->getNumVertex(); i++){

        vector<int> vt = mesh->VT(i);
        map<vector<int>, Simplex_Graph* > lower_star = compute_lower_star(i, &vt);

        vector<int> v;
        v.push_back(i);

        if(lower_star.size() == 0){
            //se il lower star è vuoto => minimo trovato
            min++;
        }
        else{
            //altrimenti eseguo l'algoritmo
            map<vector<int>, int> critical_simplexes_lower = map<vector<int>, int>();
            gradient* gradient_vector_in_lower = compute_gradient(v, &lower_star, &critical_simplexes_lower);

            for(simplices_map::iterator it = gradient_vector_in_lower->begin(vector<int>(1)); it != gradient_vector_in_lower->end(vector<int>(1)); it++){
                setVE(it->first[0], it->second);
            }

            for(simplices_map::iterator it = gradient_vector_in_lower->begin(vector<int>(2)); it != gradient_vector_in_lower->end(vector<int>(2)); it++){
                int v = 3;
                v -= mesh->getTopSimplex(it->second).vertex_index(it->first[0]);
                v -= mesh->getTopSimplex(it->second).vertex_index(it->first[1]);

                setEF(mesh->getTopSimplex(it->second).TV(v), it->second);
            }

            delete gradient_vector_in_lower;
        }


        for(map<vector<int>, Simplex_Graph* >::iterator it = lower_star.begin(); it != lower_star.end(); it++)
            delete it->second;
    }

    cout << min << " " << selle1 << " " << max << endl;
}


map<vector<int>, Simplex_Graph*> FormanGradientVector::compute_lower_star(int v, vector<int>* vt){

    //per ogni tetraedro salvo i sottosimplessi che mi interessa visitare.... potrei fare una struttura a grafo per avere easy le relazioni boundary-coboundary
    map<vector<int>, Simplex_Graph*> mappa_costruzione;

    for(int i=0; i<vt->size(); i++){
        int t = vt->at(i); //per ogni tetraedro devo vedere quali sono i suoi vertici che sono compresi nel lower star (scalarfield < di quello di v)
        vector<int> simplex;
        simplex.push_back(v); //il vertice ci va per forza

        Triangle t1 = mesh->getTopSimplex(t);


        for(int j=0; j<3; j++){
            int v1 = t1.TV(j);
         //   if(mesh->getVertex(v).getZ() > mesh->getVertex(v1).getZ()){
             if(filtration[v]>filtration[v1]){
                simplex.push_back(v1);   // All the vertices that are lower than v.
            }
//            if(mesh->getVertex(v).getZ() == mesh->getVertex(v1).getZ() && v != v1)
//                cout << v << " " << mesh->getVertex(v).getZ() << " " << v1 << " " << mesh->getVertex(v1).getZ() << endl;
            //assert(mesh->getVertex(v).getZ() != mesh->getVertex(v1).getZ() || v == v1);
        }
        if(simplex.size() == 1) continue; //v is minimum

        simplex= sort_simplex(&simplex);

        Simplex_Graph* simplex_grafo = new Simplex_Graph();
        simplex_grafo->simplex = simplex;
        simplex_grafo->coboundary = vector<Simplex_Graph*>();
        if(mappa_costruzione.find(simplex) == mappa_costruzione.end())
            mappa_costruzione[simplex] = simplex_grafo; //metto in relazione il simplesso top con la sua rappresentazione
                                                        //relate the top simplex to its representation
        //mappa costruzione: construction map                                                
        //estraggo tutti i sottosimplessi (facce) costruendo il loro coboundary
        if(simplex.size()>2){
            for(int j=1; j<simplex.size(); j++){
                //qua giro nel simplesso ed elimino i vertici uno ad uno per costruire i sottosimplessi (non elimino mai il primo vertice)
                vector<int> sub_simplex = simplex;
                sub_simplex.erase(sub_simplex.begin()+j); // cancello la posizione i-esima


                    map<vector<int>, Simplex_Graph*>::iterator it = mappa_costruzione.find(sub_simplex);
                    Simplex_Graph* sopra;
                    if(it == mappa_costruzione.end()){
                        sopra = new Simplex_Graph();
                        sopra->simplex = sub_simplex;
                        sopra->coboundary = vector<Simplex_Graph*>();
                        mappa_costruzione[sub_simplex] = sopra;
                    }
                    else sopra = it->second;

                push_in_coboundary(sopra, mappa_costruzione[simplex]);
            }
        }
    }

    return mappa_costruzione;
}

gradient* FormanGradientVector::compute_gradient(vector<int>& v, map<vector<int>, Simplex_Graph*>* lower_star_structure, map<vector<int>, int>* critical_points){

    gradient* gradient_v = new gradient(); //gradiente locale

    vector<int> minimal_edge = vector<int>(0); //estraggo l'edge più piccolo (c'è sicuramente perchè altrimenti il lower star sarebbe vuoto)
    vector<vector<int> > other_edges; //già che devo cercare gli edge me li salvo già tutto così ho pronto questo

    for(map< vector<int>, Simplex_Graph* >::iterator it = lower_star_structure->begin(); it != lower_star_structure->end(); it++){
        if(it->first.size() == 2){
            bool same=false;

            if(minimal_edge.size()==0){
                minimal_edge = it->first;
            }
            else if(is_lower(it->first, minimal_edge, &same)){
//                if(!same && minimal_edge != it->first){
                    other_edges.push_back(minimal_edge);
                    minimal_edge = it->first;
//                }
            }
            else other_edges.push_back(it->first);
        }
    }

    //inizializzo le due "code con priorità" pq_one e pq_zero, utilizzo le liste perchè poi devo elminare elementi da pq_zero (elementi anche interni)
    //e con le priority_queue viene una menata inefficente
    list<vector<int> > pq_one = list<vector<int> >();
    list<vector<int> > pq_zero = list<vector<int> >();

    gradient_v->insert(v, minimal_edge[1]); //aggiungo il primo edge
    for(int i=0; i<other_edges.size(); i++){
        push_ordered(&pq_zero, other_edges.at(i));
    }


    vector<Simplex_Graph*> coboundary_faces = (*lower_star_structure)[minimal_edge]->coboundary;
    for(int i=0; i<coboundary_faces.size(); i++){
        coboundary_faces.at(i)->simplex;
        if(num_unpared_faces(coboundary_faces.at(i)->simplex, gradient_v, critical_points) == 1){
            push_ordered(&pq_one, coboundary_faces.at(i)->simplex);

        }
    }

    while(pq_one.size() != 0 || pq_zero.size() != 0){
        while(pq_one.size() != 0){
            vector<int> popped = pq_one.front();
            pq_one.pop_front();
            if(num_unpared_faces(popped, gradient_v, critical_points) == 0)
                push_ordered(&pq_zero, popped);
            else{

                //istanzio il nuovo pezzo di gradiente e rimuovo da PQ_Zero
                vector<int> pair = unique_pairable_face(popped, gradient_v, critical_points);

                for(int i=0; i<popped.size();i++){
                    int j=0;
                    for(j=0; j<pair.size();j++){
                        if(pair.at(j) == popped.at(i)) break;
                    }
                    if(j == pair.size()) break;
                }

                if(pair.size() == 2){
                    //assert(pair.size() == 2);
                    //assert(popped.size() == 3);
                    gradient_v->insert(pair, vector_to_index(popped));
                }
                else{
                    //assert(pair.size() == 1);
                    //assert(popped.size() == 2);
                    gradient_v->insert(pair, extract_missing_id(pair,popped));
                }

                pq_zero.remove(pair);

                //aggiungo le coboundary faces del simplesso "popped" se hanno una sola faccia con cui essere paired
                vector<Simplex_Graph*> coboundary_faces_popped = (*lower_star_structure)[popped]->coboundary;
                for(int i=0; i<coboundary_faces_popped.size(); i++){

                    if(num_unpared_faces(coboundary_faces_popped.at(i)->simplex, gradient_v, critical_points) == 1){
                        push_ordered(&pq_one, coboundary_faces_popped.at(i)->simplex);
                     }
                }

                //aggiungo le coboundary faces del simplesso "pair" se hanno una sola faccia con cui essere paired
                vector<Simplex_Graph*> coboundary_faces_pair= (*lower_star_structure)[pair]->coboundary;
                for(int i=0; i<coboundary_faces_pair.size(); i++){

                    if(num_unpared_faces(coboundary_faces_pair.at(i)->simplex, gradient_v, critical_points) == 1)
                        push_ordered(&pq_one, coboundary_faces_pair.at(i)->simplex);
                }
            }

        }

        if(pq_zero.size() != 0){
            vector<int> critico = pq_zero.front();
            pq_zero.pop_front();

            //aggiungo come punto critico quello con valore inferiore nella coda pq_zero
            if(critico.size() == 2){ (*critical_points)[critico] = critico.front();
                selle1++;
            }
            if(critico.size() == 3){
                (*critical_points)[critico] = critico.front();
                max++;
            }
            //aggiungo tutti i suoi simplessi nel coboundary in pq_one

            vector<Simplex_Graph*> coboundary_faces_critico= (*lower_star_structure)[critico]->coboundary;
            for(int i=0; i<coboundary_faces_critico.size(); i++){
                if(num_unpared_faces(coboundary_faces_critico.at(i)->simplex, gradient_v, critical_points) == 1)
                    push_ordered(&pq_one, coboundary_faces_critico.at(i)->simplex);
            }

        }
    }


    return gradient_v;
}

int FormanGradientVector::num_unpared_faces(vector<int>& co_face, gradient* gradient_v, map<vector<int>, int>* critical_points){

    // preso il complesso co_face
        //estraggo le sue facce e per ogniuna di esse
            //controllo se sono state accoppiate come testa
            //controllo che non siano state accoppiate come coda cercando, come testa, uno delle sue cofacce

    //NOTA BENE: lavorando con un i-simplesso esso può essere accoppiato solo come testa di un
            //(i-1)-simplesso o
            //coda di un (i+1)-simplesso

    //ESEMPIO:
        //co_face = triangolo
        //estraggo i suoi 3 edges
        //ogni edge può essere presente come chiave di ricerca nel vector field (allora è una coda che punta ad una faccia)
        //oppure può essere presente come testa, in quel caso devo cercare i suoi vertici che saranno la coda

    int num_unpaired=0;
    vector<int> simplex = co_face;

    for(int i=1; i<simplex.size(); i++){
        vector<int> simpl_face = simplex;//prendo la sua faccia
        simpl_face.erase(simpl_face.begin()+i);


        if(gradient_v->find(simpl_face) != gradient_v->end(simpl_face)){
            continue;//se già lei è presente nel field allora niente da contare
        }

        if(critical_points->find(simpl_face) != critical_points->end()){
            continue; //se è già un critical point allora niente
        }


        int j=-1;
        if(simpl_face.size() > 1){
            for(j=1; j<simpl_face.size();j++){
                vector<int> simpl_face_second = simpl_face;//prendo la sua faccia
                simpl_face_second.erase(simpl_face_second.begin()+j);

                simplices_map::iterator it = gradient_v->find(simpl_face_second);
                if(it != gradient_v->end(simpl_face_second) /*&& it->second == simpl_face*/){
                    vector<int> simpl;

                    if(simpl_face_second.size() == 1){
                        simpl = simpl_face_second;
                        simpl.push_back(it->second);
                    }
                    else{
                        for(int pos=0; pos<3; pos++){
                            simpl.push_back(mesh->getTopSimplex(it->second).TV(pos));
                        }
                    }
                    simpl = sort_simplex(&simpl);

                    if(simpl == simpl_face)
                        break;
                }

            }
        }
        if(j==simpl_face.size()){num_unpaired++;}
    }

    return num_unpaired;
}

vector<int> FormanGradientVector::unique_pairable_face(vector<int>& s_face, gradient* gradient_v, map<vector<int>, int>* critical_points){

    vector<int> simplex = s_face;
    for(int i=1; i<simplex.size(); i++){
        vector<int> simpl_face = simplex;//prendo la sua faccia
        simpl_face.erase(simpl_face.begin()+i);

        if(gradient_v->find(simpl_face) != gradient_v->end(simpl_face)) continue;//se già lei è presente nel field allora niente da contare
        if(critical_points->find(simpl_face) != critical_points->end()) continue; //se è già un critical point allora niente

        int j=-1;
        if(simpl_face.size() > 1){
            for(j=1; j<simpl_face.size();j++){
                vector<int> simpl_face_second = simpl_face;//prendo la sua faccia
                simpl_face_second.erase(simpl_face_second.begin()+j);

                simplices_map::iterator it = gradient_v->find(simpl_face_second);
                if(it != gradient_v->end(simpl_face_second) /*&& it->second == simpl_face*/){

                    vector<int> simpl;

                    if(simpl_face_second.size() == 1){
                        simpl = simpl_face_second;
                        simpl.push_back(it->second);
                    }
                    else{
                        for(int pos=0; pos<3; pos++){
                            simpl.push_back(mesh->getTopSimplex(it->second).TV(pos));
                        }
                    }
                    simpl = sort_simplex(&simpl);

                    if(simpl == simpl_face)
                        break;
                }
            }
        }
        if(j==simpl_face.size()) return simpl_face;
    }


}



void FormanGradientVector::push_in_coboundary(Simplex_Graph* sface, Simplex_Graph* co_sface){

    int i=0;
    for(; i<sface->coboundary.size(); i++){
        if(sface->coboundary.at(i)->simplex == co_sface->simplex){
            break;
        }
    }

    if(i==sface->coboundary.size()) sface->coboundary.push_back(co_sface);
}


void FormanGradientVector::writeBoxVTK(vector<double> box){

//    vector<double> _1 = vector<double>(box.begin(), box.begin()+3);
//        vector<double> _2 = vector<double>(box.begin(), box.begin()+3); _2[0] = box[3];
//        vector<double> _3 = vector<double>(box.begin(), box.begin()+3); _3[0] = box[3]; _3[2] = box[5];
//        vector<double> _4 = vector<double>(box.begin(), box.begin()+3); _4[2] = box[5];
//        vector<double> _5 = vector<double>(box.begin(), box.begin()+3); _5[1] = box[4];
//        vector<double> _6 = vector<double>(box.begin(), box.begin()+3); _6[1] = _2[1];
//        vector<double> _7 = vector<double>(box.begin()+3, box.end());
//        vector<double> _8 = vector<double>(box.begin(), box.begin()+3); _8[1] = _4[1];

//        vector<vector<double> > vertici;
//        vertici.push_back(_1);
//        vertici.push_back(_2);
//        vertici.push_back(_3);
//        vertici.push_back(_4);
//        vertici.push_back(_5);
//        vertici.push_back(_6);
//        vertici.push_back(_7);
//        vertici.push_back(_8);


    FILE* file;
    file = fopen("box.vtk", "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS 2 float\n");


    fprintf(file, "%f %f %f\n", box.at(0), box.at(1), box.at(2));
    fprintf(file, "%f %f %f\n", box.at(3), box.at(4), box.at(5));
    fprintf(file, "\n\n");

}


bool FormanGradientVector::valid_gradient(){

    for(int i=0; i<mesh->getTopSimplexesNum(); i++){
        if(!mesh->is_alive(i)) continue;
        Triangle t = mesh->getTopSimplex(i);
        int trovato_uno=0;
        for(int j=0; j<3; j++){
            Edge* e = t.TE(j);
            int t1 = getEF(e);
            if(t1==i)trovato_uno++;
            if(t1 != -1 && t1 == i){

                if(t.TT(j) != -1){
                    Triangle t2 = mesh->getTopSimplex(t.TT(j));
                    //assert(t2.contains(e->EV(1)) && t2.contains(e->EV(1)));
                    //assert(convert_compressed_to_expand(forman_gradient[t.TT(j)]).is_edge_unpaired(t2.vertex_index(e->EV(0)), t2.vertex_index(e->EV(1))));

                }
            }
            else if(t1 == -1){
                if(!(is_edge_critical(e->EV(0), e->EV(1)) ||
                       (getVE(e->EV(0)) != NULL && getVE(e->EV(0))->EV(1) == e->EV(1)) ||
                       (getVE(e->EV(1)) != NULL && getVE(e->EV(1))->EV(1) == e->EV(0)))){

                    cout << "Edge " << e->EV(0) << " " << e->EV(1) << endl;
                    if(getVE(e->EV(0)) != NULL)
                        cout << getVE(e->EV(0))->EV(0) << " " << getVE(e->EV(0))->EV(1) << endl;
                    if(getVE(e->EV(1)) != NULL)
                        cout << getVE(e->EV(1))->EV(0) << " " << getVE(e->EV(1))->EV(1) << endl;
                    //assert(is_edge_critical(e->EV(0), e->EV(1)) || (getVE(e->EV(0)) != NULL && getVE(e->EV(0))->EV(1) == e->EV(1)) || (getVE(e->EV(1)) != NULL && getVE(e->EV(1))->EV(1) == e->EV(0)));
                }
            }
            delete e;
        }
        //assert(trovato_uno<=1);
    }

    for(int i=0; i<mesh->getNumVertex();i++){
        if(!mesh->is_v_alive(i)) continue;
        Edge* e = getVE(i);
        if(e!=NULL){
            vector<int> et = mesh->ET(*e);
            for(int j=0; j<et.size(); j++){
                Triangle t = mesh->getTopSimplex(et[j]);
                //assert(convert_compressed_to_expand(forman_gradient[et[j]]).get_vertex_pair(t.vertex_index(i)) == t.vertex_index(e->EV(1)));
            }
        }
        delete e;
    }

    return true;
}
