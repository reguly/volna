#include "volna_util.h"
#include <vector>
#include <list>
#include <set>
#include <algorithm>

// PT_SCOTCH header - use sequential functions
#ifdef HAVE_PTSCOTCH
  #include <scotch.h>
#endif
// ParMETIS header - use sequential functions
#ifdef HAVE_PARMETIS
  #include <metis.h>
#endif

//
// Create self neighboring map (eg. edges-to-edges map)
//
void op_create_self_map(op_map setToSet, op_set set, op_map cellsToEdges, op_map edgesToCells) {

}


//
// Prints the max. distance of adjacent cells
//
void printBandwidth(op_map map) {
  int ncell = map->from->size;
  int edge_id;
  int *ccell = (int*) malloc(3*ncell*sizeof(int));
  int max_dist = 0;
  int dist = 0;
  double sum_dist = 0;
  int n_dist = 0;
  double avg_dist = 0;
  for(int i=0; i<ncell; i++) {
    for(int j=0; j<3; j++) {
      // If the neighbor is a boundary cell
      if(map->map[i*3+j] == -1 ) {
        dist = 0;
      } else {
        dist = abs(i-map->map[i*3+j]);
        sum_dist += dist;
        n_dist++;
      }
      if(dist > max_dist) {
//        op_printf("%d : %d\n",i,map->map[i*3+j]);
        max_dist = dist;
      }
    }
  }
  avg_dist = sum_dist / (double)n_dist;
  op_printf("Max distance in cell-cell map = %d\n", max_dist);
  op_printf("Avg distance in cell-cell map = %g\n", avg_dist);
  free(ccell);
}

//
// Prints the adjacency list to a file
//
void printAdjList(const char* filename, op_map map) {
  FILE *fid;
  fid = fopen(filename,"w");
  for(int i=0; i < map->from->size; i++) {
    fprintf(fid,"%d %d %d\n", map->map[i*3],map->map[i*3+1],map->map[i*3+2]);
  }
  fprintf(fid,"\n");
  fclose(fid);
}

//
// Sort Q according to deg array (increasing)
//
void sortDeg(int *Q, int *deg, int degree) {
  int tmp = 0;
  int tmp_ind = 0;

  for(int k=0; k<degree; k++) {
    for(int i=0; i<degree-k; i++) {
      if(deg[i+1] > deg[i]) {
        tmp = deg[i];
        deg[i] = deg[i+1];
        deg[i+1] = tmp;
        tmp_ind = Q[i];
        Q[i] = Q[i+1];
        Q[i+1] = tmp_ind;
      }
    }
  }
}

// Adjacent vertices with their degrees
struct VrtDegree {
  int vrtid; // Vertex id
  int vrtdeg; // Vertex degree
  bool operator<(const VrtDegree &rhs) const {return vrtdeg < rhs.vrtdeg;};
  bool operator==(const VrtDegree &rhs) const {return vrtid == rhs.vrtid;};
};

bool isLess(VrtDegree v0, VrtDegree v1) {return v0.vrtdeg < v1.vrtdeg; };
bool isEq(VrtDegree v0, VrtDegree v1) {return v0.vrtid == v1.vrtid; };

//
// Get the complement, adjacent, ordered (by degree) vertices
//
void getCompOrdAdj(std::list<VrtDegree> *Q, int *xadj, int *adjncy, int vertex, std::vector<int> *R) {
  int ind = 0;
  VrtDegree tmp;
  std::vector<int>::iterator it_vec;
  std::list<VrtDegree>::iterator it_list;
//  int vertex;
//  for(int j=0; j<(*R).size(); j++) {
//    vertex = (*R)[j];
    for(int i=xadj[vertex]; i<xadj[vertex+1]; i++) {
      ind = adjncy[i];

      tmp.vrtid = ind;
      tmp.vrtdeg = xadj[ind+1] - xadj[ind];

      it_vec = find((*R).begin(), (*R).end(), ind);
//      op_printf("it_vec = %d     ind = %d     (*R).end() = %d \n", *it_vec, ind, *(*R).end());
      it_list = find((*Q).begin(), (*Q).end(), tmp);
//      op_printf("it_list = %d     ind = %d    (*Q).end() = %d\n", (*it_list).vrtid, ind, (*(*Q).end()).vrtid);
      if(it_vec == (*R).end() && it_list==(*Q).end()) {
//        tmp.vrtid = ind;
//        tmp.vrtdeg = xadj[ind+1] - xadj[ind];
        (*Q).push_back(tmp);
      }
    }
//  }
//  std::sort((*Q).begin(), (*Q).end(), isLess);
    (*Q).sort();
}

//
// Append sets: R = R U (Q.vrtid)
//
void appSets(std::vector<int> *R, std::vector<VrtDegree> *Q) {
  for(int i=0; i<(*Q).size(); i++) {
    (*R).push_back((*Q)[i].vrtid);
  }
  (*Q).clear();
}


/*
 * Get adjacency list from two maps, where map1->map2 refers back to
 * the original set.
 */
void op_get_adjlist_indirect(int **xadj, int **adjncy, int *n_xadj, int *n_adjncy, op_map map, op_set set) {
  if(map->from != set && map->to != set) {
    op_printf("Error: At least one set in the map has to refer to the"
        "given set. \n");
    exit(-1);
  }

  //
  // Create adjacency list
  //

  // Find which side of the map points to the given set: 1 = from ; 0 = to
  int set_is_from = (map->from == set) ? 1 : 0;
  int set_size = 0;
  int iset_size = 0;

  if(map->from == set) {
    set_size = map->from->size;;
    iset_size = map->to->size;
  } else {
    set_size = map->to->size;;
    iset_size = map->from->size;
  }
  int map_dim = map->dim;


  int cell_id = 0;

  //
  // Print map
  //
  for(int i=0; i<set_size; i++) {
//    op_printf("map->map %d : ", i);
    for(int j=0; j<map_dim; j++) {
      cell_id = map->map[i*map_dim+j];
//        op_printf(" %d ", cell_id);
    }
//    op_printf("\n");
  }
//  op_printf("\n");


  std::vector<std::set<int> > imap;
  imap.resize(iset_size);
  std::vector<std::set<int> > map_self;
  map_self.resize(set_size);


  // Create inverse map, eg. from edges-to-cells the cells-to-edges

  //
  // Get inverse map
  //
  for(int i=0; i<set_size; i++) {
    for(int j=0; j<map_dim; j++) {
      cell_id = map->map[i*map_dim+j];
      imap[ cell_id ].insert(i);
//      op_printf("imap %d : %d \n", cell_id, i);
    }
//    op_printf("\n");
  }
//  op_printf("\n");


  //
  // Print inverse map
  //
//  for(int i=0; i<map->to->size; i++) {
//    op_printf("imap %d : ", i);
//    std::set<int>::iterator it = imap[i].begin();
//    for(; it != imap[i].end(); it++) {
//      op_printf(" %d ", *it);
//    }
//    op_printf("\n");
//  }
//  op_printf("\n");




  //
  // Get map_self
  //
  for(int i=0; i<set_size; i++) {
    for(int j=0; j<map_dim; j++) {
//      if(j==0) op_printf("map_self %d : ", i);

      std::set<int>::iterator it = imap[map->map[i*map_dim+j]].begin();
      for(; it != imap[map->map[i*map_dim+j]].end(); it++) {
        map_self[i].insert( *it );
//        op_printf(" %d ",*it);
      }
    }
//    op_printf("\n");
  }
//  op_printf("\n");


  //
  // Print map_self
  //
//  for(int i=0; i<set_size; i++) {
//    op_printf("map_self %d : ", i);
//    std::set<int>::iterator it = map_self[i].begin();
//    for(; it != map_self[i].end(); it++) {
//      op_printf(" %d ", *it);
//    }
//    op_printf("\n");
//  }
//  op_printf("\n");




//  printBandwidth(map);
//  op_printf("Reordering...\n");
//  //#ifdef HAVE_PARMETIS
  *n_xadj = set_size+1;

  *xadj = (int*) calloc(*n_xadj, sizeof(int));
  // Temporarily allocate array enough large to store data
  std::vector<int> adjncy_vec;
  int *tmp = NULL;


  *n_adjncy = 0;
  int i = 0;
  (*xadj)[i] = 0;
   std::set<int>::iterator it = map_self[i].begin();
   for(; it != map_self[i].end(); it++) {
      (*xadj)[i+1]++;
      adjncy_vec.push_back( *it );
      (*n_adjncy)++;
  }

  for(i=1; i<set_size; i++) {
    (*xadj)[i+1] = (*xadj)[i];
    std::set<int>::iterator it = map_self[i].begin();
    for(; it != map_self[i].end(); it++) {
        (*xadj)[i+1]++;
        adjncy_vec.push_back( *it );
        (*n_adjncy)++;
      }
    }

  // Allocate new array to store exact sized adjacency list
  *adjncy = (int*) malloc((*n_adjncy)*sizeof(int));
  for(i=0; i < *n_adjncy; i++) {
    (*adjncy)[i] = adjncy_vec[i];
  }
  adjncy_vec.clear();

  return;
}



/*
 * Get adjacency list from a map on self refering sets, eg. cell-to-cell
 */
void op_get_adjlist_direct(int **xadj, int **adjncy, int *n_xadj, int *n_adjncy, op_map map) {
  //
  // Create adjacency list
  //
  int set_size = map->from->size;
  int map_dim = map->dim;
  //for(int x=0; x<10; x++) {
  printBandwidth(map);
  op_printf("Reordering...\n");
  //#ifdef HAVE_PARMETIS
//  int vtxdist = set_size;
  *n_xadj = set_size+1;
  int tmp_n_adjncy = set_size*map_dim;
  *xadj = (int*) calloc(*n_xadj, sizeof(int));
  // Temporarily allocate array enough large to store data
  *adjncy = (int*) calloc(tmp_n_adjncy, sizeof(int));
  int *tmp = NULL;


  *n_adjncy = 0;
  int i = 0;
  (*xadj)[i] = 0;
//    op_printf("%d neigh: ",i);
  for(int j=0; j<3; j++) {
    if(map->map[i*3+j] != -1) {
//            op_printf("%d ",map->map[i*3+j]);
      (*xadj)[i+1]++;
      (*adjncy)[*n_adjncy] = map->map[i*3+j];
      (*n_adjncy)++;
    }
  }
//    op_printf("\n");
  for(i=1; i<set_size; i++) {
    (*xadj)[i+1] = (*xadj)[i];
//        op_printf("%d neigh: ",i);
    for(int j=0; j<3; j++) {
      if(map->map[i*3+j] != -1) {
//                op_printf("%d ",map->map[i*3+j]);
        (*xadj)[i+1]++;
        (*adjncy)[*n_adjncy] = map->map[i*3+j];
        (*n_adjncy)++;
      }
    }
//        op_printf("\n");
  }

  tmp = *adjncy;
  // Allocate new array to store exact sized adjacency list
  *adjncy = (int*) malloc((*n_adjncy)*sizeof(int));
  for(i=0; i < *n_adjncy; i++) {
    (*adjncy)[i] = tmp[i];
  }
  free(tmp);

  return;
}

/*
 * Get permutation list based on specified reordering
 */
void op_get_permutation(int **perm, int **iperm, op_map map, op_set set) {
  if(map->from != set && map->to != set) {
    op_printf("Error: At least one set in the map has to refer to the"
        "given set. \n");
    exit(-1);
  }
  int set_size = set->size;
  int n_xadj;
  int n_adjncy;
  int *xadj = NULL;
  int *adjncy = NULL;
  *perm = (int*) malloc(set_size * sizeof(int));
  *iperm = (int*) malloc(set_size * sizeof(int));

  if(map->from == map->to) {
    op_get_adjlist_direct(&xadj, &adjncy, &n_xadj, &n_adjncy, map);
  } else {
    op_get_adjlist_indirect(&xadj, &adjncy, &n_xadj, &n_adjncy, map, set);
  }

  //
  // Print adjacency list
  //
//  op_printf("xadj: ");
//  for(int i=0; i<n_xadj; i++) {
//    op_printf("%d ",xadj[i]);
//  }
//  op_printf("\n");
//  op_printf("adjcny: ");
//  for(int i=0; i<n_adjncy; i++) {
//    op_printf("%d ",adjncy[i]);
//  }
//  op_printf("\n");


    //
    // METIS reordering - for fill-in reduction
    //
    ////  idxtype *vtxdist = (idxtype *)xmalloc(sizeof(idxtype)*(comm_size+1));
    //  int mesg = 0;
    //
    //  idx_t options[METIS_NOPTIONS];
    //  METIS_SetDefaultOptions(options);
    //
    //  options[METIS_OPTION_NUMBERING] = 0;
    //  options[METIS_OPTION_CCORDER] = 1;
    //
    //  mesg = METIS_NodeND( &ncell, xadj, adjncy, NULL, options, perm, iperm);
    //  switch(mesg) {
    //  case METIS_OK:
    //    op_printf("METIS: OK\n");
    //    break;
    //  case METIS_ERROR_INPUT:
    //    op_printf("METIS: Error in input\n");
    //    break;
    //  case METIS_ERROR_MEMORY:
    //    op_printf("METIS: Error in memory access\n");
    //    break;
    //  case METIS_ERROR:
    //    op_printf("METIS: Error\n");
    //    break;
    //  }
    //  int * tmp;
    //  tmp = perm;
    //  perm = iperm;
    //  iperm = tmp;


    //
    // Using SCOTCH for reordering
    //
    SCOTCH_Num baseval = 0; // start numbering from 0
    SCOTCH_Num vertnbr = set_size; // number of vertices in graph = number of cells in mesh
    SCOTCH_Num edgenbr = n_adjncy;
    op_printf("vertnbr = %d \n",vertnbr);
    op_printf("edgenbr = %d \n",edgenbr);

    SCOTCH_Graph *graphptr = SCOTCH_graphAlloc();
    SCOTCH_graphInit(graphptr);

    SCOTCH_Strat *stratptr = NULL;
    SCOTCH_Num *verttab = xadj;

    SCOTCH_Num *vendtab = &verttab[1]; // = NULL; // Used to calculate vertex degree = verttab[i+1] - verttab[i]
    SCOTCH_Num *velotab = NULL; // Vertex load = vertex weight
    SCOTCH_Num *vlbltab = NULL;
    SCOTCH_Num *edgetab = adjncy;

    SCOTCH_Num *edlotab = NULL; // Edge load = edge weight
    SCOTCH_Num *permtab = (SCOTCH_Num*) malloc(set_size*sizeof(SCOTCH_Num));
    SCOTCH_Num *peritab = (SCOTCH_Num*) malloc(set_size*sizeof(SCOTCH_Num));
    SCOTCH_Num *cblkptr = (SCOTCH_Num*) malloc(set_size*sizeof(SCOTCH_Num));
    SCOTCH_Num *rangtab = NULL;//(SCOTCH_Num*) malloc(1 + ncell*sizeof(SCOTCH_Num));
    SCOTCH_Num *treetab = NULL;//(SCOTCH_Num*) malloc(ncell*sizeof(SCOTCH_Num));

    int mesg = 0;
    op_printf("Running SCOTCH_graphBuild() \n");
    mesg = SCOTCH_graphBuild(graphptr, baseval, vertnbr, verttab, vendtab, velotab, vlbltab, edgenbr, edgetab, edlotab);
    if(mesg != 0){
      op_printf("Error during SCOTCH_graphBuild() \n");
      exit(-1);
    }
    op_printf("Done SCOTCH_graphBuild() \n");


    SCOTCH_Strat *straptr = SCOTCH_stratAlloc();
    SCOTCH_stratInit(straptr);

    char * strategyString = "g";
//    char * strategyString = "(g{pass=100})";
    op_printf("Running SCOTCH_stratGraphOrder() \n");
    mesg = SCOTCH_stratGraphOrder(straptr, strategyString);
    if(mesg != 0){
      op_printf("Error during setting strategy string. \n");
      exit(-1);
    }
    op_printf("Done SCOTCH_stratGraphOrder() \n");

    op_printf("Running SCOTCH_graphOrder() \n");
    mesg = SCOTCH_graphOrder(graphptr, straptr, permtab, peritab, cblkptr, rangtab, treetab);
    if(mesg != 0){
      op_printf("Error during SCOTCH_graphOrder() \n");
      exit(-1);
    }
    op_printf("Done SCOTCH_graphOrder() \n");

    //  // Reorder cell-cell map according to the SCOTCH permutation vector


    SCOTCH_graphExit(graphptr);
    SCOTCH_stratExit(straptr);


//    printAdjList("orig_mat.txt",map);






    //  std::vector<int> R;
    //  std::list<VrtDegree> Q;
    //  int Rhead = 0;
    //  int Qhead = 0;
    //  int head = 0;
    ////  int tmp = 0;
    //  int tmp_ind = 0;
    //  int ind = 0;
    //  int degree = 0;
    //  int deg_tmp = 0;
    //  int max = 0;
    //  int max_ind = 0;
    ////  int C = 0;
    //  int found = 0;
    //
    //  R.push_back(0);
    //
    //  getCompOrdAdj(&Q, xadj, adjncy, R[0], &R);
    //
    //  int C;
    //  while(Q.size() != 0) { //&& R.size()!=ncell) {
    //    C = Q.front().vrtid;
    //    Q.pop_front();
    //    R.push_back(C);
    //    getCompOrdAdj(&Q, xadj, adjncy, C, &R);
    ////    op_printf("R size = %d    Q size = %d \n", R.size(), Q.size());
    //  }
    ////  op_printf("R: ");
    //  for(int i=0; i<R.size(); i++) {
    ////    op_printf("(%d : %d)  \n", Q[i].vrtid, Q[i].vrtdeg);
    ////    op_printf("(%d : %d)  \n", i, R[i]);
    //  }
    ////  op_printf("\n");
    //
    //
    //
    //  //
    //  // Check
    //  //
    //  int *check = (int*) calloc(ncell, sizeof(int));
    //  int sum =0;
    //  for(int i=0; i<ncell; i++) check[R[i]] ++;
    //  for(int i=0; i<ncell; i++) sum += check[i];
    //  if(sum != ncell) {
    //    op_printf("vector sum not eq. to ncell\n");
    //    exit(-1);
    //  }
    //  for(int i=0; i<ncell; i++) {
    //    if(check[i] != 1) {
    //      op_printf("vector element not eq, 1\n");
    //      exit(-1);
    //    }
    //  }
    //


    // Fill up the inverse permutation and permutation arrays
    for(int i=0; i<set_size; i++) {
      (*iperm)[i] = peritab[i];
      (*perm)[peritab[i]] = i;
    }
//    op_printf("iperm: ");
//    for(int i=0; i<set_size; i++) {
//      op_printf(" %d ", (*iperm)[i]);
//    }
//    op_printf("\n");





    //
    // Check if permutation map reference every element
    //
    int *check = (int*) calloc(set_size, sizeof(int));
    int sum =0;
    for(int i=0; i<set_size; i++) check[(*iperm)[i]]++;
    for(int i=0; i<set_size; i++) sum += check[i];
    if(sum != set_size) {
      op_printf("vector sum not eq. to ncell\n");
      exit(-1);
    }
    for(int i=0; i<set_size; i++) {
      if(check[i] != 1) {
        op_printf("vector element not eq, 1\n");
        exit(-1);
      }
    }
    free(check);
  }


//void op_reorder_map_from(op_map map, ) {
//
//}

void op_reorder_map(op_map map, int *perm, int *iperm, op_set set) {
  if(map->from != set && map->to != set) {
    op_printf("Error: The map has to have at least one set on which it "
        "is referring.\n");
    exit(-1);
  }
  int from_size = map->from->size;
  int set_size = set->size;
  int map_dim = map->dim;
  if(map->to == set) {
    // Reorder cell-cell map according to the SCOTCH permutation vector
    // Frist: set the neighbor cell indices according to the map
    for(int i=0; i < from_size * map_dim; i++) {
      if(map->map[i] > -1) {
        map->map[i] = perm[map->map[i]];
      }
    }
  }
  if(map->from == set) {
    // Second: switch row according to permutations
    int *tmp_map;
    int *new_map = (int*) malloc(map_dim*set_size*sizeof(int));
    for(int i=0; i < from_size; i++) {
      for(int j=0; j < map_dim; j++) new_map[map_dim*i+j] = map->map[map_dim*iperm[i]+j];
    }
    tmp_map = map->map;
    map->map = new_map;
    free(tmp_map);
  }

//    printAdjList("reordered_mat.txt",map);
  if(map->from == map->to)
    printBandwidth(map);


//    // Reuse tmp_map for reordering cellsToNodes
//    for(int i=0; i < cellsToNodes->from->size; i++) {
//      for(int j=0; j < 3; j++) new_map[3*i+j] = cellsToNodes->map[3*iperm[i]+j];
//    }
//    tmp_map = cellsToNodes->map;
//    cellsToNodes->map = new_map;
//    new_map = tmp_map;
//
//    // Reuse tmp_map for reordering cellsToEdges
//    for(int i=0; i < cellsToEdges->from->size; i++) {
//      for(int j=0; j < 3; j++) new_map[3*i+j] = cellsToEdges->map[3*iperm[i]+j];
//    }
//    tmp_map = cellsToEdges->map;
//    cellsToEdges->map = new_map;
//    new_map = tmp_map;
//
//    // Reorder edgesToCells
//    for(int i=0; i < edgesToCells->from->size*2; i++) {
//      if(edgesToCells->map[i] != -1) {
//        edgesToCells->map[i] = perm[edgesToCells->map[i]];
//      }
//    }
}

void op_reorder_dat(op_dat dat, int *perm, int *iperm, op_set set) {
  if(dat->set != set) {
    op_printf("Error: The dataset is not associated with the given set. \n");
    exit(-1);
  }
  int set_size = set->size;
  int dat_dim = dat->dim;
    //
    // Reorder op_dats according to the reordered cells
    //
  float *tmp_dat = (float*) malloc(dat_dim * set_size *sizeof(float));
  float *data_ptr = NULL;
  // Reuse tmp_map for reordering cellCenters
  data_ptr = (float*) dat->data;
  for(int i=0; i < set_size; i++) {
    for(int j=0; j < dat_dim; j++) {
      tmp_dat[dat_dim*i+j] = data_ptr[dat_dim*iperm[i]+j];
    }
  }
  for(int i=0; i < set_size; i++) {
    for(int j=0; j < dat_dim; j++) {
      data_ptr[dat_dim*i+j] = tmp_dat[dat_dim*i+j];
    }
  }
  free(tmp_dat);
}

//
//
//
//
//
//  //
////  // Check if all the elements are present in the permutation array
////  //
////  int *check = (int*) calloc(ncell, sizeof(int));
////  for(int i=0; i<ncell; i++) {
////    check[perm[i]]++;
////  }
////  int n_one = 0; // number of elements in array that are equal to one
////  int n_neone = 0; // number of elements in array that are not equal to one
////  for(int i=0; i<ncell; i++) {
////    check[i] == 1 ? n_one++ : n_neone++;
////  }
////  if(n_one != ncell) {
////    op_printf("Error: the permutation vector does not contain reference to all the cells.");
////    exit(-1);
////  }
////
////    FILE *fid;
////    fid = fopen("orig_mat.txt","w");
//////    fprintf(fid,"%d %d\n", ncell, n_adjncy);
////    for(int i=0; i < ncell; i++) {
//////      for(int j=0; j<3; j++) {
//////      if(cellsToCells->map[i*3+j] != -1){
//////        fprintf(fid,"%d ", cellsToCells->map[i*3+j]);
//////      }
//////      }
//////      fprintf(fid,"\n");
////      fprintf(fid,"%d %d %d\n", cellsToCells->map[i*3],cellsToCells->map[i*3+1],cellsToCells->map[i*3+2]);
////    }
////    fprintf(fid,"\n");
////    fclose(fid);
////
////
////  // Reorder cell-cell map according to the permutation vector
////  for(int i=0; i < cells->size*3; i++) {
////    if(cellsToCells->map[i] != -1) {
////      cellsToCells->map[i] = iperm[cellsToCells->map[i]];
////    }
////  }
////  free(check);
//////  free(xadj);
//////  free(adjncy);
////  free(perm);
////  free(iperm);
//////}
//
//
//
//
//
//
//
////  op_printf("perm: ");
////  for(int i=ncell-100; i<ncell; i++)
////    op_printf("%d=>%d ",i,  perm[i]);
////  op_printf("\n");
//
////  // Reorder cell-cell map according to the permutation vector
////  for(int i=0; i < cells->size*3; i++) {
////    if(cellsToCells->map[i] != -1) {
////      cellsToCells->map[i] = perm[cellsToCells->map[i]];
////    }
////  }
//
////  fid = fopen("reordered_mat.txt","w");
////  for(int i=0; i < ncell; i++) {
////    fprintf(fid,"%d %d %d\n", cellsToCells->map[i*3],cellsToCells->map[i*3+1],cellsToCells->map[i*3+2]);
////  }
////  fprintf(fid,"\n");
////  fclose(fid);
////  printBandwidth(cellsToCells);
//
//
//
//
//
//
//
//
////  for(int i=0; i < ncell; i++) {
////    ((float*)cellCenters->data)[i*2] =
////  }
////
////
////  cellVolumes
//
////  op_printf("Reordering done.\n ");
//
//
//
//
////  mesg = ParMETIS_V3_NodeND( idx_t *vtxdist, idx t *xadj, idx t *adjncy, idx t *numflag, idx t *options, idx t *order,
////  idx t *sizes, MPI Comm *comm
////  )
//
////  op_printf("xadj: ");
////  for(int i=0; i<ncell; i++) {
////    op_printf("%d ",xadj[i]);
////  }
////  op_printf("\n");
////  op_printf("adjcny: ");
////  for(int i=0; i<n_adjncy; i++) {
////     op_printf("%d ",adjncy[i]);
////  }
////  op_printf("\n");


