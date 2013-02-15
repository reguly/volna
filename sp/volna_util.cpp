#include "volna_util.h"
#include "volna_common.h"
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include "op_lib_core.h" // To find inverse map in OP_map_list


// PT_SCOTCH header - use sequential functions
#ifdef HAVE_PTSCOTCH
  #include <scotch.h>
#endif
// ParMETIS header - use sequential functions
#ifdef HAVE_PARMETIS
  #include <parmetis.h>
#endif

/*
 * Prints the max. distance of adjacent cells
 */
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
      // If the neighbour is a boundary cell
      if(map->map[i*3+j] == -1 ) {
        dist = 0;
      } else {
        dist = abs(i-map->map[i*3+j]);
        sum_dist += dist;
        n_dist++;
      }
      if(dist > max_dist) {
        // op_printf("%d : %d\n",i,map->map[i*3+j]);
        max_dist = dist;
      }
    }
  }
  avg_dist = sum_dist / (double)n_dist;
  op_printf("\n--------------------- Map statistics -----------------------\n");
  op_printf("Map: %s\n", map->name);
  op_printf("Matrix bandwidth (2*MAX+1)                         = %d\n", 2*max_dist+1);
  op_printf("Avg. dist. between neighboring cells cell-cell map = %g\n", avg_dist);
  op_printf("------------------------------------------------------------\n");
  free(ccell);
}

/*
 * Prints the adjacency list to a file
 */
void printAdjListToFile(const char* filename, op_map map) {
  FILE *fid;
  fid = fopen(filename,"w");
  for(int i=0; i < map->from->size; i++) {
    fprintf(fid,"%d %d %d\n", map->map[i*3],map->map[i*3+1],map->map[i*3+2]);
  }
  fprintf(fid,"\n");
  fclose(fid);
}


/*
 * Print adjacency list to console
 */
void printAdjListToConsole(int *xadj, int n_xadj, int *adjncy, int n_adjncy) {
  op_printf("xadj: ");
  for(int i=0; i<n_xadj; i++) {
    op_printf("%d ",xadj[i]);
  }
  op_printf("\n");
  op_printf("adjcny: ");
  for(int i=0; i<n_adjncy; i++) {
    op_printf("%d ",adjncy[i]);
  }
  op_printf("\n");
}


/*
 * Sort Q according to deg array (increasing)
 */
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

/*
 * Get the complement, adjacent, ordered (by degree) vertices
 */
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

/*
 * Append sets: R = R U (Q.vrtid)
 */
void appSets(std::vector<int> *R, std::vector<VrtDegree> *Q) {
  for(int i=0; i<(*Q).size(); i++) {
    (*R).push_back((*Q)[i].vrtid);
  }
  (*Q).clear();
}


/*
 * Remove non-referenced elements
 */
void op_remove_nonref(int *cnode, int ncell, int *isRefNode, int *nnode, float **x) {
  //
  // Make an exclusive scan (all-prefix-sum or cumulative vector sum)
  //
  op_printf("Removing non-referenced elements...\n");
  //  isRefNode[0] = 0;
  //  int prev = isRefNode[0];
  int *newNum = (int*) malloc((*nnode)*sizeof(int));
  int prev = 0;
  int tmp = 0;
  int num = 0;
  //  for(int i=1; i<nnode; i++) {
  for(int i=0; i<(*nnode); i++) {
    if(isRefNode[i] == 1) {
      num++;
      newNum[i] = num; // if node is referenced give it a unique number (numbering start from 1!)
    }else{
      newNum[i] = 0; // if node is not referenced make it 0
    }
    //    tmp = isRefNode[i-1] + prev;
    //    prev = isRefNode[i];
    //    isRefNode[i] = tmp;
    //    op_printf(" %d ", newNum[i]);
  }

  op_printf("  Updating list of node coordinates, number of nodes and cell-node map...\n");
  //  int newnnode = isRefNode[nnode-1] + 1; // Number of nodes that are referenced by any cell
  int newnnode = num; // Number of nodes that are referenced by any cell
  float* newx = (float*) malloc(MESH_DIM * newnnode * sizeof(float)); // Node coordinates that are referenced by any cell
  int n, nRef;

  for(int i=0; i<N_NODESPERCELL*ncell; i++) {
    n = cnode[i];
    //    nRef  = isRefNode[n]; // isRefNode gives the new number of a referenced node
    nRef = newNum[n] - 1; // The new number of a referenced node has to start from 0
    cnode[i] = nRef;
    newx[2*nRef  ] = (*x)[2*n  ];
    newx[2*nRef+1] = (*x)[2*n+1];
  }
  free(newNum);
  *x = newx;
  *nnode = newnnode;
  op_printf("  done.\n");
  op_printf("done.\n");
}


/*
 * Get adjacency list from two maps, where map1->map2 refers back to
 * the original set.
 */
void op_get_adjlist_indirect(int **xadj, int **adjncy, int *n_xadj, int *n_adjncy, op_map map, op_set set) {
  if(map->from != set && map->to != set) {
    op_printf("\nError: At least one set in the map has to refer to the"
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
  // Create inverse map
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

void print_map_row(op_map map, int row) {
  // Cehck if row is in the range of the map
  if(row < 0 || row > map->from->size) {
    op_printf("\nERROR: %d row number is not within the range 0...%d of "
        "map %s.\n", row, map->from->size, map->name);
    exit(-1);
  // If OK print the given row of the map
  } else {
    op_printf("\n%s maps row %d: ", map->name, row);
    for(int j=0; j<map->dim; j++) {
      op_printf(" %d ",map->map[ row*map->dim + j]);
    }
    op_printf("\n");
  }
}

/*
 * Check if permutation vector references every element
 */
void check_map(op_map map) {
  if(map->from == NULL || map->to == NULL) {
    op_printf("\nERROR: map has NULL set on at least one side.\n");
    exit(-1);
  }
  int n_error = 0;
  // If map is an identity map
  if(map->from == map->to ) {
    int set_size = map->from->size;
    int map_dim = map->dim;
    int elem = 0;
    int *ref = (int*) calloc(set_size, sizeof(int));
    int n_refed = 0; // Number of referenced elements (in one row of the map)
    for(int i=0; i < set_size; i++) {
      n_refed = 0;
      for(int j=0; j < map_dim; j++) {
        elem = map->map[i*map_dim+j];
        // Check if element is in the sets range
        if(elem < 0 || elem > set_size-1) {
          op_printf("\nERROR: element %d in the identity map"
              " %s is not in the range 0...%d of the set. ",
              i, map->name, set_size-1);
          print_map_row(map, i);
          n_error++;
        // If yes increase reference number
        } else {
          ref[ elem ] ++; // increase reference number
          // if an element is not referencing other elements its value
          // is most probably set to 0 or its own index. Only increase
          // reference number if the references seem to be valid.
          if((elem != i && elem != 0)) {
            n_refed++;
          }
        }
      }
      if(n_refed == 0) {
        op_printf("\nERROR: element %d in the identity map"
            " %s is either referencing itself or 0s. ",
            i, map->name);
        print_map_row(map, i);
        n_error++;
      }
    }
    for(int i=0; i < set_size; i++) {
      if(ref[i] == 0) {
        op_printf("\nERROR: element %d in the identity map"
            " %s is not referenced by any other elements. ",
            i, map->name);
        print_map_row(map, i);
        n_error++;
      }
      if(ref[i] > map_dim) {
        op_printf("\nERROR: element %d in the identity map "
            " %s is referenced by more than %d "
            "elements. ",i ,map->name, map->dim);
        print_map_row(map, i);
        n_error++;
      }
    }
    free(ref);
  }
  // If map is a not an identity map
  else {
    // FIRST: check the referenced elements are in range and count the references of the "to" elements
    int map_dim = map->dim;
    int from_size = map->from->size;
    int to_size = map->to->size;
    int *ref = (int*) calloc(to_size, sizeof(int));
    int elem = 0;
    for(int i=0; i < from_size; i++) {
      for(int j=0; j < map_dim; j++) {
        elem = map->map[i*map_dim+j];
        // Check if element is in the sets range
        if(elem < 0 || elem > to_size-1) {
          op_printf("\nERROR: element %d in the (non-identity) map"
              " %s is not in the range 0...%d of the set. ",
              i, map->name, to_size-1);
          print_map_row(map, i);
          n_error++;
        // If yes increase reference number
        } else {
          ref[ elem ]++; // increase reference number
        }
      }
    }
    // SECOND: Check if the every "to" element is referenced at least once, i.e. theres is no outlier
    for(int i=0; i < to_size; i++) {
      if(ref[i] == 0) {
        op_printf("\nERROR: element %d in the (non-identity) map"
            " %s is not referenced by any other elements. ",
            i, map->name);
        print_map_row(map, i);
        n_error++;
      }
    }

//    //
//    // Checking inverse set
//    //
//    free(ref);
//    ref = (int*) calloc(from_size, sizeof(int));
//    op_map imap = NULL;
//    int imap_dim = 0;
//    bool found = false;
//    int i=0;
//    // Look for inverse setin map list
//    while(i<OP_map_index && !found) {
//      if(OP_map_list[i]->from == map->to && OP_map_list[i]->to == map->from) {
//        found = true;
//        op_printf("SEEEEEEEEEEEEEETTTTTTTTTTTT FOUND\n");
//        imap = OP_map_list[i];
//        int imap_dim = imap->dim;
//        for(int i=0; i < to_size; i++) {
//          for(int j=0; j< from_size; j++) {
//            elem = imap->map[i*imap_dim+j];
//            // Check if element is in the sets range
//            if(elem < 0 || elem > from_size-1) {
//              op_printf("\nERROR: [while checking inverse map] element %d in the (non-identity) map"
//                  " %s is not in the range 0...%d of the set. ",
//                  i, imap->name, from_size-1);
//              print_map_row(imap, i);
//              n_error++;
//            // If yes increase reference number
//            } else {
//              ref[ elem ]++; // increase reference number
//            }
//          }
//        }
//      }
//      i++;
//    }
//    free(ref);
  }

  if(n_error != 0) {
    op_printf("\n %d ERRORs occured during the check of map references."
        " Exiting... \n", n_error);
    exit(-1);
  }
}

/*
 * Check if permutation vector references every element
 */
void is_all_referenced_once(int set_size, int* iperm) {
  int *check = (int*) calloc(set_size, sizeof(int));
  int sum =0;
  for(int i=0; i<set_size; i++) check[iperm[i]]++;
  for(int i=0; i<set_size; i++) sum += check[i];
  if(sum != set_size) {
    op_printf("\nVector sum not eq. to ncell\n");
    exit(-1);
  }
  for(int i=0; i<set_size; i++) {

    if(check[i] != 1) {
      op_printf("\ncheck[%d] = %d \n", i, check[i]);
      op_printf("\nVector element != 1\n");
      exit(-1);
    }
  }
  free(check);
}


/*
 * Get permutation list based on the GPS reordering alg. or based on the
 * given map. If (map->from == set) a self set-to-set map is created and
 * this is used by the GPS alg. In case the map->to==set, this map is
 * used to reorder the set by visiting the elements of map->from and
 * assigning new numbers to the set elements according to this visiting
 * sequence.
 */
#ifdef HAVE_PTSCOTCH
void op_get_permutation(int **perm, int **iperm, op_map map, op_set set) {
  if(map->from != set && map->to != set) {
    op_printf("\nError: At least one set in the map has to refer to the"
        "given set. \n");
    exit(-1);
  }
  int set_size = set->size;
  if(map->from == set) {
//    op_printf("Reordering %s set with GPS...", set->name);
    int n_xadj;
    int n_adjncy;
    int *xadj = NULL;
    int *adjncy = NULL;
    //  *perm = (int*) malloc(set_size * sizeof(int));
    //  *iperm = (int*) malloc(set_size * sizeof(int));

    if(map->from == map->to) {
      op_get_adjlist_direct(&xadj, &adjncy, &n_xadj, &n_adjncy, map);
    } else {
      op_get_adjlist_indirect(&xadj, &adjncy, &n_xadj, &n_adjncy, map, set);
    }

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
    mesg = SCOTCH_graphBuild(graphptr, baseval, vertnbr, verttab, vendtab, velotab, vlbltab, edgenbr, edgetab, edlotab);
    if(mesg != 0){
      op_printf("Error during SCOTCH_graphBuild() \n");
      exit(-1);
    }

    SCOTCH_Strat *straptr = SCOTCH_stratAlloc();
    SCOTCH_stratInit(straptr);

    char * strategyString = "g";
    //    char * strategyString = "(g{pass=100})";
    mesg = SCOTCH_stratGraphOrder(straptr, strategyString);
    if(mesg != 0){
      op_printf("Error during setting strategy string. \n");
      exit(-1);
    }

    mesg = SCOTCH_graphOrder(graphptr, straptr, permtab, peritab, cblkptr, rangtab, treetab);
    if(mesg != 0){
      op_printf("Error during SCOTCH_graphOrder() \n");
      exit(-1);
    }

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

    // If the map points to the given set, order the set incrementally
    // according to the map->from sets ordering
//    op_printf("done.\n");
  }
  else if(map->to == set) {
    //
    // Obtain new ordering of nodes based on reordered cells
    //
    // Create permutation vector for nodes
//    op_printf("Reordering %s based on %s map...\n", set->name, map->name);
    int *nodes_perm  = (int*) calloc( set->size, sizeof(int));
    int *nodes_iperm = (int*) calloc( set->size, sizeof(int));
    int counter = 0;
    for(int i=0; i<map->from->size; i++) {
      for(int j=0; j<map->dim; j++) {
        // If node didn't get a new numbering so far
        if(nodes_perm[map->map[i*3+j]] == 0) {
          // label it with a the next unique number that is larger than 0
          counter++;
          nodes_perm[map->map[i*3+j]] = counter;
        }
      }
    }
    for(int i=0; i<counter; i++) nodes_perm[i]--;
    //    op_printf("Number of nodes in set                = %d \n", nnode);
    op_printf("Number of nodes in set                = %d \n", set->size);
    op_printf("Number of nodes in permutation vector = %d \n", counter);

    // Check
    if(counter != set_size) {
      op_printf("The number of nodes in the set and number of "
          "reordered nodes doesn't correspond!\n");
      exit(-1);
    }

    // Create inverse permutation vector
    for(int i=0; i<counter; i++) {
      nodes_iperm[nodes_perm[i]] = i;
    }

    // Fill up the inverse permutation and permutation arrays
    for(int i=0; i<set_size; i++) {
      (*iperm)[i] = nodes_iperm[i];
      (*perm)[i] = nodes_perm[i];
    }
  }

  // Check if every set element is referenced exactly once
  is_all_referenced_once(set_size, *iperm);
//  op_printf("done.\n");
}
#endif

/*
 * Reorder "from" and "to" sets of a map according to the permutation
 * and inverse permutation vectors
 */
void op_reorder_map(op_map map, int *perm, int *iperm, op_set set) {
  if(map->from != set && map->to != set) {
    op_printf("\nError: The map has to have at least one set on which it "
        "is referring.\n");
    exit(-1);
  }
//  op_printf("Reordering %s map elements based on new ordering of %s set.", map->name, set->name);
  int from_size = map->from->size;
  int set_size = set->size;
  int map_dim = map->dim;
  if(map->to == set) {
    // Reorder cell-cell map according to the SCOTCH permutation vector
    // First: set the neighbour cell indices according to the map
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
//  op_printf("done.\n");
  //    printAdjList("reordered_mat.txt",map);
  if(map->from == map->to) {
    printBandwidth(map);
  }
}

/*
 * Reorder set associated data according to the inverse permutation
 * vector
 */
void op_reorder_dat(op_dat dat, int *iperm, op_set set) {
  if(dat->set != set) {
    op_printf("\nError: The dataset is not associated with the given set. \n");
    exit(-1);
  }
//  op_printf("Reordering %s dat based on new ordering of %s set.", dat->name, set->name);
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
//  free(tmp_dat);
//  op_printf("done.\n");
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


