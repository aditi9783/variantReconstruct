//
//  snpgraph.h
//  variantReconstruct
//
//  Created by Aditi Gupta on 8/5/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//

#ifndef variantReconstruct_snpgraph_h
#define variantReconstruct_snpgraph_h

#include <vector>
#include <map>
#include <string>
using namespace std;

class snpgraph{
    int numnodes;
    string graphcenter; // the snp which is the growing point of the connectivity graph
    map<string, int> nodeindex; // key: snpkey, value: index in adj matrix
    vector<string> nodeindex_vec; // vec index is same as index in adj matrix for the snpkey at that index
    vector<vector<int> > adjacencyMat;

public:
    void setgraphcenter( string snpkey ) { graphcenter = snpkey; nodeindex[snpkey] = 0; nodeindex_vec.push_back(snpkey); }
    void initializeGraph( map<string, int> ); // read in shared_snps of a given snp
    void updateGraph( string, map<string, int>, vector<string> ); // read in snps shared by a snp in the initial graph
    void printGraph();
    
    string returnGraphCenter() { return graphcenter; }
    map<string, int> returnGraphNodes() { return nodeindex; }
    vector<string> returnGraphNodesVec() { return nodeindex_vec; }
    vector<vector<int> > returnAdjMat() { return adjacencyMat; }
    
};


#endif
