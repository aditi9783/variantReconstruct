//
//  snpgraph.cpp
//  variantReconstruct
//
//  Created by Aditi Gupta on 8/5/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//

#include <iostream>
#include "snpgraph.h"
#include <string>

void snpgraph::initializeGraph( map<string, int> linkedsnps ) {
    
    int index = 1; // index 0 already belongs to the graphcenter
    int nshared = linkedsnps.size(); // number of shared/linked snps
    // initialize adjacency matrix
    vector<int> temp;
    for (int i = 0; i != nshared+1; i++) { // +1 for the graphcenter node
        for (int j = 0; j != nshared+1; j++) {
            temp.push_back(0);
        }
        adjacencyMat.push_back(temp);
    }
    
    // read the shared snps, update nodeindex vector representing snps that are linked
    for (map<string, int>::iterator ss = linkedsnps.begin(); ss != linkedsnps.end(); ++ss ) {
        nodeindex[ss->first] = index;
        adjacencyMat[0][index] = ss->second; // index 0 is for the graphcenter, i.e. the snp whose shared_snp list is being used here to initialize the adjacency matrix
        adjacencyMat[index][0] = ss->second; // make the matrix symmetric
        index++;
    }
    
    numnodes = nodeindex.size();
    // print node index map
/*    cout << "in initilize graph" << endl;
    for (map<string, int>::iterator ni = nodeindex.begin(); ni != nodeindex.end(); ++ni) {
        cout << "\t" << ni->first << "\t" << ni->second << endl;
    }
 */
}

void snpgraph::updateGraph( string sharedskey, map<string, int> linkedsnps ) {
    // here, the linkedsnps are the snps that share a read with one of the snps in nodeindex, i.e., the snps that share a read with the graphcenter
    
    // get the index of the snp that matches sharedskey
    int skeyidx = nodeindex[ sharedskey ];
    
    //update the links that this snp may have with other nodes in this graph
    for (map<string, int>::iterator ss = linkedsnps.begin(); ss != linkedsnps.end(); ++ss ) {
        if (nodeindex.find(ss->first) == nodeindex.end()) { // if ss->first is not present in nodeindex map, do nothing
            continue;
        } else {    
            int ssidx = nodeindex[ss->first]; // get index for the shared snp
            adjacencyMat[skeyidx][ssidx] = ss->second;
            adjacencyMat[ssidx][skeyidx] = ss->second;
        }
    }
    // print node index map
/*    cout << "in update graph" << endl;
    for (map<string, int>::iterator ni = nodeindex.begin(); ni != nodeindex.end(); ++ni) {
        cout << "\t" << ni->first << "\t" << ni->second << endl;
    }
 */

}

void snpgraph::printGraph() {
    
    cout << "graph center: " << graphcenter << endl;
    cout << "snps linked to this graph center:" << endl;
    // print nodes with the index first
    for (map<string, int>::iterator n = nodeindex.begin(); n != nodeindex.end(); ++n) {
        cout << n->first << " " << n->second << endl;
    }
    
    // print adjacency matrix
    for (int i = 0; i != numnodes; i++) {
        for (int j =0; j != numnodes; j++) {
            cout << adjacencyMat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
