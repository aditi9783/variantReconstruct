//
//  main.cpp
//  variantReconstruct
//
//  Created by Aditi Gupta on 7/18/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "snp.h"
#include "snpgraph.h"
#include <stdio.h>
#include <vector>
#include <map>
#include <math.h>
#include <algorithm>
#include <regex.h>

using namespace std;
 
vector<string> shallow_snps; // snps at positions where depth is < 1000.
vector<string> snps_to_remove; // snpkeys of snps that are removed beacuse of low depth and fewer connectivity to other snps
map< string, snp > snpmap;   // key: string with snp pos and base, value: corresponding snp object
int numsnp = 0;
map<int, vector<int> > hapmap; // key: snp pos, val: vector of snp pos that belong to a clique centered at the snp pos which is key
map<string, int> snpkey_index; // numerical index for each snpkey, starting from 0 till numsnp-1;
map<int, vector<string> > snpbins; // key: bin identifier, value: list of snps in that bin
vector<vector<string> > cliques;
int LINKTHRES = 2;
int NUMREADTHRES = 2;
float MINSNPFREQ = 0.002;
int MINDEPTH = 1000; // coverage depth below which a snp is considered 'shallow', and not used as building blocks of variants

ofstream outfile; // output file
ofstream linkmat_out; // contains linkmatrix
ofstream truesnps; // print snps that have sufficient read depth, snp freq, and connectivity to other snps
ofstream cliquestofile; // file to print cliques to

// declare functions
void hapread_open(const char *);
void snp_open(const char *);
void linkmatrix();
void removeWeakSNPs();

void getHaps();
vector<vector<string> > disjointedSNPs ( vector<string> ); // take a list of snps, return disjointed set of these snps
vector<string> BFS ( vector<string> );
void printSNPset ( vector<vector<string> > );

void getHaps2();
vector<string> sortedSNPsByFreq();
vector<vector<int> > findAllCliques( vector<vector<int> > ); // takes in adjacency matrix, returns vector of cliques
vector<vector<string> > printCliques( vector<vector<int> >, vector<string> );
void printVec (vector<int>);
vector<vector<int> > uniqueCliques ( vector<vector<int> > );
vector<vector<string> > uniqueStringCliques ( vector<vector<string> > );


int usage() {
    cout << "USAGE:\n"
    << "-r\tHapread file: contains snp pos, snp base, space-separated list of read ids that contain the snp\n"
    << "-s\tSNP positions as obtained by first error pass: contains snp pos, snp base, snp coverage, mpileup coverage at that position\n" << endl;
    return 0;
}

int main(int argc, const char * argv[]) {
    outfile.open("out1.txt");
    linkmat_out.open("linkmatrix.txt");
    truesnps.open("true_snps.txt");
    cliquestofile.open("cliques.txt");
    
    map<int, vector<float> > bindefined;

    
    // check if enough parameters have been passed
    if (argc < 3)
        usage();
    else { 
        for(int i = 1; i < argc; i++) {
            //outfile<<argv[i]<<" "<<argv[i+1]<<endl;
            if ( strcmp(argv[i], "-s") == 0 ) {
                snp_open(argv[i+1]);
            } else if ( strcmp(argv[i], "-r") == 0 ) {
                hapread_open(argv[i+1]);
            }
        } 
        // print the snps that are shared by a given snp
   /*     for (map<string, snp >::iterator s = snpmap.begin(); s != snpmap.end(); ++s) {
            cout << s->first << endl;
            s->second.printSharedSNP();
        }
    */
        removeWeakSNPs(); // SNPs with low depth, low freq, and low connectivity to other SNPs
        // print remaining snps which may be considered as true snps
        numsnp = snpmap.size();
        cout << "Number of snps: " << numsnp << endl;
        for (map<string, snp>::iterator s=snpmap.begin(); s!= snpmap.end(); ++s) {        
            truesnps << s->second.position() << " " << s->second.base() << endl;
        }
        
        // generate and print link matrix to a file
        linkmatrix();
        
        // getHaps();
        getHaps2();
        
        cout << cliques.size() << " cliques obtained!\nRemoving redundant cliques.\n";
        // sometimes cliques completely contained in other cliques are found as well, remove those redundant cliques
        vector<vector<string> > uniqstrclq = uniqueStringCliques(cliques);
        
        cout << "Printing unique cliques to file\n";
        // print these cliques to file
        for (vector<vector<string> >::size_type u=0; u!= uniqstrclq.size(); u++) {
            for (vector<string>::size_type v=0; v!= uniqstrclq[u].size(); v++) {
                cliquestofile << uniqstrclq[u][v] << " ";
            }
            cliquestofile << endl;
        }
    }    
    return 0;
}



void snp_open(const char * snpfile) {
    FILE * file;
    int depth;
    int nreads;
    float frac;
    char base;
    int pos;
    file = fopen (snpfile , "r");
    if (file == NULL) perror ("Error opening snpfile.");
    else {
        while ( ! feof (file) ) {
            snp s;
            stringstream snpkey;  // concatenated string of snp pos and base
            depth = 0;
            nreads = 0.0;
            fscanf(file, "%d\t%c\t%d\t%d\n", &pos, &base, &nreads, &depth);
            snpkey << pos << base; 
            frac = (float)nreads / (float)depth;
            s.define_snp(pos, base, frac, snpkey.str());
            outfile << snpkey.str() << ":" << s.position() << ", " << s.base() << ", " << s.cov() << " " << depth << endl;
                        
            if (depth < MINDEPTH) { // only keep positions where depth is at least MINDEPTH to make reliable estimates of frac
                shallow_snps.push_back( snpkey.str() );
            }
            
            snpmap[snpkey.str()] = s;
        }
        fclose (file);
        //   outfile << endl;
    }
    numsnp = snpmap.size();
    outfile << "number of snps " << numsnp << endl;
}

void hapread_open( const char * hapreadfile ) {
    FILE * file;
    int rid, pos;
    int prevrid = 0; // first read id will be 0
    char base;
    int nsnps = 0; // number of snps that appear together in a read. If a read has only one snp, ignore that read
    vector<string> snpread;    // list of string representation of snp, i.e. snpkey
    snpread.clear(); // make sure it is empty
    
    file = fopen( hapreadfile, "r" );
    if (file == NULL) perror ("Error opening hapread file");
    else {
        while (! feof (file) ){
            stringstream snpkey;
            
            fscanf(file, "%d\t%d\t%c\n", &rid, &pos, &base);
            snpkey << pos << base;
            //           outfile << rid << " " << prevrid << " " << snpkey.str() << endl;
            if (rid == prevrid) { // if read id is same, add snp to the growing list of snpkeys that share this read
                snpread.push_back(snpkey.str());
                nsnps++;
            } else {
                if (nsnps > 1) { // update sharedSNP field only if there are > 1 snps in a read
                    //   outfile << "printing snpread vector" << endl;
                    // if read id is different from previous one, take all snps that shared the prev read and update their shared_snp vector
                    for (vector<string>::size_type i = 0; i != snpread.size(); i++) {
                        //       outfile << i << ": " << snpread[i] << endl;
                        // update number of unique reads that contain snp
                        snpmap[ snpread[i] ].updateNumDiffReads();
                        for (vector<string>::size_type j = 0; j != snpread.size(); j++) {
                            //    outfile << i << " " << snpread[i] << " " << j << " " << snpread[j] << endl;
                            if (i != j ) {
                                snpmap[ snpread[i] ].updateSharedSNP( snpread[j] );
                            }
                        }
                    }
                }
                snpread.clear();    // clear out the snps in prev read
                snpread.push_back(snpkey.str());  // 1st snp in new read
                nsnps = 1;
            }
            prevrid = rid;
        }
        // update sharedSNP for the last read
        if (nsnps > 1) {
            for (vector<string>::size_type i = 0; i != snpread.size(); i++) {
                for (vector<string>::size_type j = 0; j != snpread.size(); j++) {
                    if (i != j) {
                        //    outfile << i << " " << snpread[i] << " " << j << " " << snpread[j] << endl;
                        snpmap[ snpread[i] ].updateSharedSNP( snpread[j] );
                    }
                }
            }
        }    
        fclose(file);
    }
}

void linkmatrix() {
    int linkmat[numsnp][numsnp];
    
    // initialize link matrix
    for(int i=0; i <numsnp; i++) {
        for(int j=0; j <numsnp; j++) {
            linkmat[i][j] = 0;
        }
    }
    // print snpkeys
    for (map<string, snp>::iterator s=snpmap.begin(); s!= snpmap.end(); ++s) {        
        linkmat_out << s->first << ",";
    }
    linkmat_out << endl;
    
    outfile << " number of snps: inside linkmatrix " << snpmap.size() << endl;
    for (map<string, snp>::iterator s=snpmap.begin(); s!= snpmap.end(); ++s) {
        int index1 = snpkey_index[s->first];
        map<string,int> sharedsnps = s->second.returnSharedSNP(); // get snps that share read with the current snp
        
        for (map<string, int>::iterator a=sharedsnps.begin(); a != sharedsnps.end(); ++a) {
            int index2 = snpkey_index[a->first];
            linkmat[index1][index2] = a->second;
        }
        // print linkmatrix
        // linkmat_out << s->first;
        for (int i=0; i < numsnp; i++) {
            linkmat_out << linkmat[index1][i] << ",";
        }
        linkmat_out << endl;        
    }
}

void removeWeakSNPs() {
    
    vector<string> snps_to_remove(shallow_snps); // copy shallow snps to snps_to_remove
    
    for (map<string, snp >::iterator s = snpmap.begin(); s != snpmap.end(); ++s) {
        // if freq is less than thres OR
        // if number of reads that have this snp is less than thres OR
        // if number of times the snp shares read with other snps is less than thres, 
        // add snp to dump list
        if (s->second.cov() < MINSNPFREQ || s->second.returnNumReads() < NUMREADTHRES || s->second.returnNumSharedSNPs() < LINKTHRES) {
            snps_to_remove.push_back(s->first);
        }
    }
    // remove these snps from snpmap
    for (vector<string>::size_type i = 0; i < snps_to_remove.size(); i++) {
        snpmap.erase( snps_to_remove[i] );
    }
    
    // for remaining snps, check the shared snp list and remove snps that are in snps_to_remove
    for (map<string, snp >::iterator s = snpmap.begin(); s != snpmap.end(); ++s) {
        map<string, int > snps_shared = s->second.returnSharedSNP();
        vector<string> shared_snps_to_remove;
        shared_snps_to_remove.clear();
        
        for (map<string, int >::iterator ss = snps_shared.begin(); ss != snps_shared.end(); ++ss) {
            for (vector<string>::size_type i = 0; i < snps_to_remove.size(); i++) {
                if (ss->first == snps_to_remove[i]) {
                    shared_snps_to_remove.push_back(ss->first);
                }
            }
        }
        // remove the shared snps
        for (vector<string>::size_type i = 0; i < shared_snps_to_remove.size(); i++) {
            s->second.deleteSharedSNP(shared_snps_to_remove[i]);
        }
    }
    
    // create snpkey index for building link matrix
    snpkey_index.clear();
    int counter = 0;
    for (map<string, snp>::iterator s=snpmap.begin(); s!= snpmap.end(); ++s) {
        snpkey_index[ s->first ] = counter;
        counter++;
    }
    
    outfile << "\noriginal number of snps:" << numsnp << endl;
    
    numsnp = snpmap.size();
    outfile << "\nafter dumping snps:" << numsnp << endl;
}

void getHaps2 () {
    vector<string> sortedsnps;
    vector<string> snps_seen;
    
    sortedsnps = sortedSNPsByFreq();
    // read snps in the ascending order of their freq
    // build adjacency matrix for all snps that share a read with the snp in consideration
    // find cliques in the graph represented by the adjacency matrix -> these are probable variants
    // also save connected components that do not form clique
    
    for (vector<string>::size_type i = 0; i != sortedsnps.size(); i++) {
        snpgraph snplinks;
        
        snp s = snpmap[ sortedsnps[i] ];
    //    s.printSharedSNP();
        map<string, int> shared_snps = s.returnSharedSNP();
        cout << "SNP " << i << endl;
        
        // erase snps that have been looked at before, to avoid rediscovering same cliques
        for (vector<string>::size_type i = 0; i < snps_seen.size(); i++) {
            shared_snps.erase( snps_seen[i] );
        }
      
            
        snplinks.setgraphcenter(sortedsnps[i]);
        snplinks.initializeGraph( shared_snps );
        // for snps that share a read with current snp, see if those snps are linked as well
        for (map<string, int>::iterator ss = shared_snps.begin(); ss != shared_snps.end(); ++ss) {
            snp nbgsnp = snpmap[ ss->first ];
            snplinks.updateGraph(ss->first, nbgsnp.returnSharedSNP(), snps_seen);
        }
        
    //    snplinks.printGraph();
        snps_seen.push_back( sortedsnps[i] );
        vector<string> nodeidxvec = snplinks.returnGraphNodesVec();
    /*    cout << "Printing node index vector:\n";
        for (int m=0; m != nodeidxvec.size(); m++) {
            cout << m << ": " << nodeidxvec[m] << endl;
        }
     */
        
        // find all cliques in the adjacency matrix of the current snp
        vector<vector<int> > snpcliques;
        
        if (snplinks.returnAdjMat().size() > 1) { // adj matrix will have size of at least 1, i.e. the snp it is centered at. But size = 1 means no other links, so do not consider such a matrix
            snpcliques = findAllCliques( snplinks.returnAdjMat() );
            vector<vector<int> > uniq_snpcliques = uniqueCliques ( snpcliques );
           // printCliques(snpcliques);
        //    cout << "Printing unique snpcliques:\n";
            vector<vector<string> > strclq = printCliques(uniq_snpcliques, nodeidxvec); // get string rep of cliques
            cliques.insert(cliques.end(), strclq.begin(), strclq.end());
            cout << "Number of cliques: " << cliques.size() << endl;
        }
    }
}

vector<vector<int> > findAllCliques( vector<vector<int> > M ) { 
    // reads in adjacency matrix and returns all cliques.
    
    // three pointers: 
    // 1. row 0 (first row) is the center of graph, meaning it is linked to all other nodes, so no need to re-check its links
    // 2. as each row is read (from second row onwards), the node represented by the row index is going to be in all cliques that will be discovered from that row, so no need to check its links as well, for that row alone
    // 3. only need to consider right inverted triangle of the matrix as it is symmetric, i.e. elements to the right of diagonal
    //      => this also implies no need to look at last row as there are no elements to the right of diagonal
    vector<vector<int> > cliques; 
   // int m_size = M.size();
   // cout << "M size: " << m_size << endl;
    
    for (vector<vector<int> >::size_type i=1; i != M.size()-1; i++) { // i=1 because no need to look at first row. see pointer 1.
        vector<int> base_clique_i; // new clique
        base_clique_i.push_back(0);
        base_clique_i.push_back(i); // basic clique i has at least node 0 and i.
        vector<int> nbg_i; // nodes that are linked to node i
        for (vector<int>::size_type j = i+1; j != M[i].size(); j++) { // only elts to right of diagonal, see pointer 3.
            if (M[i][j] > 0) { // nodes i and j are linked: j is neighbor of i
                nbg_i.push_back( j );
            }
        }
   /*     cout << "row " << i << ". nbg_i: ";
        printVec(nbg_i);
        cout << endl;
    */
        
        int num_nbg = nbg_i.size();
        if (num_nbg == 0) { // i has no neighbors, node 0 and i thus form a clique
            cliques.push_back(base_clique_i);
       /*     cout << "num nbg = 0, pushing in clique ";
            printVec(base_clique_i);
            cout << endl;
        */
        } else if (num_nbg == 1) {
            vector<int> c = base_clique_i; // new clique
            c.push_back(nbg_i[0]); // add only neighbor of i to clique
            cliques.push_back(c);
         /*   cout << "num nbg = 1, pushing in clique ";
            printVec(c);
            cout << endl;
          */

        } else { // i has multiple neighbors, and thus clique candidates
            vector<int> separator; // among all neigbors of i, list of nodes that are unlinked, and thus would separate neighbors in different cliques
            // check if neighbors of i are linked as well: look at all neighbors of i in pairs
            for (vector<int>::size_type k=0; k != nbg_i.size(); k ++) {
                
                for (vector<int>::size_type l=k+1; l != nbg_i.size(); l++) {
                  //  cout << "checking edge between nbgs " << nbg_i[k] << " and " << nbg_i[l] << ". val: " << M[ nbg_i[k] ][ nbg_i[l] ];
                    
                    if (M[ nbg_i[k] ][ nbg_i[l] ] == 0) { // no link, these nodes are clique separators
                        separator.push_back( nbg_i[k] );
                        separator.push_back( nbg_i[l] );
                   //     cout << ": does not exist!" << endl;
                    } else {
                   //     cout << ": edge exists!" << endl;
                    }
                }
            }
            
            if ( separator.empty() ) { // no seprators, add all neighbors of i to clique.
                vector<int> c = base_clique_i; // new clique
                c.insert(c.end(), nbg_i.begin(), nbg_i.end()); // add all neighbors of i to clique
                cliques.push_back(c);
            /*    cout << "num nbg = " << num_nbg << ", separator empty, pushing in clique ";
                printVec(c);
                cout << endl;
             */
            } else { // there are separators
                // retain the separator nodes only once. remove duplicates.
                sort( separator.begin(), separator.end() );
                separator.erase( unique( separator.begin(), separator.end() ), separator.end() );
                
                vector<int> non_separators; // neighbors of i that are not separators, and thus form complete cl.
                for (int x=0; x != nbg_i.size(); x++) {
                    int flag = 0;
                    for (int y=0; y != separator.size(); y++) {
                        if (nbg_i[x] == separator[y]) {
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == 0) {
                        non_separators.push_back( nbg_i[x] );
                    }
                }
                // each clique will have one separator, all non_separators, and base_clique
                // append base clique to non separators
                non_separators.insert(non_separators.end(), base_clique_i.begin(), base_clique_i.end()); // add all neighbors of i to clique
                
            /*    cout << "num nbg = " << num_nbg << ", separator NOT empty, non separators ";
                printVec(non_separators);
                cout << endl;
             */
                

                for (int y=0; y != separator.size(); y++) {
                    vector<int> c = non_separators; // new clique
                    c.push_back(separator[y]);
                    cliques.push_back(c);
               /*     cout << "pushing in clique :";
                    printVec(c);
                    cout << endl;
                */
                }
            }
        }
    }
    return cliques;
}

void printVec ( vector<int> v ) {
    for (int i = 0; i != v.size(); i++) {
        cout << v[i] << " ";
    }
    cout << endl;
}

vector<vector<string> > printCliques( vector<vector<int> > c, vector<string> nodedesc ){
    vector<vector<string> > strcliqs;
    
    for (vector<vector<int> >::size_type i = 0; i != c.size(); i++) {
        vector<string> sc;
   //     cout << "clique " << i << " : ";
   //     for (vector<int>::size_type j = 0; j != c[i].size(); j++) {
   //         cout << " " << c[i][j];
   //     }
   //     cout << "\t";
        for (vector<int>::size_type j = 0; j != c[i].size(); j++) {
   //         cout << " " << nodedesc[ c[i][j] ];
            sc.push_back( nodedesc[ c[i][j] ] );
        }
   //     cout << endl;
        strcliqs.push_back( sc );
    }
    return strcliqs;
}

vector<vector<int> > uniqueCliques ( vector<vector<int> > c ) {
    vector<vector<int> > uniq;
    map<int, vector<int> > cliq_size; // key: cliq size, value: vector of cliq ids of that size
    vector<int> sizes;
    
    for (vector<vector<int> >::size_type i = 0; i != c.size(); i++) { // for each clique, save clique id, 'i', in correct size vector
        int cs = c[i].size(); // clique i has size cs
 //       cout << "clique " << i << " size " << cs << endl;
        cliq_size[cs].push_back(i); 
        sizes.push_back(cs);
    }    
    sort(sizes.begin(), sizes.end());
    int numsizes = sizes.size();
    // for all cliques of a given size, see if it occurs in cliques that are bigger, and save it only if it doesn't and thus is unique
    for (int i=0; i != numsizes; i++) { // sizes are sorted, which are also keys of cliq_size map
        
        vector<int> cs_ids = cliq_size[ sizes[i] ]; // ids of cliques of this size
        int numc = cs_ids.size();
        for (int a = 0; a != numc; a++) {
            int cid = cs_ids[a]; // clique id for a clique of this size
            vector<int> cliq = c[cid]; // get clique using this id
            int flag = 0; // if a bigger cliq is found that contains this current cliq, flag will be set to 1
                
            // test if this clique is present in bigger cliques
            for (int j=i+1; j != numsizes; j++) { // sizes are sorted, which are also keys of cliq_size map. only looking at sizes that are bigger
                if (flag == 1) {
                    break; // cliq has been matched, no need to look into bigger cliques
                } else {
                    vector<int> bcs_ids = cliq_size[ sizes[j] ]; // ids of cliques of this size
                    int numbc = bcs_ids.size();            
                    for (int b = 0; b != numbc; b++) {
                        int bcid = bcs_ids[b]; // clique id for the bigger cliq
                        vector<int> bcliq = c[bcid]; // get clique using this id
                            
                        // now compare cliq and bigger cliq
                        int match_size = 0;
                        for (int x=0; x != cliq.size(); x++) {
                            for (int y=0; y != bcliq.size(); y++) {
                                if (cliq[x] == bcliq[y]) {
                                    match_size++;
                                    break;
                                }
                            }
                        }
                        if (match_size == cliq.size()) {
                            flag = 1; // cliq is not unique
                            break;
                        }
                    }
                }
            }
            if (flag == 0) {
                uniq.push_back( cliq );
            } else {
                continue; // if clique is not unique, then the bigger clique that contains this cliq will be saved.
            }
        }
    }
    return uniq;
}

vector<vector<string> > uniqueStringCliques ( vector<vector<string> > c ) {
    vector<vector<string> > uniq;
    map<int, vector<int> > cliq_size; // key: cliq size, value: vector of cliq ids of that size
    vector<int> sizes;
    
    for (vector<vector<int> >::size_type i = 0; i != c.size(); i++) { // for each clique, save clique id, 'i', in correct size vector
        int cs = c[i].size(); // clique i has size cs
        //       cout << "clique " << i << " size " << cs << endl;
        cliq_size[cs].push_back(i); 
        sizes.push_back(cs);
    }    
    sort(sizes.begin(), sizes.end());
    int numsizes = sizes.size();
    // for all cliques of a given size, see if it occurs in cliques that are bigger, and save it only if it doesn't and thus is unique
    for (int i=0; i != numsizes; i++) { // sizes are sorted, which are also keys of cliq_size map
        cout << "Clique size: " << sizes[i] << endl;
        vector<int> cs_ids = cliq_size[ sizes[i] ]; // ids of cliques of this size
        int numc = cs_ids.size();
        for (int a = 0; a != numc; a++) {
            int cid = cs_ids[a]; // clique id for a clique of this size
            vector<string> cliq = c[cid]; // get clique using this id
            int flag = 0; // if a bigger cliq is found that contains this current cliq, flag will be set to 1
            
            // test if this clique is present in bigger cliques
            for (int j=i+1; j != numsizes; j++) { // sizes are sorted, which are also keys of cliq_size map. only looking at sizes that are bigger
                if (flag == 1) {
                    break; // cliq has been matched, no need to look into bigger cliques
                } else {
                    vector<int> bcs_ids = cliq_size[ sizes[j] ]; // ids of cliques of this size
                    int numbc = bcs_ids.size();            
                    for (int b = 0; b != numbc; b++) {
                        int bcid = bcs_ids[b]; // clique id for the bigger cliq
                        vector<string> bcliq = c[bcid]; // get clique using this id
                        
                        // now compare cliq and bigger cliq
                        int match_size = 0;
                        for (int x=0; x != cliq.size(); x++) {
                            for (int y=0; y != bcliq.size(); y++) {
                                if (cliq[x] == bcliq[y]) {
                                    match_size++;
                                    break;
                                }
                            }
                        }
                        if (match_size == cliq.size()) {
                            flag = 1; // cliq is not unique
                            break;
                        }
                    }
                }
            }
            if (flag == 0) {
                uniq.push_back( cliq );
            } else {
                continue; // if clique is not unique, then the bigger clique that contains this cliq will be saved.
            }
        }
    }
    return uniq;
}

vector<string> sortedSNPsByFreq() {
    vector<string> sortedsnps;
    multimap<float, string> snpfreq; // create a multimap as multiple snps can have same freq
    
    for (map<string, snp>::iterator m = snpmap.begin(); m != snpmap.end(); ++m) {
        snpfreq.insert( std::pair<float, string>(m->second.cov(), m->first) );
    }    
    // multimap sorts elements by key, so read in that order and save snpkey
    for (multimap<float, string>::iterator s=snpfreq.begin(); s!=snpfreq.end(); ++s) {
        sortedsnps.push_back(s->second);
        cout << s->first << " " << s->second << endl;
    }
    // check the snps are stored in correct order
 /*   for (vector<string>::size_type i = 0; i != sortedsnps.size(); i++) {
        cout << sortedsnps[i] << endl;;
    }
  */
    return sortedsnps;
}

void getHaps() {
    vector<vector<string> >haps;
    vector< vector<string> > snpsets;
    
    // for each snpbin, identify disjointed snpsets (subgraphs) using DFS
    for (map<int, vector<string> >::iterator b = snpbins.begin(); b != snpbins.end(); ++b) {
        snpsets.clear();
        snpsets = disjointedSNPs ( b->second );
        // printSNPset(snpsets);
        
        // add snpsets to haps
        haps.insert( haps.end(), snpsets.begin(), snpsets.end() );
    }
    outfile << "\nprinting final haps " << endl;
    printSNPset(haps);
}

vector< vector<string> > disjointedSNPs( vector<string> snpl ) {
    vector<vector<string> > subsets;
    vector<string> oneset;
    int starting_nsnps = snpl.size(); // number of snps started out with
    int snpsaccountedfor = 0;
    
    vector<string> Q;
    
    while (snpsaccountedfor < starting_nsnps) {
        Q.clear();
        Q.push_back( snpl[0] ); // add first snp to Q
        
        oneset = BFS( snpl );
        // add this connected snp set to the subsets
        subsets.push_back( oneset );
        snpsaccountedfor += oneset.size();
        
        // delete snps in the oneset from snpl
        for (int i = 0; i < oneset.size(); i++ ) {
            for (int j = 0; j < snpl.size(); j ++ ) {
                if ( oneset[i].compare(snpl[j]) == 0 ) {
                    // remove jth elt from snpl
                    snpl.erase(snpl.begin()+j);
                }
            }
        }
    }    
    return subsets;
}

vector<string> BFS ( vector<string> snpsinbin ) {
    vector<string> connectedgraph;
    vector<string> queue;
    queue.push_back( snpsinbin[0] );
    connectedgraph.push_back( snpsinbin[0] );
    
    while (!queue.empty()) {
        snp s = snpmap[ queue[0] ]; // get first element
        queue.pop_back(); // erase last elt, making it a queue (not a stack, thats why BFS not DFS).
        
        map<string, int> sharedsnps = s.returnSharedSNP(); // get all snps that share a read with s
        for ( map<string, int>::iterator n= sharedsnps.begin(); n != sharedsnps.end(); ++n ) { 
            for (vector<string>::size_type i = 0; i != snpsinbin.size(); i ++ ) { // get bin neighbors
                // if bin neighor also shares a read with the snp, then they form a connected graph
                // require snp s share at least LINKTHRES reads with bin neighbor: a stricter definition of neighbor
                if (n->first == snpsinbin[i] and n->second >= LINKTHRES) {
                    // if n does not exist in connectedgraph, add to it
                    int flag = 0;
                    for (int i = 0; i < connectedgraph.size(); i++ ) {
                        if ( n->first.compare(connectedgraph[i]) == 0 )  {
                            flag = 1; // neighbor snp is present in connected list
                            break;
                        } else {
                            continue;
                        }
                    }
                    if ( flag == 0 ) { // n does not exist in connectedgraph
                        queue.push_back( n->first );
                        connectedgraph.push_back( n->first );
                    }
                }
            }
        }
    }
    return connectedgraph;
}
 
void printSNPset ( vector<vector<string> > ss ) {
    for (vector<vector<string> >::size_type i = 0; i != ss.size(); i++ ) {
        outfile << "snpset " << i << endl;
        for (int j=0; j < ss[i].size(); j++) {
            outfile << ss[i][j] << " ";
        }
        outfile << endl;
    }
}
