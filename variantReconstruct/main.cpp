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

#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include <boost/graph/graph_concepts.hpp>


using namespace std;
using namespace boost;


struct all_clique_visitor
{
    vector< deque<unsigned long> > all_cliques;
    
    template <typename Clique, typename Graph>
    inline void clique(const Clique& p, const Graph& g)
    {
        all_cliques.push_back(p);
    }
    
    inline vector<deque<unsigned long> > returnAllCliques() 
    { 
        return all_cliques;
    }    
        
};

struct max_clique_visitor
{
    max_clique_visitor(std::size_t& max)
    : maximum(max)
    { }
    
    template <typename Clique, typename Graph>
    inline void clique(const Clique& p, const Graph& g)
    {
        BOOST_USING_STD_MAX();
        maximum = max BOOST_PREVENT_MACRO_SUBSTITUTION (maximum, p.size());
    }
    std::size_t& maximum;
};
 
vector<string> shallow_snps; // snps at positions where depth is < 1000.
vector<string> snps_to_remove; // snpkeys of snps that are removed beacuse of low depth and fewer connectivity to other snps
map< string, snp > snpmap;   // key: string with snp pos and base, value: corresponding snp object
int numsnp = 0;
map<string, int> snpkey_index; // numerical index for each snpkey, starting from 0 till numsnp-1;
map<int, vector<string> > snpbins; // key: bin identifier, value: list of snps in that bin
int LINKTHRES = 2;
int NUMREADTHRES = 2;
float MINSNPFREQ = 0.002;
int MINDEPTH = 1000; // coverage depth below which a snp is considered 'shallow', and not used as building blocks of variants

ofstream outfile; // output file
ofstream linkmat_out; // contains linkmatrix
ofstream truesnps; // print snps that have sufficient read depth, snp freq, and connectivity to other snps


// declare functions
void hapread_open(const char *);
void snp_open(const char *);
void linkmatrix();
void removeWeakSNPs();
map<int, vector<float> > identifybins();
void assignbins( map<int, vector<float> > );
void printSNPbins();

void getHaps();
vector<vector<string> > disjointedSNPs ( vector<string> ); // take a list of snps, return disjointed set of these snps
vector<string> BFS ( vector<string> );
void printSNPset ( vector<vector<string> > );

void getHaps2();
vector<string> sortedSNPsByFreq();

void getHaps3();


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
        
       // getHaps2();
        getHaps3();
        
        /*
        bindefined = identifybins();
        assignbins( bindefined );
        printSNPbins();
        getHaps();
        */
    }
    
    return 0;
}

void getHaps3 () {
    vector<string> sortedsnps;
    
    sortedsnps = sortedSNPsByFreq();
    // read snps in the ascending order of their freq
    // build adjacency matrix for all snps that share a read with the snp in consideration
    // find cliques in the graph represented by the adjacency matrix -> these are probable variants
    // also save connected components that do not form clique
    
    for (vector<string>::size_type i = 0; i != sortedsnps.size(); i++) {
        snpgraph snplinks;
        
        snp s = snpmap[ sortedsnps[i] ];
        s.printSharedSNP();
        map<string, int> shared_snps = s.returnSharedSNP();
        
        snplinks.setgraphcenter(sortedsnps[i]);
        snplinks.initializeGraph( shared_snps );
        // for snps that share a read with current snp, see if those snps are linked as well
        for (map<string, int>::iterator ss = shared_snps.begin(); ss != shared_snps.end(); ++ss) {
            snp nbgsnp = snpmap[ ss->first ];
            snplinks.updateGraph(ss->first, nbgsnp.returnSharedSNP());
        }
        
        snplinks.printGraph();
        
        map<string, int> nodes = snplinks.returnGraphNodes();
        vector<string> vertexvec;
        vertexvec.resize( nodes.size() );
        vector<vector<int> >adjmat = snplinks.returnAdjMat();
        // create vector list of nodes
        for (map<string, int>::iterator n= nodes.begin(); n != nodes.end(); ++n) {
            vertexvec[n->second] = n->first;
        }
        int num_vertices = vertexvec.size(); // number of nodes in the graph
        
        // create a typedef for the Graph type
      //  typedef adjacency_list<vecS, vecS, undirectedS> Graph;
        typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
        
        Graph g(num_vertices);

        all_clique_visitor acv;
        
        // create edge list
        typedef std::pair<int, int> Edge;
        for (int k = 0; k < num_vertices; k++) {
            for (int j = k+1; j < num_vertices; j++) {
                if (adjmat[k][j] > 0) {
                    add_edge(k, j, g);
                }
            }
        }
        
        // get the property map for vertex indices
        typedef property_map<Graph, vertex_index_t>::type IndexMap;
        IndexMap index = get(vertex_index, g);
        
        std::cout << "vertices(g) = ";
        typedef graph_traits<Graph>::vertex_iterator vertex_iter;
        std::pair<vertex_iter, vertex_iter> vp;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first)
            std::cout << index[*vp.first] <<  " ";
        std::cout << std::endl;
        
      //  std::size_t ret = 0;
      //  bron_kerbosch_all_cliques(g, find_max_clique(ret));

        bron_kerbosch_all_cliques(g, acv, 2);
        /*
        vector<deque<unsigned long> > all_cliques = acv.returnAllCliques();
        for (vector<deque<unsigned long> >::size_type c = 0; c != all_cliques.size(); c++) {
            cout << "clique " << c << endl;
            for (deque<unsigned long>::size_type cs = 0; cs != all_cliques[c].size(); cs++) {
                cout << all_cliques[c][cs] << " ";
            }
            cout << endl;
        }
         */
        
    }
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
    
    sortedsnps = sortedSNPsByFreq();
    // read snps in the ascending order of their freq
    // build adjacency matrix for all snps that share a read with the snp in consideration
    // find cliques in the graph represented by the adjacency matrix -> these are probable variants
    // also save connected components that do not form clique
    
    for (vector<string>::size_type i = 0; i != sortedsnps.size(); i++) {
        snpgraph snplinks;
        
        snp s = snpmap[ sortedsnps[i] ];
        s.printSharedSNP();
        map<string, int> shared_snps = s.returnSharedSNP();
            
        snplinks.setgraphcenter(sortedsnps[i]);
        snplinks.initializeGraph( shared_snps );
        // for snps that share a read with current snp, see if those snps are linked as well
        for (map<string, int>::iterator ss = shared_snps.begin(); ss != shared_snps.end(); ++ss) {
            snp nbgsnp = snpmap[ ss->first ];
            snplinks.updateGraph(ss->first, nbgsnp.returnSharedSNP());
        }
        
        snplinks.printGraph();
    }
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

map<int, vector<float> > identifybins() { // identify boundries of snp frequency bins
    
    // get list of frequencies for concerned snps
    vector<float> freqs;
    map<int, vector<float> > freqbins; // key: bin identifier, or rank, value: list of freqs in that bin
    map<int, float> freqbin_means; // mean value of frequencies in a bin
    int bincounter = 0;
    map<int, vector<float> > binboundries; // key: bin identifier, value: min and max freqs in this bin
    
    for (map<string, snp>::iterator m = snpmap.begin(); m != snpmap.end(); ++m) {
        freqs.push_back(m->second.cov());
   //     outfile << m->second.cov() << endl;

        if (m->second.cov() == 0.0) {
            outfile << m->first << endl;
        }
    }
    
    sort(freqs.begin(), freqs.end());

    // add first freq to first bin
    freqbins[ bincounter ].push_back(freqs[0]);
    freqbin_means[ bincounter ] = freqs[0];
    bincounter++;
    outfile << "first bin entry " << freqs[0] << endl;
    
    
    for (int i=1; i<freqs.size(); i++) {
        // check which bin is the best fit- i.e. freq within a fifth of bin mean freq
        int flag = 0; // flag = 1 if a bin is assigned, else stays 0 => create new bin
        for (map<int, float>::iterator bm=freqbin_means.begin(); bm != freqbin_means.end(); ++bm) {
            float upper = bm->second + bm->second/5.0;
            float lower = bm->second - bm->second/5.0;
            if (freqs[i] >= lower and freqs[i] <= upper) {
                freqbins[ bm->first ].push_back( freqs[i] );
                flag = 1;
                break;
            } else {
                continue;
            }
        }
        if (flag == 0) { // i.e., no suitable bin was found
            freqbins[ bincounter ].push_back(freqs[i]);
            freqbin_means[ bincounter ] = freqs[i];
            bincounter++;
        }
    }
    // print bins
    outfile << "\nPrinting snp bins\n";
    for (map<int, vector<float> >::iterator b=freqbins.begin(); b != freqbins.end(); ++b) {
        outfile << "\n" << b->first << " : " << endl;
        binboundries[ b->first ].push_back( b->second[0] ); // min freq value in bin
        binboundries[ b->first ].push_back( b->second[ b->second.size()-1 ] ); // last (max) freq in bin
        for (int i = 0; i <b->second.size(); i++) {
            outfile << b->second[i] << ", ";
        }
    }    
    return binboundries;
}

void assignbins ( map<int, vector<float> > bindef ) {
    
    // for each snp in snpmap, find the bin based on snp freq
    
    for (map<string, snp>::iterator s = snpmap.begin(); s != snpmap.end(); ++s) {
        float snpf = s->second.cov();
        cout << s->second.getsnpkey() << endl;
        s->second.printSharedSNP();
        
        for (map<int, vector<float> >::iterator b = bindef.begin(); b != bindef.end(); ++b) {
            if ( snpf >= b->second[0] and snpf <= b->second[1] ) { // snpfreq is within the upper and lower range of bin
                s->second.setbin( b->first );
                // update global variable snpbins
                snpbins[ b->first ].push_back( s->second.getsnpkey() );
                break;
            }
        }
    }
}

void printSNPbins () {
    for (map<int, vector<string> >::iterator s= snpbins.begin(); s != snpbins.end(); ++s) {
        outfile << "bin is " << s->first << endl;
        for (int i = 0; i < s->second.size(); i++ ) {
            outfile << s->second[i] << endl;
        }
    }
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
