//
//  snp.h
//  variantReconstruct
//
//  Created by Aditi Gupta on 7/18/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//

#ifndef variantReconstruct_snp_h
#define variantReconstruct_snp_h

#include <vector>
#include <map>
#include <string>

using namespace std;

class snp{
    int pos;
    char nt;
    float frac_cov;     // number of reads with snp base/coverage
    int snpbin;
    int numdiffreads = 0; // number of different reads that the snp was found on
    int numsharedsnps = 0; // number of times a snp share reads with other snps
    string snpkey;
    vector<int> hap;    //list of hap indexes that this snp is part of
    vector<float> hfreq; // frequencies of haps that this snp is part of (to check if snp frac_cov matches sum of hap freq)
    map<string, int> shared_snp; // key: snp that shares the same read, val: number of times this snp shares the read
    
public:
    void define_snp(int, char, float, string); 
    int position() { return pos; }
    char base() { return nt; }
    float cov() { return frac_cov; }
    string getsnpkey() { return snpkey; }
    void updateHfreq( float a ) { hfreq.push_back(a); }
    float hapFreqSum();
    void updateSharedSNP (string newsnp) { shared_snp[newsnp]++; numsharedsnps++; }
    void updateNumDiffReads() { numdiffreads++; };
    int returnNumReads() { return numdiffreads; };
    int returnNumSharedSNPs() { return numsharedsnps; }
    void deleteSharedSNP (string newsnp) { shared_snp.erase(newsnp); } 
    void printSharedSNP();
    map<string,int> returnSharedSNP() { return shared_snp; }
    void setbin( int a ) { snpbin = a; }
    int getbin() { return snpbin; }
    void hap_assoc(int);
    vector<int> gethaps() { return hap; }
    
};




#endif
