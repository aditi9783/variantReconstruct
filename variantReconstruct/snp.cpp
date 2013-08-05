//
//  snp.cpp
//  variantReconstruct
//
//  Created by Aditi Gupta on 7/23/13.
//  Copyright (c) 2013 Michigan State University. All rights reserved.
//

#include <iostream>
#include <string>
// #include <main>
#include "snp.h"


void snp::define_snp (int a, char b, float c, string d) {     // set values for snp
    pos = a;
    nt = b;
    frac_cov = c;
    snpkey = d;
}

float snp::hapFreqSum() {
    float sum = 0.0;
    for (vector<float>::iterator i = hfreq.begin(); i != hfreq.end(); ++i) {
        sum = sum + *i;
    }
    return sum;
}


void snp::hap_assoc (int a) {   // a is hap index, and is added to hap vector for this snp
    hap.push_back(a);
}

void snp::printSharedSNP(){    // print map that has key as snps that share a read, and value is support
    cout << "\tsnpkey\tsupport\n";
    for (map<string, int>::iterator s = shared_snp.begin(); s != shared_snp.end(); ++s) {
        cout << "\t" << s->first << "\t" << s->second << endl;
    }
    cout << "Number of different reads: " << numdiffreads << endl;
}

