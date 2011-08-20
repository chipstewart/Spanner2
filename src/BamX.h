//
//  BamX.h
//  Spanner
//
//  Created by Chip Stewart on 6/9/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
#ifndef BAMX_H
#define BAMX_H
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <string>
#include <algorithm>
#include <vector>
#include <list>
#include <map>
#include <iterator>

#include "UtilityFunctions.h"
#include "Histo.h"
#include "Pars.h"
#include "FastaFile.h"
#include "sam.h"
#include "BamUtil.h"


//using namespace re2;
using namespace std;

class BamX
{
    
public:
	
    BamX();	// default constructor
    BamX(pars & );	// optional constructor
    //int fetch_func(const bam1_t *b, void *data);
    pars  Params;    
    BamUtil Bam;
    FastaObj Reference;
    map<string,BamHeaderContainer,less<string> > outputBam;
    map<string,BamHeaderContainer,less<string> >::iterator ioutputBam;
    int Nread;
    int Npair;
    int Nproper;
    string ReadGroup;
    hists Histos;
    char  MateMode;
    // 0: FR (illumina RP short); 1: R1R2/F2F1 (454);  2: F1F2/R2R1 (SOLiD); 3: RF (illumina long)
    int Qmin;
    int LRmin;
    int LFlow;
    int LFhigh;
    bool IlluminizeBam;
    double maxmismatchPC;
    int FragLengthWindow;
    C_region region;
    double FragmentTailPercent;
    string outputDirectory;
    
    bool outFragTailBam;
    bool outInterChromBam;
    bool outUniqueMultipleBam;
    bool outUniquePartialBam;
    bool outUniqueUnmappedBam;
    
}; // end class 

#endif
