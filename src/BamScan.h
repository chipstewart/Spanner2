//
//  BamScan.h
//  SpannerScan
//
//  Created by Chip Stewart on 6/9/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
#ifndef BAMSCAN_H
#define BAMSCAN_H
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


// Technology mate modes
#define MATEMODE_ILLUMINA        0
#define MATEMODE_454             1
#define MATEMODE_SOLID           2
#define MATEMODE_ILLUMINA_LONG   3
#define MATEMODE_PACBIO          4

//using namespace re2;
using namespace std;

typedef struct {
    bool limit;
    int anchor;
    int start;
    int end;
 } region_t;

// Basic mapped read record
class BamRead
{
    friend ostream &operator<<(ostream &, const BamRead &);
public:
    int pos;           // position in contig (unpadded)
    unsigned short len;         // length of this read aligment in contig coordinates
    unsigned short anchor;      // anchor index  
    char sense;                 // forward ('F') or reverse complement ('R')
    char q;                     // mapping quality
  	char q2;                    // mapping quality of next best alignment
    unsigned short nmap;        // number of mappings for this read 
    unsigned short mm;          // number of mismatches     
  	char mob;                   // char for special contig hit
  	string cigar;               // cigar
    unsigned int flag;          // bam flag
    unsigned int isize;         // bam isize
    
    BamRead();  
    BamRead(const BamRead &);  
    BamRead(int p, unsigned short l, unsigned short a, char s, char q,char q2, unsigned short n,unsigned short m);
    ~BamRead(){};  
    BamRead&operator=(const BamRead &rhs);
    int operator==(const BamRead &rhs) const;
    int operator<(const BamRead &rhs) const;  
}; 



class BamReads
{
    public:
    BamReads();  
    map<string, BamRead, less<string> > mapping;
    map<string, BamRead, less<string> >::iterator iRead;
};

class BamPair
{
    public:
    BamPair(BamRead &, BamRead &, string &, char);  
    BamRead read[2];  
    string  queryName;
    char model; 
    bool ProperPair;
    bool convert2IlluminaShortPair(char);
    int FragmentLength();   
    int SpanLength();
};

class BamPairs
{
    public:
    list<BamPair> Pairs;
};

class BamScan
{
	friend ostream &operator<<(ostream &, const BamScan &);
    
public:
	
    BamScan();	// default constructor
    BamScan(pars & );	// optional constructor
    int fetch_func(const bam1_t *b, void *data);
    void histo_init();
    void histo_single();
    void histo_print();
    void extractMateMode();
    //bool extractChromRegion();
    double extractGCfromSeq(string &, int, int);
    pars  Params;    
    BamReads SingleReads;
    BamPairs Fragments;
    FastaObj Reference;
    int Nread;
    int Npair;
    int Nproper;
    string ReadGroup;
    hists Histos;
    char  MateMode;
    // 0: FR (illumina RP short); 1: R1R2/F2F1 (454);  2: F1F2/R2R1 (SOLiD); 3: RF (illumina long)
    int Qmin;
    int LRmin;
    double maxmismatchPC;
    int FragLengthWindow;
    region_t region;
}; // end class 

#endif
