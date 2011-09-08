//
//  BamUtil.h
//  SpannerScan
//
//  Created by Chip Stewart on 6/15/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
#ifndef BAMUTIL_H
#define BAMUTIL_H
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
#include "sam.h"


// Technology mate modes
#define MATEMODE_ILLUMINA        0
#define MATEMODE_454             1
#define MATEMODE_SOLID           2
#define MATEMODE_ILLUMINA_LONG   3
#define MATEMODE_PACBIO          4

int  calcReadLength(bam1_t *);

// container class to stick bam1_t allocated data into a string object



class BamContainer
{
    friend ostream &operator<<(ostream &, const BamContainer &);
public: 
    BamContainer();
    BamContainer(bam1_t*);
    BamContainer&operator=(const BamContainer &);
    
    bam1_t b;                   // bam record
    vector<uint8_t> packeddata; // variable length part stashed into stl vector  
    char sense;                 // forward ('+') or reverse complement ('-')
    unsigned short len;         // length of this read aligment in contig coordinates
    char q;                     // mapping quality
  	char q2;                    // mapping quality of next best alignment
    unsigned short nmap;        // number of mappings for this read 
    unsigned short mm;          // number of mismatches     
  	char mob;                   // char for special contig hit (A=alu,L=L1,S=SVA,...)
    int clip[2];                // clipped bases ends  5' [0] and 3' [1]  
};


class BamContainerPair
{
   friend ostream &operator<<(ostream &, const BamContainerPair &);
    
public: 
    BamContainerPair();
    BamContainerPair(BamContainer &, BamContainer &);
    BamContainerPair&operator=(const BamContainerPair &);
    bool Illuminize(bool);
    //bool IlluminizeSpanOnly();
    bool calcFragmentLengths();
    bool flipSense(char, bool);
    //bool flipSenseSpanOnly(char);

    vector<BamContainer> BamEnd;
    //vector<SpanRead> SpanEnd;
    int FragmentLength;
    int SpanLength;
    int ISIZE;
    int MateMode;
    int StrombergModel;
    string Orientation;
    bool SpanIlluminaShortConvention;
    bool BamIlluminaShortConvention;
    bool UniquePartial;
    bool UniqueMultiple;
    bool UniqueUnique;
    bool Inverted;  
 
};

/*
class BamOutStream 
{
    
public:
    BamOutStream();
    BamOutStream(string &, string &,bam_header_t *);
    BamOutStream&operator=(const BamOutStream &);
    bool write(bam1_t *);
    bool write(bam1_t *, bam1_t *);
    void close();
    string filename;
    string programgroup;
    samfile_t *fp;
    bam_header_t *h;
    bam1_t *b;
    int Npair;
    int Nread;    
};
*/

class BamHeaderContainer 
{
    
public:
    BamHeaderContainer();
    BamHeaderContainer(bam_header_t *,string &);
    BamHeaderContainer&operator=(const BamHeaderContainer &);
    bam_header_t* header();
    //bam_header_t h;
    vector<string> target_name;
    vector<uint32_t> target_len;
    string text;
    int32_t n_targets;
    size_t l_text, n_text;    
};

class BamFileContainer 
{
    
public:
    BamFileContainer();
    BamFileContainer(string &, string &,bam_header_t *);
    BamFileContainer&operator=(const BamFileContainer &);
    bool write(bam1_t *);
    bool write(bam1_t *, bam1_t *);
    void close();
    samfile_t* fp();
    string filename;
    string programgroup;
    samfile_t fp1;    
    BamHeaderContainer h;
    int Npair;
    int Nread;    
};

// Basic mapped read record
class BamUtil
{
    friend ostream &operator<<(ostream &, const BamUtil &);
    
public:
    
    BamUtil();
    BamUtil(string &);
    BamContainerPair  getNextBamPair();
    list<BamContainer> getRemainingBamSingleEnds();
    //vector<BamOutStream> outputBam;
    
    void extractMateMode();
    
    samfile_t* fp;
    bam1_t *b;
    
    string filename; 
    int Nread;
    int MateMode;
    string ReadGroup;
    string bamheadertext; 
    map<string, BamContainer, less<string> > bs;
    
}; 




#endif

