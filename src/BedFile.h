/*
 *  BedFile.h
 *  Spanner
 *
 *  Created by Chip Stewart on 1/28/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *
 */

#ifndef BEDFILE_H
#define BEDFILE_H

#include <iostream>
#include <ostream>
#include <istream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <limits.h>
#include <map>

// local
#include "UtilityFunctions.h"

using namespace std;

class C_BedRecord {
    friend ostream &operator<<(ostream &, const C_BedRecord &);
public:    
    C_BedRecord(); // empty constructor
    C_BedRecord(string &);
    C_BedRecord(string &, int, short, string &);
    C_BedRecord(string &, int, short, string &, int);
    C_BedRecord(string &, int, short, string &, int, char);
    C_BedRecord(string &, int, short, string &, int, char,int, int);
    C_BedRecord& operator=(const C_BedRecord &rhs);
    bool operator==(const C_BedRecord &rhs);
    bool operator<(const C_BedRecord &rhs);  
    string chr;
    int pos;
    short length;
    string name;
    int score;
    char strand;
    int thickStart;
    int thickEnd;
    char itemRgb[3];
    short blockCount;
    vector<int> blockSizes;
    vector<int> blockStarts;
    int NX;
};


// structure for Bed File in one chromosome
class C_BedChr {
  friend ostream &operator<<(ostream &, const C_BedChr &);
  public:
    C_BedChr(); // empty constructor
    C_BedChr(string &); // create bed record for chr
    C_BedChr(string &, string &); // load  bed record for chr
    // ~C_BedChr();    
    C_BedChr&operator=(const C_BedChr &rhs);
    C_BedChr SelectRegion(int, int);  // return records within chr window p1 p2
    int push_back(C_BedRecord &);
    void setChr(const string &);
    void setHeader(const string &);
    int N;
    string chr;
    string header;
    string filename;
    int plim[2];
    short NX;  // number of fields in bed
    vector<C_BedRecord> b;
};    


// structure for Bed File- map of chromosomes
class C_BedFile {
  friend ostream &operator<<(ostream &, const C_BedFile &);
  public:
    C_BedFile(); // empty constructor
    C_BedFile(string &); // load bed info from file 
    C_BedFile(string &,string &); // load bed info from file within region  
    ~C_BedFile();
    C_BedFile&operator=(const C_BedFile &rhs);
    C_BedChr SelectChr(string &);  // return records within chr
    C_BedChr SelectRegion(string &, int, int);  // return records within chr window p1 p2
    void setNames(const string &);
    void setHeader(const string &);
    int Nchr;
    map<string, C_BedChr, less<string> > c;
    string filename;
    string header;
};

#endif


 

