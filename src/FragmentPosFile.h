//
//  FragmentPosFile.h
//  SpannerX
//
//  Created by Chip Stewart on 2/20/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//


#ifndef FRAGMENTPOSFILE_H
#define FRAGMENTPOSFILE_H

// standard includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <iterator>

// local includes
#include "UtilityFunctions.h"

using namespace std;

// single fragment class Contig class 
class FragmentPosObj
{
friend ostream &operator<<(ostream &, const FragmentPosObj &);
    
public:
    FragmentPosObj();
    FragmentPosObj( int inum, int ichr1, int istd1, int ista1, int iend1,  int ichr2, int istr2, int ista2, int iend2, int iq1, int iq2, int ifle); 
    FragmentPosObj(const vector<string> &);    
    int num;
    int chr1; 
    char strand1; 
    int start1; 
    int end1; 
    int chr2; 
    char strand2; 
    int start2; 
    int end2; 
    char qual1; 
    char qual2;
    char fle;
    FragmentPosObj&operator=(const FragmentPosObj &);
    int operator==(const FragmentPosObj &) const;
    int operator<(const FragmentPosObj &) const;
    
}; // end class 

// fragment position list
class FragmentPosFileObj
{
public:
    FragmentPosFileObj();
    FragmentPosFileObj( const string &);
    vector<vector<FragmentPosObj> > fragmentPosList;
    string posfile;

    int find(const FragmentPosObj &);
}; // end class 

#endif
