//
//  ReadPosFile.cpp
//  SpannerX
//
//  Created by Chip Stewart on 2/20/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "FragmentPosFile.h"
// constructor

FragmentPosObj::FragmentPosObj()
{
    num=0;
    chr1=0; 
    strand1=0; 
    start1=0; 
    end1=0; 
    chr2=0; 
    strand2=0; 
    start2=0; 
    end2=0; 
    qual1=0; 
    qual2=0;
    fle=0; 
}

FragmentPosObj::FragmentPosObj( int inum, int ichr1, int istd1, int ista1, int iend1, int ichr2, int istr2, int ista2, int iend2, int iq1, int iq2, int ifle)
{
    num=inum;
    chr1=ichr1; 
    strand1=istd1; 
    start1=ista1; 
    end1=iend1; 
    chr2=ichr2; 
    strand2=istr2; 
    start2=ista2; 
    end2=iend2; 
    qual1=iq1; 
    qual2=iq2;
    fle=ifle; 
}

FragmentPosObj::FragmentPosObj(const vector<string> & xs)
{
    size_t ns=xs.size();
    if (ns!=12) {
        return;
    }
    num=string2Int(xs[0]);
    chr1=string2Int(xs[1]);
    strand1=string2Int(xs[2]);
    start1=string2Int(xs[3]);
    end1=string2Int(xs[4]);
    chr2=string2Int(xs[5]);
    strand2=string2Int(xs[6]);
    start2=string2Int(xs[7]);
    end2=string2Int(xs[8]);
    qual1=string2Int(xs[9]);
    qual2=string2Int(xs[10]);
    fle=string2Int(xs[11]);

}

FragmentPosObj& FragmentPosObj::operator=(const FragmentPosObj &src)  
{                             
    num=src.num;
    chr1 = src.chr1;
    strand1=src.strand1;
    start1= src.start1;
    end1= src.end1;
    chr2 = src.chr2;
    strand2=src.strand2;
    start2= src.start2;
    end2= src.end2;
    qual1= src.qual1;
    qual2= src.qual2;
    fle= src.fle;
    return *this;
}

int FragmentPosObj::operator==(const FragmentPosObj &rhs) const
{
    if( this->chr1 != rhs.chr1) return 0;
    if( this->strand1 != rhs.strand1) return 0;
    if( this->start1 != rhs.start1) return 0;
    if( this->chr2 != rhs.chr2) return 0;
    if( this->strand2 != rhs.strand2) return 0;
    if( this->start2 != rhs.start2) return 0;
    if ((rhs.end1>0)&&(rhs.end2>0)) {
        if( this->end1 != rhs.end1) return 0;
        if( this->end2 != rhs.end2) return 0;        
    }
    return 1;
}

int FragmentPosObj::operator<(const FragmentPosObj &rhs) const
{
    if( this->chr1 < rhs.chr1 ) return 1;
    if( this->chr1 > rhs.chr1 ) return 0;
    if( this->start1 < rhs.start1 ) return 1;
    if( this->start1 > rhs.start1 ) return 0;
    if( this->chr2 < rhs.chr2 ) return 1;
    if( this->chr2 > rhs.chr2 ) return 0;
    if( this->start2 < rhs.start2 ) return 1;
    if( this->start2 > rhs.start2 ) return 0;
    if( this->qual1 < rhs.qual1 ) return 1;
    if( this->qual1 > rhs.qual1 ) return 0;
    if( this->qual2 < rhs.qual2 ) return 1;
    if( this->qual2 > rhs.qual2 ) return 0;    
    if( this->strand1 < rhs.strand1 ) return 1;
    if( this->strand2 < rhs.strand2 ) return 1;
    if ((rhs.end1>0)&&(rhs.end2>0)) {
        if( this->end1 > rhs.end1) return 1;
        if( this->end1 < rhs.end1) return 0;
        if( this->end2 > rhs.end2) return 1;
        if( this->end2 < rhs.end2) return 0;
    }
    return 0;
}

FragmentPosFileObj::FragmentPosFileObj()
{
    this->posfile="";
    fragmentPosList.clear();
}
//---------------------------------------------------------------------------------
// FastaObj is the container class for a Fasta file sequence 
//---------------------------------------------------------------------------------
// constructor
FragmentPosFileObj::FragmentPosFileObj(const string & fn)
{
	this->posfile=fn;
    if (fn.length()==0) {
        return;
    }
        
    // input line
    string line;
    
    
    //----------------------------------------------------------------------------
    // input Frag Pos file
    //----------------------------------------------------------------------------
    
    // open
    ifstream fpIn(this->posfile.c_str(), ios::in);
    
    if (!fpIn) {
        cerr << "Unable to open file: " << this->posfile << endl;
        exit(1);
    }
    
    int nl=0;
    while (getline(fpIn, line)) {
        
        nl++;
        if (line.size()<2) continue;
        
        if (nl<2) { // skip header 
            size_t found;
            found =line.find_first_of("idchrstq");
            if (found!=string::npos) continue;
        }
        
        vector<string> s1;
        int ns= split(s1 , line, " \t"); 
        if (ns==12) {
            
            FragmentPosObj fp1(s1);
            /*
            int inum=string2Int(s1[0]);
            int ichr1=string2Int(s1[1]);
            int istd1=string2Int(s1[2]);
            int ista1=string2Int(s1[3]);
            int iend1=string2Int(s1[4]);
            int ichr2=string2Int(s1[5]);
            int istr2=string2Int(s1[6]);
            int ista2=string2Int(s1[7]);
            int iend2=string2Int(s1[8]);
            int iq1=string2Int(s1[9]);
            int iq2=string2Int(s1[10]);
            int ifle=string2Int(s1[11]);            
            
            FragmentPosObj fp1(inum, ichr1, istd1, ista1, iend1, ichr2,  istr2, ista2, iend2,  iq1, iq2,  ifle);
            */
            fragmentPosList.push_back (fp1);
        }
        
    }
    //fragmentPosList.sort();
    sort(fragmentPosList.begin(),fragmentPosList.end());
    
}

int FragmentPosFileObj::find(const FragmentPosObj & fp1)
{
    
    if (binary_search  (fragmentPosList.begin(), fragmentPosList.end(), fp1)) {
        cout << "found!\n"; 
        return 1;
    } 
    return 0;        
}



