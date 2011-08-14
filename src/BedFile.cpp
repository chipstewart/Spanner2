/*
 *  BedFile.cpp
 *  Spanner
 *
 *  Created by Chip Stewart on 1/28/10.
 *  Copyright 2010 Boston College. All rights reserved.
 *
 */

#include "BedFile.h"

C_BedFile::C_BedFile() { // empty constructor
    Nchr=0;
    c.clear();
    filename="";
    header="";

}

C_BedFile::C_BedFile(string & fn, string & r1 ) { // load bed info from file 

    bool sortbed=true;
    
    this->filename=fn;
    if (fn.length()==0) {
        return;
    }
    
    C_region reg(r1);
    
    // input line
    string line;
    
    // open
    ifstream bed(this->filename.c_str(), ios::in);
    
    if (!bed) {
        cerr << "Unable to open file: " << this->filename<< endl;
        exit(1);
    }
    
    this->header = "";
    
    int nl=0;
    while (getline(bed, line)) {
        
        nl++;
        
        if (line.compare(0,3,"chr") != 0) {
            this->header=this->header+"\n"+line;
        } else {            
            C_BedRecord b1(line);
            int p2=b1.pos+b1.length;
            if (!reg.overlap(b1.chr,b1.pos,p2) ) 
                continue;                        
            if (c.count(b1.chr)<1) {
                C_BedChr emptyBed;
                emptyBed.chr=b1.chr;
                emptyBed.header=header;
                emptyBed.NX=b1.NX;
                c[b1.chr] = emptyBed;
                c[b1.chr].push_back(b1);
            } else { 
                c[b1.chr].push_back(b1);
            }
        }
    }
    
    bed.close();
    
    
    if (sortbed) {
        map<string, C_BedChr, less<string> >::iterator ib;
        ib=c.begin();
        while (ib!=c.end())  
        {
            sort((*ib).second.b.begin(), (*ib).second.b.end() );
        }
    
    }
    
            
}



//------------------------------------------------------------------------------
// Bed class constructor 
//------------------------------------------------------------------------------
C_BedChr::C_BedChr() {
    N = 0;
    plim[0] =INT_MAX;
    plim[1] =0;
    chr= "";
    header= "";
    NX=0;
    b.clear();
}





C_BedChr& C_BedChr::operator=(const C_BedChr &rhs) {
    N=rhs.N ;
    chr=rhs.chr;
    filename=rhs.filename;
    header=rhs.header;
    plim[0]=rhs.plim[0];
    plim[1]=rhs.plim[1];
    b=rhs.b ;
    return *this;
}


int C_BedChr::push_back(C_BedRecord & b1) {
    N++;
    b.push_back(b1);
    if (b1.pos<plim[0])
        plim[0]=b1.pos;
    if ((b1.pos+b1.length)<plim[1])
        plim[1]=b1.pos+b1.length;  
    return N;
}


//------------------------------------------------------------------------------
// BedRecord class constructor 
//------------------------------------------------------------------------------
C_BedRecord::C_BedRecord() {
    chr = "";
    pos=0;
    length=0;
    name="";
    score=0;
    strand=' ';
    thickStart=0;
    thickEnd=0;
    itemRgb[0]=0;
    itemRgb[1]=0;
    itemRgb[2]=0;
    blockCount=0;
    blockSizes.clear();
    blockStarts.clear();
    NX = 0;
}

C_BedRecord::C_BedRecord(string & line) {
    
    // construct a stream from the string
    stringstream strstr(line);
    // use stream iterators to copy the stream to the vector as whitespace separated strings
    istream_iterator<string> it(strstr);
    istream_iterator<string> end;
    vector<string> f(it, end);
    
    if (f.size()<3) {
        cerr << " too few fields:\n " << line << endl;
        return;
    }
    
    // chrom start end name score strand 
    // thickstart thickend RGB(3v) blockcount blocksizes(v) blockstarts(v)
    chr=f[0];
    pos=atoi(f[1].c_str());
    int p2=  atoi(f[2].c_str());
    length = p2-pos;
    if (f.size()>3) name=f[3];
    if (f.size()>4) score = atoi(f[4].c_str());
    if (f.size()>5) strand = f[5][0];
    if (f.size()>6) thickStart = atoi(f[6].c_str());
    if (f.size()>7) thickEnd = atoi(f[7].c_str());
    if (f.size()>8) { 
        vector<string> rgb;
        int nrgb=split(rgb , f[8] , "," );
        if (nrgb==3) {
            for (int irgb=0; irgb<3; irgb++) {
                itemRgb[irgb]=atoi(rgb[irgb].c_str());
            }            
        }
    }
    if (f.size()>9) blockCount = atoi(f[9].c_str());    
    if (f.size()>10) {      
        vector<string> b;
        int nb=split(b , f[10] , "," );
        for (int ib=0; ib<nb; ib++) {
            blockSizes.push_back(atoi(b[ib].c_str()));
        }            
    }
    if (f.size()>11) {   
        vector<string> b;
        int nb=split(b , f[11] , "," );
        for (int ib=0; ib<nb; ib++) {
            blockStarts.push_back(atoi(b[ib].c_str()));
        }            
    }
    NX=int(f.size());
}

C_BedRecord::C_BedRecord(string & c, int p, short l, string & n) {
    NX=4;
    chr=c;
    pos=p;
    length=l;
    name=n;    
}

C_BedRecord::C_BedRecord(string & c, int p, short l, string & n, int s) {
    NX=5;
    chr=c;
    pos=p;
    length=l;
    name=n;    
    score=s;    
}

C_BedRecord::C_BedRecord(string & c, int p, short l, string & n, int s, char w) {
    NX=6;
    chr=c;
    pos=p;
    length=l;
    name=n;    
    score=s;    
    strand=w;    
}

C_BedRecord::C_BedRecord(string & c, int p, short l, string & n, int s, char w, int ts, int te) {
    NX=8;
    chr=c;
    pos=p;
    length=l;
    name=n;    
    score=s;    
    strand=w;    
    thickStart=ts;    
    thickEnd=te;    
}

    
ostream & operator<<(ostream &output, const C_BedRecord & b1) {
    int p2=b1.pos+b1.length;
    output << b1.chr << "\t" << b1.pos << "\t"  << p2 ;
    if (b1.NX>3) output << "\t" << b1.name ;
    if (b1.NX>4) output << "\t" << b1.score ;
    if (b1.NX>5) output << "\t" << b1.strand ;
    if (b1.NX>6) output << "\t" << b1.thickStart ;
    if (b1.NX>7) output << "\t" << b1.thickEnd ;
    if (b1.NX>6) output << "\t" << b1.thickStart ;
    if (b1.NX>8) output << "\t" << b1.itemRgb[0] << ","<< b1.itemRgb[1] << ","<< b1.itemRgb[2];
    if (b1.NX>9) output << "\t" << b1.blockCount ;
    if (b1.NX>10) {
        int i=1;
        output << "\t" << b1.blockSizes[0];
        while (i<b1.blockSizes.size()) {
            output << "," << b1.blockSizes[i];   
            i++;
        }
        i=1;
        output << "\t" << b1.blockStarts[0];
        while (i<b1.blockStarts.size()) {
            output << "," << b1.blockStarts[i];        
            i++;
        }
    }
    return output;
}

C_BedRecord& C_BedRecord::operator=(const C_BedRecord &rhs)
{
    chr = rhs.chr;
    pos=rhs.pos;
    length=rhs.length;
    name=rhs.name;
    score=rhs.score;
    strand=rhs.strand;
    thickStart=rhs.thickStart;
    thickEnd=rhs.thickEnd;
    itemRgb[0]=rhs.itemRgb[0];
    itemRgb[1]=rhs.itemRgb[1];
    itemRgb[2]=rhs.itemRgb[2];
    blockCount=rhs.blockCount;
    blockSizes=rhs.blockSizes;
    blockStarts=rhs.blockStarts;
    NX=rhs.NX;
    return (*this);
}

bool C_BedRecord::operator==(const C_BedRecord &rhs)
{
    if (pos!=rhs.pos) return false;
    if (length!=rhs.length) return false;
    if (NX!=rhs.NX) return false;
    if (chr.compare(rhs.chr)!=0 ) return false;
    if ( (NX>3) & (name.compare(rhs.name)!=0 ) ) return false;
    if ( (NX>4) & (score!=rhs.score) ) return false;
    if ( (NX>5) & (strand!=rhs.strand) ) return false;
    if ( (NX>6) & (thickStart!=rhs.thickStart) ) return false;
    if ( (NX>7) & (thickEnd!=rhs.thickEnd) ) return false;
    if ( (NX>8) & (itemRgb[0]!=rhs.itemRgb[0]) ) return false;
    if ( (NX>8) & (itemRgb[1]!=rhs.itemRgb[1]) ) return false;
    if ( (NX>8) & (itemRgb[2]!=rhs.itemRgb[2]) ) return false;
    if ( (NX>9) & (blockCount!=rhs.blockCount) ) return false;
    if (NX>10) {
        if (blockSizes.size()!=rhs.blockSizes.size()) return false;
        int i=0;
        while (i<blockSizes.size()) {
            if( blockSizes[i]!=rhs.blockSizes[i]) return false;        
            i++;
        }
        if (blockStarts.size()!=rhs.blockStarts.size()) return false;
        i=0;
        while (i<blockStarts.size()) {
            if( blockStarts[i]!=rhs.blockStarts[i]) return false;        
            i++;
        }
    }
    return true;
}

bool C_BedRecord::operator<(const C_BedRecord &rhs)
{
    int a=chr.compare(rhs.chr);
    if (a!=0) return a<0;
    if (pos!=rhs.pos) return (pos<rhs.pos);
    if (length!=rhs.length) return (length<rhs.length);
    if (NX!=rhs.NX) return (NX<rhs.NX);
    return true;
}


