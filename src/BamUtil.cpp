//
//  BamUtil.cpp
//  SpannerScan
//
//  Created by Chip Stewart on 6/15/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "BamUtil.h"

BamUtil::BamUtil() 
{    
}

BamUtil::BamUtil(string & fname) 
{
    //samfile_t *fp;
    //bam_header_t *bam_header;
    if ((this->fp = samopen(fname.c_str(), "rb", 0)) == 0) {
		fprintf(stderr, "samopen: Fail to open BAM file %s\n", fname.c_str());
		exit(101);
	}
    bam1_t *b0 = bam_init1();
    b = b0;
    Nread=0;
    bamheadertext = this->fp->header->text;
    ReadGroup = extractBamTag(bamheadertext,"@RG");
    extractMateMode();
    
}

void BamUtil::extractMateMode()
{
    string s=upperCase(extractBamSubTag(ReadGroup,"PL:"));
    char seqtech1 = toupper(s[0]);
    size_t found=s.find("LONG");
	switch (seqtech1) {
		case '4':
			MateMode=MATEMODE_454;
			break;
		case 'S':
			MateMode=MATEMODE_SOLID;
			break;		
        case 'P':
			MateMode=MATEMODE_PACBIO;
			break;
	    case 'I':
			MateMode=MATEMODE_ILLUMINA;            
            if (found!=string::npos) { 
                MateMode=MATEMODE_ILLUMINA_LONG;
            }
			break;            
		default:
			MateMode=MATEMODE_ILLUMINA;
			break;
	}
    
}



BamContainerPair BamUtil::getNextBamPair() 
{
    //bam1_t b;
    BamContainerPair bv1;
    BamContainerPair bv=bv1;
    bv.MateMode=MateMode;
    
    while (samread(this->fp, this->b) >= 0) {
        Nread++;
        BamContainer bc1(this->b);
        string qid = bam1_qname(this->b);
        
        //cout << "BamContainerPair\t" << calcReadLength(this->b) << endl;
        

        // mate already in MappedReads?     
        if (bs.count(qid) == 1) {
            BamContainer bc2;
            bc2=bs[qid];
            bs.erase(qid);
            bv.BamEnd.push_back(bc2);
            bv.BamEnd.push_back(bc1);
            break; //return bv;
        } else { 
            bs[qid]=bc1;        
        }

    }
    
    //bv.SpanEnd.push_back(SpanRead(bv.BamEnd[0]));
    //bv.SpanEnd.push_back(SpanRead(bv.BamEnd[1]));
    // flip strands to Illumina Short convention
    //bv.Illuminize();
    //bv.calcFragmentLengths();
    
    BamContainerPair bv2=bv;
    return bv2;
}


list<BamContainer> BamUtil::getRemainingBamSingleEnds()
{
    
    list<BamContainer> bo;    
    map<string, BamContainer, less<string> >::iterator im;
    im=bs.begin();
    while (im!=bs.end())  
    {
        BamContainer b1=(*im).second; 
        //bo.push_back(b.[(*im).first]);
        bo.push_back(b1);
        bs.erase((*im).first);
        im++;
    }
    return bo;
}
 
BamOutStream::BamOutStream()
{
    filename="";
    programgroup="";
    this->fp=0;
    Npair=0;
    Nread=0;    
};

BamOutStream::BamOutStream(string & f, string & r,bam_header_t *h0)
{
    filename=f;
    programgroup=r;
   
    h = bam_header_init();
    *h = *h0;
    h->hash = h->dict = h->rg2lib = 0;
    h->text = (char*)calloc(h->l_text + 1, 1);
    memcpy(h->text, h0->text, h->l_text);
    h->target_len = (uint32_t*)calloc(h->n_targets, 4);
    h->target_name = (char**)calloc(h->n_targets, sizeof(void*));
    for (int i = 0; i < h->n_targets; ++i) {
      h->target_len[i] = h0->target_len[i];
      h->target_name[i] = strdup(h0->target_name[i]);
    }
    
    // add programgroup line to header
    size_t len = r.size();
    size_t x = h->l_text + 1;
    size_t y = h->l_text + len + 1; // 1 byte null
    kroundup32(x); 
    kroundup32(y);
    if (x < y) h->text = (char*)realloc(h->text, y);
    strncpy(h->text + h->l_text,r.data(), len); // we cannot use strcpy() here.
    h->l_text += len;
    h->text[h->l_text] = 0;
    
    
    //append_header_text(header, programgroup.c_str(),programgroup.size())
    if ((this->fp = samopen(filename.c_str(), "wb", h)) == 0) {
		fprintf(stderr, "samopen: Fail to open output BAM file %s\n", filename.c_str());
		exit(101);
	}
    b= bam_init1();    
    Npair=0;
    Nread=0;
};

BamOutStream& BamOutStream::operator=(const BamOutStream &src)  
{                   
    filename=src.filename;
    programgroup=src.programgroup;
    
    h = bam_header_init();
    *h = *(src.h);
    h->hash = h->dict = h->rg2lib = 0;
    h->text = (char*)calloc(h->l_text + 1, 1);
    memcpy(h->text, src.h->text, h->l_text);
    h->target_len = (uint32_t*)calloc(h->n_targets, 4);
    h->target_name = (char**)calloc(h->n_targets, sizeof(void*));
    for (int i = 0; i < h->n_targets; ++i) {
        h->target_len[i] = src.h->target_len[i];
        h->target_name[i] = strdup(src.h->target_name[i]);
    }
    fp=src.fp;
    Npair=src.Npair;
    Nread=src.Nread;
    return *this;
}

bool BamOutStream::write(bam1_t *b)
{
    Nread++;
    return samwrite(fp, b)>0;
}

bool BamOutStream::write(bam1_t *b1, bam1_t *b2 )
{
    Nread++;
    Nread++;
    Npair++;
    int s1=samwrite(fp, b1);
    int s2=samwrite(fp, b2);
    return (s2>0)&(s1>0);
}

void BamOutStream::close()
{
    cout << "close " << filename << ":\t" << Npair << endl; 
    samclose(fp);
}

BamContainer::BamContainer()
{
    b.core.tid=0;
    b.core.pos=0;
    b.core.bin=0;
    b.core.flag=0;
    b.core.l_qseq=0;
    b.core.mtid=0;
    b.core.mpos=0;
    b.core.isize=0;
    b.l_aux = 0;
    b.data_len = 0;
    b.m_data = 512;    
    packeddata.resize (b.m_data,0);
    
    sense='.';     // forward ('F') or reverse complement ('R')
    len=0;         // length of this read aligment in contig coordinates
    q=0;           // mapping quality
  	q2=0;          // mapping quality of next best alignment
    nmap=0;        // number of mappings for this read 
    mm=0;          // number of mismatches     
  	mob='.';       // char for special contig hit (A=alu,L=L1,S=SVA,...)
    clip[0]=0;     // clipped bases ends  5' [0]
    clip[1]=0;     // clipped bases ends 3' [1]  
    
}

BamContainer::BamContainer(bam1_t* b1)
{
    b.core = b1->core;
    b.l_aux = b1->l_aux;
    b.data_len = b1->data_len;
    b.m_data = b1->m_data;
    packeddata=vector<uint8_t>(b.data_len);
    memcpy(&packeddata[0],b1->data,b1->data_len); 
    b.data=&packeddata[0];
    
        
    // Span parts 
    uint32_t *cigar = bam1_cigar(b1);
    const bam1_core_t *c = &b1->core;
    if (b1->core.tid < 0) return;
    //int len0=bam_cigar2qlen(c, cigar);
    len =bam_cigar2qlen(c, cigar);
    mm=0;
    clip[0]=0;
    clip[1]=0;
    for (int i = 0; i < int(c->n_cigar); ++i) {
        int op = cigar[i]&0xf;
		if (op == BAM_CDEL || op == BAM_CREF_SKIP)
			mm += cigar[i]>>4; 
        if (op == BAM_CSOFT_CLIP & i==0)
            clip[0] += cigar[i]>>4;
        if (op == BAM_CSOFT_CLIP & i>0)
            clip[1] += cigar[i]>>4;        
    }    
    
    sense= (c->flag&BAM_FREVERSE)? '-' : '+';
    q= char (c->qual);
    
    // mismatches
    uint8_t *pmd, *pnm;
    char *md = NULL;
    pnm = bam_aux_get(b1, "NM");
    if (pnm) {
        mm = bam_aux2i(pnm);
    } else {
        mm=0;
        pmd= bam_aux_get(b1,  "MD");    
        if (pmd) {
            md = bam_aux2Z(pmd);
            string smd=md;
            bool insert=false;
            for (size_t i=0; i<smd.size(); i++) {
                if (isdigit(smd[i])) {
                    insert=false;
                    continue;
                } else if (smd[i]=='^') { 
                    insert=true;
                } else if (!insert) {
                    mm++;
                }
            }
        }   
    }
    
    // extra info
    q2=0;  // next best mq
    nmap=1;  //num mappings
    
    // single end mapQ when available; 
    uint8_t *psm;
    psm = bam_aux_get(b1, "SM");
    if (psm) {
        q = bam_aux2i(psm);
    }
    
    // next mapQ: q scaled by as/xs 
    uint8_t *pas,*pxs;
    pas = bam_aux_get(b1, "AS");
    pxs = bam_aux_get(b1, "XS");
    if ((psm!=0)&(pxs!=0)) {
        q2 = bam_aux2i(pas);
        double xq=round(q*double(bam_aux2i(pxs))/double(q2));
        q2 = xq>255? 255: char(xq); 
    }
    
    
    // multiple maps; 
    uint8_t *px0, *px1;
    px0 = bam_aux_get(b1, "X0");
    if (px0) {
        nmap+= bam_aux2i(px0);
    }
    px1 = bam_aux_get(b1, "X1");
    if (px1) {
        nmap+= bam_aux2i(px1);
    }
 
    
}


BamContainer& BamContainer::operator=(const BamContainer &src)  
{                             
    b = src.b;
    packeddata=src.packeddata;
    b.data=&packeddata[0];
    sense=src.sense;  
    len=src.len;  
    q=src.q;  
    q2=src.q2;  
    nmap=src.nmap;  
    mm=src.mm;  
    mob=src.mob;  
    len=src.len;  
    len=src.len;  
    clip[0]=src.clip[0];  
    clip[1]=src.clip[1];  
    return *this;
}

ostream &operator<<(ostream &output, const BamContainer & x)
{
    output << x.b.core.tid << "\t" << x.b.core.pos << "\t" << x.b.core.isize; 
    return output;
}


BamContainerPair::BamContainerPair()
{
    BamEnd.clear();
    //SpanEnd.clear();
    FragmentLength=0;
    SpanLength=0;
    ISIZE=0;
    MateMode=0;
    StrombergModel=0;
    Orientation="";
    SpanIlluminaShortConvention=false;
    BamIlluminaShortConvention=false;
    UniquePartial=false;
    UniqueMultiple=false;
    UniqueUnique=false;
    Inverted=false;    
}

BamContainerPair& BamContainerPair::operator=(const BamContainerPair &src)  
{                             
    BamEnd = src.BamEnd;
    FragmentLength=src.FragmentLength;
    SpanLength= src.SpanLength;
    ISIZE=      src.ISIZE;
    MateMode=   src.MateMode;
    StrombergModel=src.StrombergModel;
    Orientation=src.Orientation;
    BamIlluminaShortConvention=src.BamIlluminaShortConvention;
    SpanIlluminaShortConvention=src.SpanIlluminaShortConvention;
    UniquePartial=src.UniquePartial;
    UniqueMultiple=src.UniqueMultiple;
    UniqueUnique=src.UniqueUnique;
    Inverted=   src.Inverted;    
    return *this;
}


// convert MappedRead data  to Illumina short proper pair convention 
bool BamContainerPair::Illuminize(bool doBamToo)
{
    StrombergModel = 255;
    Orientation="XXXX";
    // chromosome order
    if (BamEnd[0].b.core.tid==BamEnd[1].b.core.tid) { 
    
        Orientation[0]=BamEnd[0].sense=='+'?'F':'R';
        Orientation[1]=BamEnd[0].b.core.flag&BAM_FREAD1? '1' : '2';
        Orientation[2]=BamEnd[1].sense=='+'?'F':'R';
        Orientation[3]=BamEnd[1].b.core.flag&BAM_FREAD1? '1' : '2';
        
        if (Orientation=="F1F2") StrombergModel=1;
        if (Orientation=="F1R2") StrombergModel=2;
        if (Orientation=="R1F2") StrombergModel=3;
        if (Orientation=="R1R2") StrombergModel=4;
        if (Orientation=="F2F1") StrombergModel=5;
        if (Orientation=="F2R1") StrombergModel=6;
        if (Orientation=="R2F1") StrombergModel=7;
        if (Orientation=="R2R1") StrombergModel=8;
    }
    
    
    if (MateMode==MATEMODE_454) { // 454: flip first
        
        
        if (BamEnd[0].b.core.flag&BAM_FREAD1) {
            flipSense(0,doBamToo);
            //BamEnd[0].sense=(SpanEnd[0].sense=='+'?'-':'+');
            // bam1_strand(b) ((BamEnd[0].b.core.flag&BAM_FREVERSE) != 0)
        } else {
            flipSense(1,doBamToo);
            //SpanEnd[1].sense=(SpanEnd[1].sense=='+'?'-':'+');
        }
        
    } else if (MateMode==MATEMODE_SOLID) { // SOLiD : flip second
        if (BamEnd[0].b.core.flag&BAM_FREAD1) {
            flipSense(1,doBamToo);
            //SpanEnd[1].sense=(SpanEnd[1].sense=='+'?'-':'+');
        } else {
            flipSense(0,doBamToo);
            //SpanEnd[0].sense=(SpanEnd[0].sense=='+'?'-':'+');
        }   
        
    } else if  (MateMode==MATEMODE_ILLUMINA_LONG) { // flip both  
        flipSense(0,doBamToo);
        flipSense(1,doBamToo);
        //SpanEnd[0].sense=(SpanEnd[0].sense=='+'?'-':'+');
        //SpanEnd[1].sense=(SpanEnd[1].sense=='+'?'-':'+');
    }
    BamIlluminaShortConvention=doBamToo;
    SpanIlluminaShortConvention=true;
    return true;
}


bool BamContainerPair::calcFragmentLengths()
{
    FragmentLength=INT_MIN;
    SpanLength=INT_MIN;
    ISIZE=0;
    
    //if (SpanEnd[0].anchor!=SpanEnd[1].anchor)
    if (BamEnd[0].b.core.tid!=BamEnd[1].b.core.tid)
        return false;

    //if (SpanEnd[0].sense==SpanEnd[1].sense)
    if (BamEnd[0].sense==BamEnd[1].sense)
        return false;

    if (BamEnd[0].sense=='+') {

        int p2=BamEnd[1].b.core.pos+BamEnd[1].len;
        int p1=BamEnd[0].b.core.pos;
        FragmentLength=p2-p1;
        p2=BamEnd[1].b.core.pos;
        p1=BamEnd[0].b.core.pos+BamEnd[0].len;
        SpanLength=p2-p1;
        
    } else { 
        
        int p2=BamEnd[0].b.core.pos+BamEnd[0].len;
        int p1=BamEnd[1].b.core.pos;
        FragmentLength=p2-p1;
        p2=BamEnd[0].b.core.pos;
        p1=BamEnd[1].b.core.pos+BamEnd[1].len;
        SpanLength=p2-p1;
    }
    
    if (BamEnd[1].b.core.flag&BAM_FREVERSE) {
        
        int p2=BamEnd[1].b.core.pos+BamEnd[1].len;
        int p1=BamEnd[0].b.core.pos;
        ISIZE=p2-p1;
          
    } else { 
        
        int p2=BamEnd[0].b.core.pos+BamEnd[0].len;
        int p1=BamEnd[1].b.core.pos;
        ISIZE=p2-p1;
    }

    
    return true;

}

// flip orientation of both Bam and Span records. 
bool BamContainerPair::flipSense(char e, bool doBamToo)
{
    if (e>1) 
        return false;
    if (e<0) 
        return false;    
    char em=e+1;
    if (em>1) em=0;
    
    if (!bam1_strand(&(BamEnd[e].b))) {  // F -> R        
        //SpanEnd[e].flag = SpanEnd[e].flag | BAM_FREVERSE;        
        //SpanEnd[em].flag = SpanEnd[em].flag | BAM_FMREVERSE;        
        BamEnd[e].sense='-';
        if (doBamToo) {
            BamEnd[e].b.core.flag = BamEnd[e].b.core.flag | BAM_FREVERSE;        
            BamEnd[em].b.core.flag = BamEnd[em].b.core.flag | BAM_FMREVERSE;        
        }
    } else {  // R -> F
        //SpanEnd[e].flag = SpanEnd[e].flag &  ~BAM_FREVERSE;        
        //SpanEnd[em].flag = SpanEnd[em].flag & ~BAM_FMREVERSE;        
        BamEnd[e].sense='+';
        if (doBamToo) {
            BamEnd[e].b.core.flag = BamEnd[e].b.core.flag & ~BAM_FREVERSE;                
            BamEnd[em].b.core.flag = BamEnd[em].b.core.flag & ~BAM_FMREVERSE;                
        }
    }
                            
    return true;
}

/*
// flip orientation of both Bam and Span records. 
bool BamContainerPair::flipSenseSpanOnly(char e)
{
    if (e>1) 
        return false;
    if (e<0) 
        return false;    
    char em=e+1;
    if (em>1) em=0;
    
    if (bam1_strand(&(BamEnd[e].b))) {  // F -> R        
        SpanEnd[e].flag = SpanEnd[e].flag | BAM_FREVERSE;        
        SpanEnd[em].flag = SpanEnd[em].flag | BAM_FMREVERSE;        
        SpanEnd[e].sense='-';
    } else {  // R -> F
        SpanEnd[e].flag = SpanEnd[e].flag &  ~BAM_FREVERSE;        
        SpanEnd[em].flag = SpanEnd[em].flag & ~BAM_FMREVERSE;        
        SpanEnd[e].sense='+';
    }
    
    return true;
}
*/


ostream &operator<<(ostream &output, const BamContainerPair & x)
{
    string qid = bam1_qname(&x.BamEnd[0].b);
    output << qid << "\t"<< x.Orientation << "\t" << x.BamEnd[0] << "\t" << x.BamEnd[1] << "\t" << x.FragmentLength<< "\t" << x.ISIZE;
    return output;
}

//------------------------------------------------------------------------------
// readmap (basic info for one read alignment)
//------------------------------------------------------------------------------
/*
SpanRead::SpanRead()   // Constructor
{
    pos = 0;
    anchor = 0;
    len = 0;
    sense = 0;
    q=0;
    q2=0;
    nmap=0;
    mm=0;
    mob=0;
    cigar="";
    flag=0;
    isize=0;
    clip[0]=0;
    clip[1]=0;
    
}

SpanRead::SpanRead(int p1, unsigned short l1, unsigned short a1, char s1, char q1, 
                 char qq2, unsigned short n1,unsigned short m1)
{
    pos = p1;
    anchor= a1;
    len = l1;
    sense = s1;
    q = q1;
    q2=qq2;
    nmap=n1;
    mm = m1;    
}

SpanRead::SpanRead(bam1_t * b) {
    
    pos = 0;
    anchor = 0;
    len = 0;
    sense = 0;
    q=0;
    q2=0;
    nmap=0;
    mm=0;
    mob=0;
    cigar="";
    flag=0;
    isize=0;
    clip[0]=0;
    clip[1]=0;
    
    //SpanRead r;
    uint32_t *cigar = bam1_cigar(b);
    //cigar=*cigar1;
    //uint8_t *pseq =bam1_seq(b);
    const bam1_core_t *c = &b->core;
    int i, l, mm,c3,c5;
    if (b->core.tid < 0) return;
    len=bam_cigar2qlen(c, cigar);
    int pend =bam_cigar2qlen(c, cigar);
    for (i = l = mm = c3=c5= 0; i < c->n_cigar; ++i) {
        int op = cigar[i]&0xf;
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
            l += cigar[i]>>4;
		if (op == BAM_CDEL || op == BAM_CREF_SKIP)
			mm += cigar[i]>>4; 
        if (op == BAM_CSOFT_CLIP & i==0)
            c5 += cigar[i]>>4;
        if (op == BAM_CSOFT_CLIP & i>0)
            c3 += cigar[i]>>4;        
    }    

    
    
    pos= (unsigned int) c->pos;
    len= (unsigned short) l;
    anchor= (unsigned short) c->tid;
    sense= (c->flag&BAM_FREVERSE)? '-' : '+';
    q= char (c->qual);
    flag=c->flag;    
    isize=c->isize;    
    
    // mismatches
    uint8_t *pmd, *pnm;
    char *md = NULL;
    pnm = bam_aux_get(b, "NM");
    if (pnm) {
        mm = bam_aux2i(pnm);
    } else {
        mm=0;
        pmd= bam_aux_get(b,  "MD");    
        if (pmd) {
            md = bam_aux2Z(pmd);
            string smd=md;
            bool insert=false;
            for (int i=0; i<smd.size(); i++) {
                if (isdigit(smd[i])) {
                    insert=false;
                    continue;
                } else if (smd[i]=='^') { 
                    insert=true;
                } else if (!insert) {
                    mm++;
                }
            }
        }   
    }
    
    // extra info
    q2=0;  // next best mq
    nmap=1;  //num mappings
        
    // single end mapQ; 
    uint8_t *psm;
    psm = bam_aux_get(b, "SM");
    if (psm) {
        q = bam_aux2i(psm);
    }
    
    // next mapQ: q scaled by as/xs 
    uint8_t *pas,*pxs;
    pas = bam_aux_get(b, "AS");
    pxs = bam_aux_get(b, "XS");
    if ((psm!=0)&(pxs!=0)) {
        q2 = bam_aux2i(pas);
        q2 = q*double(bam_aux2i(pxs))/double(q2);
    }
    
    
    // multiple maps; 
    uint8_t *px0, *px1;
    px0 = bam_aux_get(b, "X0");
    if (px0) {
        nmap+= bam_aux2i(px0);
    }
    px1 = bam_aux_get(b, "X1");
    if (px1) {
        nmap+= bam_aux2i(px1);
    }
                    
}

SpanRead::SpanRead(BamContainer & b1) {
    
    bam1_t* b = &b1.b;    
    SpanRead  m1(b);
    (*this) = m1;
}

SpanRead::SpanRead(const SpanRead &copyin)   // Copy constructor to handle pass by value.
{                             
    pos = copyin.pos;
    anchor= copyin.anchor;
    len = copyin.len;
    sense = copyin.sense;
    q = copyin.q;
    q2=copyin.q2;
    nmap=copyin.nmap;
    mm = copyin.mm;
    mob=copyin.mob;
    cigar=copyin.cigar;
    flag=copyin.flag;
    isize=copyin.isize;
}

ostream &operator<<(ostream &output, const SpanRead & x)
{
    output << x.anchor << "\t" << x.pos << "\t" << x.len 
    << "\t" << x.sense <<"\t" << int(x.q) <<"\t" << int(x.mm) << int(x.flag) << int(x.isize);
    return output;
}

SpanRead& SpanRead::operator=(const SpanRead &rhs)
{
    this->pos = rhs.pos;
    this->anchor = rhs.anchor;
    this->len = rhs.len;
    this->sense = rhs.sense;
    this->q = rhs.q;
    this->mm = rhs.mm;
    this->nmap = rhs.nmap;
    this->mob = rhs.mob;
    this->cigar = rhs.cigar;
    this->flag = rhs.flag;
    this->isize = rhs.isize;
    return *this;
}

int SpanRead::operator==(const SpanRead &rhs) const
{
    if( this->pos != rhs.pos) return 0;
    if( this->anchor != rhs.anchor) return 0;
    if( this->len != rhs.len) return 0;
    if( this->sense != rhs.sense) return 0;
    if( this->q != rhs.q) return 0;
    if( this->q2 != rhs.q2) return 0;
    if( this->nmap != rhs.nmap) return 0;
    if( this->mm != rhs.mm) return 0;
    if( this->flag != rhs.flag) return 0;
    if( this->isize != rhs.isize) return 0;
    return 1;
}

// This function is required for built-in STL list functions like sort
int SpanRead::operator<(const SpanRead &rhs) const
{
    if( this->anchor < rhs.anchor ) return 1;
    if( this->anchor == rhs.anchor && this->pos < rhs.pos ) return 1;
    if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->len < rhs.len) return 1;
    if( this->anchor == rhs.anchor && this->pos == rhs.pos 
       && this->len == rhs.len && this->sense < rhs.sense) return 1;
    if( this->anchor == rhs.anchor && this->pos == rhs.pos && this->len == rhs.len 
       && this->sense == rhs.sense && this->q < rhs.q) return 1;
	if( this->anchor == rhs.anchor && this->pos == rhs.pos && this->len == rhs.len 
       && this->sense == rhs.sense && this->q == rhs.q &&  this->mm < rhs.mm) return 1;
	if( this->anchor == rhs.anchor && this->pos == rhs.pos && this->len == rhs.len 
       && this->sense == rhs.sense && this->q == rhs.q &&  this->mm == rhs.mm &&  this->nmap < rhs.nmap) return 1;
    return 0;
}

*/

int calcReadLength(bam1_t * b)
{
    int i, l, mm,c3,c5;
    if (b->core.tid < 0) return 0;
    uint32_t *cigar = bam1_cigar(b);
    for (i = l = mm = c3=c5= 0; i < int(b->core.n_cigar); ++i) {
        int op = cigar[i]&0xf;
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
            l += cigar[i]>>4;
        if (op == BAM_CDEL || op == BAM_CREF_SKIP)
            mm += cigar[i]>>4; 
        if (op == BAM_CSOFT_CLIP & i==0)
            c5 += cigar[i]>>4;
        if (op == BAM_CSOFT_CLIP & i>0)
            c3 += cigar[i]>>4;        
    }   
    return l;
}

