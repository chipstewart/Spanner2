//
//  BamScan.cpp
//  SpannerScan
//
//  Created by Chip Stewart on 6/9/11.
//  Copyright 2011 Broad. All rights reserved.
//
#include "BamScan.h"

BamScan::BamScan(pars & Params1)	// optional constructor
{
    // parameters
    Params=Params1;
    Nread=0;
    Npair=0;
    Nproper=0;
    region.limit=false;
 
	samfile_t *fp;
    bam_header_t *bam_header;
    
    string s = Params.getInput();
    //char *c=(char*) s.c_str();
    
	if ((fp = samopen(s.c_str(), "rb", 0)) == 0) {
		fprintf(stderr, "samopen: Fail to open BAM file %s\n", s.c_str());
		exit(101);
	}
    
    // parameters
    string r = Params.getString("ChromRegion");
    int maxReads = Params.getInt("MaxReads");
    Qmin = Params.getInt("Qmin");
    LRmin = Params.getInt("MinReadLength");
    maxmismatchPC=Params.getDouble("FractionMaxMisMatches");
    FragLengthWindow=Params.getInt("FragmentLengthWindow");
    int cmd_MateMode=Params.getInt("ReadPairSenseConfig");
    string ReferenceFastaFile=Params.getString("ReferenceFastaFile");

    if (ReferenceFastaFile.size()>0) {
        FastaObj RF1(ReferenceFastaFile, "");
        Reference=RF1;
        RF1.seq.clear();  // free some memory 
    }
    
    bam_header= fp->header;
    string bamheadertext = bam_header->text;
    ReadGroup = extractBamTag(bamheadertext,"@RG");
    
    //region
    if (r.size()>0) {
        int r1,r2,r3;
        if ( bam_parse_region(bam_header, r.c_str(), &r1, &r2, &r3)==0) {
            region.limit=true;
            region.anchor=r1;
            region.start=r2;
            region.end=r3;
        } else {
            cerr << "region not found\t" << r << endl;
            exit(111);
        }
        
    }
    
    cout << ReadGroup << endl << endl;
    
    extractMateMode();
    
    if (cmd_MateMode>=0)   MateMode=cmd_MateMode;          
    
    histo_init();
    
    bam1_t *b = bam_init1();

    while (samread(fp, b) >= 0) {
        Nread++;
        fetch_func(b, fp);   
        if (Nread>=maxReads) break; 
    }
    
    
    histo_print();
    
    bam_destroy1(b);
	samclose(fp);
    
}

int BamScan::fetch_func(const bam1_t *b, void *data)
{
	samfile_t *fp = (samfile_t*)data;
	uint32_t *cigar = bam1_cigar(b);
    uint8_t *pseq =bam1_seq(b);
	const bam1_core_t *c = &b->core;
	int i, l, mm,c3,c5,ql;
	if (b->core.tid < 0) return 0;
	for (i = l = mm = c3=c5=ql= 0; i < int(c->n_cigar); ++i) {
		int op = cigar[i]&0xf;
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
			l += cigar[i]>>4;
        if (op == BAM_CMATCH || op == BAM_CDEL )
			ql += cigar[i]>>4;
		if (op == BAM_CDEL || op == BAM_CREF_SKIP)
			mm += cigar[i]>>4; 
        if (op == BAM_CSOFT_CLIP & i==0)
			c5 += cigar[i]>>4;
        if (op == BAM_CSOFT_CLIP & i>0)
			c3 += cigar[i]>>4;        
	}    
    
    /*
     printf("%s\t%d\t%d\t%s\t%d\t%c\n", fp->header->target_name[c->tid],
		   c->pos, c->pos + l, bam1_qname(b), c->qual, (c->flag&BAM_FREVERSE)? '-' : '+');
    */
    
    unsigned int p= (unsigned int) c->pos;
    unsigned short len= (unsigned short) l;
    //unsigned short querylen= (unsigned short) ql;
    unsigned short a= (unsigned short) c->tid;
    char s= (c->flag&BAM_FREVERSE)? '-' : '+';
    char q= char (c->qual);
    
    if (region.limit) {
        if (int(a)!=region.anchor) return 1;
        if (int(p)<region.start) return 1;
        if (int(p)>region.end) return 1;        
    }
    
    if (len<LRmin)  return 1;
   
    
    // mismatches
 	uint8_t *pmd, *pnm;
    char *md = NULL;
    pnm = bam_aux_get(b, "NM");
    if (pnm) {
        mm = bam_aux2i(pnm);
    } else {
        pmd= bam_aux_get(b,  "MD");    
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
    char  q2=0;  // next best mq
    unsigned short n=1;  //num mappings
    unsigned short m=mm; //mismatches
 
    
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
        double q0 = (q*double(bam_aux2i(pxs))/double(q2));
        if (q0>255) q0=255;
        q2 = char(q0);
    }

    
    // multiple maps; 
    uint8_t *px0, *px1;
    px0 = bam_aux_get(b, "X0");
    if (px0) {
        n+= bam_aux2i(px0);
    }
    px1 = bam_aux_get(b, "X1");
    if (px1) {
        n+= bam_aux2i(px1);
    }
    
    double mismatchPC=100.0*double(mm)/c->l_qseq;
    if (mismatchPC>maxmismatchPC)  return 1;

    string chrom=fp->header->target_name[c->tid];
    
    double gc=0;
    //for (int i  = 0; i < c->l_qseq; ++i) {
    for (int i=0; i<c->l_qseq; i++) {
        //cout<<i<< endl;
        int b4 = bam1_seqi(pseq,i);
		if ( (b4 == 2) || (b4 == 4) ) {
            gc+=1;
        }
    }
    gc=gc/c->l_qseq;
    
         
    BamRead r1(p, len, a, s,  q, q2,  n, m);
    r1.flag=c->flag;
    
    string qid = bam1_qname(b);
 
    Histos.h["RQ"].Fill1(r1.q);
    
    // cycle histos
    uint8_t*  bq=bam1_qual(b);
    for (int ib=0; ib<c->l_qseq; ib++) {
        int ic=ib;
        if (s=='-') {
            ic=int(c->l_qseq)-ib-1;
        } 
        Histos.h["C0"].Fill1(ic+1);
        double bq1=double(bq[ib]);
        Histos.h["C1"].FillW(ic+1,bq1);
        if (bq1>2.0) {
            Histos.h["C2"].Fill1(ic+1);
        }        
        if ( (ib>=c5) & (ib<=(int(c->l_qseq)-c3) ) ) {
            Histos.h["CM"].Fill1(ic+1);            
            Histos.h["CQ"].FillW(ic+1,bq1);
        }
        
    }

    
    if( r1.q>Qmin) {
        Histos.h["LR"].Fill1(r1.len);
        Histos.h["RS"].Fill1(r1.anchor+1);
        Histos.h["NA"].Fill1(r1.nmap);
        Histos.h["MM"].Fill1(r1.mm);
        Histos.h["GC"].Fill1(gc);
        Histos.h["CL"].Fill1(-c3);
        Histos.h["CL"].Fill1(c5);       
                        
    }
    
    // mate already in Bamreads?     
    if (SingleReads.mapping.count(qid) == 1) {
        
        Npair++;
        //cout<<"Key exists once"<<endl;
        SingleReads.iRead = SingleReads.mapping.find(qid);
        BamRead r0(SingleReads.iRead->second);
        SingleReads.mapping.erase(qid);
        BamPair fragment(r0,r1,qid, MateMode);
        bool sameA=(fragment.read[0].anchor==fragment.read[1].anchor);
        bool diffS=(fragment.read[0].sense!=fragment.read[1].sense);
        bool Qok=(fragment.read[0].q>=Qmin)&(fragment.read[1].q>=Qmin);
        int LF=fragment.FragmentLength();
        bool RefGC=Reference.seq.size()>0;

        if ((sameA)&(Qok)) {
            int PM = fragment.model;
            Histos.h["OC"].Fill1(PM);
        }
        if ((sameA)&(diffS)&(Qok)) {
            Histos.h["LF"].Fill1(LF);
            int SL=fragment.SpanLength();
            Histos.h["SL"].Fill1(SL);
            if ((LF>0)&(LF<FragLengthWindow)&RefGC) {
                gc =  extractGCfromSeq(chrom,fragment.read[0].pos,fragment.read[1].pos+fragment.read[1].len);
                if (gc>=0) {
                    Histos.h["GF"].Fill1(gc);
                    Histos.h["G0"].Fill1(LF);
                    Histos.h["G1"].FillW(LF,gc);
                    Histos.h["G2"].FillW(LF,gc*gc);
                }
            }
        }
        
        int RC0= fragment.read[0].nmap<2? fragment.read[0].nmap: 2;
        int RC1= fragment.read[1].nmap<2? fragment.read[1].nmap: 2;            
        
        // weird pair histo
        if (Qok) {    
            
             Histos.h["RC"].Fill1(RC0+RC1);
            
            if (!sameA) {
                Histos.h["WP"].FillW(4,2);
            } else if (!diffS) {
                Histos.h["WP"].FillW(3,2);
            } else {
                if (abs(LF)<=FragLengthWindow) {
                    Histos.h["WP"].FillW(1,2);
                } else {
                    Histos.h["WP"].FillW(2,1);                
                }
            }
        } else {
            Histos.h["WP"].FillW(6,2);      
            
            RC0= fragment.read[0].q<Qmin? RC0 : 0;
            RC1= fragment.read[1].q<Qmin? RC1 : 0;
            Histos.h["RC"].Fill1(RC0+RC1);
            
        }    
                
    }
                

    if (SingleReads.mapping.count(qid) == 0) {
        
        //cout<<"Key does not exist"<<endl;
        SingleReads.mapping[qid]=r1;
    } 
    
    if (SingleReads.mapping.count(qid) > 1) {
        cout<<"Key too many "<<endl;
        exit(110);
    } 
        
	return 0;
}

void BamScan::extractMateMode()
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

double BamScan::extractGCfromSeq(string & chrom, int pos0, int pos1)
{
    string s= "";
    for (size_t i=0; i<Reference.seqNames.size(); i++) {
        string s1=Reference.seqNames[i];
        size_t pos=s1.find(" ");
        if (string::npos != pos ) { 
           s1=s1.substr(0,pos);
        }
        if (chrom==s1) {
            s=Reference.seqNames[i];
            break;
        }
    }
    if (s.size()<1) return -1.0;
    //string seq=Reference.seq[s].substr(pos0,pos1-pos0);
    double gc=-1.0;
    if (pos0==pos1) return gc;
    if (pos0<0) return gc;
    if (pos1<1) return gc;
    if (pos0>int(Reference.seq[s].size())) return gc;
    if (pos1>int(Reference.seq[s].size())) return gc;
    gc=0;
    for (int b=pos0; b<pos1; b++) {
        if(Reference.seq[s][b]=='G'||Reference.seq[s][b]=='C') gc++;        
    }
    return gc/(pos1-pos0);
    
}
void BamScan::histo_init()
{
    
    // labels: 0=no map; 1=one-map; N=multi-map;
    string pairCountLab1[] = {"0-0","0-1","1-1","1-N","N-N"};
    vector<string> pairCountLab(pairCountLab1, pairCountLab1 + 5);
    
    // labels: F1R2=(first in pos order)(next); F=forward sense, R=reverse;  1=mate 1, 2=mate 2;
	string pairModelLab1[] = {"F1F2","F1R2","R1F2","R1R2","F2F1","F2R1","R2F1","R2R1"};
    vector<string> pairModelLab(pairModelLab1, pairModelLab1 + 8);
    
    string weirdPairLab1[] = {"PropPair","WidePair","FlipPair","XPair","SingleEnd","lowQ"};
    vector<string> weirdpairLab(weirdPairLab1, weirdPairLab1 + 6);
  
    
    // fragment length distribution
    int W = Params.getInt("FragmentLengthWindow");
    double FL0=-W-0.5;
    double FL1=W+0.5;
    int NW=2*W+1;
    hist fragStats;
    fragStats.Initialize(NW,FL0,FL1);  
    fragStats.setTitle("LF fragment mapping length");
    Histos.h["LF"]=fragStats;

    hist readlenStats;
    readlenStats.Initialize(1001,-0.5,1000.5);  
    readlenStats.setTitle("LR read length");
    Histos.h["LR"]=readlenStats;
    
    hist mapMultStats;
    mapMultStats.Initialize(1001,-0.5,1000.5);  
    mapMultStats.setTitle("NA read mapping multiplicity");
    Histos.h["NA"]=mapMultStats;
            
    hist SpanLenStats;
    SpanLenStats.Initialize(NW,FL0,FL1);  
    SpanLenStats.setTitle("SL span mapping length");  
    Histos.h["SL"]=SpanLenStats;
	
    hist mapQStats;
    mapQStats.Initialize(101,-0.5,100.5); 
    mapQStats.setTitle("RQ read map quality");
    Histos.h["RQ"]=mapQStats;
    
    hist pairComboStats;
    pairComboStats.Initialize(6,-0.5,5.5); 
    pairComboStats.setTitle("RC pair multiplicity combinations");
    pairComboStats.setXlabel("combo");
    pairComboStats.setBinLabels(pairCountLab);
    Histos.h["RC"]=pairComboStats;
    
    hist pairModelStats;
    pairModelStats.Initialize(8,0.5,8.5); 
    pairModelStats.setTitle("OC pair orienation model combinations");
    pairModelStats.setXlabel("model");
    pairModelStats.setBinLabels(pairModelLab);
    Histos.h["OC"]=pairModelStats;
    
    hist refStats;
    refStats.Initialize(101,-0.5,100.5); 
    refStats.setTitle("RS mapped reference");
    Histos.h["RS"]=refStats;
    
    hist weirdPairStats;
    weirdPairStats.Initialize(6,0.5,6.5); 
    weirdPairStats.setTitle("XP weird pair combos");
    weirdPairStats.setXlabel("combo");
    weirdPairStats.setBinLabels(weirdpairLab);
    Histos.h["WP"]=weirdPairStats;
    
    hist mismatchStats;
    mismatchStats.Initialize(101,-0.5,100.5); 
    mismatchStats.setTitle("MM mismatch in reads");
    Histos.h["MM"]=mismatchStats;
    
    hist GCStats;
    GCStats.Initialize(101,-0.005,1.005); 
    GCStats.setTitle("GC content in reads");
    Histos.h["GC"]=GCStats;

    hist TrimStats;
    TrimStats.Initialize(401,-200.5,200.5); 
    TrimStats.setTitle("CL soft clip trimmed ends of reads");
    TrimStats.setXlabel("-3':+5'");
    Histos.h["CL"]=TrimStats;
    
    hist GFStats;
    GFStats.Initialize(101,-0.005,1.005); 
    GFStats.setTitle("GF <GC> content averaged on fragments");
    Histos.h["GF"]=GFStats;

    hist GLStats0;
    GLStats0.Initialize(NW,FL0,FL1); 
    GLStats0.setTitle("G0 GC normalization for fragment length ");
    Histos.h["G0"]=GLStats0;    
    
    hist GLStats1;
    GLStats1.Initialize(NW,FL0,FL1); 
    GLStats1.setTitle("G1 <GC> dependence on fragment length");
    Histos.h["G1"]=GLStats1;
    
    hist GLStats2;
    GLStats2.Initialize(NW,FL0,FL1); 
    GLStats2.setTitle("G2 sigma_GC  dependence on fragment length");
    Histos.h["G2"]=GLStats2;
    
    hist CycleStats0;
    CycleStats0.Initialize(1001,-0.5,1000.5);   
    CycleStats0.setTitle("C0 Cycle base count ");
    Histos.h["C0"]=CycleStats0;

    hist CycleStats1;
    CycleStats1.Initialize(1001,-0.5,1000.5);   
    CycleStats1.setTitle("C1 Cycle base quality");;
    Histos.h["C1"]=CycleStats1;

    hist CycleStats2;
    CycleStats2.Initialize(1001,-0.5,1000.5);   
    CycleStats2.setTitle("C2 Cycle base quality>2 count");;
    Histos.h["C2"]=CycleStats2;

    hist CycleStatsM;
    CycleStatsM.Initialize(1001,-0.5,1000.5);   
    CycleStatsM.setTitle("CM Cycle base mapped ");;
    Histos.h["CM"]=CycleStatsM;

    hist CycleStatsQ;
    CycleStatsQ.Initialize(1001,-0.5,1000.5);   
    CycleStatsQ.setTitle("CQ Cycle base mapped quality ");;
    Histos.h["CQ"]=CycleStatsQ;

    
}

void BamScan::histo_single()
{
     for ( SingleReads.iRead=SingleReads.mapping.begin() ; SingleReads.iRead != SingleReads.mapping.end(); SingleReads.iRead++ ) {
         if ((*SingleReads.iRead).second.q>=Qmin) {
             Histos.h["WP"].Fill1(5);                
             Histos.h["RC"].Fill1(1);                
         } else {
             Histos.h["WP"].Fill1(6);       
             Histos.h["RC"].Fill1(0);   
         }
    }
}

void BamScan::histo_print()
{
    
    if (Histos.h.count("G1")>0) {
        for (size_t b=0; b<Histos.h["G1"].n.size(); b++) {
            double x0=Histos.h["G0"].n[b];
            double x1=Histos.h["G1"].n[b];
            double x2=Histos.h["G2"].n[b];
            if (x0>0) {
                x1=x1/x0;
                Histos.h["G1"].n[b]=x1;
                x2=sqrt(fabs(x2/x0-x1*x1));
                Histos.h["G2"].n[b]=x2;
            }
        }
        Histos.h["G1"].Finalize();
    }
    
    map<string, hist, less<string> >::iterator ih;         
    for ( ih=Histos.h.begin() ; ih != Histos.h.end(); ih++ ) {
        (*ih).second.Finalize(); 
        cout << (*ih).second << endl;
    }
}

//------------------------------------------------------------------------------
// readmap (basic info for one read alignment)
//------------------------------------------------------------------------------
BamRead::BamRead()   // Constructor
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
    string cigar1("");
    cigar=cigar1;
    flag=0;
    isize=0;

}

BamRead::BamRead(int p1, unsigned short l1, unsigned short a1, char s1, char q1, 
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


BamRead::BamRead(const BamRead &copyin)   // Copy constructor to handle pass by value.
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

ostream &operator<<(ostream &output, const BamRead & x)
{
    output << x.anchor << "\t" << x.pos << "\t" << x.len 
    << "\t" << x.sense <<"\t" << int(x.q) <<"\t" << int(x.mm) << int(x.flag) << int(x.isize);
    return output;
}

BamRead& BamRead::operator=(const BamRead &rhs)
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

int BamRead::operator==(const BamRead &rhs) const
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
int BamRead::operator<(const BamRead &rhs) const
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

//------------------------------------------------------------------------------
// BamReads  (map of BamRead objects: collection of un-mated pairs )
//------------------------------------------------------------------------------

BamReads::BamReads() {
    mapping.clear();
}

BamPair::BamPair(BamRead & r0, BamRead & r1, string & qid, char MateMode1) {
    if (r0.pos<r1.pos) {
        read[0]=r0;  
        read[1]=r1;
    } else {
        read[1]=r0;  
        read[0]=r1;        
    }
        
    convert2IlluminaShortPair(MateMode1);
    queryName=qid;
    ProperPair = (read[0].sense=='+'&read[1].sense=='-')&(r0.anchor==r1.anchor)&(read[0].pos<read[1].pos);
}

// convert platforms to Illumina short proper pair convention 
bool BamPair::convert2IlluminaShortPair(char MateMode1)
{
    model = 255;
    
    // chromosome order
    if (read[0].anchor==read[1].anchor) { 
        if (read[0].pos>read[1].pos) 
        {
            BamRead readtmp=read[1];
            read[1]=read[0];
            read[0]=readtmp;
        }
        string M="XXXX";
        M[0]=read[0].sense=='+'?'F':'R';
        M[1]=read[0].flag&BAM_FREAD1? '1' : '2';
        M[2]=read[1].sense=='+'?'F':'R';
        M[3]=read[1].flag&BAM_FREAD1? '1' : '2';
        
        if (M=="F1F2") model=1;
        if (M=="F1R2") model=2;
        if (M=="R1F2") model=3;
        if (M=="R1R2") model=4;
        if (M=="F2F1") model=5;
        if (M=="F2R1") model=6;
        if (M=="R2F1") model=7;
        if (M=="R2R1") model=8;
    }

    
    if (MateMode1==MATEMODE_454) { // 454: flip first
       
        if (read[0].flag&BAM_FREAD1) {
            read[0].sense=(read[0].sense=='+'?'-':'+');
        } else {
            read[1].sense=(read[1].sense=='+'?'-':'+');
        }
        
    } else if (MateMode1==MATEMODE_SOLID) { // SOLiD : flip second
        if (read[0].flag&BAM_FREAD1) {
            read[1].sense=(read[1].sense=='+'?'-':'+');
        } else {
            read[0].sense=(read[0].sense=='+'?'-':'+');
        }   
        
    } else if  (MateMode1==MATEMODE_ILLUMINA_LONG) { // flip both        
        read[0].sense=(read[0].sense=='+'?'-':'+');
        read[1].sense=(read[1].sense=='+'?'-':'+');
    }
    
    return true;
}

int BamPair::FragmentLength()
{
    int LF=INT_MIN;
    if (read[0].sense==read[1].sense) return LF;
    if (read[0].sense=='+') {
        LF=int(read[1].pos+read[1].len)-read[0].pos;
    } else  {
        LF=int(read[0].pos+read[0].len)-read[1].pos;
    }    
    return LF;
}

int BamPair::SpanLength()
{
    int LS=INT_MIN;
    if (read[0].sense==read[1].sense) return LS;
    if (read[0].sense=='+') {
        LS=read[1].pos-(read[0].pos+read[0].len);
    } else  {
        LS=read[0].pos-(read[1].pos+read[1].len);
    }    
    return LS;
}
