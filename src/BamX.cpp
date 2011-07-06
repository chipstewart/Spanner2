//
//  BamX.cpp
//  Spanner
//
//  Created by Chip Stewart on 6/9/11.
//  Copyright 2011 Broad. All rights reserved.
//
#include "BamX.h"

BamX::BamX(pars & Params1)	// optional constructor
{
    // parameters
    Params=Params1;
    Nread=0;
    Npair=0;
    Nproper=0;
    LFlow=INT_MIN;
    LFhigh=INT_MAX;
    region.limit=false;
    IlluminizeBam=0;

    outFragTailBam=false;
    outInterChromBam=false;
    outUniqueMultipleBam=false;
    outUniquePartialBam=false;
    outUniqueUnmappedBam=false;
 
    //output file
	//samfile_t *fp;
    bam_header_t *bam_header;
    
    string s = Params.getInput();
    BamUtil bam1(s);
    Bam = bam1;

    string filename=extractfilename(s);
    
    // parameters
    string r = Params.getString("ChromRegion");
    int maxReads = Params.getInt("MaxReads");
    Qmin = Params.getInt("Qmin");
    LRmin = Params.getInt("MinReadLength");
    maxmismatchPC=Params.getDouble("FractionMaxMisMatches");
    FragLengthWindow=Params.getInt("FragmentLengthWindow");
    int cmd_MateMode=Params.getInt("ReadPairSenseConfig");
    string ReferenceFastaFile=Params.getString("ReferenceFastaFile");
    FragmentTailPercent=Params.getDouble("FragmentTailPercent");
    IlluminizeBam=Params.getInt("Illuminize")>0;
    outputDirectory=Params.getString("OutputDirectory");
    int minLR=Params.getInt("MinReadLength"); 
    
    
    string StatFile=Params.getString("StatFile");
    if (StatFile.size()>0) {
        hists H1(StatFile);
        hist HLF=H1.h["LF"];
        hist HLR=H1.h["LR"];
        Params.setHist("LF",HLF);
        Params.setHist("LR",HLR);        
        H1.h.clear();  // free some memory 
        if (FragmentTailPercent>0) {
            LFlow=int(HLF.p2xTrim(FragmentTailPercent/100.));   
            LFhigh=int(HLF.p2xTrim(1-FragmentTailPercent/100.));   
        }
    }

    if (ReferenceFastaFile.size()>0) {
        FastaObj RF1(ReferenceFastaFile, "");
        Reference=RF1;
        RF1.seq.clear();  // free some memory 
    }
    
    bam_header= Bam.fp->header;
    string bamheadertext = bam_header->text;
    ReadGroup = extractBamTag(bamheadertext,"@RG");
    
    outFragTailBam=FragmentTailPercent>0;
    outInterChromBam=true;
    outUniqueMultipleBam=true;
    outUniquePartialBam=true;
    outUniqueUnmappedBam=true;
    // output Bams
    outputBam.clear();
    
    if (outFragTailBam) {
        string outfile=outputDirectory+"/"+filename+".fragtail.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:FragmentTail\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["FT"]=BamOutStream(outfile, q , bam_header);  
    }
     
    if (outInterChromBam) {
        string outfile=outputDirectory+"/"+filename+".interchrom.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:InterChromPairs\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["IC"]=BamOutStream(outfile, q , bam_header);  
    }
    
    if (outUniqueMultipleBam) {
        string outfile=outputDirectory+"/"+filename+".uMult.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:uniqMultiplyMappedPairs\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["UM"]=BamOutStream(outfile, q , bam_header);  
    }
    
    if (outUniquePartialBam) {        
        string outfile=outputDirectory+"/"+filename+".uPart.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:uniqPartiallyMappedPairs\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["UP"]=BamOutStream(outfile, q , bam_header);  
    }

    if (outUniqueUnmappedBam) {        
        string outfile=outputDirectory+"/"+filename+".uUnmapped.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:uniqUnMappedPairs\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["UZ"]=BamOutStream(outfile, q , bam_header);  
    }

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
    
    //extractMateMode();
    
    if (cmd_MateMode>=0)   MateMode=cmd_MateMode;          
     
    BamContainerPair bampair;
    
    bool more = true;
    while (more)
    {
        bampair=Bam.getNextBamPair();
        bampair.Illuminize(IlluminizeBam);  
        bampair.calcFragmentLengths();
        Npair++;
        if (Npair>=maxReads) break; 
        more=(bampair.BamEnd[1].packeddata.size()>1);
        //if (bampair.BamEnd[0].b.core.tid==bampair.BamEnd[1].b.core.tid) 
        //    cout<< bampair << endl;
        
        if (outFragTailBam) { 
            if ((bampair.BamEnd[0].q<Qmin)|(bampair.BamEnd[1].q<Qmin)) continue;
            bool FT=(bampair.FragmentLength>LFhigh)|((bampair.FragmentLength<LFlow)&(bampair.FragmentLength>INT_MIN));
            if (FT) {
                outputBam["FT"].write(&(bampair.BamEnd[0].b),&(bampair.BamEnd[1].b));
            }
        }
        if (outInterChromBam) { 
            if ((bampair.BamEnd[0].q<Qmin)|(bampair.BamEnd[1].q<Qmin)) continue;
            bool IC=(bampair.BamEnd[0].b.core.tid!=bampair.BamEnd[1].b.core.tid);
            if (IC) {
                outputBam["IC"].write(&(bampair.BamEnd[0].b),&(bampair.BamEnd[1].b));
            }
        }
        if (outUniqueMultipleBam) { 
            if ((bampair.BamEnd[0].q<Qmin)&(bampair.BamEnd[1].q<Qmin)) continue;            
            bool UM=(bampair.BamEnd[0].nmap!=bampair.BamEnd[1].nmap);            
            if (UM) {
                outputBam["UM"].write(&(bampair.BamEnd[0].b),&(bampair.BamEnd[1].b));
            }
        }
        if (outUniquePartialBam) { 
            if ((bampair.BamEnd[0].q<Qmin)|(bampair.BamEnd[1].q<Qmin)) continue;
            int c0=bampair.BamEnd[0].clip[0]+bampair.BamEnd[0].clip[1];
            int c1=bampair.BamEnd[0].clip[0]+bampair.BamEnd[0].clip[1];
            if  (((c1>minLR)&(c0>minLR))|((c0<minLR)&(c1<minLR))) continue; 
            bool UP=((c1+c0)>minLR);
            if (UP) {
                outputBam["UP"].write(&(bampair.BamEnd[0].b),&(bampair.BamEnd[1].b));
            }
        }
        if (outUniqueUnmappedBam) { 
            if ((bampair.BamEnd[0].q>Qmin)|(bampair.BamEnd[1].q>Qmin)) continue;
            bool z0=((bampair.BamEnd[0].b.core.flag&BAM_FUNMAP)>0);
            bool z1=((bampair.BamEnd[1].b.core.flag&BAM_FUNMAP)>0);            
            bool UZ=(z0|z1)&(!(z1&z0));
            if (UZ) {
                outputBam["UZ"].write(&(bampair.BamEnd[0].b),&(bampair.BamEnd[1].b));
            }
        }
      
        
        //cout<< bampair.Orientation << "\t"<< bampair.FragmentLength << "\t" <<bampair.BamEnd[1].b.core.pos << endl;
        
    }
    
    
    for (ioutputBam=outputBam.begin(); ioutputBam!=outputBam.end(); ioutputBam++) {
        (*ioutputBam).second.close();
    }
    /*
    if (FragmentTailPercent>0) 
        outputBam["FT"].close();
    */
    
    samclose(Bam.fp);
    
    
}



