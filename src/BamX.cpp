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
    Nout=0;
    LFlow=INT_MIN;
    LFhigh=INT_MAX;
    region.limit=false;
    IlluminizeBam=0;

    outFragTailBam=false;
    outInterChromBam=false;
    outUniqueMultipleBam=false;
    outUniquePartialBam=false;
    outUniqueUnmappedBam=false;
    outAllPairsBam=false;
    outReadPairPosBam=false;
    
    //output file
	//samfile_t *fp;
    bam_header_t *bam_header;
    
    string s = Params.getInput();
    BamUtil bam1(s);
    Bam = bam1;

    string filename=extractfilename(s);
    
    // parameters
    string fragPosFile = Params.getString("ReadPairPosFile");
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
    int SplitBracketMin=Params.getInt("SplitBracketMin"); 
    int SplitBaseQmin=Params.getInt("SplitBaseQmin"); 
    
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
    
    int dbg = Params.getInt("Dbg");
    time(&tprev);
    
    if (ReferenceFastaFile.size()>0) {
        FastaObj RF1(ReferenceFastaFile, "");
        Reference=RF1;
        RF1.seq.clear();  // free some memory 
    }
    
    bam_header= Bam.fp->header;
    string bamheadertext = bam_header->text;
    ReadGroup = extractBamTag(bamheadertext,"@RG");
    
    outAllPairsBam=(r.size()>0);
    if (!outAllPairsBam) { 
        outFragTailBam=true; //FragmentTailPercent>=0;
        outInterChromBam=true;
        outUniqueMultipleBam=true;
        outUniquePartialBam=true;
        outUniqueUnmappedBam=true;
    }
    // output Bams
    outputBam.clear();
    
    /*
    // test BamHeaderContainer
    vector<BamHeaderContainer> x;
    string sv=SpannerVersion;    
    string q="@PG\tID:FragmentTail\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
    while (true) {
        string outfile=outputDirectory+"/"+filename+".fragtail.bam";
        q=q+"\n@PG\tID:FragmentTail\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        BamHeaderContainer x1( bam_header, q); 
        x.push_back(x1);
        bam_header_t* h1=x[x.size()-1].header();
        cout<< h1->text << endl;
    }
    cout << x.size() << endl;
    */
    
    samfile_t *fpFT=0;
    samfile_t *fpIC=0;
    samfile_t *fpUM=0;
    samfile_t *fpUP=0;
    samfile_t *fpUZ=0;
    samfile_t *fpAP=0;
    samfile_t *fpWP=0;
    
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
    
    
    //fragPosFile
    if (fragPosFile.size()>0) {
        
        FragmentPosFileObj fp(fragPosFile);
        if (fp.fragmentPosList.size()>0) {
            FragPos=fp;
        } else {
            cerr << "Read Pair Pos file not found\t" <<  fragPosFile << endl;
            exit(112);
        }
        outFragTailBam=false; 
        outInterChromBam=false;
        outUniqueMultipleBam=false;
        outUniquePartialBam=false;
        outUniqueUnmappedBam=false;
        outReadPairPosBam=true;
        
    }

    
    if (outAllPairsBam) {
        string outfile=outputDirectory+"/"+filename+"."+r+".bam";
        string sv=SpannerVersion;
        string q="@PG\tID:Region\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["AP"]=BamHeaderContainer(bam_header,q); 
        bam_header_t* h1=outputBam["AP"].header();
        if ((fpAP = samopen(outfile.c_str(), "wb", h1)) == 0) {
            fprintf(stderr, "samopen: Fail to open output BAM file %s\n", filename.c_str());
            exit(160);
        }
    }

    
    if (outFragTailBam) {
        string outfile=outputDirectory+"/"+filename+".fragtail.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:FragmentTail\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["FT"]=BamHeaderContainer(bam_header,q); 
        bam_header_t* h1=outputBam["FT"].header();
        if ((fpFT = samopen(outfile.c_str(), "wb", h1)) == 0) {
            fprintf(stderr, "samopen: Fail to open output BAM file %s\n", filename.c_str());
            exit(161);
        }
    }
     
    if (outInterChromBam) {
        string outfile=outputDirectory+"/"+filename+".interchrom.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:InterChromPairs\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["IC"]=BamHeaderContainer(bam_header,q);   
        bam_header_t* h1=outputBam["IC"].header();
        if ((fpIC = samopen(outfile.c_str(), "wb", h1)) == 0) {
            fprintf(stderr, "samopen: Fail to open output BAM file %s\n", filename.c_str());
            exit(162);
        }
    }
    
    if (outUniqueMultipleBam) {
        string outfile=outputDirectory+"/"+filename+".uMult.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:uniqMultiplyMappedPairs\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["UM"]=BamHeaderContainer(bam_header,q); 
        bam_header_t* h1=outputBam["IUM"].header();
        if ((fpUM = samopen(outfile.c_str(), "wb", h1)) == 0) {
            fprintf(stderr, "samopen: Fail to open output BAM file %s\n", filename.c_str());
            exit(163);
        }
    }
    
    if (outUniquePartialBam) {        
        string outfile=outputDirectory+"/"+filename+".uPart.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:uniqPartiallyMappedPairs\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["UP"]=BamHeaderContainer(bam_header,q);  
        bam_header_t* h1=outputBam["UP"].header();
        if ((fpUP = samopen(outfile.c_str(), "wb", h1)) == 0) {
            fprintf(stderr, "samopen: Fail to open output BAM file %s\n", filename.c_str());
            exit(164);
        }
    }

    if (outUniqueUnmappedBam) {        
        string outfile=outputDirectory+"/"+filename+".uUnmapped.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:uniqUnMappedPairs\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["UZ"]=BamHeaderContainer(bam_header,q); 
        bam_header_t* h1=outputBam["UZ"].header();
        if ((fpUZ = samopen(outfile.c_str(), "wb", h1)) == 0) {
            fprintf(stderr, "samopen: Fail to open output BAM file %s\n", filename.c_str());
            exit(165);
        }

    }

    if (outReadPairPosBam) {        
        string outfile=outputDirectory+"/"+filename+".weirdpairs.bam";
        string sv=SpannerVersion;
        string q="@PG\tID:weirdpairs\tPN:SpannerX\tVN"+sv+"\tCL:"+Params.getCmdLine();
        outputBam["WP"]=BamHeaderContainer(bam_header,q); 
        bam_header_t* h1=outputBam["WP"].header();
        if ((fpWP = samopen(outfile.c_str(), "wb", h1)) == 0) {
            fprintf(stderr, "samopen: Fail to open output BAM file %s\n", filename.c_str());
            exit(165);
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
        // skip if neither end within region
        more=(bampair.BamEnd.size()>1);
        
        Npair++;
        if (Npair>=maxReads) break; 

        //
        if ( (dbg!=0)&&(elapsedtime()>float(dbg))) {
			time(&tprev);
			cout << " pairs:" << Npair << "\toutput:" << Nout;
			cout << "\tchr:" << bampair.BamEnd[0].b.core.tid+1;
            cout << "\tpos:" << bampair.BamEnd[0].b.core.pos;
            cout << endl;			
		}  
        
        if (!more) continue; 
        if (region.limit) {
            bool overlap = false;
            for (int e=0; e<=1; e++) { 
                int a1=bampair.BamEnd[e].b.core.tid;
                int p1=bampair.BamEnd[e].b.core.pos;
                int p2=p1+bampair.BamEnd[e].len;
                overlap=region.overlap(a1,p1,p2);
                if (overlap) break; 
            }        
            if (!overlap) continue;
        }
    
        
        bampair.Illuminize(IlluminizeBam);  
        bampair.calcFragmentLengths();
        more=(bampair.BamEnd[1].packeddata.size()>1);
        //if (bampair.BamEnd[0].b.core.tid==bampair.BamEnd[1].b.core.tid) 
        //    cout<< bampair << endl;
        
        bool bothmap = ((bampair.BamEnd[0].b.core.flag&BAM_FUNMAP)==0)&&((bampair.BamEnd[0].b.core.flag&BAM_FMUNMAP)==0);
            
        
        if (outAllPairsBam) {
            Nout++;
            int s1=samwrite(fpAP, &(bampair.BamEnd[0].b));
            int s2=samwrite(fpAP, &(bampair.BamEnd[1].b));
            if ((s1*s2)>0) {
                continue;
            } else {
                cerr << "bad write to pairs.bam" << endl;
                exit(150);
            }
        }

        
        if (outReadPairPosBam) {
            int ichr1=bampair.BamEnd[0].b.core.tid+1;
            int istd1=bampair.BamEnd[0].sense=='+'? 0: 1;
            int ista1=bampair.BamEnd[0].b.core.pos;
            int iq1=bampair.BamEnd[0].q;
            int ichr2=bampair.BamEnd[1].b.core.tid+1;
            int istd2=bampair.BamEnd[1].sense=='+'? 0: 1;
            int ista2=bampair.BamEnd[1].b.core.pos;
            int iq2=bampair.BamEnd[1].q;

            FragmentPosObj  fp1(0,ichr1,istd1,ista1,0,ichr2,istd2,ista2,0,iq1, iq2,0);
            if (FragPos.find(fp1)) {
                Nout++;
                int s1=samwrite(fpWP, &(bampair.BamEnd[0].b));
                int s2=samwrite(fpWP, &(bampair.BamEnd[1].b));
                if ((s1*s2)>0) {
                    continue;
                } else {
                    cerr << "bad write to weirdpairs.bam" << endl;
                    exit(156);
                }
            }
        }        
        bool ok[2];
        for (int e=0; e<2; e++) {
            uint8_t*  bq=bam1_qual(&(bampair.BamEnd[e].b));
            int LR=bampair.BamEnd[0].b.core.l_qseq;
            double bok=0;
            for (int ib=0; ib<LR; ib++) {
                if (bq[ib]>SplitBaseQmin) {
                    bok++;
                }
            }
            ok[e]=(bok>LRmin);
        }
        
        if (! (ok[0]&ok[1]) )
            continue;
        
        if ( (outFragTailBam) & ((bampair.BamEnd[0].q>=Qmin)|(bampair.BamEnd[1].q>=Qmin)) ) {            
            bool FT=(bampair.FragmentLength>LFhigh)|((bampair.FragmentLength<LFlow)&(bampair.FragmentLength>INT_MIN))&bothmap;
            if (FT && (fpFT!=0)) {
                Nout++;
                int s1=samwrite(fpFT, &(bampair.BamEnd[0].b));
                int s2=samwrite(fpFT, &(bampair.BamEnd[1].b));
                //if (outputBam["FT"].write(&(bampair.BamEnd[0].b),&(bampair.BamEnd[1].b))) {
                if ((s1*s2)>0) {
                    continue;
                } else {
                    cerr << "bad write to fragtail.bam" << endl;
                    exit(151);
                }
            }
        }
        
        if ((outInterChromBam) & ((bampair.BamEnd[0].q>=Qmin)&(bampair.BamEnd[1].q>=Qmin))) { 
            bool IC=(bampair.BamEnd[0].b.core.tid!=bampair.BamEnd[1].b.core.tid)&&bothmap;
            if (IC && (fpIC!=0)) {
                Nout++;
                int s1=samwrite(fpIC, &(bampair.BamEnd[0].b));
                int s2=samwrite(fpIC, &(bampair.BamEnd[1].b));
                if ((s1*s2)>0) {
                    continue;
                } else {
                    cerr << "bad write to interchrom.bam" << endl;
                    exit(152);
                }
            }
        }
        if ((outUniqueMultipleBam) & ((bampair.BamEnd[0].q>=Qmin)|(bampair.BamEnd[1].q>=Qmin))){
            int im=bampair.BamEnd[0].nmap>1? 0: 1;
            int iu=bampair.BamEnd[0].q>=Qmin? 0: 1;
            bool UM=(bampair.BamEnd[iu].nmap>1)&&(iu!=im)&&bothmap;            
            if (UM && (fpUM!=0)) {
                Nout++;
                int s1=samwrite(fpUM, &(bampair.BamEnd[0].b));
                int s2=samwrite(fpUM, &(bampair.BamEnd[1].b));
                if ((s1*s2)>0) {
                    continue;
                } else {
                    cerr << "bad write to uMult.bam" << endl;
                    exit(153);
                }
            }
        }
        if ( (outUniquePartialBam) && ((bampair.BamEnd[0].q>=Qmin)|(bampair.BamEnd[1].q>=Qmin)) && bothmap) {            
            int c0=bampair.BamEnd[0].clip[0]+bampair.BamEnd[0].clip[1];
            int LR=bampair.BamEnd[0].b.core.l_qseq;
            bool split0=((LR-c0)>SplitBracketMin)&(c0>SplitBracketMin);
            int ib0=0;
            if ((split0)&(bampair.BamEnd[0].clip[0]>SplitBracketMin)) {
                ib0=bampair.BamEnd[0].clip[0];
            } else if ((split0)&(bampair.BamEnd[0].clip[1]>SplitBracketMin) ) {
                ib0=LR-bampair.BamEnd[0].clip[1];
            }
            split0=split0&(ib0>0);
            if (split0) {
                uint8_t*  bq=bam1_qual(&(bampair.BamEnd[0].b));
                for (int ib=(ib0-SplitBracketMin); ib<(ib0+SplitBracketMin); ib++) {
                    if (bq[ib]<SplitBaseQmin) {
                        split0=false;
                        break;
                    }
                }
            }
            
            int c1=bampair.BamEnd[1].clip[0]+bampair.BamEnd[1].clip[1];
            LR=bampair.BamEnd[1].b.core.l_qseq;
            bool split1=((LR-c0)>SplitBracketMin)&(c1>SplitBracketMin);;
            int ib1=0;
            if ((split1)&(bampair.BamEnd[1].clip[0]>SplitBracketMin)) {
                ib1=bampair.BamEnd[1].clip[0];
            } else if ((split1)&(bampair.BamEnd[1].clip[1]>SplitBracketMin) ) {
                ib1=LR-bampair.BamEnd[1].clip[1];
            }
            split1=split1&(ib1>0);
            if (split1) {
                uint8_t*  bq=bam1_qual(&(bampair.BamEnd[1].b));
                for (int ib=(ib1-SplitBracketMin); ib<(ib1+SplitBracketMin); ib++) {
                    if (bq[ib]<SplitBaseQmin) {
                        split1=false;
                        break;
                    }
                }
            }
            bool UP=(split0|split1)&((c1+c0)>minLR);
            if (UP && (fpUP!=0)) {
                Nout++;
                int s1=samwrite(fpUP, &(bampair.BamEnd[0].b));
                int s2=samwrite(fpUP, &(bampair.BamEnd[1].b));
                if ((s1*s2)>0) {
                    continue;
                } else {
                    cerr << "bad write to uPart.bam" << endl;
                    exit(154);
                }
            }
        }
        if ( (outUniqueUnmappedBam) & ((bampair.BamEnd[0].q>=Qmin)|(bampair.BamEnd[1].q>=Qmin)) ) {
            bool z0=((bampair.BamEnd[0].b.core.flag&BAM_FUNMAP)>0);
            bool z1=((bampair.BamEnd[1].b.core.flag&BAM_FUNMAP)>0);  
            
            
            uint8_t*  bq=bam1_qual(&(bampair.BamEnd[0].b));
            for (int nb,ib=0; ib<bampair.BamEnd[0].b.core.l_qseq; ib++) {
                if (bq[ib]<SplitBaseQmin) {
                    nb++;
                }
            }

            
            bool UZ=(z0|z1)&(!(z1&z0));
            if (UZ && (fpUZ!=0)) {
                Nout++;
                int s1=samwrite(fpUZ, &(bampair.BamEnd[0].b));
                int s2=samwrite(fpUZ, &(bampair.BamEnd[1].b));
                if ((s1*s2)>0) {
                    continue;
                } else {
                    cerr << "bad write to uUnmapped.bam" << endl;
                    exit(155);
                }
            }
        }
        
        
        //cout<< bampair.Orientation << "\t"<< bampair.FragmentLength << "\t" <<bampair.BamEnd[1].b.core.pos << endl;
        
    }
    
    if (outReadPairPosBam) {
        samclose(fpWP);
    } else {        
        if (outAllPairsBam) {
            samclose(fpAP);
        } else {       
            samclose(fpFT);
            samclose(fpIC);    
            samclose(fpUP);    
            samclose(fpUM);
            samclose(fpUZ);
        }
    }
    
    /*
     for (ioutputBam=outputBam.begin(); ioutputBam!=outputBam.end(); ioutputBam++) {
        (*ioutputBam).second.close();
    }
     
    if (FragmentTailPercent>0) 
        outputBam["FT"].close();
    */
    
    samclose(Bam.fp);
    
    
}

float BamX::elapsedtime() {
    time_t t;
    time(&t);
    float et = difftime (t,tprev);
    //float et=float(t-tprev)*1.0/float(CLOCKS_PER_SEC);
    return et;
}


