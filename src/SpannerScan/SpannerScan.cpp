//
//  SpannerScan.cpp
//  SpannerScan
//
//  Created by Donald Stewart on 5/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Pars.h"
#include "BamScan.h"
#include "SpannerVersion.h"

//int main (int argc, const char * argv[])
int main (int argc, char * argv[])
{        
    string cmdOpt="mf:x:c:d:r:q:p:l:t:h";
    string inputoutput =" SpannerScan -options <input bam file>  >  <output stat file> ";
    pars Params(argc,argv,cmdOpt,inputoutput);
    //pars Params(argc,argv);

    int dbg = Params.getInt("Dbg");
    if (dbg>0) {
     for(int i = 0; i < argc; i++) {
        cerr <<  argv[i] << " ";
     }
     cerr << endl;
    }
    
    BamScan Scan(Params);
    
    return 0;
}

