//
//  SpannerScan.cpp
//  SpannerScan
//
//  Created by Donald Stewart on 5/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Pars.h"
#include "BamX.h"
#include "SpannerVersion.h"

//int main (int argc, const char * argv[])
int main (int argc, char * argv[])
{        
    // options for BamX
    string cmdOpt="mo:P:S:x:c:d:r:q:p:l:t:h:I";
    string inputoutput ="SpannerX -o <outputdir> -options <input bam file>  >  stdout ";
    pars Params(argc,argv,cmdOpt,inputoutput);

    int dbg = Params.getInt("Dbg");
    if (dbg>0) {
     for(int i = 0; i < argc; i++) {
        cerr <<  argv[i] << " ";
     }
     cerr << endl;
    }
    
    BamX Extract(Params);
    
    return 0;
}
