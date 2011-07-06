//
//  SpannerBam.cpp
//  SpannerBam
//
//  Created by Chip Stewart on 6/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "ParamBam.h"
#include "BamScan.h"
#include "SpannerScanVersion.h"

//int main (int argc, const char * argv[])
int main (int argc, char * argv[])
{        
    ParamBam Params(argc,argv);
    
    int dbg = Params.getDbg();
    if (dbg>0) {
        for(int i = 0; i < argc; i++) {
            cerr <<  argv[i] << " ";
        }
        cerr << endl;
    }
    
    BamScan Scan(Params);
    
    return 0;
}

