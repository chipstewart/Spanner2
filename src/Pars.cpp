//
//  Pars.cpp
//  Spanner2scan
//
//  Created by Chip Stewart on 6/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "Pars.h"
#include "SpannerVersion.h"

using std::string;
using namespace std; 

pars::pars()
{
    cmdOpt="mP:S:f:x:c:d:r:q:p:l:t:h";
    program="Spanner";
    inputoutputcmdline=" <input bam file>  >  < output file> ";
}

pars::pars(int argc,  char * argv[], string & cmdOpt1, string & inout1 )
{
    cmdOpt=cmdOpt1;
    program=argv[0];
    inputoutputcmdline=inout1;
    
    doubles.clear();
    ints.clear();
    strings.clear();
    histos.clear();
    
    // defaults
    ints["FragmentLengthWindow"]=10000;     
    ints["ReadPairSenseConfig"]=-1;
    ints["MaxReads"]=INT_MAX;
    ints["MinReadLength"]=30;   
    ints["Masking"]=0;        
    ints["Qmin"]=20;
    ints["Dbg"]=0;
    ints["Illuminize"]=1;
    
    doubles["FragmentTailPercent"]=0.1;
    doubles["FractionMaxMisMatches"]=double(100.0);      
    
    strings["StatFile"]="";
    strings["ReferenceFastaFile"]="";
    strings["ChromRegion"]="";
 
    histos["LF"]=hist();
    histos["LR"]=hist();
    

    getCommandLineParameters(argc, argv);
    
} 

// constructor with default parameters
void pars::getCommandLineParameters(int argc, char * argv[])
{
    int c;
    
    /*
     if (argc == 1) {
     fprintf(stderr, "Usage: SpannerScan [options] <in.bam> \n");
     fprintf(stderr, "SpannerScan Version %s \n",SpannerScanVersion);
     exit(100);
     }
     */
    
    bool help=(argc==1);
    
    // master list of all Spanner options
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"masking",                no_argument,       0, 'm'},
        {"illuminze",              no_argument,       0, 'I'},
        {"fragmenttailpercent",    optional_argument, 0, 'P'},
        {"statfile",               optional_argument, 0, 'S'},
        {"outputdirectory",        optional_argument, 0, 'o'},
        {"fragmentlengthwindow",   optional_argument, 0, 'f'},
        {"maxreads",               optional_argument, 0, 'x'},
        {"chromregion",            optional_argument, 0, 'c'},
        {"debug",                  optional_argument, 0, 'd'},
        {"referencefasta",         optional_argument, 0, 'r'},
        {"qmin",                   optional_argument, 0, 'q'},
        {"maxmismatchpercent",     optional_argument, 0, 'p'},
        {"minreadlength",          optional_argument, 0, 'l'},
        {"technology_matemode",    optional_argument, 0, 't'},
        {"help",                   no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    int xi;
    double xd;
    string xs,s;
    
    cmdline="";
    for (int i=0; i<argc; i++) {
        cmdline = cmdline + " " + argv[i];
    }
    
    while (1)
    {
        
        /* getopt_long stores the option index here. */
        int option_index = 0;
        
        
        c = getopt_long (argc, argv, cmdOpt.c_str(),
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
        
        size_t found=cmdOpt.find(c);
        if (found==string::npos)
            break;
        
        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                /*
                 fprintf(stderr,"option %s", long_options[option_index].name);
                 if (optarg)
                 fprintf(stderr," with arg %s\n", optarg);
                 */
                break;
                
            case 'm':
                //fprintf(stderr,"option -m\n");
                setInt("Masking",1);
                break;
                
            case 'I':
                //fprintf(stderr,"option -m\n");
                setInt("Illuminize",1);
                break;
                
            case 'P':
                //fprintf(stderr,"option -P with value `%s'\n", optarg);
                xs=optarg;
                xd=string2Int(xs);
                setDouble("FragmentTailPercent",xd);
                break;
                
            case 'S':
                //fprintf(stderr,"option -S with value `%s'\n", optarg);
                xs=optarg;
                setString("StatFile",xs);
                break;

            case 'o':
                //fprintf(stderr,"option -S with value `%s'\n", optarg);
                xs=optarg;
                setString("OutputDirectory",xs);
                break;
                                
            case 'f':
                //fprintf(stderr,"option -w with value `%s'\n", optarg);
                xs=optarg;
                xi=string2Int(xs);
                setInt("FragmentLengthWindow",xi);
                break;
                
            case 'x':
                //fprintf(stderr,"option -x with value `%s'\n", optarg);
                xs=optarg;
                xi=string2Int(xs);
                setInt("MaxReads",xi);
                break;
                
            case 'c':
                //fprintf(stderr,"option -c with value `%s'\n", optarg);
                xs=optarg;
                setString("ChromRegion",xs);
                break;
                
            case 'd':
                //fprintf(stderr,"option -d with value `%s'\n", optarg);
                xs=optarg;
                xi=string2Int(xs);
                setInt("Dbg",xi);
                break;
                
            case 'r':
                //fprintf(stderr,"option -r with value `%s'\n", optarg);
                xs=optarg;
                setString("ReferenceFastaFile",xs);
                break;
                
            case 'q':
                //fprintf(stderr,"option -q with value `%s'\n", optarg);
                xs=optarg;
                xi=string2Int(xs);
                setInt("Qmin",xi);
                break;
                
                
            case 'p':
                //fprintf(stderr,"option -p with value `%s'\n", optarg);
                xs=optarg;
                xd =string2Double(xs);
                setDouble("FractionMaxMisMatches",xd);
                break;
                
                
            case 'l':
                //fprintf(stderr,"option -l with value `%s'\n", optarg);
                xs=optarg;
                xi=string2Int(xs);
                setInt("MinReadLength",xi);
                break;
                
            case 't':
                //fprintf(stderr,"option -m with value `%s'\n", optarg);
                xs=optarg;
                xi=string2Int(xs);
                setInt("ReadPairSenseConfig",xi);
                break;
                
            case 'h':
                help=true;
                break;
                
            case '?':
                /* getopt_long already printed an error message. */
                break;
                
            default:
                abort ();
        }
    }
    
    /* Print any remaining command line arguments (not options). */
    /*
     if (optind < argc)
     {
     fprintf(stderr,"Input Bam : ");
     while (optind < argc)
     fprintf(stderr,"%s ", argv[optind++]);
     fprintf(stderr,"\n");
     }
     */
    
    // input file
    input=argv[argc-1];
    
    if (help) {
        fprintf(stderr,"SpannerX -options <input bam file>  >  stdout  \n");
        fprintf(stderr,"\tVersion: %s \n",SpannerVersion);
        fprintf(stderr,"\toptions:\n");
        
        int o=0;
        while (long_options[o].val>0) {
            
            /*
             if (long_options[o].has_arg==no_argument) {
             fprintf(stderr,"\t-%c: --%s : flag \n",long_options[o].val,long_options[o].name);
             
             } else if (long_options[o].has_arg==optional_argument) {
             */   
            {
                char v=char(long_options[o].val);
                
                size_t found=cmdOpt.find(v);
                if (found==string::npos)
                    break;

                
                switch (v) {
                        
                    case 'I':
                        fprintf(stderr,"\t-%c: --%s  :\t Illuminze pair orientation in output bam\n",long_options[o].val,long_options[o].name);
                        break;
                    
                    case 'S':
                        fprintf(stderr,"\t-%c: --%s <statfile> :\t library stat file\n",long_options[o].val,long_options[o].name);
                        break;
                        
                    case 'P':
                        fprintf(stderr,"\t-%c: --%s <FragmentTailPercent> :\t p-value on tails of fragment length \n",long_options[o].val,long_options[o].name);
                        break;
                        
                        
                    case 'o':
                        fprintf(stderr,"\t-%c: --%s <output area> :\t output file directory\n",long_options[o].val,long_options[o].name);
                        break;

                    case 'l':
                        xi=getInt("MinReadLength");
                        fprintf(stderr,"\t-%c: --%s <%d> :\t\t\t min allowed read length \n",long_options[o].val,long_options[o].name,xi);
                        break;
                        
                    case 'p':
                        xd=getDouble("FractionMaxMisMatches");
                        fprintf(stderr,"\t-%c: --%s <%.1f> :\t\t max allowed mismatch percentage \n",long_options[o].val,long_options[o].name,xd);
                        break;
                        
                    case 'q':
                        xi=getInt("Qmin");
                        fprintf(stderr,"\t-%c: --%s <%d> :\t\t\t\t min allowed mapQ \n",long_options[o].val,long_options[o].name,xi);
                        break;
                        
                    case 'r':
                        fprintf(stderr,"\t-%c: --%s <reference fasta> :\t used for masking + gc\n",long_options[o].val,long_options[o].name);
                        break;
                        
                    case 'd':
                        fprintf(stderr,"\t-%c: --%s  <0> :\t\t\t\t debug level (0=none)\n",long_options[o].val,long_options[o].name);
                        break;
                        
                    case 'c':
                        fprintf(stderr,"\t-%c: --%s <chr:start-end> :\t\t region to process \n",long_options[o].val,long_options[o].name);
                        break;
                        
                    case 'x':
                        fprintf(stderr,"\t-%c: --%s <%d> :\t\t\t max reads to process \n",long_options[o].val,long_options[o].name,INT_MAX);
                        break;
                        
                    case 'f':
                        xi=getInt("FragmentLengthWindow");
                        fprintf(stderr,"\t-%c: --%s <%d> :\t\t frag length histogram range\n",long_options[o].val,long_options[o].name,xi);
                        break;
                        
                    case 'm':
                        fprintf(stderr,"\t-%c: --%s :\t\t\t\t\t impose masking from Reference Fasta\n",long_options[o].val,long_options[o].name);
                        break;
                        
                    case 't':
                        xi=getInt("ReadPairSenseConfig");
                        fprintf(stderr,"\t-%c: --%s : \t\t\t <overide bam ReadGroup PL tag> \n",long_options[o].val,long_options[o].name);
                        fprintf(stderr,"\t\t0: FR (IL short); 1: R1R2/F2F1 (454);  2: F1F2/R2R1 (SOLiD); 3: RF (IL long)\n");
                        
                        break;
                        
                    case 'h':
                        fprintf(stderr,"\t-%c: --%s :\t\t\t\t\t this message \n",long_options[o].val,long_options[o].name);
                        break;
                        
                    case '?':
                        /* getopt_long already printed an error message. */
                        break;
                        
                    default:
                        exit(0);
                        
                }
            }            
            o++;
        }
        exit(0);
    }
}


pars& pars::operator=(const pars &rhs)
{
    input=rhs.input;
    output=rhs.output;
    ints=rhs.ints;
    doubles=rhs.doubles;
    strings=rhs.strings;
    histos=rhs.histos;
    program=rhs.program;
    inputoutputcmdline=rhs.inputoutputcmdline;
    cmdOpt=rhs.cmdOpt; 
    cmdline=rhs.cmdline; 
    
    return *this;
}

int pars::getInt(const string & s) {
    if(ints.count(s)>0)
        return ints[s];
    else
        return INT_MIN;    
}
bool pars::setInt(const string & s, const int xi) {
    bool b = ints.count(s)>0;
    ints[s]=xi;
    return b;
}

double  pars::getDouble(const string & s) {
    if(doubles.count(s)>0)
        return doubles[s];
    else
        return 1e-100;    
}
bool pars::setDouble(const string & s, const double  xd) {
    bool b = ints.count(s)>0;
    doubles[s]=xd;
    return b;
}

string  pars::getString(const string & s) {
    if(strings.count(s)>0)
        return strings[s];
    else
        return "";    
}
bool pars::setString(const string & s, const string & xs) {
    bool b = ints.count(s)>0;
    strings[s]=xs;
    return b;
}

hist  pars::getHist(const string & s) {
    if(histos.count(s)>0)
        return histos[s];
    else
        return hist();    
}

bool pars::setHist(const string & s, const hist & xh) {
    bool b = ints.count(s)>0;
    histos[s]=xh;
    return b;
}


string  pars::getInput() {
    return input;
}

string  pars::getOutput() {
    return output;
}


string pars::getPartype(const string & s) {
    string o="";
    if(ints.count(s)>0)
        o="Int";
    if(doubles.count(s)>0)
        o="Double";    
    if(strings.count(s)>0)
        o="String";
    if(histos.count(s)>0)
        o="Histogram";
    return o;
}

string  pars::getCmdOpt() {
    return cmdOpt;
}

string  pars::getProgram() {
    return program;
}

string  pars::getInputOutput() {
    return inputoutputcmdline;    
}

string  pars::getCmdLine() {
    return cmdline;
}