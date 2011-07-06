//
//  Pars.h
//  Spanner2scan
//
//  Created by Chip Stewart on 6/24/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//
#ifndef PARS_H
#define PARS_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>

// uses

#include "UtilityFunctions.h"
#include "Histo.h"
#include "SpannerVersion.h"

class pars
{
	friend ostream &operator<<(ostream &, const pars &);
    
public:
	
    pars();	// default constructor
    pars(int argc,  char * argv[],string &, string &);	// optional constructor
    
	pars&operator=(const pars &rhs);
    
    void getCommandLineParameters(int argc, char * argv[]);
     
	int getInt(const string &);
	bool setInt(const string & , const int);

	double getDouble(const string &);
	bool setDouble(const string &, const double);

	string getString(const string &);
	bool setString(const string &, const string &);
    
	hist  getHist(const string &);
	bool setHist(const string &, const hist &);
    
    string getPartype(const string &);
	
    vector<string> partypes;

    string getInput();
	string getOutput();
    string getCmdOpt();    
    string getProgram();
    string getCmdLine();
    string getInputOutput();   
    
private:
    
    map<string, int, less<string> > ints;
    map<string, int, less<string> >::iterator iti;

    map<string, double, less<string> > doubles;
    map<string, double, less<string> >::iterator itd;
    
    map<string, string, less<string> > strings;
    map<string, string, less<string> >::iterator its;

    map<string, hist, less<string> > histos;
    map<string, hist, less<string> >::iterator ith;
    
    map<string, string, less<string> > partype;
    map<string, string, less<string> >::iterator itp;
    
    string input;
    string output;
    string program;
    string cmdline;
    string inputoutputcmdline;
 
    string cmdOpt;  // list of short options allowed for this instance
    
}; 

#endif