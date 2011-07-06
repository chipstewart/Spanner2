/*
 *  Histo.h
 *  Spanner
 *
 *  Created by Chip Stewart on 11/3/07.
 *  Copyright 2007 Boston College. All rights reserved.
 *
 */
#ifndef HISTO_H
#define HISTO_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <iterator>

using namespace std;


// structure for histogram 
class hist {
  friend ostream &operator<<(ostream &, const hist &);
  public:
    hist(const vector<double> & , const int , const double , const double); // constructor (do it all)
    hist(); // empty constructor 
    hist(string &); // load histogram from file 
    ~hist();
    hist&operator=(const hist &rhs);
    void Initialize(const int , const double , const double);  // constructor with bins only 
    void Fill(const vector<double> & , const int , const double , const double);  
    void Fill(const vector<int> & , const int , const double , const double);  
    void Fill(const vector<short> & , const int , const double , const double);
    void Fill1(const double);  
    void Fill1(const int);  
    void Fill1(const short);  
    void FillW(const double,const double );  
    void FillW(const int,const double);  
    void FillW(const short,const double);  
    void setTitle(const string &);
    void setXlabel(const string &);
    void setBinLabels(vector<string> &);
    void  Finalize();  
    double x2p(double); 
    double p2x(double); 
    double p2xTrim(double); 
    double x2pTrim(double); 
  	hist collapse( int);
    hist expand();
    int Nbin;
    double xlow;
    double xhigh;
    double dx;
    double Ntot;
    double Nin;
    double Nover;
    double Nunder;
    double mode;
    double median;
    vector<double> n;
    vector<double> xc;
    vector<double> c;
    double mean;
    double std;
    double sumx;
    double sumxx;
    string title;
    string xlabel;
    vector<string> binlabels;    
    bool normalize;
    double mode1;
    bool collapsed;
	bool expanded;
};

class hists {
	//friend ostream &operator<<(ostream &, const hists &);
public:
    hists();
    hists(string &); // constructor 
    hists(ifstream &); // constructor 
    map<string, hist, less<string> > h;
    vector<string> ReadGroupTag;	
};


class C_HistoGroups {
	friend ostream &operator<<(ostream &, const C_HistoGroups &);
public:
	C_HistoGroups(); // constructor 
	C_HistoGroups(string &);
    hists getHistos(string &);
 	vector<hists> Groups;	
    map<string,int,less<string> > ReadGroupIndex; 
};



class StatObj {
  friend ostream &operator<<(ostream &, const StatObj &);
  public:
    StatObj(const vector<double> & , const int , const double , const double ); // constructor 
    StatObj();      // constructor 
    StatObj(string &); // load stats from file 
    void Fill(const vector<double> & , const int , const double , const double );  
    void Fill(const vector<int> & , const int , const double , const double );  
    void Fill(const vector<short> & , const int , const double , const double );  
    void Initialize(const int , const double , const double );  
    void Fill1(const double);  
    void Fill1(const int);  
    void Fill1(const short);  
    void Finalize();  
    int N;
    double mean;
    double std;
    double sumx;
    double sumxx;
    double threshold(const double);
    hist h;
};

#endif


