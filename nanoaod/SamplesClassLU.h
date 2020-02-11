//I am going to define a class for the Samples for my analisis. 
//I am going to consider SAMPLES class. There may be many samples with different names and characteristics but all of them will share some common properties like all of the ID, number of events, xsection, number of files, etc. 
//& means it's a reference to a pointer. I allow the function to change the variable

#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <iomanip>

using namespace std;

struct SAMPLES
{
	static double TotalAnalysisLumi;

	//Data members
	TString  SampleId;
	TString  SampleFile;
	TString	 SampleFlag;
	int	 SampleNFiles; 
	//int	 SampleNEvents;
	double	 SampleXSection;
	double 	 SampleWeight = 1.;
	double   SampleGenWeight;
	double   IndividualLuminosity = 0.0;
	double   EquivalentLuminosity = 0.0;

	double 	 Inverted_Luminosity;

	long long NeventsSumFiles;

	// Constructor
	explicit SAMPLES (const TString& id,  const TString& dir, const TString& file, int nfiles, double xsec, double luminosity, const TString Flag, long long TotalNumberOfEvents);
	//OTRO CONSTRUCTOR PARA LOS FICHEROS
	explicit SAMPLES(const char* samplefichtxtfile);
	
	// Destructor
	~SAMPLES () {}	

	//Bolean variables to classify samples depending on if they are: MCSignal, MCBckg, DATA.
	bool isSignalSample () {return ((SampleFlag=="isMCSignal")?true:false);} 

	bool isMCBckgSample () {return ((SampleFlag=="isMCBckg")?true:false);}

	bool isMCSample     () {return ((SampleFlag=="isMCBckg" || SampleFlag=="isMCSignal")?true:false);}

	bool isData     () {return ((SampleFlag=="isData")?true:false);}

	//Member Functions()
	void Definitions(const TString& id,  const TString& pathfile, int nfiles, double xsec, double luminosity, const TString& Flag, long long TotalNumberOfEvents);

	const TString &GetSampleId () {return SampleId;}
	const int &GetNumberOfFilesInSample(){return SampleNFiles;}
	
	const double &GetSampleXSection () {return SampleXSection;}
	const double &GetSampleGenWeight() {return SampleGenWeight;}
	const double &GetInvertedLuminosity ()
	{
		Inverted_Luminosity = 1.0/EquivalentLuminosity;
		return Inverted_Luminosity;
	}

	void Readingtxtfiles(const char* name);

};

/*struct FICH
{
	explicit FICH(const char* name, const TString& SampleID, const TString& WhatIsIt, const TString& Path, double xsec, double luminosity, const int Skimmed_nano_trees, long  NumberOfEvents, double TotalEventWeight);

	~FICH() {}

	void Readingtxtfiles(const char* name);
};*/


