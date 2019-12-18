#include "SamplesClassLU.h"

// Initialize static member of class SAMPLES
double SAMPLES::TotalAnalysisLumi = 0;

//meter metodos aqui


	SAMPLES::SAMPLES (const TString& id,  const TString& dir, const TString& file, int nfiles, double xsec, double luminosity, const TString Flag, long long TotalNumberOfEvents)
	{

		Definitions(id, dir, file, nfiles, xsec, luminosity, Flag, TotalNumberOfEvents);

	}

	void SAMPLES::Definitions(const TString& id,  const TString& dir, const TString& file, int nfiles, double xsec, double luminosity, const TString& Flag, long long TotalNumberOfEvents)
	{
		SampleId = id;
		SampleFile = dir+file;
		SampleNFiles = nfiles;
		SampleXSection = xsec;
		SampleFlag = Flag;

		if(isMCSample())
		{
			double EquivalentLuminosity = TotalNumberOfEvents/xsec;
			SampleWeight = SAMPLES::TotalAnalysisLumi/EquivalentLuminosity;
		}
	 
		else if(Flag == "isData")
		{
			IndividualLuminosity = luminosity;
			TotalAnalysisLumi += luminosity;
		}
		
		//SampleNEvents=getNumberOfEntries(filename);
	}

//LEER FICHEROS //PASAR LOS NOMBRES DE LOS FICHEROS.TXT
	FICH::FICH(const char* name, const TString& SampleID, const TString& WhatIsIt, const TString& Path, double xsec, double luminosity, const int Skimmed_nano_trees, long  NumberOfEvents, double TotalEventWeight)
	{
		Readingtxtfiles(name);
	}

	void FICH::Readingtxtfiles(const char* name) 
	{
		string line;
		TString line_tstring;
		ifstream inputcard;
		inputcard.open(name);
		string items[7];

		TString SampleID;
		TString WhatIsIt;
		TString Path;
		double luminosity;
		double xsection;
		int Skimmed_nano_trees;
		long NumberOfEvents;
		double TotalEventWeight;

		if(inputcard.is_open())
		{
			unsigned int i = 0;
			while( getline (inputcard, line) )
			{
				line_tstring = (TString)line;
				if(line_tstring.Contains("#")) continue;
				items[i] = line;
				//std::cout << line << "\n";
				i = i+1;
			}	
			inputcard.close();	
		}
		else std::cout << "Unable to open file";

		SampleID = items[0];
		WhatIsIt = items[1];
		Path 	 = items[2];

		if(SampleID.Contains("Data")) luminosity = std::stod(items[3]);
		else xsection = std::stod(items[3]);

		Skimmed_nano_trees = std::stod(items[4]);
		NumberOfEvents = std::stod(items[5]);
		TotalEventWeight = std::stod(items[6]);

		//printf("%f\n", TotalEventWeight);
		/*std::cout << std::setprecision(20);
		if (SampleID.Contains("Data")) 
		std::cout << SampleID <<"\t" << WhatIsIt << "\t" << Path << "\t" << luminosity << "\t" << Skimmed_nano_trees << "\t" << NumberOfEvents << "\t" << TotalEventWeight << std::endl;
		else 
		std::cout << SampleID <<"\t" << WhatIsIt << "\t" << Path << "\t" << xsection << "\t" << Skimmed_nano_trees << "\t" << NumberOfEvents << "\t" << TotalEventWeight << std::endl;
		return 0;*/

	}
