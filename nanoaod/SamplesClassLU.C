#include "SamplesClassLU.h"
#include "TString.h"
// Initialize static member of class SAMPLES
double SAMPLES::TotalAnalysisLumi = 0;

//meter metodos aqui


	SAMPLES::SAMPLES (const TString& id,  const TString& dir, const TString& file, int nfiles, double xsec, double luminosity, const TString Flag, long long TotalNumberOfEvents)
	{
		SampleWeight = 1.;
		IndividualLuminosity = 0.0;
		EquivalentLuminosity = 0.0;
		NeventsSumFiles = 0;


		Definitions(id, dir+file, nfiles, xsec, luminosity, Flag, TotalNumberOfEvents);
	}

	SAMPLES::SAMPLES(const char* samplefichtxtfile)
	{
		SampleWeight = 1.;
		IndividualLuminosity = 0.0;
		EquivalentLuminosity = 0.0;
		NeventsSumFiles = 0;

		Readingtxtfiles(samplefichtxtfile);
	}

	void SAMPLES::Definitions(const TString& id,  const TString& pathfile, int nfiles, double xsec, double luminosity, const TString& Flag, long long TotalNumberOfEvents)
	{
		SampleId = id;
		SampleFile = pathfile;
		SampleNFiles = nfiles;
		SampleXSection = xsec;
		SampleFlag = Flag;
		SampleGenWeight = TotalNumberOfEvents;

		if(isMCSample() && !id.Contains("WJets"))
		{
			EquivalentLuminosity = TotalNumberOfEvents/xsec; //pb-1
			SampleWeight = SAMPLES::TotalAnalysisLumi/EquivalentLuminosity;

			std::cout << "id: " << id << "\t xsec: " << xsec << "\t TotalNumberOfEvents: " << TotalNumberOfEvents << "\t EquivalentLuminosity: " << EquivalentLuminosity <<"\t TotalAnalysisLumi: " << TotalAnalysisLumi << "\t SampleWeight: " << SampleWeight << std::endl;
		}
		//else SampleWeight = 1.;
	 
		if(Flag == "isData")
		{
			//SampleWeight = 1.;
			IndividualLuminosity = luminosity;
			TotalAnalysisLumi += luminosity;
		}
		
		//SampleNEvents=getNumberOfEntries(filename);
	}

	

//LEER FICHEROS //PASAR LOS NOMBRES DE LOS FICHEROS.TXT
	/*FICH::FICH(const char* name, const TString& SampleID, const TString& WhatIsIt, const TString& Path, double xsec, double luminosity, const int Skimmed_nano_trees, long  NumberOfEvents, double TotalEventWeight)
	{
		Readingtxtfiles(name);
	}*/

	void SAMPLES::Readingtxtfiles(const char* samplefichtxtfile) 
	{
		//string line;
		//TString line_tstring;
		ifstream inputcard(samplefichtxtfile);//input card is the name of the txt file with the information of the sample
		//inputcard.open(samplefichtxtfile);
		string items[7];

		if(inputcard.is_open())
		{
			unsigned int i = 0;
			string line;
			while( getline (inputcard, line ) )
			{
				//line_tstring = (TString)line;
				if(line.find('#')!=string::npos) continue;
				items[i] = line;
				//std::cout << line <<" " << i<< "\n";
				i = i+1;
			}	
			inputcard.close();	
		}
		else std::cout << "Unable to open file";

		TString SampleID = items[0];
		TString SampleFlag = items[1];
		TString Path 	 = items[2];

		double luminosity = 0.0;
		double xsection = 0.0;

		if(SampleID.Contains("Data")) luminosity = std::stod(items[3]);
		else xsection = std::stod(items[3]);

		int Skimmed_nano_trees = std::stod(items[4]);
		long NumberOfEvents = std::stod(items[5]);
		double TotalEventWeight = std::stod(items[6]);

		Definitions(SampleID, Path, Skimmed_nano_trees, xsection, luminosity, SampleFlag, TotalEventWeight);

		//printf("%f\n", TotalEventWeight);
		/*std::cout << std::setprecision(20);
		if (SampleID.Contains("Data")) 
		std::cout << SampleID <<"\t" << WhatIsIt << "\t" << Path << "\t" << luminosity << "\t" << Skimmed_nano_trees << "\t" << NumberOfEvents << "\t" << TotalEventWeight << std::endl;
		else 
		std::cout << SampleID <<"\t" << WhatIsIt << "\t" << Path << "\t" << xsection << "\t" << Skimmed_nano_trees << "\t" << NumberOfEvents << "\t" << TotalEventWeight << std::endl;
		return 0;*/

	}
