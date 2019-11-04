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
	
