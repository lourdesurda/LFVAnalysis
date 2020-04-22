//Macro to read the sample cards 

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TString.h>
#include <stdio.h>
#include <iomanip>

void Readingtxtfiles(const char* name) 
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
	std::cout << std::setprecision(20);
	if (SampleID.Contains("Data")) 
	std::cout << SampleID <<"\t" << WhatIsIt << "\t" << Path << "\t" << luminosity << "\t" << Skimmed_nano_trees << "\t" << NumberOfEvents << "\t" << TotalEventWeight << std::endl;
	else 
	std::cout << SampleID <<"\t" << WhatIsIt << "\t" << Path << "\t" << xsection << "\t" << Skimmed_nano_trees << "\t" << NumberOfEvents << "\t" << TotalEventWeight << std::endl;


	return 0;

}
