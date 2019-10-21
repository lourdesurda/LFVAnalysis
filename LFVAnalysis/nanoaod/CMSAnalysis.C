//Code for LFV Analysis for 2018 CMS data
//Written by Juan Alcaraz and it had been modified by Lourdes Urda
//August 2019

#include "TROOT.h"
#include "TSystem.h"
#include "TRint.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include "TMath.h"
#include "TF1.h"
#include "TLine.h"
#include "tdrstyle.C"
#include "CMSAnalysis.h"

TTree* CMSAnalysis::_NANOTREE = NULL;
std::vector<TBranch*> CMSAnalysis::_BUFFERBRANCHES = {NULL};
std::vector<Char_t> CMSAnalysis::_BUFFER = {0};

CMSAnalysis::CMSAnalysis() 
{
	_file = 0;
      	_currentIndex = -1; 
      	_dataIndex = -1;

     	_BUFFERBRANCHES.clear();
     	_BUFFER.clear();
};

CMSAnalysis::~CMSAnalysis() 
{ 
      if (_file) _file->Close();
};

void CMSAnalysis::AddDataSample(const TString& id, const TString& file, double luminosity, int nfiles) 
{	
	_lumi = luminosity;
  	_sampleId.push_back(id);
  	_sampleFile.push_back(file);

 	TFile file_tmp(file,"READONLY");
  	TTree* tree_tmp = 0;
  	file_tmp.GetObject("Events",tree_tmp);
 	if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", file.Data());
  	_sampleNevents.push_back(tree_tmp->GetEntriesFast());
  	int nDump = 1; while (_sampleNevents.back()/nDump>20) nDump *= 10;
  	_nDump10.push_back(nDump);

  	_sampleXsection.push_back(_sampleNevents.back()/_lumi);
  	_sampleWeight.push_back(1.);
  	_signalFlag.push_back(false); // set to false for data

        _sampleNFiles.push_back(1);
  	_dataIndex = _sampleId.size()-1;

  	// Printout summary for data files
  	printf("DATA> %s - %s: luminosity: %f /pb, events %d\n", id.Data(), _sampleFile.back().Data(), _lumi, _sampleNevents.back());
}

long long int CMSAnalysis::NumberOfEntries(const TString& path, int nfiles)
{
	long long int suma=0;
	for(int itree=1; itree<=nfiles; itree++)
	{
		TFile *file =TFile::Open(path+Form("tree_%i.root", itree));
		if(file==nullptr) continue;
		TH1F* h = (TH1F*)file->Get("eventcount_Skim");
		suma =  suma + h->GetBinContent(1);
		file->Close();
	}
	return suma;
}


void CMSAnalysis::AddDataSampleFiles(const TString& id,  const TString& dir,const TString& file, double luminosity, int nfiles) {
        _lumi = luminosity;
        _totalLuminosity+=luminosity;
        _sampleId.push_back(id);
        _sampleFile.push_back(dir+file);

        TFile file_tmp(dir+file+"_1.root","READONLY");
        TTree* tree_tmp = 0;
        file_tmp.GetObject("Events",tree_tmp);

        TString filenameexample=dir+file+"_1.root";
        if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", filenameexample.Data());

        _sampleNevents.push_back(tree_tmp->GetEntriesFast());
        _sampleNeventsTree.push_back(tree_tmp->GetEntriesFast());
        _sampleNeventsSumFiles.push_back(0);

        int nDump = 1; while (_sampleNevents.back()/nDump>20) nDump *= 10;
        _nDump10.push_back(nDump);

        _sampleXsection.push_back(_sampleNevents.back()/_lumi);
        _sampleWeight.push_back(1.);
        _signalFlag.push_back(false); // set to false for data
        _dataIndex = _sampleId.size()-1;

        _sampleNFiles.push_back(nfiles);


        // Printout summary for data files
        printf("DATA> %s - %s: luminosity: %f /pb (%f /pb), events %d\n", id.Data(), _sampleFile.back().Data(), _lumi, _totalLuminosity,_sampleNevents.back());
}

void CMSAnalysis::AddMCSampleFiles(const TString& id, const TString& dir, const TString& file, int maxevents, double xsec, int total_events_for_xsection, int nfiles) {
        if (_dataIndex<0) {
                printf(">>> Warning: you should call AddDataSample first, to define the luminosity!!!\n");
                printf(">>> SAMPLE NOT ADDED!\n");
                return;
        }
        _sampleId.push_back(id);
        _sampleFile.push_back(dir+file);

        _sampleNFiles.push_back(nfiles);

        TString filenameexample=dir+file+"_1.root";
        TFile file_tmp(filenameexample,"READONLY");
        TTree* tree_tmp = 0;
        file_tmp.GetObject("Events",tree_tmp);

        if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", filenameexample.Data());

        _sampleXsection.push_back(xsec);
        _signalFlag.push_back(false); // set to false to start with

        //        for (int file=1; file<nfiles+1; file++){
        //              TFile filecount(dir+file+"_"+".root";,"READONLY"); 
        //              TH1F *histoEvents=(TH1F*)file_tmp.Get("eventcount");
        //        if(histoEvents) std::cout<<histoEvents->GetBinContent(1);

        double equivalent_lumi = -1;
        if (total_events_for_xsection>0 && xsec>0) equivalent_lumi = total_events_for_xsection/xsec;
        else {std::cout<<"Please introduce a valid number for cross section and event count before cuts for sample"<<dir<<std::endl;}

        _sampleNevents.push_back(total_events_for_xsection); //
        _sampleWeight.push_back(_totalLuminosity/equivalent_lumi);
        _sampleNeventsSumFiles.push_back(0);//
        _sampleNeventsTree.push_back(0);

        int nDump = 1; while (_sampleNevents.back()/nDump>20) nDump *= 10;
        _nDump10.push_back(nDump);

        // Printout summary for MC samples
        printf("MC> %s - %s: xsec %f pb, events %d, original events %d, weight %f\n", id.Data(), _sampleFile.back().Data(), _sampleXsection.back(),_sampleNevents.back(),total_events_for_xsection,_sampleWeight.back());
	std::cout << "HOLA HOLA HOLA HOLA HOLA HOLA HOLA HOLA " << total_events_for_xsection << " " << std::endl;
}

bool CMSAnalysis::SetTreeFile(int i, int fileJ) {
        if (_file) {
                _file->Close();
                _file = NULL;
        }

        std::cout<<"Setting tree File"<<std::endl;

        TString thisfile= _sampleFile[i]+"_"+std::to_string(fileJ)+".root";
	//TString thisfile= _sampleFile[i];
        printf("Processing sample '%s'...\n", thisfile.Data());

        _currentIndex = i;
        _file = new TFile(thisfile,"READONLY");

        std::cout<<"Sample processed"<<std::endl;
        if (!_file) return 0;

        _NANOTREE = NULL;
        _file->GetObject("Events",_NANOTREE);

        //TFile _file(thisfile,"READONLY");
        //_file.GetObject("Events",_NANOTREE);
        if(_NANOTREE==NULL) return 0;

        TH1F *histoEvents=(TH1F*)_file->Get("eventcount_Skim");
        if(histoEvents) {
                int nEventsFromHisto=_sampleNeventsSumFiles[i]+histoEvents->GetBinContent(1);
                _sampleNeventsSumFiles[i]=nEventsFromHisto;
        }

        std::cout<<_sampleNeventsSumFiles[i]<<std::endl;      

        // Reset buffer
        _BUFFERBRANCHES.clear();
        _BUFFER.clear();

        // Determine the maximum size of the buffer for a tree entry
        // and reserve that size for it
        // Then set addresses in an internal buffer for all variables
        // Later read only the requested branches
        // This saves processing time
        unsigned int bufferSize_max = 0;
        TObjArray* branches = _NANOTREE->GetListOfBranches();
        for (Int_t ib=0; ib<branches->GetEntries(); ib++) {

                TBranch* branch = (TBranch*)branches->At(ib);
                TClass* clptr; EDataType this_type;
                branch->GetExpectedType(clptr,this_type);
                unsigned int content_size = TDataType::GetDataType(this_type)->Size();
                if (branch->GetEntryOffsetLen()>0) content_size *= branch->GetMaxBaskets();
                bufferSize_max += content_size;
        }

        _BUFFER.resize(bufferSize_max);

        unsigned int bufferIndex = 0;
        for (Int_t ib=0; ib<branches->GetEntries(); ib++) {

                TBranch* branch = (TBranch*)branches->At(ib);
                branch->SetStatus(1);
                branch->SetAddress(&_BUFFER[bufferIndex]);

                TClass* clptr; EDataType this_type;
                branch->GetExpectedType(clptr,this_type);
                unsigned int content_size = TDataType::GetDataType(this_type)->Size();
                if (branch->GetEntryOffsetLen()>0) content_size *= branch->GetMaxBaskets();
                bufferIndex += content_size;
        }

        int nentriesInTree = _NANOTREE->GetEntriesFast();
        printf("\tReading %d entries in this file \n",nentriesInTree);
        //printf("\tMaximum size of an entry in this tree is %d bytes\n", bufferSize_max);

        _sampleNeventsTree[i]=nentriesInTree;


        return 1;

}

void CMSAnalysis::AddMCSample(const TString& id, const TString& file, int maxevents, double xsec, int total_events_for_xsection) 
{
	if (_dataIndex<0) 
	{
      		printf(">>> Warning: you should call AddDataSample first, to define the luminosity!!!\n");
     		printf(">>> SAMPLE NOT ADDED!\n");
      		return;
  	}
  
	_sampleId.push_back(id);
  	_sampleFile.push_back(file);

  	TFile file_tmp(file,"READONLY");
  	TTree* tree_tmp = 0;
  	file_tmp.GetObject("Events",tree_tmp);
  	if (!tree_tmp) printf("Error reading tree for file %s; crash expected!!\n", file.Data());
  	int maxeventsInTree = tree_tmp->GetEntriesFast();
  	_sampleXsection.push_back(xsec);
  	_signalFlag.push_back(false); // set to false to start with
        _sampleNFiles.push_back(1);
  	double equivalent_lumi = maxeventsInTree/xsec;
  	if (total_events_for_xsection>0 && xsec>0) equivalent_lumi = total_events_for_xsection/xsec;
  	if (maxevents<0 || maxevents>maxeventsInTree) 
	{
      		_sampleNevents.push_back(maxeventsInTree);
  	} 
	else 
	{
      		_sampleNevents.push_back(maxevents);
      		equivalent_lumi *= double(maxevents)/maxeventsInTree;
  	}
  
	_sampleWeight.push_back(_lumi/equivalent_lumi);

  	int nDump = 1; while (_sampleNevents.back()/nDump>20) nDump *= 10;
  	_nDump10.push_back(nDump);

  	// Printout summary for MC samples
  	printf("MC> %s - %s: xsec %f pb, events %d, weight %f\n", id.Data(), _sampleFile.back().Data(), _sampleXsection.back(),_sampleNevents.back(),_sampleWeight.back());
}

void CMSAnalysis::AddMCSignalSample(const TString& id, const TString& file, int maxevents, double xsec, int total_events_for_xsection) 
{
	if (_dataIndex<0) 
	{
      		printf(">>> Warning: you should call AddDataSample first, to define the luminosity!!!\n");
     		printf(">>> SAMPLE NOT ADDED!\n");
      		return;
  	}
  
	AddMCSample(id, file, maxevents, xsec, total_events_for_xsection);
  	_signalFlag[_sampleId.size()-1] = true; // set signal flag to true
      //  _sampleNFiles.push_back(0);
}

void CMSAnalysis::AddMCSignalSampleFiles(const TString& id, const TString& dir, const TString& file, int maxevents, double xsec, int total_events_for_xsection, int nfiles) 
{
  	std::cout << "AddMCSignalSampleFiles starts " << id << std::endl << std::endl;

	if (_dataIndex<0) 
	{
      		printf(">>> Warning: you should call AddDataSample first, to define the luminosity!!!\n");
     		printf(">>> SAMPLE NOT ADDED!\n");
      		return;
  	}


	AddMCSampleFiles(id, dir, file,  maxevents,  xsec, total_events_for_xsection, nfiles);
  	_signalFlag[_sampleId.size()-1] = true; // set signal flag to true
      	//_sampleNFiles.push_back(0);
}

void CMSAnalysis::AddPlot1D(const TString& name, const TString& title, int nbins, double xmin, double xmax) 
{
      	for (unsigned int i=0; i<_sampleId.size(); ++i) 
	{
        	bool existing = false;
            	for (unsigned int j=0; j<hists_1D.size(); ++j) 
		{
                  	TString thisname = hists_1D[j]->GetName();
                  	if (thisname == _sampleId[i]+"_"+name) 
			{
                        	existing = true; 
                        	break;
                  	}
            	}
            	if (existing) continue;

            	hists_1D.push_back(new TH1D(_sampleId[i]+"_"+name, title, nbins, xmin, xmax));
            	hists_1D[hists_1D.size()-1]->Sumw2();
		//hists_1D[hists_1D.size()-1]->GetXaxis()->SetTitle(title);
      }
}
void CMSAnalysis::ScalePlots(){

        for (unsigned int i=0; i<_sampleFile.size(); ++i) {
                if (_sampleId[i].Contains("Data")) continue;
                if (_sampleNeventsSumFiles[i]==0) continue;
                if (_sampleNeventsSumFiles[i]==_sampleNevents[i]) continue;
                std::cout<<_sampleId[i]<<"   "<<_sampleNeventsSumFiles[i]<<"  "<<_sampleNevents[i]<<std::endl;
                TString prefix = _sampleId[i] + "_";
                for (unsigned int j=0; j<hists_1D.size(); ++j) {
                        TString histname = hists_1D[j]->GetName(); 
                        if (!histname.BeginsWith(prefix)) continue;
                        double correctScale=1.00*_sampleNevents[i]/_sampleNeventsSumFiles[i];
                        hists_1D[j]->Scale(correctScale);
                } 

        }

}
void CMSAnalysis::FillPlot1D(const TString& name, int isample, double value, double weight) //PileUp
{
      	for (unsigned int j=0; j<hists_1D.size(); ++j) 
	{
		//std::cout << hists_1D[j]->GetName() << " " << _sampleId[isample]+"_"+name << std::endl;
        	if (hists_1D[j]->GetName()==_sampleId[isample]+"_"+name) 
		{
                	hists_1D[j]->Fill(value,_sampleWeight[isample]*weight);

                  	return;
            	}
      	}
	//std::cout << "**************************Sample Weight " << _sampleWeight[isample]  << std::endl;
}


void CMSAnalysis::SavingHistograms(int isample, const TString& name, const TString& option) 
{
	TFile *f = TFile::Open(_sampleId[isample]+"_ToPlot.root", option);

	for (unsigned int j=0; j<hists_1D.size(); ++j) 
	{
        	if (hists_1D[j]->GetName()==_sampleId[isample]+"_"+name) 
		{
                	hists_1D[j]->Write();

            	}
      	}
	//f->ls();
	f->Close();
	//delete f;
	return;

} 
void CMSAnalysis::DrawPlot1D(const TString& name, const TString& suffix, const TString& dir) //name of histogram
{
      	setTDRStyle();
	
      	int CanvasDivision = 2;

      	gROOT->SetStyle("Plain"); 

      	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0); 

      	TLegend* leg = new TLegend(0.52,0.6,0.88,0.88);
      	//TLegend* leg = new TLegend(0.80,0.75,1.00,1.00); // to see the whole histogram
      	leg->SetFillColor(4001);
 	leg->SetBorderSize(0);

        THStack* hMCStack = new THStack(name,name+"_histograms");
      	TH1D* hData = nullptr;
	std::vector<TH1D*> hSignal;
	TH1D* hBckgTotal;
	TH1D* RatioDataMC = nullptr;

      	// Get data; stack all MC except signal MC
      	unsigned int nhists = hists_1D.size();

	int colors[13] = {kBlue+1,  kYellow+1, kPink+1, kGreen+1, kOrange+1, kAzure+1, kMagenta+1, kCyan+1, kSpring+1, kTeal+1,   kViolet+1, kGray, kRed+1, };
	//int colors[10] = {kRed-9, kMagenta-9, kBlue-9, kPink-9, kViolet-9, kAzure-9, kGreen-9, kCyan-9, KSpring-9, kTeal-9], 
	
      	int mcindex = -1;

	TH1D* suma = nullptr;	
	TH1D* suma_signal = nullptr;
	//std::cout << "TEST 1 "  << nhists << " " << hists_1D.size() << " "  << hists_1D[99] << std::endl;

	int alreadyFilledData=0;

      	for (unsigned int j=0; j<nhists; ++j) 
	{
		//std::cout<<"Histograma: "<< hists_1D[j]->GetName()<<std::endl;
	
        	TString histname = hists_1D[j]->GetName();
            	TString suffix_hist = "_" + name;

		//std::cout << "Estoy aqui " << name << " " << j << " " << histname << " " << suffix_hist <<  std::endl;

            	if (!histname.EndsWith(suffix_hist)) continue;

            	for (unsigned int i=0; i<_sampleFile.size(); ++i) 
		{
			//std::cout << "TEST 1.1" << std::endl;
                  	TString prefix = _sampleId[i] + "_";

                  	if (!histname.BeginsWith(prefix)) continue;
			if(_sampleId[i].Contains("Data"))
			{
				if(alreadyFilledData==0) 
				{
					hData = hists_1D[j];
					alreadyFilledData =1;
					hData->SetName(_sampleId[i]+"_"+histname);
					hData->GetXaxis()->SetTitle(hists_1D[j]->GetTitle());
					hData->SetTitle(hists_1D[j]->GetTitle());
					hData->GetYaxis()->SetTitle("Events");
					leg->AddEntry(hData, "Data", "P");
					//std::cout << "TEST 1.2" << std::endl;
				}
				else hData->Add(hists_1D[j]);
				hData->SetMarkerStyle(20);
				hData->SetMarkerSize(1.0);
				hData->GetXaxis()->SetTitle(hists_1D[j]->GetTitle());
				hData->SetTitle(hists_1D[j]->GetTitle());
				//hData->GetYaxis()->SetTitle("Events");
				break;
			}
			else 
			{
                        	mcindex++;
                        	int color = colors[mcindex%10];
                        	hists_1D[j]->SetLineWidth(3);
                        	hists_1D[j]->SetLineColor(TColor::GetColorDark(color));
				//hists_1D[j]->SetLineColor(TColor::GetColorBright(color));
 				hists_1D[j]->SetMarkerStyle(7);
                        	hists_1D[j]->SetMarkerSize(0.5);
				hists_1D[j]->GetYaxis()->SetTitle("Events");
				//hists_1D[j]->GetXaxis()->SetTitle(hists_1D[j]->GetTitle());
				hists_1D[j]->SetTitle(hists_1D[j]->GetTitle());
                       	 	leg->AddEntry(hists_1D[j],_sampleId[i].Data(),"F");
				
                        	// Do not add the signal component to the stack yet
                        	if (!_signalFlag[i]) 
				{
					hMCStack->Add(hists_1D[j]); 
                        		hists_1D[j]->SetFillColor(color);
					hists_1D[j]->SetTitle(hists_1D[j]->GetTitle());

					//hMCStack->GetXaxis()->SetTitle(hists_1D[j]->GetTitle());
					//std::cout << "TEST 1.4" << std::endl;
					if(suma==nullptr) {suma = (TH1D*)hists_1D[j]->Clone(); }
					else suma->Add(hists_1D[j]);
				}
                        	break;
                  	}

            	}
				
      	}

	if(suma==nullptr) 
	{
		std::cout << "ERROR : Histogram " << name << " not found"  << std::endl;
		return;
	}


	RatioDataMC = (TH1D*)hData->Clone("hRatio");
	RatioDataMC->SetName("RatioDataMC_"+name);
	std::cout << "TEST 2 " << suma << " RatioDMC"  << RatioDataMC << std::endl;
	
	RatioDataMC->Divide(suma);

     	// Add the signals components to the stack now
	std::cout << "TEST 3" << std::endl;

	for(unsigned int i=0; i<_sampleFile.size(); ++i)
	{
		if(!_signalFlag[i]) continue;
		for (unsigned int j=0; j<nhists; ++j) 
		{
                  	if (hists_1D[j]->GetName()==_sampleId[i]+"_"+name) 
			{
				TH1D* hSignal_tmp = hists_1D[j];
                        	hSignal_tmp->SetLineWidth(3);
                        	//hSignal_tmp->SetLineColor(kBlue);
				//hSignal_tmp->SetLineColor(kRed);
				hSignal.push_back(hSignal_tmp);
				//hSignal[0]->SetLineColor(kBlue);
				//hSignal[1]->SetLineColor(kRed);
                        	break;
                  	}
            	}
	}
	suma_signal = (TH1D*) hSignal[0]->Clone();
	suma_signal->SetName("sumasignal_clone");
	//for(auto &h : hSignal)
	//{
	suma_signal->Add(hSignal[1]);
	//}

	std::cout << "SUMA " << suma_signal->Integral() << std::endl;
	TString c1name = "c1_" + name;
      	TCanvas* c1 = new TCanvas(c1name.Data(),c1name.Data(),10,10,1000,1000);
	//TCanvas* c1 = new TCanvas(c1name.Data(),c1name.Data(),-1);
std::cout << "TEST 3" << std::endl;
	TPad* pad1 = new TPad("pad1", "The pad 70% of the height",0.0,0.3,1.0,1.0,31);
	TPad* pad2 = new TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.3,32);

	pad1->SetFillColor(kWhite);
	pad2->SetFillColor(kWhite);
	pad1->Draw();
	pad2->Draw();
      	std::cout << "TEST 3" << std::endl;
	pad1->cd();
	pad1->SetBottomMargin(0.1);
	pad1->SetLeftMargin(0.2);
      	//hData->SetTitle("");
      	hData->GetYaxis()->SetTitleOffset(1.8);
      	hData->GetXaxis()->SetTitleOffset(1.8);
	hData->GetYaxis()->SetTitle("Events");
	hData->SetLabelSize(0.05, "XYZ");
	
	hData->SetMaximum(hData->GetMaximum()*1.5);
 	//if(hMCStack->GetMaximum()<hData->GetMaximum()) hMCStack->SetMaximum(hData->GetMaximum()*1.5); 
	
	for(int bin=0; bin<hData->GetNbinsX(); bin++)
	{
		if(suma->GetBinContent(bin)!=0 && suma_signal->GetBinContent(bin)/suma->GetBinContent(bin) > 0.5)
		hData->SetBinContent(bin, -1000);
	}

      	hData->Draw("e");
	hMCStack->Draw("same hist F");

      	TLatex* preliminary = new TLatex(0.2,0.95,"CMS preliminary");
      	preliminary->SetNDC();
      	preliminary->SetTextFont(42);
      	preliminary->Draw();

        TLatex* luminosity = new TLatex(0.6,0.95,Form("L = %.1f fb^{-1} (13 TeV)",_totalLuminosity/1000.));
	luminosity->SetNDC();
      	luminosity->SetTextFont(42);
      	luminosity->Draw();
	std::cout << "TEST 3" << std::endl;
      	hData->Draw("esame F");
	
	for(auto &h : hSignal) h->Draw("hist same");
	
	leg->SetNColumns(2);
      	leg->Draw();

      	gPad->SetTicks(1,1);
      	gPad->RedrawAxis();
	
	pad2->cd();
      	setTDRStyle();
	pad2->SetTopMargin(0.04);
	pad2->SetLeftMargin(0.2);
	pad2->SetBottomMargin(0.3);
	//gStyle->SetOptStat(0);

	//int binmax=RatioDataMC->GetMaximumBin(); double x = RatioDataMC->GetBinCenter(binmax);

	RatioDataMC->Draw("P");

	RatioDataMC->SetMarkerColor(kRed+1);//roja (signal is not included in the Stack
	RatioDataMC->SetMaximum(2);
	RatioDataMC->SetMinimum(0);
	RatioDataMC->SetMarkerSize(0.8);
	RatioDataMC->GetYaxis()->SetTitle("Data/MC");
      	RatioDataMC->GetYaxis()->SetTitleOffset(0.45);
	RatioDataMC->SetLabelSize(0.08, "xy");
	RatioDataMC->SetTitleSize(0.08, "xy");
	//RatioDataMC->SetTitleSize(0.08, "y");
	//RatioDataMC->GetYaxis()->SetTitle(0.1);
	//RatioDataMC->GetXCD axis()->SetTitle(0.1);
	std::cout << "TEST 3" << std::endl;
	//RatioDataMC->SetLabelSize(0.02, "y");

	RatioDataMC->GetYaxis()->SetNdivisions(4);
	//RatioDataMC->GetXaxis()->SetTextFont(22);
	pad2->SetGrid();
	
      	TString myDir = dir;
      	if ((myDir!="") && (myDir!=".")) 
	{
           	gSystem->MakeDirectory(myDir);
      	}
	else 
	{
            myDir = ".";
      	} 

      	if (suffix!="") 
	{  
            	c1->SaveAs(myDir+"/"+name+"_"+suffix+".png");

      	} 
	else 
	{
            	c1->SaveAs(myDir+"/"+name+".png");
		c1->SaveAs(myDir+"/"+name+".C");

      	}
	std::cout << "TEST 3" << std::endl;
}

void CMSAnalysis::AddPlot1D_bare(const TString& name, const TString& title, int nbins, double xmin, double xmax) {
      bool existing = false;
      for (unsigned int j=0; j<hists_1D.size(); ++j) {
            TString thisname = hists_1D[j]->GetName();
            if (thisname == name) {
                  existing = true; 
                  break;
            }
      }
      if (existing) return;

      hists_1D.push_back(new TH1D(name, title, nbins, xmin, xmax));
      hists_1D[hists_1D.size()-1]->Sumw2();
}
  
void CMSAnalysis::FillPlot1D_bare(const TString& name, double value, double weight) {

      for (unsigned int j=0; j<hists_1D.size(); ++j) {
            if (hists_1D[j]->GetName()==name) {
                  hists_1D[j]->Fill(value,weight);
                  return;
            }
      }
}
  
void CMSAnalysis::DrawPlot1D_bare(const TString& name, const TString& suffix, const TString& dir) {
      setTDRStyle();

      //gROOT->SetStyle("Pub"); 
      gROOT->SetStyle("Plain"); 

     // gStyle->SetPadGridX(true); 
     // gStyle->SetPadGridY(true);
      //gStyle->SetOptStat(111122222);
      //gStyle->SetOptStat(1111);
      gStyle->SetOptStat(0);

      TLegend* leg = new TLegend(0.52,0.6,0.85,0.85);
      //TLegend* leg = new TLegend(0.80,0.75,1.00,1.00); // to see the whole histogram
      leg->SetFillColor(kWhite);
      leg->SetBorderSize(0);
      THStack* hMCStack = new THStack(name,name+" histograms");
      TH1D* hData;

      // Get data; stack all MC except signal MC
      unsigned int nhists = hists_1D.size();
      //int colors[10] = {46, 47, 48, 49, 50, 51, 52, 53, 54, 55};
      //int colors[10] = {6, 7, 8, 9, 11, 44, 46, 2, 3, 4}; 
      //int colors[10] = {TColor::GetColor("#ff55ff"), TColor::GetColor("#00ffff"), TColor::GetColor("#ff5500"), 9, 11, 44, 46, 2, 3, 4}; 
      int colors[10] = {TColor::GetColor("#00ddff"), TColor::GetColor("#0055ff"), kMagenta, kRed, 11, 44, 46, 2, 3, 4}; 
      int mcindex = -1;
      for (unsigned int j=0; j<nhists; ++j) {
            TString histname = hists_1D[j]->GetName();
            if (histname==name) {
                  hData = hists_1D[j];
                  hData->SetMarkerStyle(20);
                  hData->SetMarkerSize(0.8);
                  leg->AddEntry(hData,"Data","P");
                  break;
            }
      }
      
      TString c1name = "c1_" + name;
      TCanvas* c1 = new TCanvas(c1name.Data(),c1name.Data(),10,10,600,600);
      hData->SetXTitle(hData->GetTitle());
      hData->SetTitle("");
      hData->SetTitleOffset(1.3);
      if (hData->GetMinimum()>0.) hData->SetMinimum(0.);
      hData->Draw("e");
      //TLatex* preliminary = new TLatex(0.56,0.92,"CMS preliminary");
      TLatex* preliminary = new TLatex(0.05,0.92,"CMS preliminary");

      preliminary->SetNDC();
      preliminary->SetTextFont(42);
      preliminary->Draw();
      hData->Draw("esame");
      leg->Draw();

      gPad->SetTicks(1,1);
      gPad->RedrawAxis();

      TString myDir = dir;
      if ((myDir!="") && (myDir!=".")) {
            gSystem->MakeDirectory(myDir);
      } else {
            myDir = ".";
      } 

      if (suffix!="") {
            c1->SaveAs(myDir+"/"+name+"_"+suffix+".root");
            c1->SaveAs(myDir+"/"+name+"_"+suffix+".jpg");
            c1->SaveAs(myDir+"/"+name+"_"+suffix+".pdf");
      } else {
            c1->SaveAs(myDir+"/"+name+".root");
            c1->SaveAs(myDir+"/"+name+".jpg");
            c1->SaveAs(myDir+"/"+name+".pdf");
      }

	
}
  
void CMSAnalysis::Normalize_Histogram(const TString& name) {
      unsigned int nhists = hists_1D.size();

      double dataInt = 0.;
      double mcInt = 0.;
      for (unsigned int j=0; j<nhists; ++j) {
            TString histname = hists_1D[j]->GetName();
            TString suffix_hist = "_" + name;
            if (!histname.EndsWith(suffix_hist)) continue;

            for (unsigned int i=0; i<_sampleFile.size(); ++i) {
                  TString prefix = _sampleId[i] + "_";
                  if (!histname.BeginsWith(prefix)) continue;
                  if (histname==_sampleId[_dataIndex]+"_"+name) {
                        dataInt = hists_1D[j]->Integral();
                        break;
                  } else {
                        mcInt += hists_1D[j]->Integral();
                        break;
                  }
            }
      }

      for (unsigned int j=0; j<nhists; ++j) {
            TString histname = hists_1D[j]->GetName();
            TString suffix_hist = "_" + name;
            if (!histname.EndsWith(suffix_hist)) continue;

            for (unsigned int i=0; i<_sampleFile.size(); ++i) {
                  TString prefix = _sampleId[i] + "_";
                  if (!histname.BeginsWith(prefix)) continue;
                  if (histname==_sampleId[_dataIndex]+"_"+name) {
                        hists_1D[j]->Scale(1./dataInt);
                        break;
                  } else {
                        hists_1D[j]->Scale(1./mcInt);
                        break;
                  }
            }
      }
}
  
void CMSAnalysis::Normalize_Histogram_bare(const TString& name) {
      unsigned int nhists = hists_1D.size();

      double dataInt = 0.;
      for (unsigned int j=0; j<nhists; ++j) {
            TString histname = hists_1D[j]->GetName();
            if (histname==name) {
                  dataInt = hists_1D[j]->Integral();
                  hists_1D[j]->Scale(1./dataInt);
                  return;
            }
      }
}
  
void CMSAnalysis::Divide_Histogram_bare(const TString& name1, const TString& name2) {
      unsigned int nhists = hists_1D.size();

      int i1 = -1;
      int i2 = -1;
      for (unsigned int j=0; j<nhists; ++j) {
            TString histname = hists_1D[j]->GetName();
            if (histname==name1) i1 = j;
            if (histname==name2) i2 = j;
            if (i1>=0 && i2>=0) break;
      }

      if (i1>=0 && i2>=0) {
            hists_1D[i1]->Divide(hists_1D[i2]);
      } else {
            printf("CMSANalysis::Divide_Histogram_bare: histograms missing, i1=%d, i2=%d\n", i1, i2);
      }
}
  
TH1D* CMSAnalysis::GetMCHistogram(const TString& name) {

      TH1D* newhist = NULL;
      
      unsigned int nhists = hists_1D.size();
      for (unsigned int j=0; j<nhists; ++j) {
            TString histname = hists_1D[j]->GetName();
            TString suffix_hist = "_" + name;
            if (!histname.EndsWith(suffix_hist)) continue;

            for (unsigned int i=0; i<_sampleFile.size(); ++i) {
                  TString prefix = _sampleId[i] + "_";
                  if (!histname.BeginsWith(prefix)) continue;
                  if (histname!=_sampleId[_dataIndex]+"_"+name) {
                        if (newhist==NULL) {
                              newhist = (TH1D*)hists_1D[j]->Clone(name+"_allMC");
                        } else {
                              newhist->Add(hists_1D[j]);
                        }
                        break;
                  }
            }
      }

      return newhist;
}
  
TH1D* CMSAnalysis::GetDataHistogram(const TString& name) {

      unsigned int nhists = hists_1D.size();
      for (unsigned int j=0; j<nhists; ++j) {
            TString histname = hists_1D[j]->GetName();
            TString suffix_hist = "_" + name;
            if (!histname.EndsWith(suffix_hist)) continue;

            for (unsigned int i=0; i<_sampleFile.size(); ++i) {
                  TString prefix = _sampleId[i] + "_";
                  if (!histname.BeginsWith(prefix)) continue;
                  if (histname==_sampleId[_dataIndex]+"_"+name) {
                        return hists_1D[j];
                  }
            }
      }

      // IF we are here the histogram was not found
      return NULL;
}
  
void CMSAnalysis::SetTree(int i) {
      if (_file) {
            _file->Close();
            _file = NULL;
      }

      printf("Processing sample '%s'...\n", _sampleFile[i].Data());
      _currentIndex = i;
      _file = new TFile(_sampleFile[i],"READONLY");

      _NANOTREE = NULL;
      _file->GetObject("Events",_NANOTREE);

      // Reset buffer
      _BUFFERBRANCHES.clear();
      _BUFFER.clear();

      // Determine the maximum size of the buffer for a tree entry
      // and reserve that size for it
      // Then set addresses in an internal buffer for all variables
      // Later read only the requested branches
      // This saves processing time
      unsigned int bufferSize_max = 0;
      TObjArray* branches = _NANOTREE->GetListOfBranches();
      for (Int_t ib=0; ib<branches->GetEntries(); ib++) {

            TBranch* branch = (TBranch*)branches->At(ib);
            TClass* clptr; EDataType this_type;
            branch->GetExpectedType(clptr,this_type);
            unsigned int content_size = TDataType::GetDataType(this_type)->Size();
            if (branch->GetEntryOffsetLen()>0) content_size *= branch->GetMaxBaskets();
            bufferSize_max += content_size;
      }

      _BUFFER.resize(bufferSize_max);

      unsigned int bufferIndex = 0;
      for (Int_t ib=0; ib<branches->GetEntries(); ib++) {

            TBranch* branch = (TBranch*)branches->At(ib);
            branch->SetStatus(1);
            branch->SetAddress(&_BUFFER[bufferIndex]);

            TClass* clptr; EDataType this_type;
            branch->GetExpectedType(clptr,this_type);
            unsigned int content_size = TDataType::GetDataType(this_type)->Size();
            if (branch->GetEntryOffsetLen()>0) content_size *= branch->GetMaxBaskets();
            bufferIndex += content_size;
      }

      int nentriesInTree = _NANOTREE->GetEntriesFast();
      printf("\tReading %d entries from a total of %d\n", _sampleNevents[i], nentriesInTree);
      //printf("\tMaximum size of an entry in this tree is %d bytes\n", bufferSize_max);

}

Int_t CMSAnalysis::SetEntry(int entryNumber) {

      // Set the entry
      if (_NANOTREE->LoadTree(entryNumber)<0) return -1;

      // Periodic printout, filesize dependent
      if (entryNumber%_nDump10[_currentIndex]==0) {
            printf("... event index %d\n", entryNumber);
      }

      return 0;
}

Int_t CMSAnalysis::GetEntry(int entryNumber) {

      // Set the entry
      if (SetEntry(entryNumber)!=0) return -1;

      // Read all active branches
      _NANOTREE->GetEntry(entryNumber); 

      return 0;
}
