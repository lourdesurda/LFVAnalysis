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
#include "TMath.h"
#include "TRandom2.h"
#include "CMSAnalysis.h"
#include <vector>
#include <TLorentzVector.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <TBranch.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TList.h>
#include <TAxisEditor.h>
#include <TAxis.h>
#include <TGraphErrors.h>

#include <TClonesArray.h>
#include <TVectorF.h>
#include <TTimeStamp.h>
#include <TH1F.h>
#include <TH2F.h>
#include <new>
#include <map>
#include <TPad.h>

const double electronMass = 0.0005109989461; //GeV
const double muonMass = 0.1056583745;
const double METMass = 0.0;

struct TopAnalysis : public CMSAnalysis 
{
	TopAnalysis()
	{
		init(); TFile("LFV_Analysis.root", "RECREATE").Close();
	};
	
	virtual ~TopAnalysis(){};

	void init () 
	{
		nbLooseSel=0; nbMediumSel=0;			
	}	

	//Preselection
	bool Preselect();   
	
	//Selection
	bool Select();

//Variables analisis

	//unsigned int njtsel;
 	unsigned int nbLooseSel;
      	unsigned int nbMediumSel;

	int Indexmusel;
	int Indexelsel;

	int nmusel;
	int nelsel;
	int njtsel;


// Select if there are at least two leptons wit some phase space cuts
	
  	const double ptmuCut = 29.;//25
  	const double etamuCut = 2.4;//2.5
  	const double isomuCut = 0.15;
  	const double ptelCut = 10.;//25
  	const double etaelCut = 2.5;
  	const double isoelCut = 0.06;//0.10
  	const double ptjtCut = 30.; //take it into account for Bjets//30
  	const double etajtCut = 4.;//2.65
  	const double btagLooseCut = 0.1241;//0.5426; // loose btag point pfCSV
  	const double btagMediumCut = 0.4184;//0.8484 ; // medium btag point pfCSV
  	const double btagTightCut = 0.7527;//0.9535 ; // tight btag point pfCSV
};
struct ASCII
{
		std::map<std::string, double> Samples;
		std::vector<std::string> HistogramsList;

		void txtCounter(const TString& name)
		{
			std::ifstream txtfile(name);
			std::string input = " ";
			do
			{
				txtfile >> input;
				if(txtfile.eof()) break;
				if (input!= " ") HistogramsList.push_back (input);
	
			}while (true);
			txtfile.close();
		}

		void mapsfilling(const TString& name)
		{
			std::ifstream txtfile(name);
			std::string input_sample;
			double input_xsec;

			while(true)
			{
				txtfile >> input_sample >> input_xsec;
				if(txtfile.eof()) break;
				Samples.insert({input_sample, input_xsec});			
			}
			txtfile.close();
		}
		
};

struct FOURVECTORS
{
	
	static double InvariantMass(const TLorentzVector &Particle1, const TLorentzVector &Particle2)
	{
		double _Invariant_Mass = (Particle1+Particle2).M(); return _Invariant_Mass;
	}

	static double TranverseMass(const TLorentzVector &Particle1, const TLorentzVector &Particle2)
	{
		double TotalEnergy = Particle1.Et()+Particle2.Et();
		double TotalPT = (Particle1+Particle2).Pt();
		return sqrt(TotalEnergy*TotalEnergy-TotalPT*TotalPT);
	}

	static double CollinearMass(const TLorentzVector &Particle1, const TLorentzVector &Particle2, const TLorentzVector &Particle3, const float Particle3_Particle2_phi) //Muon, Electron, MET
	{
		double ptnu = abs(Particle3.Et()*cos(Particle3_Particle2_phi));
		double visfrac = Particle2.Pt()/(Particle2.Pt()+ptnu);
		double m_e_Mass = (Particle1+Particle2).M();
		return m_e_Mass/sqrt(visfrac);
	}

};

struct ANGLES
{
	static double DeltaPhi(double phi1, double phi2)
	{         
		double dphi = phi1-phi2;         
		if ( dphi > M_PI )
		{                 
			dphi -= 2.0*M_PI;         
		} 
		else if ( dphi <= -M_PI ) 
		{                 
			dphi += 2.0*M_PI;         
		}         
		return dphi;
	}

	static double DeltaEta(double eta1, double eta2)
	{
		double deta = eta1-eta2;
		return deta;
	}

	static double DeltaR(double phi1, double phi2, double eta1, double eta2)
	{
		double deltar = sqrt( DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) + DeltaEta(eta1, eta2)*DeltaEta(eta1, eta2) );
		return deltar;
	}
};

struct HISTS : public CMSAnalysis 
{
	static void Add(CMSAnalysis &cms)
	{
		cms.AddPlot1D("hMet", "Missing E_{T} [GeV]", 100, -10., 200.);
  		cms.AddPlot1D("hMuon", "Muon pT [GeV]", 50, 25.0, 100.);
  		cms.AddPlot1D("hElectron", "Electron pT [GeV]", 50, 10.0, 150.);
  		cms.AddPlot1D("hJet", "Jet pT [GeV]", 50, 0.0, 250.);
        	//cms.AddPlot1D("hMetJets30", "Missing E_{T} with 2 Jets of 30GeV [GeV]", 50, 0., 500.);

  		cms.AddPlot1D("hMuonPhi", "Muon Phi", 100, -3.15, 3.15);
  		cms.AddPlot1D("hMuonEta", "Muon Eta", 100, -3.0, 3.0);

 		cms.AddPlot1D("hElectronPhi", "Electron Phi", 100, -3.15, 3.15);
  		cms.AddPlot1D("hElectronEta", "Electron Eta", 100, -3.0, 3.0);

  		cms.AddPlot1D("hJetEta", "Most energetic jet Eta", 100, -3.0, 3.0);
		/*cms.AddPlot1D("hPVs","",50,0,50);         
		cms.AddPlot1D("hPVs2","",50,0,50);

		cms.AddPlot1D("hNJets","",10,0,10);         
		cms.AddPlot1D("hNJets2","",10,0,10);*/

		cms.AddPlot1D("hnVertex", "Number of Vertex", 100, 0, 100);
  		cms.AddPlot1D("hnMuon", "Number of Muons", 6, 0, 6);
  		cms.AddPlot1D("hnElectron", "Number of Electrons", 7, 0, 7);
  		cms.AddPlot1D("hnJet", "Number of Jets", 25, 0, 25);	

  		cms.AddPlot1D("hmetPhi", "MET Phi", 100, -3.15, 3.15);

		cms.AddPlot1D("hMuonElectronPhi", "Muon-Electron Phi", 100, -3.15, 3.15);
		cms.AddPlot1D("hMuonMetPhi", "Muon-MET Phi", 100, -3.15, 3.15);
		cms.AddPlot1D("hElectronMetPhi", "Electron-Met Phi", 100, -3.15, 3.15);

		cms.AddPlot1D("hInvariantMassMuonElectron", "Invariant Mass Muon-Electron [GeV]", 100, 0, 200);
		cms.AddPlot1D("hTransverseMassMuonMET", "Transverse Mass Muon-MET [GeV]", 100, 0, 300);
		cms.AddPlot1D("hTransverseMassElectronMET", "Transverse Mass Electron-MET [GeV]", 100, 0, 300);

		cms.AddPlot1D("hnbLooseSel", "Number of BJets Loose", 5, -0.5, 4.5);
		cms.AddPlot1D("hnbMediumSel", "Number BJets Medium", 5, -0.5, 4.5);
		
		cms.AddPlot1D("hdeltaR", "#DeltaR Muon-Electron", 13, -2.4, 10.0);
		cms.AddPlot1D("hCollinearMass", "Collinear Mass Muon-Electron [GeV]", 100, 0, 350);
		
		cms.AddPlot1D("hnMuoSel", "Number of selected Muons", 10, 0, 10);
		cms.AddPlot1D("hnEleSel", "Number of selected Electrons", 10, 0, 10);
		cms.AddPlot1D("hnJetSel", "Number of selected jets", 11, -1, 10);

		//cms.AddPlot1D("hCortes", "", 10, 0, 10);    
        	/*cms.AddPlot1D("hJetPt", "Jet Pt [GeV]", 50, 0., 500.);         
		cms.AddPlot1D("hJetEta", "Jet Eta", 50, -5., 5.);         
		cms.AddPlot1D("hEleIDMVAFallV280","",2,0,2);         
		cms.AddPlot1D("hEleIDMVAFallV190","",2,0,2);         
		cms.AddPlot1D("hEleIso03","",100,0,1);         

     
		cms.AddPlot1D("hVisibleMassNoVeto", "Visible Mass (all) [GeV]", 300, 0., 300.);         
		cms.AddPlot1D("hVisibleMassBTAGSel", "Visible Mass (Btag Sel) [GeV]", 300, 0., 300.);        
 		cms.AddPlot1D("hVisibleMass", "Visible Mass (Btag Veto) [GeV]", 300, 0., 300.);         
		cms.AddPlot1D("hVisibleMass0Jets", "Visible Mass 0Jets (BTag Veto) [GeV]", 300, 0., 300.);         
		cms.AddPlot1D("hVisibleMass1Jets", "Visible Mass 1Jets (BTag Veto) [GeV]", 300, 0., 300.);         
		cms.AddPlot1D("hVisibleMass2Jets", "Visible Mass 2Jets (BTag Veto) [GeV]", 300, 0., 300.);*/
	}
	static void Fill(TopAnalysis &cms, unsigned int iSample)
	{
		Float_b(MET_pt); 		cms.FillPlot1D("hMet", iSample, MET_pt);	    
	    	VFloat_b(Electron_pt); 		cms.FillPlot1D("hElectron", iSample, Electron_pt[cms.Indexelsel]);	  
            	VFloat_b(Muon_pt); 		cms.FillPlot1D("hMuon", iSample, Muon_pt[cms.Indexmusel]);
            	VFloat_b(Jet_pt); 		cms.FillPlot1D("hJet", iSample, Jet_pt[0]);		
   	
	    	Int_b(PV_npvs);			cms.FillPlot1D("hnVertex", iSample, PV_npvs);
    	    	Int_b(nMuon);			cms.FillPlot1D("hnMuon", iSample, nMuon);
    	    	Int_b(nElectron); 		cms.FillPlot1D("hnElectron", iSample, nElectron);
    	    	Int_b(nJet); 			cms.FillPlot1D("hnJet", iSample, nJet);//nJet

           	Float_b(MET_phi);		cms.FillPlot1D("hmetPhi", iSample, MET_phi);
            	VFloat_b(Electron_phi);		cms.FillPlot1D("hElectronPhi", iSample, Electron_phi[cms.Indexelsel]);
            	VFloat_b(Muon_phi);		cms.FillPlot1D("hMuonPhi", iSample, Muon_phi[cms.Indexmusel]);

            	VFloat_b(Electron_eta);		cms.FillPlot1D("hElectronEta", iSample, Electron_eta[cms.Indexelsel]);
            	VFloat_b(Muon_eta);		cms.FillPlot1D("hMuonEta", iSample, Muon_eta[cms.Indexmusel]);
	        VFloat_b(Jet_eta); 		cms.FillPlot1D("hJetEta", iSample, Jet_eta[0]);		

		cms.FillPlot1D("hMuonElectronPhi", iSample, ANGLES::DeltaPhi(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel]));
		cms.FillPlot1D("hMuonMetPhi", iSample, ANGLES::DeltaPhi(Muon_phi[cms.Indexmusel], MET_phi));
		cms.FillPlot1D("hElectronMetPhi", iSample, ANGLES::DeltaPhi(Electron_phi[cms.Indexelsel], MET_phi));

		TLorentzVector LorentzElec;	LorentzElec.SetPtEtaPhiM((double)Electron_pt[cms.Indexelsel], (double) Electron_eta[cms.Indexelsel], (double) Electron_phi[cms.Indexelsel], electronMass);
		TLorentzVector LorentzMuon;	LorentzMuon.SetPtEtaPhiM(Muon_pt[cms.Indexmusel], Muon_eta[cms.Indexmusel], Muon_phi[cms.Indexmusel], muonMass);
		TLorentzVector LorentzMET;	LorentzMET.SetPtEtaPhiM(MET_pt, 0.0, MET_phi, METMass);

		cms.FillPlot1D("hInvariantMassMuonElectron", iSample, FOURVECTORS::InvariantMass(LorentzMuon, LorentzElec));
		cms.FillPlot1D("hTransverseMassMuonMET", iSample, FOURVECTORS::TranverseMass(LorentzMuon, LorentzMET));
		cms.FillPlot1D("hTransverseMassElectronMET", iSample, FOURVECTORS::TranverseMass(LorentzElec, LorentzMET));

		cms.FillPlot1D("hnbLooseSel", iSample, cms.nbLooseSel);
		cms.FillPlot1D("hnbMediumSel", iSample, cms.nbMediumSel);

		//cms.FillPlot1D("hdeltaR", iSample, sqrt( (Muon_phi[cms.Indexmusel]-Electron_phi[cms.Indexelsel])*(Muon_phi[cms.Indexmusel]-Electron_phi[cms.Indexelsel]) 
							//+ (Muon_eta[cms.Indexmusel]-Electron_eta[cms.Indexelsel])*(Muon_eta[cms.Indexmusel]-Electron_eta[cms.Indexelsel]) ) );
		//cms.FillPlot1D("hdeltaR", iSample, sqrt( (DeltaPhi(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel])*DeltaPhi(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel]) 
							//+ (Muon_eta[cms.Indexmusel]-Electron_eta[cms.Indexelsel])*(Muon_eta[cms.Indexmusel]-Electron_eta[cms.Indexelsel]) ) );

		cms.FillPlot1D("hdeltaR", iSample, ANGLES::DeltaR(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel], Muon_eta[cms.Indexmusel], Electron_eta[cms.Indexelsel]));

		//cms.FillPlot1D("hCollinearMass", iSample, FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET, MET_phi-Electron_phi[cms.Indexelsel]) );

		cms.FillPlot1D("hCollinearMass", iSample, FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET, ANGLES::DeltaPhi(Electron_phi[cms.Indexelsel], MET_phi)) );

		cms.FillPlot1D("hnMuoSel", iSample, cms.nmusel);
		cms.FillPlot1D("hnEleSel", iSample, cms.nelsel);
		cms.FillPlot1D("hnJetSel", iSample, cms.njtsel);
	
		/*std::cout << " N JET sel when the loop is over " << cms.njtsel << " nJets iniciales "  << nJet << std::endl;
		std::cout << " N ELE sel when the loop is over " << cms.nelsel << " nElectrons iniciales "  << nElectron << std::endl;
		std::cout << " N MUO sel when the loop is over " << cms.nmusel << " nMuons iniciales "  << nMuon << std::endl;*/
	}

	static void Save(TopAnalysis &cms, unsigned int iSample, ASCII &txtFiles)
	{
		TString option = "RECREATE"; 

		for(auto &name: txtFiles.HistogramsList)
		{
			cms.SavingHistograms(iSample, name, option);
			if(option =="RECREATE") option = "UPDATE";
		}
	}

	static void Draw(CMSAnalysis &cms, ASCII &txtFiles)
	{
		for(auto &name: txtFiles.HistogramsList)
		{
			std::cout << "Name of Histograms " << name << std::endl;
		 	cms.DrawPlot1D(name);
			std::cout << "Name of Histograms que se ha ploteado" << name << std::endl;
			
		}
	}

};

int main(int argc, char** argv)
{
  	// Directory where files are sitting
	TString dir = "/eos/user/c/cepeda/LFVNANO/Skimmed2/"; //TString dir = "/eos/user/c/cepeda/LFVNANO/MUESKIM/"; 
	
	// Initialize analysis structure
  	TopAnalysis cms; 
	ANGLES angulos;
  	// Add ****DATA**** sample. Input is: (titleId, file, luminosity) // One should add this sample first, to define the luminosity for normalizations
	
	double lumi_invpbA =13.5*1000.;        
	cms.AddDataSampleFiles("DataA", dir+"SingleMuRun2018A/", "tree", lumi_invpbA, 177); 

	double lumi_invpbB =6.8*1000.;        
	cms.AddDataSampleFiles("DataB", dir+"SingleMuRun2018B/", "tree", lumi_invpbB, 85); 
       
	double lumi_invpbC =6.6*1000.;        
	cms.AddDataSampleFiles("DataC", dir+"SingleMuRun2018C/", "tree", lumi_invpbC, 71);
        
	double lumi_invpbD = 32.*1000.; //RUND
  	cms.AddDataSampleFiles("DataD", dir+"SingleMuRun2018D/", "tree", lumi_invpbD, 1);

  	int maxevents = -1;// Set maximum number of events in MC//int maxevents = 3000000; for fast tests but reasonabe MC statistics int maxevents = 500000; for very fast tests with MC int maxevents = 50000; for very fast tests with MC 	

  	// Add ****SIGNAL**** sample. Input is: (titleId, file, maxevents, xsection in pb, total_events_for_xsection) 	

	TString inputSignalFile = "GluGlu_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8/";	
	double xsec_signal = 48.58*0.1;

	//cms.AddMCSignalSample("GG", inputSignalFile, maxevents, xsec_signal, 8000);
	cms.AddMCSignalSampleFiles("GG (B = 10%)", inputSignalFile, "tree", maxevents, xsec_signal, 8000, 1);

	TString inputSignalFile_VBF = "VBF_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8/";	
	double xsec_signal_VBF = 3.782*0.1;
	//cms.AddMCSignalSample("VBF", inputSignalFile_VBF, maxevents, xsec_signal_VBF, 114800);
	cms.AddMCSignalSampleFiles("VBF (B = 10%)", inputSignalFile_VBF, "tree", maxevents, xsec_signal_VBF, 114800, 1);

	// Add ****MC SAMPLES****. Input is: (titleId, file, maxevents, xsection in pb)

	ASCII txtFiles;

	txtFiles.mapsfilling("NewSamples_Sept2019.txt");
	TString legendname[5] = {"WW", "W+Jets", "TW", "#bar{T}W", "T#bar{T}"};
	int ss = 0;
	for(auto &s : txtFiles.Samples)
	{
		TString inputBkgdFile = s.first;
		double xsec_bkgd = s.second; 
		cms.AddMCSampleFiles(legendname[ss], dir+inputBkgdFile, "tree", maxevents, xsec_bkgd, cms.NumberOfEntries(dir+inputBkgdFile, 1), 1);
		ss++;
	}

	double xsec_DY = 2075.14*3; // DY 13 TeV, Mll>50 GeV         
    	TString inputBkgdFile = "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/";
	cms.AddMCSampleFiles("DY", dir+inputBkgdFile, "tree", maxevents, xsec_DY, cms.NumberOfEntries(dir+inputBkgdFile, 129), 129);

	//cms.AddMCSampleFiles("DY", dir+"DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/", "tree", maxevents, xsec_DY, maxevents_forXsection_DY, 129);

	//cms.AddMCSample("WZ_background_DAS", "WZ_background_DAS.root", maxevents, 27.6, 370000); //XSDB from DAS 
	//cms.AddMCSample("ZZ_background_DAS", "ZZ_background_DAS.root", maxevents, 12.14, 637000);

  	// Initialize 1D histograms

	HISTS::Add(cms);
  	
  	// Loop on samples
  	unsigned int nsamples = cms.GetNumberOfSamples();

	//txtFiles.txtCounter("NameOfHistograms.txt");

  	for (unsigned int iSample=0; iSample<nsamples; ++iSample) 
	{
		std::cout<<"Procesando muestra " << iSample << " " << cms.GetSampleId(iSample) << " " <<cms.GetNumberOfFilesInSample(iSample)<<std::endl;

            	for (unsigned int iFileInSample=1; iFileInSample<=cms.GetNumberOfFilesInSample(iSample); iFileInSample++)
		{
			// Set tree                 
			bool foundTree=true;                 
	                      
			if(!cms.SetTreeFile(iSample,iFileInSample)) foundTree=false;                 
		
			if (!foundTree) continue;
      			
			// Loop on events for the current sample
      			unsigned int nevents = 4000000;//cms.GetNumberOfEvents(iSample);//coger entradas por tree


      			for (unsigned int iEvent=0; iEvent<nevents; iEvent++)
			{
          	 		// Set next entry
           	 		if (cms.SetEntry(iEvent)<0) break;

           	 		// The event number has typeId "unsigned long long" because it can be a big number
           	 		//auto eventnumber = cms.Get<unsigned long long>("event");
           	 		//printf("Event number: %llu\n", eventnumber);

           	 		// Preselect (summary branch must be read before)
           	 		if (!cms.Preselect()) continue;
          
           	 		// Select (event branch must be read before)
           	 		if (!cms.Select()) continue;

				/*double myweight = 1. ;
				if(!cms.GetSampleId(iSample).Contains("Data") && !cms.GetSampleId(iSample).Contains("GG") && !cms.GetSampleId(iSample).Contains("VBF"))
				{
					Float_b(puWeight);
					myweight= puWeight;
				}*/
            			// Fill histograms			
				HISTS::Fill(cms, iSample);
			
      			}			
		}
			
		//HISTS::Save(cms, iSample, txtFiles);

  	}  	//To see things interactively (commnent otherwise) 

	//cms.ScalePlots();
  	auto app = new TRint("Top Analysis", &argc, argv); 

	txtFiles.txtCounter("NameOfHistograms.txt");

  	HISTS::Draw(cms, txtFiles);

  	// To see things interactively (commnent otherwise) 

  	if (!gROOT->IsBatch()) app->Run();

  	return 0;
}

bool TopAnalysis::Preselect()
{
	// Cut on summary information to save processing time
  	UInt_b(nElectron);
  	UInt_b(nMuon);
	Bool_b(HLT_IsoMu27);
	Bool_b(HLT_IsoMu24);
  	//printf("nmu: %d, nel: %d\n", nMuon, nElectron);

  	if (nMuon+nElectron<2) 	return false; // look for events with >= 2 leptons
	if (!HLT_IsoMu27 && !HLT_IsoMu24) return false;

  	UInt_b(nJet);

  	//printf("njt: %d\n", nJet);
 	//if (nJet>=2) 	return false; // look for events with at least 2 jets
  	return true;    
}

bool TopAnalysis::Select() 
{
  	//init();
  	nmusel = 0;
  	nelsel = 0;
	njtsel = 0;
	Indexmusel = -1;
	Indexelsel = -1;
	nbLooseSel=0; 
	nbMediumSel=0;
	//int Indexmusel = -1;

  	UInt_b(nMuon);
	VInt_b(Muon_charge);
  	VFloat_b(Muon_pt);
	VFloat_b(Muon_phi);
  	VFloat_b(Muon_eta);
  	VBool_b(Muon_tightId);
	VFloat_b(Muon_pfRelIso04_all);

  	for (unsigned int im=0; im<nMuon; ++im) 
	{
     		if (!Muon_tightId[im]) 			continue;
      		if (Muon_pfRelIso04_all[im]>isomuCut) 	continue;
      		if (Muon_pt[im]<ptmuCut) 		continue;
      		if (fabs(Muon_eta[im])>etamuCut) 	continue;

      		//printf("Muon_pt[%d]: %f, Muon_eta[%d]: %f\n", im, Muon_pt[im], im, Muon_eta[im]);
      		nmusel++;
		if (Indexmusel == -1) Indexmusel = im;

		//std::cout << " Item: " << im << " nmusel in the loop " << nmusel << " ID " << im << std::endl;
  	}
  	  //std::cout << " NmuoSEL when the loop is over " << nmusel << " ID " << Indexmusel << std::endl;
	if (nmusel!=1) return false; 
	//std::cout << " 1.NmuoSEL when the loop is over " << nmusel << " muon iniciales " << nMuon << std::endl;

	//int Indexelsel = -1;

  	UInt_b(nElectron);
	VInt_b(Electron_charge);
  	VFloat_b(Electron_pt);
	VFloat_b(Electron_phi);
  	VFloat_b(Electron_eta);
  	//VInt_b(Electron_cutBased); //old cut
	VBool_b(Electron_mvaFall17V1Iso_WP90); //I replace cutBased for this new cut
  	VFloat_b(Electron_pfRelIso03_all);

  	for (unsigned int ie=0; ie<nElectron; ++ie) 
	{
      		//if (Electron_cutBased[ie]<3) 		continue; // 0:fail, 1:veto, 2:loose, 3:medium, 4:tight

      		if (Electron_pfRelIso03_all[ie]>isoelCut) continue;
		if (!Electron_mvaFall17V1Iso_WP90[ie])	continue;
      		if (Electron_pt[ie]<ptelCut) 		continue;
      		if (fabs(Electron_eta[ie])>etaelCut) 	continue;
      		nelsel++;
		if (Indexelsel == -1) Indexelsel = ie;
		//std::cout << " Item: " << ie << " nelesel in the loop " << nelsel << " ID " << ie << std::endl;
  	}

	if (nelsel!=1) return false; 
	//std::cout << " NeleSEL when the loop is over " << nelsel << " electrones iniciales " << nElectron << std::endl;
	//Once I have selected both muon and electron, I want to check if the sum of their charges is neutral

	if(Electron_charge[Indexelsel]+Muon_charge[Indexmusel]!=0) return false;

	//unsigned int njtsel = 0;
  	UInt_b(nJet);
  	VFloat_b(Jet_pt);
  	VFloat_b(Jet_eta);
	VFloat_b(Jet_phi);
  	//VFloat_b(Jet_btagCSVV2);
	VFloat_b(Jet_btagDeepB);
  	for (unsigned int ij=0; ij<nJet; ++ij) 
	{
		//std::cout << "Jet pt[" << ij <<"] "<< Jet_pt[ij] <<  std::endl;
      		if (Jet_pt[ij]<ptjtCut) 	continue;
		//std::cout << "	Jet pt[" << ij <<"] "<< Jet_pt[ij] <<  std::endl;
      		if (fabs(Jet_eta[ij])>etajtCut) continue;

		double deltaR_MuoJet = ANGLES::DeltaR(Muon_phi[Indexmusel], Jet_phi[ij], Muon_eta[Indexmusel], Jet_eta[ij]);
		double deltaR_EleJet = ANGLES::DeltaR(Electron_phi[Indexelsel], Jet_phi[ij], Electron_eta[Indexelsel], Jet_eta[ij]);
		
		if(deltaR_MuoJet<0.4 || deltaR_EleJet<0.4) continue; 

      		njtsel++;

		//std::cout << "		nJet " << nJet << " nJetSel "  << njtsel << std::endl;
      		if (Jet_btagDeepB[ij]>btagMediumCut) 
		{
           		nbMediumSel++;
            		nbLooseSel++;
      		} 
		else if (Jet_btagDeepB[ij]>btagLooseCut) 
		{
           		 nbLooseSel++;
      		}
		
		//std::cout << " Item: " << ij << " njtsel in the loop " << nbMediumSel << " " << nbLooseSel << std::endl;
  	}

	//std::cout << " NJTSEL when the loop is over " << njtsel << " nJets iniciales "  << nJet << std::endl;
  	
  	//if (njtsel>=2) 		return false; 
  	//if (nbLooseSel<1) return false; 
  	//if (nbMediumSel>0) 	return false; 
  	//printf("nmu: %d, nel: %d\n", nmu, nel);

  	return true;

}
