//Code for LFV Analysis for 2018 CMS data
//Written by Juan Alcaraz and it had been modified by Lourdes Urda
//August 2019

#include "TString.h"
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
#include <typeinfo>

//#include <RooWorkspace.h>

#include "RooGlobalFunc.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"

#include "RooEffProd.h"
#include "RooNLLVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

#include "CMSAnalysis.h"

const double electronMass = 0.0005109989461; //GeV
const double muonMass = 0.1056583745;//GeV
const double METMass = 0.0;
const double tauonMass = 1.77686; //GeV

//using namespace RooFit;

struct TopAnalysis : public CMSAnalysis 
{
	TopAnalysis()
	{
		init(); 
	};
	
	virtual ~TopAnalysis(){};

	void init () 
	{
		nbLooseSel=0; nbMediumSel=0;			
	}	

	//Preselection
	bool Preselect();   
	
	//Selection
	bool Select(char* argv[]);

	//Full Selection
	bool FullSelection(char *argv[], double event_weight, SAMPLES& sample);

//Variables analisis

	//unsigned int njtsel;
 	unsigned int nbLooseSel;
      	unsigned int nbMediumSel;

	int Indexmusel;
	int Indexelsel;

	//char* argv[];

	std::vector<int> Indexjetsel = {0};
	std::vector<int> IndexMedjetsel = {0};
	std::vector<int> Indexbtagsel = {0};

	int nmusel;
	int nelsel;
	int njtsel;
	int nbtagsel;
	int nmuonveto;
	int nelectronveto;
	int ntauonveto;

	TString mydir = "/afs/cern.ch/work/l/lurda/CMS/May_2019/ExoticHiggsDecay/codigojuan/GitLab_NANOAOD/LFVAnalysis/nanoaod/";
	//TString dir = "/eos/user/l/lurda/CMS/SamplesAfterSkimmingNovember2019/"; 
	//TString dir2 = "/eos/user/c/cepeda/LFVNANO/Skimming3/";
	//TString dirmadgraph = "/eos/user/l/lurda/CMS/SamplesToMergeJanuary2020_madgraphMLM/";
	//TString diramcatnlo = "/eos/user/l/lurda/CMS/SamplesToMergeJanuary2020_amcatnloFXFX/";

// Select if there are at least two leptons wit some phase space cuts
	
  	const double ptmuCut = 24.;//25 //WE HAVE TO SET IT TO 26
  	const double etamuCut = 2.4;//2.5
  	const double isomuCut = 0.15;
  	const double ptelCut = 13.;//10
  	const double etaelCut = 2.4;
  	const double isoelCut = 0.06;//0.10
  	const double ptjtCut = 30.; //take it into account for Bjets//30
  	const double etajtCut = 4.7;//2.65
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
		cms.AddPlot1D("hMetPt",				"Missing E_{T} [GeV]", 			 100, 0., 200.);
  		cms.AddPlot1D("hMuonPt", 			"Muon pT [GeV]", 			 80, 20.0, 100.);
  		cms.AddPlot1D("hElectronPt", 			"Electron pT [GeV]", 			 150, 0.0, 150.);
  		cms.AddPlot1D("hJetPt1", 			"First Jet pT [GeV]", 			 250, 0.0, 250.);
  		cms.AddPlot1D("hJetPt2", 			"Second Jet pT [GeV]", 			 250, 0.0, 250.);
  		cms.AddPlot1D("hMuonPhi", 			"Muon Phi", 				 100, -3.15, 3.15);
  		cms.AddPlot1D("hMuonEta", 			"Muon Eta", 				 100, -3.0, 3.0);
 		cms.AddPlot1D("hElectronPhi", 			"Electron Phi", 			 100, -3.15, 3.15);
  		cms.AddPlot1D("hElectronEta", 			"Electron Eta", 			 100, -3.0, 3.0);
  		cms.AddPlot1D("hJetEta1", 			"First jet Eta", 		 	100, -5.0, 5.0);
  		cms.AddPlot1D("hJetEta2", 			"Second jet Eta", 		 	100, -5.0, 5.0);
		cms.AddPlot1D("hnVertex", 			"Number of Vertex", 			 100, 0, 100);
  		cms.AddPlot1D("hnMuon", 			"Number of Muons", 			 1, -0.5, 5.5);
  		cms.AddPlot1D("hnElectron", 			"Number of Electrons", 			 1, -0.5, 6.5);
  		cms.AddPlot1D("hnJet", 				"Number of Jets", 			 1, 0.5, 25.5);	
  		cms.AddPlot1D("hMetPhi", 			"MET Phi", 				 100, -3.15, 3.15);
		cms.AddPlot1D("hMuonElectronPhi", 		"Muon-Electron Phi", 			 100, -3.15, 3.15);
		cms.AddPlot1D("hMuonMetPhi", 			"Muon-MET Phi", 			 100, -3.15, 3.15);
		cms.AddPlot1D("hElectronMetPhi", 		"Electron-Met Phi", 			 100, -3.15, 3.15);
		//cms.AddPlot1D("hbtagDeepB",			"b-Tag jets Deep B",			 50, 0, 1);
		cms.AddPlot1D("hInvariantMassMuonElectron", 	"Invariant Mass Muon-Electron [GeV]", 	 100, 0, 200);
		cms.AddPlot1D("hTransverseMassMuonMET", 	"Transverse Mass Muon-MET [GeV]", 	 100, 0, 200);
		cms.AddPlot1D("hTransverseMassElectronMET", 	"Transverse Mass Electron-MET [GeV]", 	 100, 0, 200);
		cms.AddPlot1D("hnbLooseSel", 			"Number of BJets Loose", 		 1, -0.5, 4.5);
		cms.AddPlot1D("hnbMediumSel", 			"Number BJets Medium", 			 1, -0.5, 4.5);
		cms.AddPlot1D("hdeltaR", 			"#DeltaR Muon-Electron", 		 50, -2.4, 10.);
		cms.AddPlot1D("hCollinearMass", 		"Collinear Mass Muon-Electron [GeV]", 	 360, 0, 360);

		cms.AddPlot1D("hCollinearMass_0jets", 		"Collinear Mass Muon-Electron [GeV]", 	 360, 0, 360);
		cms.AddPlot1D("hCollinearMass_1jet", 		"Collinear Mass Muon-Electron [GeV]", 	 360, 0, 360);
		cms.AddPlot1D("hCollinearMass_2jetsVBF", 	"Collinear Mass Muon-Electron [GeV]", 	 360, 0, 360);
		cms.AddPlot1D("hCollinearMass_2jetsGGH", 	"Collinear Mass Muon-Electron [GeV]", 	 360, 0, 360);

		cms.AddPlot1D("mjj_beforecuts", 		"Jets invariant mass (before cuts)", 	 600, 0, 600);
		cms.AddPlot1D("mjj_GGHcuts", 			"Jets invariant mass (GGH)", 	 	 600, 0, 600);
		cms.AddPlot1D("mjj_VBFcuts", 			"Jets invariant mass (VBF)", 	 	 600, 0, 600);

		cms.AddPlot1D("hnMuoSel", 			"Number of selected Muons", 		 10, 0, 10);
		cms.AddPlot1D("hnEleSel", 			"Number of selected Electrons", 	 10, 0, 10);
		cms.AddPlot1D("hnJetSel", 			"Number of selected jets", 		 1, -1, 10);
		cms.AddPlot1D("hSumMuonElectronPt", 		"Sum Muon Electorn pt", 		 200, -0.5, 200.5);

/*		cms.AddPlot1D("hgenMZFound",			"genM when Z boson is found", 		 1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZFoundElectronPair",	"genM (electron) Z boson found", 	 1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZFoundMuonPair",		"genM (muon) Z boson found", 	 	 1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZFoundTauonPair",		"genM (tauon) Z boson found", 	 	 1501, -0.5, 1500.5);

		cms.AddPlot1D("hgenMZNOTFoundElectronPairFound", 	"genM when electron pair is found",      1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZNOTFoundMuonPairFound", 		"genM when muon pair is found",      	 1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZNOTFoundTauonPairFound", 		"genM when tauon pair is found",      	 1501, -0.5, 1500.5);
	
		cms.AddPlot1D("hgenMZFoundCorrected",		 		"genM when Z boson is found (Corrected)", 	 1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZFoundElectronPairCorrected",		"genM (electron) Z boson found (Corrected)", 	 1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZFoundMuonPairCorrected",	 		"genM (muon) Z boson found (Corrected)", 	 1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZFoundTauonPairCorrected",	 		"genM (tauon) Z boson found (Corrected)", 	 1501, -0.5, 1500.5);

		cms.AddPlot1D("hgenMZNOTFoundElectronPairFoundCorrected", 	"genM when electron pair is found (Corrected)", 1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZNOTFoundMuonPairFoundCorrected", 	 	"genM when muon pair is found (Corrected)",     1501, -0.5, 1500.5);
		cms.AddPlot1D("hgenMZNOTFoundTauonPairFoundCorrected", 	 	"genM when tauon pair is found (Corrected)",    1501, -0.5, 1500.5);

		cms.AddPlot1D("hgenpTZFound",					"genpT when Z boson is found", 		 	101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZFoundElectronPair",			"genpT (electron) Z boson found", 	 	101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZFoundMuonPair",				"genpT (muon) Z boson found", 	 	 	101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZFoundTauonPair",				"genpT (tauon) Z boson found", 	 	 	101, -0.5, 100.5);

		cms.AddPlot1D("hgenpTZNOTFoundElectronPairFound", 		"genpT when electron pair is found",      	101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZNOTFoundMuonPairFound", 			"genpT when muon pair is found",      	 	101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZNOTFoundTauonPairFound", 			"genpT when tauon pair is found",      	 	101, -0.5, 100.5);
	
		cms.AddPlot1D("hgenpTZFoundCorrected",		 		"genpT when Z boson is found (Corrected)", 	 101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZFoundElectronPairCorrected",		"genpT (electron) Z boson found (Corrected)", 	 101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZFoundMuonPairCorrected",	 		"genpT (muon) Z boson found (Corrected)", 	 101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZFoundTauonPairCorrected",	 		"genpT (tauon) Z boson found (Corrected)", 	 101, -0.5, 100.5);

		cms.AddPlot1D("hgenpTZNOTFoundElectronPairFoundCorrected", 	"genpT when electron pair is found (Corrected)", 101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZNOTFoundMuonPairFoundCorrected", 	 	"genpT when muon pair is found (Corrected)",     101, -0.5, 100.5);
		cms.AddPlot1D("hgenpTZNOTFoundTauonPairFoundCorrected", 	"genpT when tauon pair is found (Corrected)",    101, -0.5, 100.5);

*/	

		//cms.AddPlot1D("hLHE_Njets", 			"LHE_Njets variable", 		 	 11, -0.5, 10.5);
		
	}
	static void Fill(TopAnalysis &cms, const SAMPLES &sample, const double event_weight)
	{
		Float_b(MET_pt); 
	    	VFloat_b(Electron_pt); 
            	VFloat_b(Muon_pt); 	
            	VFloat_b(Jet_pt);
	    	Int_b(PV_npvs);
    	    	Int_b(nMuon);	
    	    	Int_b(nElectron);
    	    	Int_b(nJet); 	
           	Float_b(MET_phi);
            	VFloat_b(Electron_phi);	
            	VFloat_b(Muon_phi);	
            	VFloat_b(Electron_eta);	
            	VFloat_b(Muon_eta);
	        VFloat_b(Jet_eta); 
		VFloat_b(Jet_btagDeepB);											 							 			VFloat_b(Muon_pfRelIso04_all);
		//UChar_b(LHE_Njets);

		//cms.FillPlot1D("hLHE_Njets", sample, LHE_Njets, event_weight);
		cms.FillPlot1D("hMetPt", sample ,  MET_pt, event_weight);
		cms.FillPlot1D("hElectronPt", sample, Electron_pt[cms.Indexelsel], event_weight);
		cms.FillPlot1D("hMuonPt", sample,  Muon_pt[cms.Indexmusel], event_weight);

		if(cms.Indexjetsel.size()==1) 
		{
			cms.FillPlot1D("hJetPt1", sample, Jet_pt[cms.Indexjetsel.at(0)], event_weight);
			cms.FillPlot1D("hJetEta1", sample, Jet_eta[cms.Indexjetsel.at(0)], event_weight);
		}
		if(cms.Indexjetsel.size()==2)
		{
			cms.FillPlot1D("hJetEta1", sample, Jet_eta[cms.Indexjetsel.at(0)], event_weight);
			cms.FillPlot1D("hJetEta2", sample, Jet_eta[cms.Indexjetsel.at(1)], event_weight);				
			cms.FillPlot1D("hJetPt1", sample, Jet_pt[cms.Indexjetsel.at(0)], event_weight);	
			cms.FillPlot1D("hJetPt2", sample, Jet_pt[cms.Indexjetsel.at(1)], event_weight);		  		  	
		}

		cms.FillPlot1D("hnVertex", sample, PV_npvs, event_weight);
		cms.FillPlot1D("hnMuon", sample, nMuon, event_weight);
		cms.FillPlot1D("hnElectron", sample, nElectron, event_weight);
		cms.FillPlot1D("hnJet", sample, nJet, event_weight);//nJet
		cms.FillPlot1D("hMetPhi", sample, MET_phi, event_weight);
		cms.FillPlot1D("hElectronPhi", sample, Electron_phi[cms.Indexelsel], event_weight);
		cms.FillPlot1D("hMuonPhi", sample, Muon_phi[cms.Indexmusel], event_weight);
		cms.FillPlot1D("hElectronEta", sample, Electron_eta[cms.Indexelsel], event_weight);
		cms.FillPlot1D("hMuonEta", sample, Muon_eta[cms.Indexmusel], event_weight);

		cms.FillPlot1D("hMuonElectronPhi", sample, ANGLES::DeltaPhi(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel]), event_weight);
		cms.FillPlot1D("hMuonMetPhi", sample, ANGLES::DeltaPhi(Muon_phi[cms.Indexmusel], MET_phi), event_weight);
		cms.FillPlot1D("hElectronMetPhi", sample, ANGLES::DeltaPhi(Electron_phi[cms.Indexelsel], MET_phi), event_weight);
		//cms.FillPlot1D("hbtagDeepB", sample, Jet_btagDeepB[cms.Indexjetsel.at(0)], event_weight);

		TLorentzVector LorentzElec;
		TLorentzVector LorentzMuon;
		TLorentzVector LorentzMET;

		LorentzElec.SetPtEtaPhiM(Electron_pt[cms.Indexelsel],Electron_eta[cms.Indexelsel],Electron_phi[cms.Indexelsel], electronMass);
		LorentzMuon.SetPtEtaPhiM(Muon_pt[cms.Indexmusel], Muon_eta[cms.Indexmusel], Muon_phi[cms.Indexmusel], muonMass);
		LorentzMET.SetPtEtaPhiM(MET_pt, 0.0, MET_phi, METMass);

		cms.FillPlot1D("hInvariantMassMuonElectron", sample, FOURVECTORS::InvariantMass(LorentzMuon, LorentzElec), event_weight);
		cms.FillPlot1D("hTransverseMassMuonMET", sample, FOURVECTORS::TranverseMass(LorentzMuon, LorentzMET), event_weight);
		cms.FillPlot1D("hTransverseMassElectronMET", sample, FOURVECTORS::TranverseMass(LorentzElec, LorentzMET), event_weight);
		cms.FillPlot1D("hnbLooseSel", sample, cms.nbLooseSel,  event_weight);
		cms.FillPlot1D("hnbMediumSel", sample, cms.nbMediumSel, event_weight);

		cms.FillPlot1D("hdeltaR", sample, ANGLES::DeltaR(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel], Muon_eta[cms.Indexmusel], Electron_eta[cms.Indexelsel]), event_weight);

		cms.FillPlot1D("hCollinearMass", sample, FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[cms.Indexelsel], MET_phi)), event_weight);

		cms.FillPlot1D("hnMuoSel", sample, cms.nmusel, event_weight);
		cms.FillPlot1D("hnEleSel", sample, cms.nelsel, event_weight);
		cms.FillPlot1D("hnJetSel", sample, cms.njtsel, event_weight);
		cms.FillPlot1D("hSumMuonElectronPt", sample, (LorentzMuon+LorentzElec).Pt() , event_weight);

	}

	static void Save(TopAnalysis &cms, const SAMPLES &sample, ASCII &txtFiles, const string& charge, const string& muon, const string& jets, const string& bjets)
	{
		TString option = "RECREATE"; 

		for(auto &name: txtFiles.HistogramsList)
		{
			std::cout << charge.c_str() << " " << muon.c_str()<< " " << jets.c_str()<< " " << bjets.c_str() << std::endl;
			cms.SavingHistograms(sample, name, option, charge, muon, jets, bjets);
			if(option =="RECREATE") option = "UPDATE";
		}
	}

	static void Draw(CMSAnalysis &cms, ASCII &txtFiles, const string& charge, const string& muon, const string& jets, const string& bjets)
	{
		for(auto &name: txtFiles.HistogramsList)
		{
			std::cout << "Name of Histograms " << name << std::endl;

		 	cms.DrawPlot1D(name, "", Form("png_%s_%s_%s_%s", charge.c_str(), muon.c_str(), jets.c_str(), bjets.c_str() ));
			std::cout << "Name of Histograms que se ha ploteado" << name << std::endl;

		}
	}

};

struct SCALEFACTORS : public CMSAnalysis
{

	static TH1D* GettingTH1Histograms(TopAnalysis &cms, TH1D* Histogram, TString& pathroot, TString& name)
	{
		Histogram = nullptr;
		Histogram = cms.ReadingFileAndGettingTH1Histogram(cms.mydir+pathroot, name);
		return Histogram;
	}
	static TH2D* GettingTH2Histograms(TopAnalysis &cms, TH2D* Histogram, TString& pathroot, TString& name)
	{
		Histogram = nullptr;
		Histogram = cms.ReadingFileAndGettingTH2Histogram(cms.mydir+pathroot, name);
		return Histogram;
	}
	static double MuonSF(TopAnalysis &cms, TH2D* MuonEfficiencyHist)
	{
		VFloat_b(Muon_pt);
  		VFloat_b(Muon_eta);

		double Muon_weight;

		if(Muon_pt[cms.Indexmusel] > 119.) 
		{
			Muon_weight = cms.ScaleFactors(MuonEfficiencyHist, 119., fabs(Muon_eta[cms.Indexmusel]));
		}
		else 
		{
			Muon_weight = cms.ScaleFactors(MuonEfficiencyHist, Muon_pt[cms.Indexmusel], fabs(Muon_eta[cms.Indexmusel]));
		}

		if(Muon_weight == 0) std::cout <<  Muon_pt[cms.Indexmusel] << " eta " << fabs(Muon_eta[cms.Indexmusel]) << std::endl;
		return Muon_weight;

	}
	static double ElectronSF(TopAnalysis &cms, TH2D* ElectronEfficiencyHist)
	{
		VFloat_b(Electron_pt);
		VFloat_b(Electron_eta);

		double Electron_weight = cms.ScaleFactors(ElectronEfficiencyHist, Electron_pt[cms.Indexelsel], Electron_eta[cms.Indexelsel]);

		return Electron_weight;

	}
	static double TriggerSF(TopAnalysis &cms, TH2D* TriggerHist)
	{
		VFloat_b(Muon_pt);
  		VFloat_b(Muon_eta);

		double Trigger_weight	= cms.ScaleFactors(TriggerHist, Muon_pt[cms.Indexmusel], fabs(Muon_eta[cms.Indexmusel]));
		if(cms.ptmuCut == 24 && Trigger_weight == 0) Trigger_weight = 1.0;
	 //	std::cout << "Muon_pt " << Muon_pt[cms.Indexmusel] << " Muon_eta " << Muon_eta[cms.Indexmusel] << std::endl;
		return Trigger_weight;
	}
	static double PileupSF(TopAnalysis &cms, TH1D* Ratio)
	{
		Float_b(Pileup_nTrueInt);

		double pileup_weight = cms.PileupReweighting(Ratio, Pileup_nTrueInt);
		
		return pileup_weight;
	}
	static double GeneratorWeightSF(TopAnalysis &cms)
	{
		Float_b(Generator_weight);

		double generator_weight = Generator_weight;	
		
		return generator_weight;
	}
	static double BtaggingSF(TopAnalysis &cms, char nbjets)
	{
		VFloat_b(Jet_pt) ;
		VInt_b(Jet_hadronFlavour);

		int inbjets = (int)nbjets - 48; //in ASCII code, the numbers (digits) start from 48

		double bTag_weight=1.;

		if (cms.IndexMedjetsel.size()==0){
			bTag_weight = cms.bTagEventWeight(cms.IndexMedjetsel.size(),-1,-1,-1,-1, "comb", "central" , inbjets , "medium");
		}

		else if (cms.IndexMedjetsel.size()==1){
			bTag_weight = cms.bTagEventWeight(cms.IndexMedjetsel.size(), Jet_pt[cms.IndexMedjetsel.at(0)],Jet_hadronFlavour[cms.IndexMedjetsel.at(0)],-1,-1, "comb", "central" , inbjets , "medium");
		}
		else if (cms.IndexMedjetsel.size()>=2){
			bTag_weight = cms.bTagEventWeight(cms.IndexMedjetsel.size(), Jet_pt[cms.IndexMedjetsel.at(0)],Jet_hadronFlavour[cms.IndexMedjetsel.at(0)],Jet_pt[cms.IndexMedjetsel.at(1)], 						Jet_hadronFlavour[cms.IndexMedjetsel.at(1)], "comb", "central" , inbjets , "medium");
		}

		return bTag_weight;
	}

	static double MergingSF(TopAnalysis &cms, TString sample, std::map<TString,double> merge)
	{
		UChar_b(LHE_Njets);
		double merge_weight;
		int numberofjets;
		TString samplename, sufix;

		if(sample.Contains("DY")) samplename = "DY";
		else if(sample.Contains("WJets")) samplename = "WJets";

		if(sample.Contains("ZTauTau")) sufix = "_ZTauTau";
                else if(sample.Contains("ZMuoMuo")) sufix = "_ZMuoMuo";
                else if(sample.Contains("ZEleEle")) sufix = "_ZEleEle";
                else sufix = "";
		
//		std::cout << sample << std::endl;

		if(sample.Contains("Inclusive"))
		{
			numberofjets = (int)LHE_Njets;
			//std::cout << "Sample " << sample << " Numberofjets " << numberofjets << std::endl;
			if(numberofjets==0) merge_weight = merge[sample];
			if(numberofjets==1) merge_weight = merge[samplename+"1Merge"+sufix];
			if(numberofjets==2) merge_weight = merge[samplename+"2Merge"+sufix];
			if(numberofjets==3) merge_weight = merge[samplename+"3Merge"+sufix];
			if(numberofjets==4) merge_weight = merge[samplename+"4Merge"+sufix];
			//else merge_weight = 1.0;
		}
		else merge_weight=merge[sample];

		if(merge_weight==0) std::cout << sample << " STOP WEIGHT IS 0 BECAUSE YOU DIDN'T ADD THE REST OF DY SAMPLES! " << merge_weight << std::endl;

		return merge_weight;

	}

        typedef struct
        {
               	bool TauFlag = false;
                bool EleFlag = false;
                bool MuoFlag = false;
        	double weight = 1.0;

        } drellyan_info;

	static drellyan_info ZptReweighting(TopAnalysis &cms, TH2D* ZptHisto ,const SAMPLES &sample)
	{

		drellyan_info x;
		Int_b(nGenPart);
		VFloat_b(GenPart_eta);
		VFloat_b(GenPart_mass);
		VInt_b(GenPart_genPartIdxMother);
		VInt_b(GenPart_pdgId);
		VInt_b(GenPart_status);
		VFloat_b(GenPart_phi);
		VFloat_b(GenPart_pt);

		double genM;
		double genpT;
		double genMlepton;
		double genpTlepton;

		int indexZ = -1;
		int nZFound = 0;

		bool ZFound = false;

		std::vector<int> selectedlepton;

		//std::cout << "new event" << std::endl;
		//std::cout << " Event, nGenPart, particle, pdgID, status, IdxMother, pt, mass, eta, phi " << std::endl;

		for(int particle = 0; particle<nGenPart; particle++)
		{
			//std::cout <<  nGenPart << " " << particle << " " << GenPart_pdgId[particle] << " " << GenPart_status[particle] << " " <<GenPart_genPartIdxMother[particle] << " " << GenPart_pt[particle] << " " << GenPart_mass[particle]  << " " << GenPart_eta[particle] << " " << GenPart_phi[particle] <<  std::endl;

			if(GenPart_pdgId[particle] == 23 && GenPart_status[particle] == 62)
			{
				ZFound = true; nZFound ++;
				if(indexZ == -1) indexZ = particle;
			}

			if(ZFound)
			{
				if(abs(GenPart_pdgId[particle]) != 13 && abs(GenPart_pdgId[particle]) != 15 && abs(GenPart_pdgId[particle]) != 11) continue;
				if(GenPart_genPartIdxMother[particle] == -1) continue;
				if(GenPart_pdgId[GenPart_genPartIdxMother[particle]] == 23 && GenPart_status[GenPart_genPartIdxMother[particle]] == 62)
				{
					//std::cout << "Lepton found after ZFound " << GenPart_pdgId[particle] << std::endl;
					selectedlepton.push_back(particle);
				}
			}

			//if(ZFound && selectedleptons.size()<2) std::cout << "ZFound but not leptons found" <<std::endl;
			//if(ZFound && selectedlepton.size()==2) break;
			else if(!ZFound)
			{
				//Cuts when it doesn't find a Z boson and has to find pair of leptons 
				if(abs(GenPart_pdgId[particle]) != 13 && abs(GenPart_pdgId[particle]) != 15 && abs(GenPart_pdgId[particle]) != 11) continue;
				if(GenPart_status[particle]<0) continue;//| GenPart_status[particle]==1 || GenPart_status[particle]==2
                        	if(GenPart_status[particle]!=23) continue;
				if(GenPart_status[particle]==1) continue; //final state in pythia
				if(GenPart_status[particle]==2) continue; //intermediate state in pythia
				if(abs(GenPart_pdgId[particle]) > 20) continue;
		        	if(GenPart_genPartIdxMother[particle] != 0) continue;

                        	//std::cout << "it is still finding particles " << GenPart_pdgId[particle] << std::endl;

				selectedlepton.push_back(particle);
			}

  		}
		//std::cout << " ZFound " << ZFound << " number of leptons found " << selectedlepton.size() << std::endl;
		//if(ZFound && selectedlepton.size()<2) std::cout << "ZFound but not a pair of leptons " << selectedlepton.size() << std::endl;
		if (nZFound>1) std::cout << "More than one Z found in the event with 62 as final status " << nZFound << std::endl;
		//if(!ZFound && selectedlepton.size()==2) std::cout << "lepton pair FOUND "<< selectedlepton.size()<< " selectedlepton[0] " << GenPart_pdgId[selectedlepton[0]] <<" selectedlepton[1] " << GenPart_pdgId[selectedlepton[1]]<< std::endl;

		int ID = abs(GenPart_pdgId[selectedlepton.at(0)]);
                double leptonMass = 0;
                if(ID==15) {leptonMass = tauonMass; x.TauFlag = true;}
                else if(ID==13) {leptonMass = muonMass; x.MuoFlag = true;}
                else if(ID==11) {leptonMass = electronMass; x.EleFlag = true;}

		if(ZFound)
		{
			genM = GenPart_mass[indexZ];
			genpT = GenPart_pt[indexZ];
		}
		else if(!ZFound)
		{
			double leptonpt = GenPart_pt[selectedlepton.at(0)];
			double antileptonpt = GenPart_pt[selectedlepton.at(1)];

			double leptoneta = GenPart_eta[selectedlepton.at(0)];
			double antileptoneta = GenPart_eta[selectedlepton.at(1)];

			double leptonphi = GenPart_phi[selectedlepton.at(0)];
			double antileptonphi = GenPart_phi[selectedlepton.at(1)];

			TLorentzVector LorentzLepton; LorentzLepton.SetPtEtaPhiM(leptonpt, leptoneta , leptonphi, leptonMass);
			TLorentzVector LorentzAntilepton; LorentzAntilepton.SetPtEtaPhiM(antileptonpt, antileptoneta, antileptonphi, leptonMass);

			genM = (LorentzLepton+LorentzAntilepton).M();
			genpT = (LorentzLepton+LorentzAntilepton).Pt();
		}

		if(genM>1000 || genpT >1000) x.weight == 1.0;
		else x.weight = cms.ScaleFactors(ZptHisto, genM, genpT);

		return x;
	}

	static double QCDEstimation(TopAnalysis &cms, RooWorkspace *w)
	{
		VFloat_b(Muon_pt);
  		VFloat_b(Muon_eta);
		VFloat_b(Electron_pt);
		VFloat_b(Electron_eta);
		VFloat_b(Electron_phi);
		VFloat_b(Muon_phi);
    	    	Int_b(nJet);

		double osss_weight;

		double dR = ANGLES::DeltaR(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel], Muon_eta[cms.Indexmusel], Electron_eta[cms.Indexelsel]);

		osss_weight = cms.QCDEstimationFunction(w, nJet, dR, Electron_pt[cms.Indexelsel], Muon_pt[cms.Indexmusel]); 

		return osss_weight;
	}
};

int main(int argc, char* argv[])//"OS", "IM", "0jets", "0bjets"
{

	// Initialize analysis structure
  	TopAnalysis cms; 

  	int maxevents = -1;
	SAMPLES::TotalAnalysisLumi = 59.266425103E3; //pb^-1

	//----- ADDING DATA SAMPLES
	cms.AddSample("SamplesDirectoryV6/Data/SingleMuRun2018A_V6.txt");
	cms.AddSample("SamplesDirectoryV6/Data/SingleMuRun2018B_V6.txt");
	cms.AddSample("SamplesDirectoryV6/Data/SingleMuRun2018C_V6.txt");
	cms.AddSample("SamplesDirectoryV6/Data/SingleMuRun2018D_V6.txt");

	//----- ADDING SIGNAL SAMPLES
	cms.AddSample("SamplesDirectoryV6/GluGlu/GluGlu_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunI0-v1_NANOAODSIM_V6.txt");
	cms.AddSample("SamplesDirectoryV6/VBF/VBF_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAu0-v1_NANOAODSIM_V6.txt");

	cms.AddSample("SamplesDirectoryV6/GluGlu/GluGluHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nan0-v1_NANOAODSIM.txt");
        cms.AddSample("SamplesDirectoryV6/GluGlu/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn0-v1_NANOAODSIM.txt");
        cms.AddSample("SamplesDirectoryV6/VBF/VBFHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nano251-v1_NANOAODSIM.txt");
        cms.AddSample("SamplesDirectoryV6/VBF/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn18N0-v1_NANOAODSIM.txt");

	//------ADDING BCKGRND SAMPLES
	cms.AddSample("SamplesDirectoryV6/W/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv0-v1_NANOAODSIM_V6.txt");
	cms.AddSample("SamplesDirectoryV6/W/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv0-v1_NANOAODSIM.txt");

		//Samples to Merge MADGRAPH

		//-WJETS
		cms.AddSample("SamplesDirectoryV6/WJets/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv6-0-v1_NANOAODSIM_V6.txt");
		cms.AddSample("SamplesDirectoryV6/WJets/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/WJets/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/WJets/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/WJets/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");

		//-DY
/*		cms.AddSample("SamplesDirectoryV6/DYmadgraph/DYJetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_RunIIAutum2019.txt");
		cms.AddSample("SamplesDirectoryV6/DYmadgraph/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/DYmadgraph/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/DYmadgraph/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/DYmadgraph/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/DYamcatnlo/DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019.txt");
*/

                cms.AddSample("SamplesDirectoryV6/ZTauTauDY/ZTauTau_DYJetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_RunIIAutum2019.txt");
                cms.AddSample("SamplesDirectoryV6/ZTauTauDY/ZTauTau_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZTauTauDY/ZTauTau_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZTauTauDY/ZTauTau_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZTauTauDY/ZTauTau_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZTauTauDY/ZTauTau_DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019.txt");

                cms.AddSample("SamplesDirectoryV6/ZMuoMuoDY/ZMuoMuo_DYJetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_RunIIAutum2019.txt");
                cms.AddSample("SamplesDirectoryV6/ZMuoMuoDY/ZMuoMuo_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZMuoMuoDY/ZMuoMuo_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZMuoMuoDY/ZMuoMuo_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZMuoMuoDY/ZMuoMuo_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZMuoMuoDY/ZMuoMuo_DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019.txt");

                cms.AddSample("SamplesDirectoryV6/ZEleEleDY/ZEleEle_DYJetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_RunIIAutum2019.txt");
                cms.AddSample("SamplesDirectoryV6/ZEleEleDY/ZEleEle_DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZEleEleDY/ZEleEle_DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZEleEleDY/ZEleEle_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZEleEleDY/ZEleEle_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
                cms.AddSample("SamplesDirectoryV6/ZEleEleDY/ZEleEle_DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019.txt");

		//Samples to Merge AMCATNLO

		//-WJETS
		/*cms.AddSample("SamplesDirectoryV6/INCLUSIVEISMISSING");*/
		/*cms.AddSample("SamplesDirectoryV6/WJets/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/WJets/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/WJets/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		//-DY
		cms.AddSample("SamplesDirectoryV6/DYamcatnlo/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAODv6_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/DYamcatnlo/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/DYamcatnlo/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/DYamcatnlo/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");*/

	cms.AddSample("SamplesDirectoryV6/tW/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIAutum1-v1_NANOAODSIM_V6.txt");
	cms.AddSample("SamplesDirectoryV6/tW/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIA1-v1_NANOAODSIM_V6.txt");

	cms.AddSample("SamplesDirectoryV6/TTbar/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/TTbar/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019.txt");

        cms.AddSample("SamplesDirectoryV6/Z/ZZ_TuneCP5_13TeV-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_u0-v1_NANOAODSIM_V6.txt");
       	cms.AddSample("SamplesDirectoryV6/W/WZ_TuneCP5_13TeV-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_u0-v1_NANOAODSIM_V6.txt");

	cms.AddSample("SamplesDirectoryV6/W/WminusHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nan0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/W/WplusHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nano0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/Z/ZHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nano25Oc0-v1_NANOAODSIM.txt");

	cms.AddSample("SamplesDirectoryV6/W/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv6-Na0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/EWK/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_RunIIAutum0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/EWK/EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/EWK/EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8_RunIIA0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/EWK/EWKZ2Jets_ZToNuNu_TuneCP5_PSweights_13TeV-madgraph-pythia8_RunIIAutu0-v1_NANOAODSIM.txt");

	//OLD WAY OF ADDING SAMPLES
  	/*cms.AddFiles("Data", "DataA", 14.027047499*1000., -1, cms.dir+"SingleMuon_Run2018A-Nano25Oct2019-v1_NANOAOD/", 555, -1, 227489240);///555//   Units = [1/picobarn]
  	cms.AddFiles("Data", "DataB", 7.060622497*1000.,  -1, cms.dir+"SingleMuon_Run2018B-Nano25Oct2019-v1_NANOAOD/", 264, -1, 110446445);//264
  	cms.AddFiles("Data", "DataC", 6.894770971*1000.,  -1, cms.dir+"SingleMuon_Run2018C-Nano25Oct2019-v1_NANOAOD/", 264, -1, 107972995);//264
  	cms.AddFiles("Data", "DataD", 31.283984136*1000., -1, cms.dir+"DeMaria_SingleMuon_Run2018D-Nano25Oct2019-v1_NANOAOD/", 1238, -1, -1);//1238
  	cms.AddFiles("MCSignal", "GG",  -1, 48.58*0.1,  cms.dir+"GluGlu_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunI0-v1_NANOAODSIM/", 1,  maxevents, 42945741.6072);
 	cms.AddFiles("MCSignal", "VBF", -1, 3.782*0.1,  cms.dir+"VBF_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAu0-v1_NANOAODSIM/", 1,  maxevents, 7631474.65134);	
  	cms.AddFiles("MCBckgr", "WW",    -1, 12.15488,  cms.dir+"WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv0-v1_NANOAODSIM/",  5,  maxevents, 85917643.789);
  	cms.AddFiles("MCBckgr", "WJets", -1, 52940.0,   cms.dirmadgraph+"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv6-0-v1_NANOAODSIM/", 11,  maxevents, 70389866.8084);//61526//11
	cms.AddFiles("MCBckgr", "WJets1", -1, 8104.0,   cms.dirmadgraph+"W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM/", 9, maxevents, 51041402.2428);//9
	cms.AddFiles("MCBckgr", "WJets2", -1, 2793.0,   cms.dirmadgraph+"W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM/", 8, maxevents, 23268796.8198);//8
	cms.AddFiles("MCBckgr", "WJets3", -1, 992.5 ,   cms.dirmadgraph+"W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM/", 4, maxevents, 14489898.9086);//4
	cms.AddFiles("MCBckgr", "WJets4", -1, 544.3 ,   cms.dirmadgraph+"W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM/", 3, maxevents, 10059700.497);//3
	cms.AddFiles("MCBckgr", "WJets0", -1, 544.3 ,   cms.diramcatnlo+"WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM/", 3, maxevents, 10059700.497);
	cms.AddFiles("MCBckgr", "WJets1", -1, 544.3 ,   cms.diramcatnlo+"WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM/", 3, maxevents, 10059700.497);
	cms.AddFiles("MCBckgr", "WJets2", -1, 544.3 ,   cms.diramcatnlo+"WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM/", 3, maxevents, 10059700.497);
  	cms.AddFiles("MCBckgr", "TW",    -1, 35.85,     cms.dir+"ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIAutum1-v1_NANOAODSIM/",  3,  maxevents, 334874732.0);
  	cms.AddFiles("MCBckgr", "TbarW", -1, 35.85,     cms.dir+"ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIA1-v1_NANOAODSIM/",  4,  maxevents, 266470418.054);
	NOUSARcms.AddFiles("MCBckgr", "DY",   -1, 2075.14*3,  cms.dir2+"DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019/",          29,  maxevents, 3.298665625309500e+12);
 	cms.AddFiles("MCBckgr", "DY",   -1,  6435.0,    cms.dir+"DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019/",    30,  maxevents, 3.44609946801E+12);//30			
 	cms.AddFiles("MCBckgr", "DY",   -1,  6435.0,    cms.dir+"DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019/",    30,  maxevents, 8.94044e+06);//30			
  	cms.AddFiles("MCBckgr", "DY0",  -1, 2075.14*3,  cms.diramcatnlo+"DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM/", 13,  maxevents, 5.68367119571e+11);
  	cms.AddFiles("MCBckgr", "DY1",  -1, 2075.14*3,  cms.diramcatnlo+"DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM/", 15,  maxevents, 4.08074701113e+11);
  	cms.AddFiles("MCBckgr", "DY2",  -1, 2075.14*3,  cms.diramcatnlo+"DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM/", 13,  maxevents, 1.62847517628e+11);
	cms.AddFiles("MCBckgr", "TTbar", -1, 88.29,    cms.dir+"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019/", 13,  maxevents, 4622080234.06);//13
	cms.AddFiles("MCBckgr", "WZ",    -1, 51.11,      cms.dir+"WZ_TuneCP5_13TeV-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_u0-v1_NANOAODSIM/",  4,  maxevents, 3884167.00546);//27.6
  	cms.AddFiles("MCBckgr", "ZZ",    -1, 16.91,     cms.dir+"ZZ_TuneCP5_13TeV-pythia */

	//Defining the lists of samples to apply the merging method
	std::vector<TString> ListWJetsToMerge {"WJetsInclusiveMerge", "WJets1Merge", "WJets2Merge", "WJets3Merge", "WJets4Merge"};
	std::vector<TString> ListDYToMerge {"DYInclusiveMerge", "DY1Merge", "DY2Merge", "DY3Merge", "DY4Merge"};
        std::vector<TString> ListDYToMerge_ZTauTau {"DYInclusiveMerge_ZTauTau", "DY1Merge_ZTauTau", "DY2Merge_ZTauTau", "DY3Merge_ZTauTau", "DY4Merge_ZTauTau"};
        std::vector<TString> ListDYToMerge_ZMuoMuo {"DYInclusiveMerge_ZMuoMuo", "DY1Merge_ZMuoMuo", "DY2Merge_ZMuoMuo", "DY3Merge_ZMuoMuo", "DY4Merge_ZMuoMuo"};
        std::vector<TString> ListDYToMerge_ZEleEle {"DYInclusiveMerge_ZEleEle", "DY1Merge_ZEleEle", "DY2Merge_ZEleEle", "DY3Merge_ZEleEle", "DY4Merge_ZEleEle"};

	//Defining the maps for merged samples that contain the information about the name of the sample and their corresponding merging weight depending on the number of jets
	std::map<TString,double> WJetsMergingMap = cms.MergingMCSamples(ListWJetsToMerge);
	std::map<TString,double> DYMergingMap = cms.MergingMCSamples(ListDYToMerge);
	std::map<TString,double> DYMergingMap_ZTauTau = cms.MergingMCSamples(ListDYToMerge_ZTauTau);
        std::map<TString,double> DYMergingMap_ZMuoMuo = cms.MergingMCSamples(ListDYToMerge_ZMuoMuo);
        std::map<TString,double> DYMergingMap_ZEleEle = cms.MergingMCSamples(ListDYToMerge_ZEleEle);
  	// Initialize 1D histograms
	HISTS::Add(cms);

  	// Loop on samples
  	unsigned int nsamples = cms.GetNumberOfSamples();

	ASCII txtFiles;

	txtFiles.txtCounter("NameOfHistograms.txt");
	//txtFiles.txtCounter("zpt_checkingplots.txt");

	//---PILE UP REWEIGHTING HISTOGRAM---
	TH1D* Ratio = nullptr;
	Ratio = cms.ReadingFileAndGettingTH1Histogram("python/pyroot/pileupweights.root", "RatioPU");

	//---MUON TIGHT ID EFFICIENCY HISTOGRAM---
	TH2D* MuonTightIDEfficiencyHist = nullptr;
	MuonTightIDEfficiencyHist = cms.ReadingFileAndGettingTH2Histogram("MuonEff_corr/RunABCD_SF_ID.root", "NUM_TightID_DEN_TrackerMuons_pt_abseta");

	//---MUON TIGHT ISO EFFICIENCY HISTOGRAM---
	TH2D* MuonTightISOEfficiencyHist = nullptr;
	MuonTightISOEfficiencyHist = cms.ReadingFileAndGettingTH2Histogram("MuonEff_corr/RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");
 
	//---ELECTRON ID SCALE FACTOR ---
	TH2D* ElectroScaleFactorHistogram = nullptr;
	ElectroScaleFactorHistogram = cms.ReadingFileAndGettingTH2Histogram("ElecEff_corr/ElectronSF_passingMVA102Xwp90isoHWWiso0p06_2018runABCD.root", "ElectronSF");

	//---TRIGGER HLT_IsoMu24 SCALE FACTOR---
	TH2D* TriggerHLTIsoMu24ScaleFactorHistogram = nullptr;
	TriggerHLTIsoMu24ScaleFactorHistogram = cms.ReadingFileAndGettingTH2Histogram("python/pyroot/Trigger_HLT_IsoMu24_weights.root", "HLT_IsoMu24_SFHist");

	//--RooWorkspace for QCD estimation
	RooWorkspace *w = nullptr;
	w = cms.ReadingFileAndGettingRooWorkspace("QCD/htt_scalefactors_legacy_2018.root", "w");

	//--Zpt Reweighting
	TH2D* ZptReweightingHistogram = nullptr;
	ZptReweightingHistogram = cms.ReadingFileAndGettingTH2Histogram("python/zptm_weights_2018_kit.root", "zptmass_histo");

	int samplecounter = -1;

	for (auto &xsmp : cms._SampleInfo)
	{
		samplecounter = samplecounter+1;
		std::cout<<"Processing sample: " << xsmp.GetSampleId() << " " << xsmp.GetNumberOfFilesInSample() <<std::endl;

		for(unsigned int iFileInSample=1; iFileInSample<=xsmp.GetNumberOfFilesInSample(); iFileInSample++)
		{
			std::cout << " LINE 381: iFileInSample " << iFileInSample << std::endl;

			// Set tree
			bool foundTree = true;

			if(!cms.SetTreeFile(xsmp,iFileInSample)) foundTree = false;

			std::cout << "NEW CHECK " << xsmp.GetSampleId() << std::endl;

			if (!foundTree) continue;

			// Loop on events for the current sample
      			long long nevents = cms._NANOTREE->GetEntriesFast();

//			if (nevents>1000) nevents=1000;

			std::cout << "Reading entries in this file: " << nevents << std::endl;

      			for (long long iEvent=0; iEvent<nevents; iEvent++)
			{
          	 		// Set next entry

           	 		if (cms.SetEntry(iEvent)<0) {std::cout << "WARNING: PROBLEMS READING NTUPLES" << std::endl; break;}

           	 		// Preselect (summary branch must be read before)
           	 		if (!cms.Preselect()) continue;

           	 		// Select (event branch must be read before)
           	 		if (!cms.Select(argv)) continue;

            			// Fill histograms
				double Pileup_weight = 1.;

				double MuonTightIDEff_weight = 1.;

				double MuonTightISOEff_weight = 1.;

				double ElectronEff_weight = 1.0;

				double Trigger_weight = 1.0;

				double Btag_weight = 1.0;

				double Generator_weight = 1.0;

				double osss_weight = 1.0;

				double merge = 1.0;

				double zpt_weight = 1.0;

				SCALEFACTORS::drellyan_info x;

				if( xsmp.SampleFlag!="isMCSignal")
				{
					osss_weight = SCALEFACTORS::QCDEstimation(cms, w);
				}

				if(!xsmp.GetSampleId().Contains("Data"))
				{
					//std::cout << "CALCULANDO PESOS" << std::endl;
					Pileup_weight = SCALEFACTORS::PileupSF(cms, Ratio);

					Generator_weight = SCALEFACTORS::GeneratorWeightSF(cms);

					MuonTightIDEff_weight = SCALEFACTORS::MuonSF(cms, MuonTightIDEfficiencyHist);

					MuonTightISOEff_weight = SCALEFACTORS::MuonSF(cms, MuonTightISOEfficiencyHist);

					ElectronEff_weight = SCALEFACTORS::ElectronSF(cms, ElectroScaleFactorHistogram);

					Trigger_weight	= SCALEFACTORS::TriggerSF(cms, TriggerHLTIsoMu24ScaleFactorHistogram);

					Btag_weight = SCALEFACTORS::BtaggingSF(cms, char(argv[4][0]));

					if( xsmp.GetSampleId().Contains("Merge") && !xsmp.GetSampleId().Contains("amcatnlo"))
					{
						if(xsmp.GetSampleId().Contains("WJets")) 
							merge = SCALEFACTORS::MergingSF(cms, xsmp.GetSampleId(), WJetsMergingMap);
						else if(xsmp.GetSampleId().Contains("DY")) 
							{
								if(xsmp.GetSampleId().Contains("ZTauTau"))merge = SCALEFACTORS::MergingSF(cms, xsmp.GetSampleId(), DYMergingMap_ZTauTau);
								else if(xsmp.GetSampleId().Contains("ZMuoMuo"))merge = SCALEFACTORS::MergingSF(cms, xsmp.GetSampleId(), DYMergingMap_ZMuoMuo);
								else if(xsmp.GetSampleId().Contains("ZEleEle"))merge = SCALEFACTORS::MergingSF(cms, xsmp.GetSampleId(), DYMergingMap_ZEleEle);
							}
					}
					else merge = 1.0;

					if(xsmp.GetSampleId().Contains("DY"))
					{
						x= SCALEFACTORS::ZptReweighting(cms, ZptReweightingHistogram, xsmp);
						if(!xsmp.GetSampleId().Contains("amcatnlo")) zpt_weight = x.weight;
						else zpt_weight= 1.0;
					}
					else zpt_weight= 1.0;

				}

				// For the data we need a fully exclusive selection.
				else if( string(argv[4]) == "0bjets" && cms.IndexMedjetsel.size()!=0 ) continue; 
				else if( string(argv[4]) == "1bjet"  && cms.IndexMedjetsel.size()!=1 ) continue; 
				else if( string(argv[4]) == "2bjets" && cms.IndexMedjetsel.size()!=2 ) continue; 

				double event_weight = Pileup_weight*MuonTightIDEff_weight*MuonTightISOEff_weight*ElectronEff_weight*Trigger_weight*Btag_weight*Generator_weight*merge*zpt_weight;
//                              double event_weight = Pileup_weight*MuonTightIDEff_weight*MuonTightISOEff_weight*ElectronEff_weight*Trigger_weight*Btag_weight*Generator_weight*merge;

				if(event_weight  ==0.)
				{
				std::cout << "SAMPLE " << xsmp.GetSampleId() << std::endl;
				std::cout << "event " << iEvent << std::endl;
				std::cout << "Pileup_weight " << Pileup_weight << std::endl;
				std::cout << "		MuonTightIDEff_weight " << MuonTightIDEff_weight << std::endl;
				std::cout << "			MuonTightISOEff_weight " << MuonTightISOEff_weight << std::endl;
				std::cout << "				ElectronEff_weight " << ElectronEff_weight << std::endl;
				std::cout << "					Trigger_weight " << Trigger_weight << std::endl;
				std::cout << "						bTag_weight " << Btag_weight << std::endl;
				std::cout << "							generator_weight " << Generator_weight << std::endl;
				std::cout << "								osss_weight " << osss_weight << std::endl;
				std::cout << "									merging " << merge << std::endl;
				std::cout << "										zpt_weight " << zpt_weight << std::endl;
				std::cout << " 											Total calculated Weight " << event_weight << std::endl << std::endl;

				}


//				std::cout << "hola" << x.TauFlag << " " << x.MuoFlag << " " << 

				if(x.TauFlag && !xsmp.GetSampleId().Contains("ZTauTau")) continue;
				else if(x.MuoFlag && !xsmp.GetSampleId().Contains("ZMuoMuo")) continue;
				else if(x.EleFlag && !xsmp.GetSampleId().Contains("ZEleEle")) continue;
				else
				{
                	                // Filling histograms
					if(string(argv[1])=="OS")
					{
//						std::cout << "FILLING" << std::endl;
						 HISTS::Fill(cms, xsmp, event_weight);
                                                 cms.FullSelection(argv, event_weight, xsmp);
					}
					else if(string(argv[1])=="SS") 
					{
//						std::cout << "FILLING" << std::endl;
						HISTS::Fill(cms, xsmp, event_weight*osss_weight);
						cms.FullSelection(argv, event_weight*osss_weight, xsmp);
					}
				}


      			}//End of events loop

		}//End of trees loop

		//For every sample I am saving the histograms written down in the txt file NameOfHistograms
		HISTS::Save(cms, xsmp, txtFiles, string(argv[1]), string(argv[2]), string(argv[3]), string(argv[4]));


  	}

	//if(nsamples>1) HISTS::Draw(cms, txtFiles, string(argv[1]), string(argv[2]), string(argv[3]), string(argv[4])); //he puesto que si pones mas de una sample, pinte, sino solo guarda el histograma
	//HISTS::Save(cms,  , txtFiles, string(argv[1]), string(argv[2]), string(argv[3]), string(argv[4]));
	//HISTS::QCDHistogram()
  	// To see things interactively (commnent otherwise) 
  	//if (!gROOT->IsBatch()) app->Run();

  	return 0;
}

bool TopAnalysis::Preselect()
{
	// Cut on summary information to save processing time
  	UInt_b(nElectron);
  	UInt_b(nMuon);
	Bool_b(HLT_IsoMu24);

	if(nElectron==0)         return false;
  	if(nMuon+nElectron<2) 	return false; // look for events with >= 2 leptons
	if(!HLT_IsoMu24)        return false;

  	return true;    
}
bool TopAnalysis::FullSelection(char *argv[], double eventweight, SAMPLES& sample)
{

	VFloat_b(Electron_pt); 
        VFloat_b(Muon_pt); 	
        VFloat_b(Electron_phi);	
       	VFloat_b(Muon_phi);	
        VFloat_b(Electron_eta);	
        VFloat_b(Muon_eta);
	Float_b(MET_pt); 
        Float_b(MET_phi);
        VFloat_b(Jet_pt);  
	VFloat_b(Jet_eta); 
	VFloat_b(Jet_phi);
	VFloat_b(Jet_mass);
	

	TLorentzVector LorentzElec;		
	TLorentzVector LorentzMuon;	
	TLorentzVector LorentzMET;


	LorentzElec.SetPtEtaPhiM(Electron_pt[Indexelsel],Electron_eta[Indexelsel],Electron_phi[Indexelsel], electronMass);
	LorentzMuon.SetPtEtaPhiM(Muon_pt[Indexmusel], Muon_eta[Indexmusel], Muon_phi[Indexmusel], muonMass);
	LorentzMET.SetPtEtaPhiM(MET_pt, 0.0, MET_phi, METMass);	


	double mt  = FOURVECTORS::TranverseMass(LorentzMuon, LorentzMET);
	double angleemu = ANGLES::DeltaPhi(Muon_phi[Indexmusel], Electron_phi[Indexelsel]);
	double angleemet = ANGLES::DeltaPhi(MET_phi, Electron_phi[Indexelsel]);

	if(string(argv[3])=="2jets")
	{
		TLorentzVector LorentzJet1;
		TLorentzVector LorentzJet2;

		LorentzJet1.SetPtEtaPhiM(Jet_pt[Indexjetsel.at(0)], Jet_eta[Indexjetsel.at(0)], Jet_phi[Indexjetsel.at(0)], Jet_mass[Indexjetsel.at(0)]);
		LorentzJet2.SetPtEtaPhiM(Jet_pt[Indexjetsel.at(1)], Jet_eta[Indexjetsel.at(1)], Jet_phi[Indexjetsel.at(1)], Jet_mass[Indexjetsel.at(1)]);

		double mjj = FOURVECTORS::InvariantMass(LorentzJet1, LorentzJet2);
		FillPlot1D("mjj_beforecuts", sample , mjj , eventweight);
		//std::cout << mjj << std::endl;
		//here we are boosting ggH signal
		if (mjj<550. && mt>15. && angleemet<0.5)
		{
		FillPlot1D("hCollinearMass_2jetsGGH",sample,FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[Indexelsel], MET_phi)),  eventweight);
		FillPlot1D("mjj_GGHcuts",sample,mjj, eventweight);
		}
		//here we are boosting VBF signal
		if (mjj>=550. && mt>15. && angleemet<0.3)   
		{
		FillPlot1D("hCollinearMass_2jetsVBF",sample,FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[Indexelsel], MET_phi)),  eventweight);
		FillPlot1D("mjj_VBFcuts", sample, mjj, eventweight);
		}
	}

	else if(string(argv[3])=="0jets")
	{
		double newptmuCut = 30.;

		if(Muon_pt[Indexmusel] < newptmuCut) return false;
		if(mt < 60.) return false;
		if(angleemet > 0.7) return false;
		if(angleemu < 2.5) return false;

		FillPlot1D("hCollinearMass_0jets",sample,FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[Indexelsel], MET_phi)), eventweight);

	}

	else if(string(argv[3])=="1jet")
	{
		if(mt < 40.) return false;
		if(angleemet > 0.7) return false;
		if(angleemu < 1.0) return false;

		FillPlot1D("hCollinearMass_1jet",sample,FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[Indexelsel], MET_phi)), eventweight);
	}


}
bool TopAnalysis::Select(char* argv[]) 
{
	////----- I am going to select one muon and one electron and veto any third lepton in the event
	// I start the counters for the interested objets to be selected. 
  	nmusel = 0;
  	nelsel = 0;
	njtsel = 0;
	nbtagsel = 0;
	Indexmusel = -1;
	Indexelsel = -1;
	nbLooseSel = 0; 
	nbMediumSel = 0;
	//These are the counters to count the leptons that we are going to veto
	nmuonveto = 0;
	nelectronveto = 0;
	ntauonveto = 0;
	//I need to empty the arrays I am working with from the event before to avoid any confusion
	Indexjetsel.clear();
	IndexMedjetsel.clear();
	Indexbtagsel.clear();

	////----- I declare the variables from the ntuple I am gonna work with
	// Muon variables
  	UInt_b(nMuon);
	VInt_b(Muon_charge);
  	VFloat_b(Muon_pt);
	VFloat_b(Muon_phi);
  	VFloat_b(Muon_eta);
  	VBool_b(Muon_tightId);
	VFloat_b(Muon_pfRelIso04_all);
	VFloat_b(Muon_dxy);
	VFloat_b(Muon_dz);
	VBool_b(Muon_mediumId);
	VUChar_b(Muon_highPtId);
	// Electron variables
	UInt_b(nElectron);
	VInt_b(Electron_charge);
  	VFloat_b(Electron_pt);
	VFloat_b(Electron_phi);
  	VFloat_b(Electron_eta);
	VBool_b(Electron_mvaFall17V1Iso_WP90); 
  	VFloat_b(Electron_pfRelIso03_all);
	VFloat_b(Electron_dxy);
	VFloat_b(Electron_dz);
	// Tauons variables
  	VFloat_b(Tau_pt);
  	UInt_b(nTau);
	VInt_b(Tau_decayMode);
	VUChar_b(Tau_idMVAoldDM);
	VFloat_b(Tau_dz);
	//Jet variables
  	UInt_b(nJet);
	VInt_b(Jet_jetId);
	VInt_b(Jet_puId);
  	VFloat_b(Jet_pt);
  	VFloat_b(Jet_eta);
	VFloat_b(Jet_phi);
	VFloat_b(Jet_btagDeepB);

	//[#selectingElectrons]
	//I start to loop over the total electron number in the event
  	for (unsigned int ie=0; ie<nElectron; ++ie) 
	{
		//Isolation: we want our selected electron to be isolated using the following cut
      		if (Electron_pfRelIso03_all[ie]>0.3)	continue;
		//Vertex geometry: we want electrons coming from the center of the colision
		if (abs(Electron_dxy[ie])>0.045)	continue;
		if (abs(Electron_dz[ie])>0.2)		continue;
		//Electron angle: we ask for a certain angle
		if (fabs(Electron_eta[ie])>etaelCut)    continue;
		//Isolation cut taking from HWW
		if (!Electron_mvaFall17V1Iso_WP90[ie])  continue;
		//Electron pt cut
		if (Electron_pt[ie]<ptelCut) 		continue;
		//The electron that we want to select should be isolated between 0.3 y 0.06, if the isolation is bigger, we should count it as a veto
		//The smaller the Electron_pfRelIso03_all variable is, the smaller is the energy surrounding the electron
      		if (Electron_pfRelIso03_all[ie]>isoelCut)
		{
			nelectronveto++;
			continue;
		}
		//Finally we count the good candidates as our selected electron
      		else nelsel++;
		//We pick the position of the electron in the array to use it later
		if (Indexelsel == -1) Indexelsel = ie;

  	}

	//We just select one electron
	if (nelsel!=1) return false;
	//We do not want events with one or more than one electron classified as a veto
	if (nelectronveto>=1) return false;
	//[#selectingElectrons]

	//[#selectingMuons]
	//I loop now over the total number of muons in the event to know which one of them I am going to pick
  	for (unsigned int im=0; im<nMuon; ++im) 
	{
		//Muon pt first pt cut because we are dealing also with veto muons
		if (Muon_pt[im]<10.) 						continue;
		//Vertex geometry: we want muons coming from the center of the colision
		if (abs(Muon_dxy[im])>0.045) 					continue;
		if (abs(Muon_dz[im]>0.2)) 					continue;
      		//Muon angle: we ask for a certain angle
		if (fabs(Muon_eta[im])>etamuCut) 				continue;
		//We want select muons with at least a medium Id
		if (!Muon_mediumId[im]) 					continue;
		//We play the same role we did with the electrons
      		if (string(argv[2]) == "IM"  && Muon_pfRelIso04_all[im]>0.3)	continue; //less Isolated muon Muon_pfRelIso04_all[im]>isomuCut
      		//if (string(argv[2]) == "NIM" && Muon_pfRelIso04_all[im]<0.3) 		continue; //more Isolated muon Muon_pfRelIso04_all[im]<isomuCut

		//I want the muon veto to be less tight than the muons I want to take as the good one
		if (!Muon_tightId[im]
      		    || (string(argv[2]) == "IM"  && Muon_pfRelIso04_all[im]>isomuCut) //less Isolated muon Muon_pfRelIso04_all[im]>isomuCut
      		 //   || (string(argv[2]) == "NIM" && Muon_pfRelIso04_all[im]<isomuCut) //more Isolated muon Muon_pfRelIso04_all[im]<isomuCut
      		    || Muon_pt[im]<ptmuCut)
		{
			nmuonveto++;
			continue;
		}
		//Finally we count the good candidates as our selected electron
		else nmusel++;
		//We pick the position of the muon in the array to use it later
		if (Indexmusel == -1) Indexmusel = im;

  	}
	//We want events with one good muon
	if (nmusel!=1) return false;
	//We want to veto the events with one or more than one muon veto
	if (nmuonveto>=1) return false;
	//[#selectingMuonsEnds]

	//[#vetoingThirdLeptonStart]
	//The missing leptons I have to loop over are the tauons, so I do it
	for(unsigned int it = 0; it <nTau; ++it)
	{
		//Tau pt cut
		if(Tau_pt[it]<20) continue;
		//Tau decaying mode is asked
		if(!Tau_decayMode[it]) continue;
		//Tau Id is set
		if(int(Tau_idMVAoldDM[it])<2) continue;
		//Some geometrical issues
		if(Tau_dz[it]>0.2) continue;
		//I count these tauons as vetos
		ntauonveto++;
	}
	//If I find one or more than one tau in the event I skip it
	if(ntauonveto>=1) return false;
	//[#vetoingThirdLeptonEnds]

	//[#DecidingTheChargeOfTheSelectedLeptonsStarts]
	//Once I have selected both muon and electron, I want to check if the sum of their charges is neutral
	if (string(argv[1]) == "OS" && Electron_charge[Indexelsel]+Muon_charge[Indexmusel] != 0) return false; // Opposite charge
	if (string(argv[1]) == "SS" && Electron_charge[Indexelsel]+Muon_charge[Indexmusel] == 0) return false; //same charge 
	//[#DecidingTheChargeOfTheSelectedLeptonsEnds]

	//[#JetsSelectionStarts]
	//I loop over the total number of jets in the event I am analyzing
  	for (unsigned int ij=0; ij<nJet; ++ij) 
	{
		//Jet pt cut because i am dealing with normal jets and bjets
		if (Jet_pt[ij]<25.) 		continue;//This cut has changed from 30 to 20 and now it is finally 25
		//I ask for a certain angle for the jets and bjets in general
      		if (fabs(Jet_eta[ij])>etajtCut) 	continue;
		//I am also care about the jet id. it should be bigger than 2
		if (Jet_jetId[ij]<2) 			continue;
		//I put again a bigger cut in pt and also I care about the pileup Id
		if (Jet_pt[ij]<50 && Jet_puId[ij]<4) 	continue;

		//I calculate the value of DeltaR between the jets and the selected muon and electron.
		double deltaR_MuoJet = ANGLES::DeltaR(Muon_phi[Indexmusel], Jet_phi[ij], Muon_eta[Indexmusel], Jet_eta[ij]);
		double deltaR_EleJet = ANGLES::DeltaR(Electron_phi[Indexelsel], Jet_phi[ij], Electron_eta[Indexelsel], Jet_eta[ij]);

		//I asked for these numbers to be bigger than 0.5
		if(deltaR_MuoJet<0.5 || deltaR_EleJet<0.5) continue;// this cut was 0.3 before
		//If I arrive to this point I got a btag jet and I selected it and store it in an array
		nbtagsel++;
		Indexbtagsel.push_back(ij);
		//Finally I distinguish between btags and jets by just using the pt value
		if (Jet_pt[ij]<ptjtCut) 		continue;
		//I count the selected jets
      		njtsel++;
		//I store them in an array
		Indexjetsel.push_back(ij);

	}
	//[#JetsSelectionEnds]
	//I skip events with more than 2 jets or bjets equally
	if (Indexjetsel.size()>2) return false; 
	if (Indexbtagsel.size()>2) return false;

	//I play with the arguments I gave to the code to run 
	//I need the exact numbers of jets that i ask for
	if (string(argv[3])=="0jets") {
		if (Indexjetsel.size()!=0) return false;
        }
	else if (string(argv[3])=="1jet") {
		if (Indexjetsel.size()!=1) return false;
        }
	else if (string(argv[3])=="2jets") {
		if (Indexjetsel.size()!=2) return false;
        }
	else {
		std::cout << " **************ERROR: NOT OPTION FOUND****************** SELECTION ROUTINE - jets " << string(argv[3]) << std::endl;
	}

	//Regarding the btag jets selection I count if the are Medium type, so that I loop over the btag selected jets
	for(auto &ij : Indexbtagsel)
	{
      		//if (string(argv[5]) == "Medium" && Jet_btagDeepB[ij]>btagMediumCut) 
		//I asked for a certain geometry
		if (fabs(Jet_eta[ij])>2.4) continue;
		//I use the btagDeep cut to decide if it is medium type or not
      		if (Jet_btagDeepB[ij]>btagMediumCut) 
		{
			//I count them if that's the case
           		nbMediumSel++;
			IndexMedjetsel.push_back(ij);
            		//nbLooseSel++;
      		}
  	}

	//Finally I check again the arguments I gave to run the code. If I want a certain number of bjets I should have at least the same quantity or less
	if (string(argv[4])=="1bjet") {
		if (IndexMedjetsel.size()==0) return false;
        }
	else if (string(argv[4])=="2bjets") {
		if (IndexMedjetsel.size()<2) return false;
        }
	else if (string(argv[4])!="0bjets") {
		std::cout << " **************ERROR: NOT OPTION FOUND****************** SELECTION ROUTINE - bjets " << string(argv[4]) << std::endl;
		return false;
	}

	return true;
}
