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
	bool FullSelection(char *argv[], double event_weight, SAMPLES& sample, TString& prefixsamplename);

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
  	const double etamuCut = 2.5;//2.5
  	const double isomuCut = 0.15;
  	const double ptelCut = 13.;//10
  	const double etaelCut = 2.5;
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
		cms.AddPlot1D("hMetPt", 			"Missing E_{T} [GeV]", 			 100, 0., 200.);
  		cms.AddPlot1D("hMuonPt", 			"Muon pT [GeV]", 			 80, 20.0, 100.);
  		cms.AddPlot1D("hElectronPt", 			"Electron pT [GeV]", 			 10, 10.0, 150.);
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
		cms.AddPlot1D("hdeltaR", 			"#DeltaR Muon-Electron", 		 50, -2.4, 10.0);
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
	static void Fill(TopAnalysis &cms, const SAMPLES &sample, const TString& prefixsamplename, const double event_weight)
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
		cms.FillPlot1D("hMetPt", sample ,  MET_pt, prefixsamplename, event_weight);
		cms.FillPlot1D("hElectronPt", sample, Electron_pt[cms.Indexelsel], prefixsamplename, event_weight);
		cms.FillPlot1D("hMuonPt", sample,  Muon_pt[cms.Indexmusel], prefixsamplename, event_weight);

		if(cms.Indexjetsel.size()==1) 
		{
			cms.FillPlot1D("hJetPt1", sample, Jet_pt[cms.Indexjetsel.at(0)], prefixsamplename, event_weight);
			cms.FillPlot1D("hJetEta1", sample, Jet_eta[cms.Indexjetsel.at(0)], prefixsamplename, event_weight);
		}
		if(cms.Indexjetsel.size()==2)
		{
			cms.FillPlot1D("hJetEta1", sample, Jet_eta[cms.Indexjetsel.at(0)], prefixsamplename, event_weight);
			cms.FillPlot1D("hJetEta2", sample, Jet_eta[cms.Indexjetsel.at(1)], prefixsamplename, event_weight);				
			cms.FillPlot1D("hJetPt1", sample, Jet_pt[cms.Indexjetsel.at(0)], prefixsamplename, event_weight);	
			cms.FillPlot1D("hJetPt2", sample, Jet_pt[cms.Indexjetsel.at(1)], prefixsamplename, event_weight);		  		  	
		}

		cms.FillPlot1D("hnVertex", sample, PV_npvs, prefixsamplename, event_weight);
		cms.FillPlot1D("hnMuon", sample, nMuon, prefixsamplename, event_weight);
		cms.FillPlot1D("hnElectron", sample, nElectron, prefixsamplename, event_weight);
		cms.FillPlot1D("hnJet", sample, nJet, prefixsamplename, event_weight);//nJet
		cms.FillPlot1D("hMetPhi", sample, MET_phi, prefixsamplename, event_weight);
		cms.FillPlot1D("hElectronPhi", sample, Electron_phi[cms.Indexelsel], prefixsamplename, event_weight);
		cms.FillPlot1D("hMuonPhi", sample, Muon_phi[cms.Indexmusel], prefixsamplename, event_weight);
		cms.FillPlot1D("hElectronEta", sample, Electron_eta[cms.Indexelsel], prefixsamplename, event_weight);
		cms.FillPlot1D("hMuonEta", sample, Muon_eta[cms.Indexmusel],prefixsamplename, event_weight);

		cms.FillPlot1D("hMuonElectronPhi", sample, ANGLES::DeltaPhi(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel]), prefixsamplename, event_weight);
		cms.FillPlot1D("hMuonMetPhi", sample, ANGLES::DeltaPhi(Muon_phi[cms.Indexmusel], MET_phi), prefixsamplename, event_weight);
		cms.FillPlot1D("hElectronMetPhi", sample, ANGLES::DeltaPhi(Electron_phi[cms.Indexelsel], MET_phi), prefixsamplename, event_weight);
		//cms.FillPlot1D("hbtagDeepB", sample, Jet_btagDeepB[cms.Indexjetsel.at(0)], event_weight);

		TLorentzVector LorentzElec;
		TLorentzVector LorentzMuon;
		TLorentzVector LorentzMET;

		LorentzElec.SetPtEtaPhiM(Electron_pt[cms.Indexelsel],Electron_eta[cms.Indexelsel],Electron_phi[cms.Indexelsel], electronMass);
		LorentzMuon.SetPtEtaPhiM(Muon_pt[cms.Indexmusel], Muon_eta[cms.Indexmusel], Muon_phi[cms.Indexmusel], muonMass);
		LorentzMET.SetPtEtaPhiM(MET_pt, 0.0, MET_phi, METMass);

		cms.FillPlot1D("hInvariantMassMuonElectron", sample, FOURVECTORS::InvariantMass(LorentzMuon, LorentzElec), prefixsamplename, event_weight);
		cms.FillPlot1D("hTransverseMassMuonMET", sample, FOURVECTORS::TranverseMass(LorentzMuon, LorentzMET), prefixsamplename, event_weight);
		cms.FillPlot1D("hTransverseMassElectronMET", sample, FOURVECTORS::TranverseMass(LorentzElec, LorentzMET), prefixsamplename, event_weight);
		cms.FillPlot1D("hnbLooseSel", sample, cms.nbLooseSel, prefixsamplename, event_weight);
		cms.FillPlot1D("hnbMediumSel", sample, cms.nbMediumSel, prefixsamplename, event_weight);

		cms.FillPlot1D("hdeltaR", sample, ANGLES::DeltaR(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel], Muon_eta[cms.Indexmusel], Electron_eta[cms.Indexelsel]), prefixsamplename, event_weight);

		cms.FillPlot1D("hCollinearMass", sample, FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[cms.Indexelsel], MET_phi)), prefixsamplename, event_weight);

		cms.FillPlot1D("hnMuoSel", sample, cms.nmusel, prefixsamplename, event_weight);
		cms.FillPlot1D("hnEleSel", sample, cms.nelsel, prefixsamplename, event_weight);
		cms.FillPlot1D("hnJetSel", sample, cms.njtsel, prefixsamplename, event_weight);
		cms.FillPlot1D("hSumMuonElectronPt", sample, (LorentzMuon+LorentzElec).Pt() , prefixsamplename, event_weight);

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
		TString samplename;

		if(sample.Contains("DY")) samplename = "DY";
		else if(sample.Contains("WJets")) samplename = "WJets";

		if(sample.Contains("Inclusive"))
		{
			numberofjets = (int)LHE_Njets;
			//std::cout << "Sample " << sample << " Numberofjets " << numberofjets << std::endl;
			if(numberofjets==0) merge_weight = merge[sample];
			if(numberofjets==1) merge_weight = merge[samplename+"1Merge"];
			if(numberofjets==2) merge_weight = merge[samplename+"2Merge"];
			if(numberofjets==3) merge_weight = merge[samplename+"3Merge"];
			if(numberofjets==4) merge_weight = merge[samplename+"4Merge"];
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
//		double zpt_weight;
//		double zpt_weightlepton;
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

		if(ZFound)
		{
			genM = GenPart_mass[indexZ];
			genpT = GenPart_pt[indexZ];
			int ID = abs(GenPart_pdgId[selectedlepton.at(0)]);
			double leptonMass = 0;
	                if(ID==15) {leptonMass = tauonMass; x.TauFlag = true;}
                        else if(ID==13) {leptonMass = muonMass; x.MuoFlag = true;}
                        else if(ID==11) {leptonMass = electronMass; x.EleFlag = true;}
		}
		else if(!ZFound)
		{
			int ID = abs(GenPart_pdgId[selectedlepton.at(0)]);
			double leptonMass = 0;
			if(ID==15) {leptonMass = tauonMass; x.TauFlag = true;}
			else if(ID==13) {leptonMass = muonMass; x.MuoFlag = true;}
			else if(ID==11) {leptonMass = electronMass; x.EleFlag = true;}

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
/*	cms.AddSample("SamplesDirectoryV6/SingleMuRun2018A_V6.txt");
	cms.AddSample("SamplesDirectoryV6/SingleMuRun2018B_V6.txt");
	cms.AddSample("SamplesDirectoryV6/SingleMuRun2018C_V6.txt");
	cms.AddSample("SamplesDirectoryV6/SingleMuRun2018D_V6.txt");

	//----- ADDING SIGNAL SAMPLES
	cms.AddSample("SamplesDirectoryV6/GluGlu_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunI0-v1_NANOAODSIM_V6.txt");
	cms.AddSample("SamplesDirectoryV6/VBF_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAu0-v1_NANOAODSIM_V6.txt");

	cms.AddSample("SamplesDirectoryV6/GluGluHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nan0-v1_NANOAODSIM.txt");
        cms.AddSample("SamplesDirectoryV6/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn0-v1_NANOAODSIM.txt");
        cms.AddSample("SamplesDirectoryV6/VBFHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nano251-v1_NANOAODSIM.txt");
        cms.AddSample("SamplesDirectoryV6/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_RunIIAutumn18N0-v1_NANOAODSIM.txt");

	//------ADDING BCKGRND SAMPLES
	cms.AddSample("SamplesDirectoryV6/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv0-v1_NANOAODSIM_V6.txt");
	cms.AddSample("SamplesDirectoryV6/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv0-v1_NANOAODSIM.txt");

		//Samples to Merge MADGRAPH

		//-WJETS
		cms.AddSample("SamplesDirectoryV6/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv6-0-v1_NANOAODSIM_V6.txt");
		cms.AddSample("SamplesDirectoryV6/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");
*/
		//-DY
		cms.AddSample("SamplesDirectoryV6/DYJetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_RunIIAutum2019.txt");
		cms.AddSample("SamplesDirectoryV6/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM.txt");
		cms.AddSample("SamplesDirectoryV6/DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019.txt");

		//Samples to Merge AMCATNLO

		//-WJETS
		/*cms.AddSample("SamplesDirectoryV6/INCLUSIVEISMISSING");*/
		/*cms.AddSample("SamplesDirectoryV6/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		//-DY
		cms.AddSample("SamplesDirectoryV6/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAODv6_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");
		cms.AddSample("SamplesDirectoryV6/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18NanoAO0-v1_NANOAODSIM");*/

/*	cms.AddSample("SamplesDirectoryV6/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIAutum1-v1_NANOAODSIM_V6.txt");
	cms.AddSample("SamplesDirectoryV6/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIA1-v1_NANOAODSIM_V6.txt");

	cms.AddSample("SamplesDirectoryV6/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv60-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019.txt");

        cms.AddSample("SamplesDirectoryV6/ZZ_TuneCP5_13TeV-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_u0-v1_NANOAODSIM_V6.txt");
       	cms.AddSample("SamplesDirectoryV6/WZ_TuneCP5_13TeV-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_u0-v1_NANOAODSIM_V6.txt");

	cms.AddSample("SamplesDirectoryV6/WminusHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nan0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/WplusHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nano0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/ZHToTauTau_M125_13TeV_powheg_pythia8_RunIIAutumn18NanoAODv6-Nano25Oc0-v1_NANOAODSIM.txt");

	cms.AddSample("SamplesDirectoryV6/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv6-Na0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_RunIIAutum0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8_RunIIA0-v1_NANOAODSIM.txt");
	cms.AddSample("SamplesDirectoryV6/EWKZ2Jets_ZToNuNu_TuneCP5_PSweights_13TeV-madgraph-pythia8_RunIIAutu0-v1_NANOAODSIM.txt");
*/
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

	//Defining the maps for merged samples that contain the information about the name of the sample and their corresponding merging weight depending on the number of jets
	std::map<TString,double> WJetsMergingMap = cms.MergingMCSamples(ListWJetsToMerge);
	std::map<TString,double> DYMergingMap = cms.MergingMCSamples(ListDYToMerge);

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
						if(xsmp.GetSampleId().Contains("WJets")) merge = SCALEFACTORS::MergingSF(cms, xsmp.GetSampleId(), WJetsMergingMap);
						else if(xsmp.GetSampleId().Contains("DY")) merge = SCALEFACTORS::MergingSF(cms, xsmp.GetSampleId(), DYMergingMap);
					}
					else merge = 1.0;

					if(xsmp.GetSampleId().Contains("DY") && !xsmp.GetSampleId().Contains("amcatnlo"))
					{
						x= SCALEFACTORS::ZptReweighting(cms, ZptReweightingHistogram, xsmp);
						zpt_weight = x.weight;
					}
					else zpt_weight= 1.0;

				}

				// For the data we need a fully exclusive selection.
				else if( string(argv[4]) == "0bjets" && cms.IndexMedjetsel.size()!=0 ) continue; 
				else if( string(argv[4]) == "1bjet"  && cms.IndexMedjetsel.size()!=1 ) continue; 
				else if( string(argv[4]) == "2bjets" && cms.IndexMedjetsel.size()!=2 ) continue; 

				double event_weight = Pileup_weight*MuonTightIDEff_weight*MuonTightISOEff_weight*ElectronEff_weight*Trigger_weight*Btag_weight*Generator_weight*merge*zpt_weight;

				/*if(event_weight  >120.)
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
				}*/

				// Full selection for mutaue channel
//				cms.FullSelection(argv, event_weight, xsmp);

				if(!xsmp.GetSampleId().Contains("DY"))
				{
					TString prefixsamplename="";
					cms.FullSelection(argv, event_weight, xsmp, prefixsamplename);
					if(string(argv[1])=="OS") HISTS::Fill(cms, xsmp, prefixsamplename, event_weight);
					else if(string(argv[1])=="SS") HISTS::Fill(cms, xsmp, prefixsamplename, event_weight*osss_weight);
				}
				else
				{
					TString prefixsamplename;
					if(x.TauFlag) prefixsamplename="ZTauTau";
					else if(x.MuoFlag) prefixsamplename="ZMuMu";
                                        else if(x.EleFlag) prefixsamplename="ZElEl";
			//		std::cout <<"PREFIX  " << prefixsamplename << " " <<x.TauFlag << " " << x.MuoFlag << " " << x.EleFlag << std::endl;
					if(string(argv[1])=="OS")
					{
						HISTS::Fill(cms, xsmp, prefixsamplename, event_weight);
						cms.FullSelection(argv, event_weight, xsmp, prefixsamplename);
					}

					else if (string(argv[1])=="SS")
					{
                                                HISTS::Fill(cms, xsmp, prefixsamplename, event_weight*osss_weight);
                                                cms.FullSelection(argv, event_weight*osss_weight, xsmp, prefixsamplename);
					}
				}
      			}//End of events loop

		}//End of trees loop
		//For every sample I am saving the histograms written down in the txt file NameOfHistograms
//		if(!xsmp.GetSampleId().Contains("DY")) 
		HISTS::Save(cms, xsmp, txtFiles, string(argv[1]), string(argv[2]), string(argv[3]), string(argv[4]));
//		else HISTS::Save(cms, xsmp, txtFiles, string(argv[1]), string(argv[2]), string(

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
bool TopAnalysis::FullSelection(char *argv[], double eventweight, SAMPLES& sample, TString& prefixsamplename)
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
		FillPlot1D("mjj_beforecuts", sample , mjj , prefixsamplename, eventweight);
		//std::cout << mjj << std::endl;
		//here we are boosting ggH signal
		if (mjj<550. && mt>15. && angleemet<0.5)
		{
		FillPlot1D("hCollinearMass_2jetsGGH",sample,FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[Indexelsel], MET_phi)), prefixsamplename, eventweight);
		FillPlot1D("mjj_GGHcuts",sample,mjj, prefixsamplename, eventweight);
		}
		//here we are boosting VBF signal
		if (mjj>=550. && mt>15. && angleemet<0.3)   
		{
		FillPlot1D("hCollinearMass_2jetsVBF",sample,FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[Indexelsel], MET_phi)), prefixsamplename, eventweight);
		FillPlot1D("mjj_VBFcuts", sample, mjj, prefixsamplename, eventweight);
		}
	}

	else if(string(argv[3])=="0jets")
	{
		double newptmuCut = 30.;

		if(Muon_pt[Indexmusel] < newptmuCut) return false;
		if(mt < 60.) return false;
		if(angleemet > 0.7) return false;
		if(angleemu < 2.5) return false;

		FillPlot1D("hCollinearMass_0jets",sample,FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[Indexelsel], MET_phi)), prefixsamplename, eventweight);

	}

	else if(string(argv[3])=="1jet")
	{
		if(mt < 40.) return false;
		if(angleemet > 0.7) return false;
		if(angleemu < 1.0) return false;

		FillPlot1D("hCollinearMass_1jet",sample,FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[Indexelsel], MET_phi)), prefixsamplename, eventweight);
	}


}
bool TopAnalysis::Select(char* argv[]) 
{
  	nmusel = 0;
  	nelsel = 0;
	njtsel = 0;
	nbtagsel = 0;
	Indexmusel = -1;
	Indexelsel = -1;
	nbLooseSel = 0; 
	nbMediumSel = 0;

	nmuonveto = 0;
	nelectronveto = 0;
	ntauonveto = 0;
	
	Indexjetsel.clear();
	IndexMedjetsel.clear();
	Indexbtagsel.clear();

	//std::cout << " FILLING JETS in selection " << Indexjetsel.size() <<" " << Indexjetsel.at(0) << " " << Indexjetsel.at(1)<< std::endl;

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

	UInt_b(nElectron);
	VInt_b(Electron_charge);
  	VFloat_b(Electron_pt);
	VFloat_b(Electron_phi);
  	VFloat_b(Electron_eta);
	VBool_b(Electron_mvaFall17V1Iso_WP90); 
  	VFloat_b(Electron_pfRelIso03_all);
	VFloat_b(Electron_dxy);
	VFloat_b(Electron_dz);

  	VFloat_b(Tau_pt);
  	UInt_b(nTau);
	VInt_b(Tau_decayMode);
	VUChar_b(Tau_idMVAoldDM);
	VFloat_b(Tau_dz);

  	UInt_b(nJet);
	VInt_b(Jet_jetId);
	VInt_b(Jet_puId);
  	VFloat_b(Jet_pt);
  	VFloat_b(Jet_eta);
	VFloat_b(Jet_phi);
	VFloat_b(Jet_btagDeepB);
	//std::cout <<" NUMBER OF LEPTONS " << nElectron << " " << nMuon << " " << nTau << std::endl;
	// ELECTRON SELECTION 
	//std::cout << "NEW EVENT" << std::endl;
  	for (unsigned int ie=0; ie<nElectron; ++ie) 
	{
      		if (Electron_pfRelIso03_all[ie]>0.3)	continue;
		if (abs(Electron_dxy[ie])>0.045)	continue;
		if (abs(Electron_dz[ie])>0.2)		continue;
		if (fabs(Electron_eta[ie])>etaelCut)    continue;
		if (!Electron_mvaFall17V1Iso_WP90[ie])  continue;
		if (Electron_pt[ie]<ptelCut) 		continue;
		//nelectronveto++;

      		if (Electron_pfRelIso03_all[ie]>isoelCut)
		{
			nelectronveto++;
			continue;
		}	
      		nelsel++;
		if (Indexelsel == -1) Indexelsel = ie;

  	}

	//std::cout << " AFTER ELECTRON SELECTION " << nelsel << " vetos " << nelectronveto << std::endl;
	if (nelsel!=1) return false; 
	//std::cout << " AFTER ELECTRON SELECTION " << nelsel << " vetos " << nelectronveto << std::endl;
	if (nelectronveto>=1) return false;

	// MUON SELECTION 
  	for (unsigned int im=0; im<nMuon; ++im) 
	{
		//std::cout << "Muon INFORMATION " << Muon_pt[im] << " " << int(Muon_highPtId[im]) << std::endl; exit(0);
     		//CUTS FOR MUON VETOS
		if (Muon_pt[im]<10.) 						continue;
		if (abs(Muon_dxy[im])>0.045) 					continue;
		if (abs(Muon_dz[im]>0.2)) 					continue;
      		if (fabs(Muon_eta[im])>etamuCut) 				continue;
		if (!Muon_mediumId[im]) 					continue;
      		if (string(argv[2]) == "IM"  && Muon_pfRelIso04_all[im]>0.3)	continue; //less Isolated muon Muon_pfRelIso04_all[im]>isomuCut
      		//if (string(argv[2]) == "NIM" && Muon_pfRelIso04_all[im]<0.3) 		continue; //more Isolated muon Muon_pfRelIso04_all[im]<isomuCut

		//nmuonveto++;

		if (!Muon_tightId[im]
      		    || (string(argv[2]) == "IM"  && Muon_pfRelIso04_all[im]>isomuCut) //less Isolated muon Muon_pfRelIso04_all[im]>isomuCut
      		    || (string(argv[2]) == "NIM" && Muon_pfRelIso04_all[im]<isomuCut) //more Isolated muon Muon_pfRelIso04_all[im]<isomuCut
      		    || Muon_pt[im]<ptmuCut)
		{
			nmuonveto++;
			continue;
		} 						
      		nmusel++;
		if (Indexmusel == -1) Indexmusel = im;

  	}

	if (nmusel!=1) return false; 
	//std::cout << " AFTER MUON SELECTION " << nmusel << " vetos " << nmuonveto << std::endl;
	if (nmuonveto>=1) return false;

	for(unsigned int it = 0; it <nTau; ++it)
	{
		//std::cout << "TAU INFORMATION " << Tau_pt[it] << " " << Tau_idMVAoldDM[it] << std::endl; exit(0);
		if(Tau_pt[it]<20) continue;
		if(!Tau_decayMode[it]) continue;
		if(int(Tau_idMVAoldDM[it])<2) continue;
		if(Tau_dz[it]>0.2) continue;
		
		ntauonveto++;
		
	}
	//std::cout << " AFTER Tau veto " << ntauonveto << std::endl;
	if(ntauonveto>=1) return false;

	//return true;
	//Once I have selected both muon and electron, I want to check if the sum of their charges is neutral

	if (string(argv[1]) == "OS" && Electron_charge[Indexelsel]+Muon_charge[Indexmusel] != 0) return false; // Opposite charge

	if (string(argv[1]) == "SS" && Electron_charge[Indexelsel]+Muon_charge[Indexmusel] == 0) return false; //same charge 

	//JETS SELECTION 
  	for (unsigned int ij=0; ij<nJet; ++ij) 
	{
		if (Jet_pt[ij]<25.) 		continue;//This cut has changed from 30 to 20 and now it is finally 25
      		if (fabs(Jet_eta[ij])>etajtCut) 	continue;
		if (Jet_jetId[ij]<2) 			continue;
		if (Jet_pt[ij]<50 && Jet_puId[ij]<4) 	continue;

		double deltaR_MuoJet = ANGLES::DeltaR(Muon_phi[Indexmusel], Jet_phi[ij], Muon_eta[Indexmusel], Jet_eta[ij]);
		double deltaR_EleJet = ANGLES::DeltaR(Electron_phi[Indexelsel], Jet_phi[ij], Electron_eta[Indexelsel], Jet_eta[ij]);
		
		if(deltaR_MuoJet<0.5 || deltaR_EleJet<0.5) continue;// this cut was 0.3 before 

		nbtagsel++;
		Indexbtagsel.push_back(ij);

		if (Jet_pt[ij]<ptjtCut) 		continue;

      		njtsel++;

		Indexjetsel.push_back(ij);

	}
	//std::cout << " AFTER JETS SELECTION " << njtsel << std::endl;

	if (Indexjetsel.size()>2) return false; 
	if (Indexbtagsel.size()>2) return false;

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
	//std::cout << Indexjetsel.size() << std::endl;

	for(auto &ij : Indexbtagsel)
	{
      		//if (string(argv[5]) == "Medium" && Jet_btagDeepB[ij]>btagMediumCut) 

		if (fabs(Jet_eta[ij])>2.4) continue;

      		if (Jet_btagDeepB[ij]>btagMediumCut) 
		{
			
           		nbMediumSel++;

			IndexMedjetsel.push_back(ij);

            		//nbLooseSel++;
      		} 
  	}

	//std::cout << " Number of bjets selected " << IndexMedjetsel.size() << " value " << Jet_btagDeepB[Indexjetsel.at(0)] << std::endl; 

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

/*	if(string(argv[3])=="0jets" && string(argv[4])=="0bjets" && (Indexjetsel.size()==0 && IndexMedjetsel.size()==0) )  return true; 

	else if(string(argv[3])=="1jet"  && string(argv[4])=="1bjet"  && (Indexjetsel.size()==1 && IndexMedjetsel.size()==1) )	return true;
	
	else if(string(argv[3])=="2jets" && string(argv[4])=="2bjets" && (Indexjetsel.size()==2 && IndexMedjetsel.size()==2) ) 	return true;

	else if(string(argv[3])=="1jet"  && string(argv[4])=="0bjets" && (Indexjetsel.size()==1 && (IndexMedjetsel.size()==0 || IndexMedjetsel.size()==1) ) ) return true;

	else if(string(argv[3])=="2jets" && string(argv[4])=="0bjets" && (Indexjetsel.size()==2 && (IndexMedjetsel.size()==0 || IndexMedjetsel.size()==1 || IndexMedjetsel.size()==2) ) ) return true;

	//std::cout << " **************ERROR: NOT OPTION FOUND****************** SELECTION ROUTINE - jets " << string(argv[3]) << " bjets " << string(argv[4]) << std::endl;
	return false;
*/
	
	//else std::cout << " You didn't give a correct number of jets and/or bjets" << std::endl;

}
