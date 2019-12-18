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
//#include <RooWorkspace.h>

const double electronMass = 0.0005109989461; //GeV
const double muonMass = 0.1056583745;
const double METMass = 0.0;

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

//Variables analisis

	//unsigned int njtsel;
 	unsigned int nbLooseSel;
      	unsigned int nbMediumSel;

	int Indexmusel;
	int Indexelsel;

	//char* argv[];

	std::vector<int> Indexjetsel;
	std::vector<int> IndexMedjetsel;

	int nmusel;
	int nelsel;
	int njtsel;

	TString mydir = "/afs/cern.ch/work/l/lurda/CMS/May_2019/ExoticHiggsDecay/codigojuan/GitLab_NANOAOD/LFVAnalysis/nanoaod/";
	TString dir = "/eos/user/l/lurda/CMS/SamplesAfterSkimmingNovember2019/"; 
	TString dir2 = "/eos/user/c/cepeda/LFVNANO/Skimming3/";

// Select if there are at least two leptons wit some phase space cuts
	
  	const double ptmuCut = 29.;//25
  	const double etamuCut = 2.4;//2.5
  	const double isomuCut = 0.15;
  	const double ptelCut = 10.;//25
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
		cms.AddPlot1D("hMet", 				"Missing E_{T} [GeV]", 			 100, 0., 200.);
  		cms.AddPlot1D("hMuon", 				"Muon pT [GeV]", 			 75, 25.0, 100.);
  		cms.AddPlot1D("hElectron", 			"Electron pT [GeV]", 			 15, 10.0, 150.);
  		cms.AddPlot1D("hJet", 				"Jet pT [GeV]", 			 250, 0.0, 250.);
  		cms.AddPlot1D("hMuonPhi", 			"Muon Phi", 				 100, -3.15, 3.15);
  		cms.AddPlot1D("hMuonEta", 			"Muon Eta", 				 100, -3.0, 3.0);
 		cms.AddPlot1D("hElectronPhi", 			"Electron Phi", 			 100, -3.15, 3.15);
  		cms.AddPlot1D("hElectronEta", 			"Electron Eta", 			 100, -3.0, 3.0);
  		cms.AddPlot1D("hJetEta", 			"Most energetic jet Eta", 		 100, -5.0, 5.0);
		cms.AddPlot1D("hnVertex", 			"Number of Vertex", 			 100, 0, 100);
  		cms.AddPlot1D("hnMuon", 			"Number of Muons", 			 6, -0.5, 5.5);
  		cms.AddPlot1D("hnElectron", 			"Number of Electrons", 			 7, -0.5, 6.5);
  		cms.AddPlot1D("hnJet", 				"Number of Jets", 			 25, 0.5, 25.5);	
  		cms.AddPlot1D("hmetPhi", 			"MET Phi", 				 100, -3.15, 3.15);
		cms.AddPlot1D("hMuonElectronPhi", 		"Muon-Electron Phi", 			 100, -3.15, 3.15);
		cms.AddPlot1D("hMuonMetPhi", 			"Muon-MET Phi", 			 100, -3.15, 3.15);
		cms.AddPlot1D("hElectronMetPhi", 		"Electron-Met Phi", 			 100, -3.15, 3.15);
		cms.AddPlot1D("hbtagDeepB",			"b-Tag jets Deep B",			 50, 0, 1);
		cms.AddPlot1D("hInvariantMassMuonElectron", 	"Invariant Mass Muon-Electron [GeV]", 	 100, 0, 200);
		cms.AddPlot1D("hTransverseMassMuonMET", 	"Transverse Mass Muon-MET [GeV]", 	 100, 0, 200);
		cms.AddPlot1D("hTransverseMassElectronMET", 	"Transverse Mass Electron-MET [GeV]", 	 100, 0, 200);
		cms.AddPlot1D("hnbLooseSel", 			"Number of BJets Loose", 		 5, -0.5, 4.5);
		cms.AddPlot1D("hnbMediumSel", 			"Number BJets Medium", 			 5, -0.5, 4.5);		
		cms.AddPlot1D("hdeltaR", 			"#DeltaR Muon-Electron", 		 50, -2.4, 10.0);
		cms.AddPlot1D("hCollinearMass", 		"Collinear Mass Muon-Electron [GeV]", 	 360, 0, 360);		
		cms.AddPlot1D("hnMuoSel", 			"Number of selected Muons", 		 10, 0, 10);
		cms.AddPlot1D("hnEleSel", 			"Number of selected Electrons", 	 10, 0, 10);
		cms.AddPlot1D("hnJetSel", 			"Number of selected jets", 		 11, -1, 10);
		
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
			
		cms.FillPlot1D("hMet", sample , MET_pt, event_weight);	    
		cms.FillPlot1D("hElectron", sample, Electron_pt[cms.Indexelsel], event_weight);	  
		cms.FillPlot1D("hMuon", sample, Muon_pt[cms.Indexmusel], event_weight);
		cms.FillPlot1D("hJet", sample, Jet_pt[cms.Indexjetsel[0]], event_weight);		  	
		cms.FillPlot1D("hnVertex", sample, PV_npvs, event_weight);
		cms.FillPlot1D("hnMuon", sample, nMuon, event_weight);
		cms.FillPlot1D("hnElectron", sample, nElectron, event_weight);
		cms.FillPlot1D("hnJet", sample, nJet, event_weight);//nJet
		cms.FillPlot1D("hmetPhi", sample, MET_phi, event_weight);
		cms.FillPlot1D("hElectronPhi", sample, Electron_phi[cms.Indexelsel], event_weight);
		cms.FillPlot1D("hMuonPhi", sample, Muon_phi[cms.Indexmusel], event_weight);
		cms.FillPlot1D("hElectronEta", sample, Electron_eta[cms.Indexelsel], event_weight);
		cms.FillPlot1D("hMuonEta", sample, Muon_eta[cms.Indexmusel],event_weight);
		cms.FillPlot1D("hJetEta", sample, Jet_eta[cms.Indexjetsel[0]], event_weight);		
		cms.FillPlot1D("hMuonElectronPhi", sample, ANGLES::DeltaPhi(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel]), event_weight);
		cms.FillPlot1D("hMuonMetPhi", sample, ANGLES::DeltaPhi(Muon_phi[cms.Indexmusel], MET_phi), event_weight);
		cms.FillPlot1D("hElectronMetPhi", sample, ANGLES::DeltaPhi(Electron_phi[cms.Indexelsel], MET_phi), event_weight);
		cms.FillPlot1D("hbtagDeepB", sample, Jet_btagDeepB[cms.Indexjetsel[0]], event_weight);

		TLorentzVector LorentzElec;		
		TLorentzVector LorentzMuon;	
		TLorentzVector LorentzMET;

		LorentzElec.SetPtEtaPhiM(Electron_pt[cms.Indexelsel],Electron_eta[cms.Indexelsel],Electron_phi[cms.Indexelsel], electronMass);
		LorentzMuon.SetPtEtaPhiM(Muon_pt[cms.Indexmusel], Muon_eta[cms.Indexmusel], Muon_phi[cms.Indexmusel], muonMass);	
		LorentzMET.SetPtEtaPhiM(MET_pt, 0.0, MET_phi, METMass);

		cms.FillPlot1D("hInvariantMassMuonElectron", sample, FOURVECTORS::InvariantMass(LorentzMuon, LorentzElec), event_weight);
		cms.FillPlot1D("hTransverseMassMuonMET", sample, FOURVECTORS::TranverseMass(LorentzMuon, LorentzMET), event_weight);
		cms.FillPlot1D("hTransverseMassElectronMET", sample, FOURVECTORS::TranverseMass(LorentzElec, LorentzMET), event_weight);
		cms.FillPlot1D("hnbLooseSel", sample, cms.nbLooseSel, event_weight);
		cms.FillPlot1D("hnbMediumSel", sample, cms.nbMediumSel, event_weight);

		cms.FillPlot1D("hdeltaR", sample, ANGLES::DeltaR(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel], Muon_eta[cms.Indexmusel], Electron_eta[cms.Indexelsel]), event_weight);

		cms.FillPlot1D("hCollinearMass", sample, FOURVECTORS::CollinearMass(LorentzMuon, LorentzElec, LorentzMET,ANGLES::DeltaPhi(Electron_phi[cms.Indexelsel], MET_phi)), event_weight);

		cms.FillPlot1D("hnMuoSel", sample, cms.nmusel, event_weight);
		cms.FillPlot1D("hnEleSel", sample, cms.nelsel, event_weight);
		cms.FillPlot1D("hnJetSel", sample, cms.njtsel, event_weight);
		
	}

	static void Save(TopAnalysis &cms, const SAMPLES &sample, ASCII &txtFiles, const string& charge, const string& muon, const string& jets, const string& bjets)
	{
		TString option = "RECREATE"; 

		for(auto &name: txtFiles.HistogramsList)
		{
			cms.SavingHistograms(sample, name, option, charge, muon, jets, bjets);
			if(option =="RECREATE") option = "UPDATE";
		}
	}

	static void Draw(CMSAnalysis &cms, ASCII &txtFiles, const string& charge, const string& muon, const string& jets, const string& bjets)
	{
	
		for(auto &name: txtFiles.HistogramsList)
		{
			std::cout << "Name of Histograms " << name << std::endl;
		 	cms.DrawPlot1D(name, "", Form("png_V6_%s_%s_%s_%s", charge.c_str(), muon.c_str(), jets.c_str(), bjets.c_str() ));
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
		VFloat_b(Jet_pt);
		VInt_b(Jet_hadronFlavour);

		int inbjets = (int)nbjets - 48; //in ASCII code, the numbers (digits) start from 48

		double bTag_weight = cms.bTagEventWeight(cms.IndexMedjetsel.size(), Jet_pt[cms.IndexMedjetsel[0]],Jet_hadronFlavour[cms.IndexMedjetsel[0]],Jet_pt[cms.IndexMedjetsel[1]], 						Jet_hadronFlavour[cms.IndexMedjetsel[1]], "comb", "central" , inbjets , "medium");
		return bTag_weight;
	}

	/*static double QCDEstimation(TopAnalysis &cms)
	{
		VFloat_b(Muon_pt);
  		VFloat_b(Muon_eta);
		VFloat_b(Electron_pt);
		VFloat_b(Electron_eta);
		VFloat_b(Electron_phi);
		VFloat_b(Muon_phi);
    	    	Int_b(nJet); 	

		double osss_weight;

	        TString file = "QCD/htt_scalefactors_legacy_2018.root";
				
		double dR = ANGLES::DeltaR(Muon_phi[cms.Indexmusel], Electron_phi[cms.Indexelsel], Muon_eta[cms.Indexmusel], Electron_eta[cms.Indexelsel]);

		osss_weight = cms.QCDEstimation(file, nJet, dR, Electron_pt, Muon_pt); 

		return osss_weight;
	}*/
	
};

int main(int argc, char* argv[])//"OS", "IM", "0jets", "0bjets"
{

	// Initialize analysis structure
  	TopAnalysis cms; 

  	int maxevents = -1; 	

	//*****USING FUNCTION AddFiles ******	
  	cms.AddFiles("Data", "DataA", 14.027047499*1000., -1, cms.dir+"SingleMuon_Run2018A-Nano25Oct2019-v1_NANOAOD/", 555, -1, 227489240);
  	cms.AddFiles("Data", "DataB", 7.060622497*1000.,  -1, cms.dir+"SingleMuon_Run2018B-Nano25Oct2019-v1_NANOAOD/", 264, -1, 110446445);
  	cms.AddFiles("Data", "DataC", 6.894770971*1000.,  -1, cms.dir+"SingleMuon_Run2018C-Nano25Oct2019-v1_NANOAOD/", 264, -1, 107972995);
  	cms.AddFiles("Data", "DataD", 31.283984136*1000., -1, cms.dir+"DeMaria_SingleMuon_Run2018D-Nano25Oct2019-v1_NANOAOD/", 1238, -1, -1);
	
  	cms.AddFiles("MCSignal", "GG",  -1, 48.58*0.1,  cms.dir+"GluGlu_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunI0-v1_NANOAODSIM/", 1,  maxevents, 42945741.6072);
 	cms.AddFiles("MCSignal", "VBF", -1, 3.782*0.1,  cms.dir+"VBF_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAu0-v1_NANOAODSIM/", 1,  maxevents, 7631474.65134);
	
  	cms.AddFiles("MCBckgr", "WW",    -1, 12.15488,  cms.dir+"WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv0-v1_NANOAODSIM/",  5,  maxevents, 85917643.789);
  	cms.AddFiles("MCBckgr", "WJets", -1, 61526.7,   cms.dir+"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv6-0-v1_NANOAODSIM/", 11,  maxevents, 70389866.8084);
  	cms.AddFiles("MCBckgr", "TW",    -1, 35.85,     cms.dir+"ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIAutum1-v1_NANOAODSIM/",  3,  maxevents, 334874732.0);
  	cms.AddFiles("MCBckgr", "TbarW", -1, 35.85,     cms.dir+"ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIA1-v1_NANOAODSIM/",  4,  maxevents, 266470418.054);
  	cms.AddFiles("MCBckgr", "TTbar", -1, 85.172,    cms.dir2+"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019/", 13,  maxevents, 4.622139340935600e+09);
  	cms.AddFiles("MCBckgr", "DY",    -1, 2075.14*3,  cms.dir2+"DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019/",          29,  maxevents, 3.298665625309500e+12);
  	cms.AddFiles("MCBckgr", "TTbar", -1, 88.29,    cms.dir+"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019/", 13,  maxevents, 4622080234.06);
 	cms.AddFiles("MCBckgr", "DY",    -1, 6435.0,   cms.dir+"DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019/",          30,  maxevents, 3.44609946801e+12);
  	cms.AddFiles("MCBckgr", "WZ",    -1, 27.6,      cms.dir+"WZ_TuneCP5_13TeV-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_u0-v1_NANOAODSIM/",  4,  maxevents, 3884167.00546);
  	cms.AddFiles("MCBckgr", "ZZ",    -1, 12.14,     cms.dir+"ZZ_TuneCP5_13TeV-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_u0-v1_NANOAODSIM/",  4,  maxevents, 1978776.75193);  	

  	// Initialize 1D histograms

	HISTS::Add(cms);
  	
  	// Loop on samples
  	unsigned int nsamples = cms.GetNumberOfSamples();

	ASCII txtFiles;

	txtFiles.txtCounter("NameOfHistograms.txt");

	//****PILE UP REWEIGHTING HISTOGRAM****
	TH1D* Ratio = nullptr; 
	Ratio = cms.ReadingFileAndGettingTH1Histogram(cms.mydir+"python/pyroot/pileupweights.root", "RatioPU");

	//****MUON TIGHT ID EFFICIENCY HISTOGRAM*****
	TH2D* MuonTightIDEfficiencyHist = nullptr;
	MuonTightIDEfficiencyHist = cms.ReadingFileAndGettingTH2Histogram(cms.mydir+"MuonEff_corr/RunABCD_SF_ID.root", "NUM_TightID_DEN_TrackerMuons_pt_abseta");

	//****MUON TIGHT ISO EFFICIENCY HISTOGRAM*****
	TH2D* MuonTightISOEfficiencyHist = nullptr;
	MuonTightISOEfficiencyHist = cms.ReadingFileAndGettingTH2Histogram(cms.mydir+"MuonEff_corr/RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");
 
	//****ELECTRON ID SCALE FACTOR *****
	TH2D* ElectroScaleFactorHistogram = nullptr;
	ElectroScaleFactorHistogram = cms.ReadingFileAndGettingTH2Histogram(cms.mydir+"ElecEff_corr/ElectronSF_passingMVA102Xwp90isoHWWiso0p06_2018runABCD.root", "ElectronSF");

	//****TRIGGER HLT_IsoMu24 SCALE FACTOR****
	TH2D* TriggerHLTIsoMu24ScaleFactorHistogram = nullptr;
	TriggerHLTIsoMu24ScaleFactorHistogram = cms.ReadingFileAndGettingTH2Histogram(cms.mydir+"python/pyroot/Trigger_HLT_IsoMu24_weights.root", "HLT_IsoMu24_SFHist");
	
	////los comentarios que empiezan por //// son los correspondientes a la clase de SAMPLES que he creado para aligerar el codigo (usar cuando est'e paralelizado)
	////for (auto &xsmp : cms._SampleInfo)
  	for (unsigned int iSample=0; iSample<nsamples; ++iSample) 
	{
		 ////std::cout<<"Procesando muestra " << xsmp.GetSampleId() << " " << xsmp.GetNumberOfFilesInSample() <<std::endl;
		std::cout<<"Procesando muestra " << iSample << " " << cms.GetSampleId(iSample) << " " <<cms.GetNumberOfFilesInSample(iSample)<<std::endl;
		
		////for(unsigned int iFileInSample=1; iFileInSample<=xsmp.GetNumberOfFilesInSample(); iFileInSample++)
		for(unsigned int iFileInSample=1; iFileInSample<=cms.GetNumberOfFilesInSample(iSample); iFileInSample++)
		{
			std::cout << " LINE 381: iFileInSample " << iFileInSample << std::endl;
			
			// Set tree                 
			bool foundTree = true;                 
	                      
			if(!cms.SetTreeFile(iSample,iFileInSample)) foundTree = false; 
			////if(!cms.SetTreeFile(1,iFileInSample)) foundTree = false;     //he puesto un 1 de momento pq para paralelizarlo solo quiero hacer una sample cada vez que corra           
			
			////std::cout << "NEW CHECK " << xsmp.GetSampleId() << std::endl;			
			
			if (!foundTree) continue;
      			
			// Loop on events for the current sample
      			long long nevents = cms._NANOTREE->GetEntriesFast();

			//if (nevents>1000) nevents=1000;

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
				
				//osss_weight = SCALEFACTORS::QCDEstimation(cms);

				////if(!xsmp.GetSampleId().Contains("Data"))	
				if(!cms.GetSampleId(iSample).Contains("Data"))
				{	
					Pileup_weight = SCALEFACTORS::PileupSF(cms, Ratio);

					Generator_weight = SCALEFACTORS::GeneratorWeightSF(cms);

					MuonTightIDEff_weight = SCALEFACTORS::MuonSF(cms, MuonTightIDEfficiencyHist);

					MuonTightISOEff_weight = SCALEFACTORS::MuonSF(cms, MuonTightISOEfficiencyHist);

					ElectronEff_weight = SCALEFACTORS::ElectronSF(cms, ElectroScaleFactorHistogram);

					Trigger_weight	= SCALEFACTORS::TriggerSF(cms, TriggerHLTIsoMu24ScaleFactorHistogram);

					Btag_weight = SCALEFACTORS::BtaggingSF(cms, char(argv[4][0]));
					
				}
				else if( (string(argv[3]) != "0jets" && string(argv[4]) == "0bjets")  && cms.IndexMedjetsel.size()>0 ) continue; 

				double event_weight = Pileup_weight * MuonTightIDEff_weight * MuonTightISOEff_weight * ElectronEff_weight * Trigger_weight * Btag_weight * Generator_weight;

				/*std::cout << "Pileup_weight " << Pileup_weight << std::endl;
				std::cout << "		MuonTightIDEff_weight " << MuonTightIDEff_weight << std::endl;
				std::cout << "			MuonTightISOEff_weight " << MuonTightISOEff_weight << std::endl;
				std::cout << "				ElectronEff_weight " << ElectronEff_weight << std::endl;
				std::cout << "					Trigger_weight " << Trigger_weight << std::endl;
				std::cout << "						bTag_weight " << Btag_weight << std::endl;
				std::cout << "							generator_weight " << Generator_weight << std::endl;
				std::cout << "								OSSS " << osss_weight << std::endl;
				std::cout << " 									Total calculated Weight " << event_weight << std::endl << std::endl;*/

				HISTS::Fill(cms, cms._SampleInfo[iSample], event_weight);
				////HISTS::Fill(cms, xsmp, event_weight);						
      			}			
		}
		//if(string(argv[6])== "SaveHistograms")
		HISTS::Save(cms, cms._SampleInfo[iSample], txtFiles, string(argv[1]), string(argv[2]), string(argv[3]), string(argv[4]));
		////HISTS::Save(cms, xsmp, txtFiles, string(argv[1]), string(argv[2]), string(argv[3]), string(argv[4]));

  	}  	

	//To see things interactively (commnent otherwise) 

  	auto app = new TRint("Top Analysis", &argc, argv); 

  	//if(string(argv[7]) == "Draw") 
	if(nsamples>1) HISTS::Draw(cms, txtFiles, string(argv[1]), string(argv[2]), string(argv[3]), string(argv[4])); //he puesto que si pones mas de una sample, pinte, sino solo guarda el histograma

  	// To see things interactively (commnent otherwise) 
  	if (!gROOT->IsBatch()) app->Run();

  	return 0;
}

bool TopAnalysis::Preselect()
{
	// Cut on summary information to save processing time
  	UInt_b(nElectron);
  	UInt_b(nMuon);
	Bool_b(HLT_IsoMu24);

	if(nElectron<0)         return false;
  	if(nMuon+nElectron<2) 	return false; // look for events with >= 2 leptons
	if(!HLT_IsoMu24)        return false;

  	return true;    
}

bool TopAnalysis::Select(char* argv[]) 
{
  	nmusel = 0;
  	nelsel = 0;
	njtsel = 0;
	Indexmusel = -1;
	Indexelsel = -1;
	nbLooseSel = 0; 
	nbMediumSel = 0;
	
	Indexjetsel.clear();
	IndexMedjetsel.clear();

  	UInt_b(nMuon);
	VInt_b(Muon_charge);
  	VFloat_b(Muon_pt);
	VFloat_b(Muon_phi);
  	VFloat_b(Muon_eta);
  	VBool_b(Muon_tightId);
	VFloat_b(Muon_pfRelIso04_all);

	UInt_b(nElectron);
	VInt_b(Electron_charge);
  	VFloat_b(Electron_pt);
	VFloat_b(Electron_phi);
  	VFloat_b(Electron_eta);
	VBool_b(Electron_mvaFall17V1Iso_WP90); 
  	VFloat_b(Electron_pfRelIso03_all);

  	UInt_b(nJet);
  	VFloat_b(Jet_pt);
  	VFloat_b(Jet_eta);
	VFloat_b(Jet_phi);
	VFloat_b(Jet_btagDeepB);

  	for (unsigned int ie=0; ie<nElectron; ++ie) 
	{
      		if (Electron_pfRelIso03_all[ie]>isoelCut) continue;
		if (!Electron_mvaFall17V1Iso_WP90[ie])	  continue;
      		if (Electron_pt[ie]<ptelCut) 		  continue;
      		if (fabs(Electron_eta[ie])>etaelCut) 	  continue;

      		nelsel++;

		if (Indexelsel == -1) Indexelsel = ie;
  	}

	if (nelsel!=1) return false; 

  	for (unsigned int im=0; im<nMuon; ++im) 
	{
     		if (!Muon_tightId[im]) 							continue;
      		if (string(argv[2]) == "IM"  && Muon_pfRelIso04_all[im]>isomuCut)       continue; //Isolated muon 
      		if (string(argv[2]) == "NIM" && Muon_pfRelIso04_all[im]<isomuCut) 	continue; //no Isolated muon 
      		if (Muon_pt[im]<ptmuCut) 						continue;
      		if (fabs(Muon_eta[im])>etamuCut) 					continue;

      		nmusel++;

		if (Indexmusel == -1) Indexmusel = im;

  	}
	if (nmusel!=1) return false; 
	
	//Once I have selected both muon and electron, I want to check if the sum of their charges is neutral

	if (string(argv[1]) == "OS" && Electron_charge[Indexelsel]+Muon_charge[Indexmusel] != 0) return false; // Opposite charge

	if (string(argv[1]) == "SS" && Electron_charge[Indexelsel]+Muon_charge[Indexmusel] == 0) return false; //same charge 

  	for (unsigned int ij=0; ij<nJet; ++ij) 
	{
      		if (Jet_pt[ij]<ptjtCut) 	continue;
      		if (fabs(Jet_eta[ij])>etajtCut) continue;

		double deltaR_MuoJet = ANGLES::DeltaR(Muon_phi[Indexmusel], Jet_phi[ij], Muon_eta[Indexmusel], Jet_eta[ij]);
		double deltaR_EleJet = ANGLES::DeltaR(Electron_phi[Indexelsel], Jet_phi[ij], Electron_eta[Indexelsel], Jet_eta[ij]);
		
		if(deltaR_MuoJet<0.4 || deltaR_EleJet<0.4) continue; 

      		njtsel++;

		Indexjetsel.push_back(ij);

	}
	
	for(auto &ij : Indexjetsel)
	{
      		//if (string(argv[5]) == "Medium" && Jet_btagDeepB[ij]>btagMediumCut) 
      		if (Jet_btagDeepB[ij]>btagMediumCut) 
		{
           		nbMediumSel++;

			IndexMedjetsel.push_back(ij);

            		//nbLooseSel++;
      		} 
		
		/*else if (string(argv[5]) == "Loose" && Jet_btagDeepB[ij]>btagLooseCut) 
		{
           		 nbLooseSel++;
      		}*/

  	}
  
  	if (Indexjetsel.size()>2) return false; 

	     if(string(argv[3])=="0jets" && string(argv[4])=="0bjets" && (Indexjetsel.size()==0 && IndexMedjetsel.size()==0) )  return true; 

	else if(string(argv[3])=="1jet"  && string(argv[4])=="1bjet"  && (Indexjetsel.size()==1 && IndexMedjetsel.size()==1) )	return true;
	
	else if(string(argv[3])=="2jets" && string(argv[4])=="2bjets" && (Indexjetsel.size()==2 && IndexMedjetsel.size()==2) ) 	return true;

	else if(string(argv[3])=="1jet"  && string(argv[4])=="0bjets" && (Indexjetsel.size()==1 && (IndexMedjetsel.size()==0 || IndexMedjetsel.size()==1) ) ) return true;

	else if(string(argv[3])=="2jets" && string(argv[4])=="0bjets" && (Indexjetsel.size()==2 && (IndexMedjetsel.size()==0 || IndexMedjetsel.size()==1 || IndexMedjetsel.size()==2) ) ) return true;

	//else std::cout << " You didn't give a correct number of jets and/or bjets" << std::endl;

}
