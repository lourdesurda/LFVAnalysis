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
		init(); //TFile("LFV_Analysis.root", "RECREATE").Close();
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

	std::vector<int> Indexjetsel;
	std::vector<int> IndexMedjetsel;

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
		cms.AddPlot1D("hMet", 				"Missing E_{T} [GeV]", 			 100, -10., 200.);
  		cms.AddPlot1D("hMuon", 				"Muon pT [GeV]", 			50, 25.0, 100.);
  		cms.AddPlot1D("hElectron", 			"Electron pT [GeV]", 			 50, 10.0, 150.);
  		cms.AddPlot1D("hJet", 				"Jet pT [GeV]", 			 50, 0.0, 250.);
  		cms.AddPlot1D("hMuonPhi", 			"Muon Phi", 				 100, -3.15, 3.15);
  		cms.AddPlot1D("hMuonEta", 			"Muon Eta", 				 100, -3.0, 3.0);
 		cms.AddPlot1D("hElectronPhi", 			"Electron Phi", 			 100, -3.15, 3.15);
  		cms.AddPlot1D("hElectronEta", 			"Electron Eta", 			 100, -3.0, 3.0);
  		cms.AddPlot1D("hJetEta", 			"Most energetic jet Eta", 		 100, -3.0, 3.0);
		cms.AddPlot1D("hnVertex", 			"Number of Vertex", 			 100, 0, 100);
  		cms.AddPlot1D("hnMuon", 			"Number of Muons", 			 6, 0, 6);
  		cms.AddPlot1D("hnElectron", 			"Number of Electrons", 			 7, 0, 7);
  		cms.AddPlot1D("hnJet", 				"Number of Jets", 			 25, 0, 25);	
  		cms.AddPlot1D("hmetPhi", 			"MET Phi", 				 100, -3.15, 3.15);
		cms.AddPlot1D("hMuonElectronPhi", 		"Muon-Electron Phi", 			 100, -3.15, 3.15);
		cms.AddPlot1D("hMuonMetPhi", 			"Muon-MET Phi", 			 100, -3.15, 3.15);
		cms.AddPlot1D("hElectronMetPhi", 		"Electron-Met Phi", 			 100, -3.15, 3.15);
		cms.AddPlot1D("hbtagDeepB",			"b-Tag jets Deep B",			 50, 0, 1);
		cms.AddPlot1D("hInvariantMassMuonElectron", 	"Invariant Mass Muon-Electron [GeV]", 	 100, 0, 200);
		cms.AddPlot1D("hTransverseMassMuonMET", 	"Transverse Mass Muon-MET [GeV]", 	 100, 0, 300);
		cms.AddPlot1D("hTransverseMassElectronMET", 	"Transverse Mass Electron-MET [GeV]", 	 100, 0, 300);
		cms.AddPlot1D("hnbLooseSel", 			"Number of BJets Loose", 		 5, -0.5, 4.5);
		cms.AddPlot1D("hnbMediumSel", 			"Number BJets Medium", 			 5, -0.5, 4.5);		
		cms.AddPlot1D("hdeltaR", 			"#DeltaR Muon-Electron", 		 50, -2.4, 10.0);
		cms.AddPlot1D("hCollinearMass", 		"Collinear Mass Muon-Electron [GeV]", 	 100, 0, 350);		
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
		VFloat_b(Jet_btagDeepB);											 							 							
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

	static void Save(TopAnalysis &cms, const SAMPLES &sample, ASCII &txtFiles)
	{
		TString option = "RECREATE"; 

		for(auto &name: txtFiles.HistogramsList)
		{
			cms.SavingHistograms(sample, name, option);
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

/*struct SCALEFACTORS 
{
	TString mydir = "/afs/cern.ch/work/l/lurda/CMS/May_2019/ExoticHiggsDecay/codigojuan/GitLab_NANOAOD/LFVAnalysis/nanoaod/";

	static TH1D* GettingTH1Histograms(CMSAnalysis &cms, TH1D* Histogram, TString& pathroot, TString& name)
	{
		Histogram = nullptr;
		Histogram = cms.ReadingFileAndGettingTH1Histogram(mydir+pathroot, name);
		return Histogram;
	}
	static TH1D* GettingTH2Histograms(CMSAnalysis &cms, TH2D* Histogram, TString& pathroot, TString& name)
	{
		Histogram = nullptr;
		Histogram = cms.ReadingFileAndGettingTH2Histogram(mydir+pathroot, name);
		return Histogram;
	}
	static double Muon(CMSAnalysis &cms, TH2D* MuonEfficiencyHist)
	{
		VFloat_b(Muon_pt);
  		VFloat_b(Muon_eta);

		if(Muon_pt[cms.Indexmusel] > 119.) 
		{
			double Muon_weight = cms.ScaleFactors(MuonEfficiencyHist, 119., fabs(Muon_eta[cms.Indexmusel]));
		}
		else 
		{
			double Muon_weight = cms.ScaleFactors(MuonEfficiencyHist, Muon_pt[cms.Indexmusel], fabs(Muon_eta[cms.Indexmusel]));
		}

		return Muon_weight;

	}
	static double ElectronSF(CMSAnalysis &cms, TH2D* ElectronEfficiencyHist)
	{
		VFloat_b(Electron_pt);
		VFloat_b(Electron_eta);

		double Electron_weight = cms.ScaleFactors(ElectronEfficiencyHist, Electron_pt[cms.Indexelsel], Electron_eta[cms.Indexelsel]);

		return Electron_weight;

	}
	static double TriggerSF(CMSAnalysis &cms, TH2D* TriggerHist)
	{
		VFloat_b(Muon_pt);
  		VFloat_b(Muon_eta);

		double Trigger_weight	= cms.ScaleFactors(TriggerHist, Muon_pt[cms.Indexmusel], fabs(Muon_eta[cms.Indexmusel]));

		return Trigger_weight;
	}
	static double PileUp(CMSAnalysis &cms, TH2D* Ratio)
	{
		Float_b(Pileup_nTrueInt);
		double pileup_weight = cms.PileupReweighting(Ratio, Pileup_nTrueInt);
		
		return pileup_weight;
	}
	static double GeneratorWeight(CMSAnalysis &cms)
	{
		Float_b(Generator_weight);
		double generator_weight = Generator_weight;	
		
		return generator_weight;
	}
};*/

int main(int argc, char** argv)
{
  	// Directory where files are sitting
	TString dir = "/eos/user/l/lurda/CMS/SamplesAfterSkimmingNovember2019/"; 
	
	// Initialize analysis structure
  	TopAnalysis cms; 

  	int maxevents = -1; 	

	//*****USING FUNCTION AddFiles ******

	//void CMSAnalysis::AddFiles(TString tipo, TString name, double luminosity, double xsection, TString path, int ntrees, double maxevents, double numberofevents)

	TString tipo = "Data";

  	cms.AddFiles(tipo, "DataA", 13.5*1000., -1, dir+"SingleMuon_Run2018A-Nano25Oct2019-v1_NANOAOD/", 555, -1, 227489240);
  	cms.AddFiles(tipo, "DataB", 6.8*1000.,  -1, dir+"SingleMuon_Run2018B-Nano25Oct2019-v1_NANOAOD/", 264, -1, 110446445);
  	cms.AddFiles(tipo, "DataC", 6.6*1000.,  -1, dir+"SingleMuon_Run2018C-Nano25Oct2019-v1_NANOAOD/", 264, -1, 107972995);
  	cms.AddFiles(tipo, "DataD", 32*1000.,   -1, dir+"SingleMuon_Run2018D-Nano25Oct2019-v1_NANOAOD/", 1238, -1, -1);
		
	tipo = "MCSignal";
	
  	cms.AddFiles(tipo,  "GG", -1, 48.58*0.1,  dir+"GluGlu_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunI0-v1_NANOAODSIM/", 1,  maxevents, 2000000);
 	cms.AddFiles(tipo, "VBF", -1, 3.782*0.1,  dir+"VBF_LFV_HToMuTau_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAu0-v1_NANOAODSIM/", 1,  maxevents, 1977300);
	
	tipo = "MCBckgr";

  	/*cms.AddFiles(tipo,"WW",    -1, 12.15488,  dir+"WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv0-v1_NANOAODSIM/",  5,  maxevents, 7758900);
  	cms.AddFiles(tipo,"WJets", -1, 61526.7,   dir+"WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv6-0-v1_NANOAODSIM/", 11,  maxevents, 70454125);
  	cms.AddFiles(tipo,"TW",    -1, 35.85,     dir+"ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIAutum1-v1_NANOAODSIM/",  3,  maxevents, 9598000);
  	cms.AddFiles(tipo,"TbarW", -1, 35.85,     dir+"ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_RunIIA1-v1_NANOAODSIM/",  4,  maxevents, 7623000);*/

  	cms.AddFiles(tipo,"TTbar", -1, 85.172,    dir+"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019/",        13,  maxevents, 4.622139340935600e+09);
  	cms.AddFiles(tipo,"DY",    -1, 2075.14*3, dir+"DYJetsToLL_M_50_TuneCP5_13TeV_amcatnloFXFX_pythia8_RunIIAutum2019/",                  29,  maxevents, 3.298665625309500e+12);
 
  	/*cms.AddFiles(tipo,"WZ",    -1, 27.6,      dir+"WZ_TuneCP5_13TeV-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_u0-v1_NANOAODSIM/",  4,  maxevents, 3885000);
  	cms.AddFiles(tipo,"ZZ",    -1, 12.14,     dir+"ZZ_TuneCP5_13TeV-pythia8_RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_u0-v1_NANOAODSIM/",  4,  maxevents, 1979000);  */	

  	// Initialize 1D histograms

	HISTS::Add(cms);
  	
  	// Loop on samples
  	unsigned int nsamples = cms.GetNumberOfSamples();

	ASCII txtFiles;
	txtFiles.txtCounter("NameOfHistograms.txt");

	TString mydir = "/afs/cern.ch/work/l/lurda/CMS/May_2019/ExoticHiggsDecay/codigojuan/GitLab_NANOAOD/LFVAnalysis/nanoaod/";

	//****PILE UP REWEIGHTING HISTOGRAM****
	TH1D* Ratio = nullptr; 
	Ratio = cms.ReadingFileAndGettingTH1Histogram(mydir+"python/pyroot/pileupweights.root", "RatioPU");

	//****MUON TIGHT ID EFFICIENCY HISTOGRAM*****
	TH2D* MuonTightIDEfficiencyHist = nullptr;
	MuonTightIDEfficiencyHist = cms.ReadingFileAndGettingTH2Histogram(mydir+"MuonEff_corr/RunABCD_SF_ID.root", "NUM_TightID_DEN_TrackerMuons_pt_abseta");

	//****MUON TIGHT ISO EFFICIENCY HISTOGRAM*****
	TH2D* MuonTightISOEfficiencyHist = nullptr;
	MuonTightISOEfficiencyHist = cms.ReadingFileAndGettingTH2Histogram(mydir+"MuonEff_corr/RunABCD_SF_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");
 
	//****ELECTRON ID SCALE FACTOR *****
	TH2D* ElectroScaleFactorHistogram = nullptr;
	ElectroScaleFactorHistogram = cms.ReadingFileAndGettingTH2Histogram(mydir+"ElecEff_corr/ElectronSF_passingMVA102Xwp90isoHWWiso0p06_2018runABCD.root", "ElectronSF");

	//****TRIGGER HLT_IsoMu24 SCALE FACTOR****
	TH2D* TriggerHLTIsoMu24ScaleFactorHistogram = nullptr;
	TriggerHLTIsoMu24ScaleFactorHistogram = cms.ReadingFileAndGettingTH2Histogram(mydir+"python/pyroot/Trigger_HLT_IsoMu24_weights.root", "HLT_IsoMu24_SFHist");
	

	//for (auto &xsmp : cms._SampleInfo)
  	for (unsigned int iSample=0; iSample<nsamples; ++iSample) 
	{
		// std::cout<<"Procesando muestra " << xsmp.GetSampleId() << " " << xsmp.NumberOfFilesInSample <<std::endl;
		std::cout<<"Procesando muestra " << iSample << " " << cms.GetSampleId(iSample) << " " <<cms.GetNumberOfFilesInSample(iSample)<<std::endl;
		
		for(unsigned int iFileInSample=1; iFileInSample<=cms.GetNumberOfFilesInSample(iSample); iFileInSample++)
		{
			std::cout << " LINE 381: iFileInSample " << iFileInSample << std::endl;
			
			// Set tree                 
			bool foundTree=true;                 
	                      
			if(!cms.SetTreeFile(iSample,iFileInSample)) foundTree=false;                 
			
			std::cout << "NEW CHECK " << cms.GetSampleId(iSample) << std::endl;			
			
			if (!foundTree) continue;
      			
			// Loop on events for the current sample
      			long long nevents = cms._NANOTREE->GetEntriesFast();

			//if (nevents>1000) nevents=1000;

			std::cout<<" LINE 395:  Reading entries in this file: "<<nevents<<std::endl;

      			for (long long iEvent=0; iEvent<nevents; iEvent++)
			{
          	 		// Set next entry
           	 		if (cms.SetEntry(iEvent)<0) {std::cout << "WARNING: PROBLEMS READING NTUPLES" << std::endl; break;}
 
           	 		// Preselect (summary branch must be read before)
           	 		if (!cms.Preselect()) continue;
          
           	 		// Select (event branch must be read before)
           	 		if (!cms.Select()) continue;

            			// Fill histograms
				double pileup_weight = 1.;

				double MuonTightIDEff_weight = 1.;
				double MuonTightISOEff_weight = 1.;

				double ElectronEff_weight = 1.0;
 
				double Trigger_weight = 1.0;

				double bTag_weight = 1.0;
		
				double generator_weight = 1.0;

				if(!cms.GetSampleId(iSample).Contains("Data"))
				{
					Float_b(Pileup_nTrueInt);
					VFloat_b(Muon_pt);
  					VFloat_b(Muon_eta);
  					VFloat_b(Jet_pt);
					VInt_b(Jet_hadronFlavour);
					Float_b(Generator_weight);
					VFloat_b(Electron_pt);
					VFloat_b(Electron_eta);


					
					pileup_weight = cms.PileupReweighting(Ratio, Pileup_nTrueInt);
					if(cms.GetSampleId(iSample).Contains("DY") )
					{
						generator_weight = Generator_weight;
					}
					else generator_weight = 1.0;
				
					if(Muon_pt[cms.Indexmusel] > 120.) 
					{
						MuonTightIDEff_weight = cms.ScaleFactors(MuonTightIDEfficiencyHist, 119., fabs(Muon_eta[cms.Indexmusel]));
						MuonTightISOEff_weight = cms.ScaleFactors(MuonTightISOEfficiencyHist, 119., fabs(Muon_eta[cms.Indexmusel]));
					}
					else 
					{
						MuonTightIDEff_weight = cms.ScaleFactors(MuonTightIDEfficiencyHist, Muon_pt[cms.Indexmusel], fabs(Muon_eta[cms.Indexmusel]));
						MuonTightISOEff_weight = cms.ScaleFactors(MuonTightISOEfficiencyHist, Muon_pt[cms.Indexmusel], fabs(Muon_eta[cms.Indexmusel]));
					}

					ElectronEff_weight = cms.ScaleFactors(ElectroScaleFactorHistogram, Electron_pt[cms.Indexelsel], Electron_eta[cms.Indexelsel]);

					Trigger_weight	= cms.ScaleFactors(TriggerHLTIsoMu24ScaleFactorHistogram, Muon_pt[cms.Indexmusel], fabs(Muon_eta[cms.Indexmusel]));
	
					bTag_weight = cms.bTagEventWeight(cms.IndexMedjetsel.size(), Jet_pt[cms.IndexMedjetsel[0]],Jet_hadronFlavour[cms.IndexMedjetsel[0]],Jet_pt[cms.IndexMedjetsel[1]], 						Jet_hadronFlavour[cms.IndexMedjetsel[1]], "comb","central" ,cms.Indexjetsel.size(), "medium");
					
				}
				//else if(cms.IndexMedjetsel.size() > 0) continue;


				//std::cout << "GENERATOR WEIGHT " << generator_weight << std::endl;
				double event_weight = pileup_weight * MuonTightIDEff_weight * MuonTightISOEff_weight * ElectronEff_weight * Trigger_weight * bTag_weight * generator_weight;
				//std::cout << "	EVENT WEIGHT " << event_weight << std::endl;
				//double event_weight_maria = pileup_weight_maria_up * MuonTightIDEff_weight * MuonTightISOEff_weight * ElectronEff_weight * Trigger_weight * bTag_weight;

				/*std::cout << "Pileup_weight " << pileup_weight << std::endl;
				std::cout << "		MuonTightIDEff_weight " << MuonTightIDEff_weight << std::endl;
				std::cout << "			MuonTightISOEff_weight " << MuonTightISOEff_weight << std::endl;
				std::cout << "				ElectronEff_weight " << ElectronEff_weight << std::endl;
				std::cout << "					Trigger_weight " << Trigger_weight << std::endl;
				std::cout << "						bTag_weight " << bTag_weight << std::endl;*/

				//std::cout << " 							PESO CALCULADO HASTA EL MOMENTO " << event_weight << std::endl << std::endl;

				HISTS::Fill(cms, cms._SampleInfo[iSample], event_weight);			
      			}			
		}
			
		HISTS::Save(cms, cms._SampleInfo[iSample], txtFiles);

  	}  	

	//To see things interactively (commnent otherwise) 

  	auto app = new TRint("Top Analysis", &argc, argv); 

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
	Bool_b(HLT_IsoMu24);
  	//printf("nmu: %d, nel: %d\n", nMuon, nElectron);

  	if (nMuon+nElectron<2) 	return false; // look for events with >= 2 leptons
	if(!HLT_IsoMu24)return false;
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

	Indexjetsel.clear();
	IndexMedjetsel.clear();

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
		Indexjetsel.push_back(ij);

	}
	
	for(auto &ij : Indexjetsel)
	{
		//std::cout << "		nJet " << nJet << " nJetSel "  << njtsel << std::endl;
      		if (Jet_btagDeepB[ij]>btagMediumCut) 
		{
           		nbMediumSel++;
			IndexMedjetsel.push_back(ij);
            		//nbLooseSel++;
      		} 

		/*else if (Jet_btagDeepB[ij]>btagLooseCut) 
		{
           		 //nbLooseSel++;
      		}*/
		
		//std::cout << " Item: " << ij << " njtsel in the loop " << nbMediumSel << " " << nbLooseSel << std::endl;
  	}

	//std::cout << " NJTSEL when the loop is over " << njtsel << " nJets iniciales "  << nJet << std::endl;
  
  	if (Indexjetsel.size()>2) 		return false; 
		
	//if((Indexjetsel.size() == 0) && (IndexMedjetsel.size() == 0 )) return true;

	//if((Indexjetsel.size() == 1) && (IndexMedjetsel.size() == 0 || IndexMedjetsel.size() == 1) ) return true;
	//if((Indexjetsel.size() == 2) && (IndexMedjetsel.size() == 0 || IndexMedjetsel.size() == 1 || IndexMedjetsel.size() == 2) ) return true;

	//if(Indexjetsel.size() == 1 && (IndexMedjetsel.size() == 0 || IndexMedjetsel.size() == 1)) return true; //std::cout << " 1 selected jet and it is btagged " << std::endl; //ttbar control region
	if((Indexjetsel.size() == 2) && (IndexMedjetsel.size() == 2)) return true;
	//else return false;
	
	//std::cout << " 2 selected jets and both of them are btagged " << std::endl; //ttbar control region
	//if(njtsel == 2 && nbMediumSel == 1) std::cout << "2 selected jets and one of them is btagged" << std::endl; //ttbar control region

  	//printf("nmu: %d, nel: %d\n", nmu, nel);
	return false;
}
