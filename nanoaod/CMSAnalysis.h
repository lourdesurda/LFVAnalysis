#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>
#include "SamplesClassLU.h"

/// Complete analysis class to do CMS analysis using CMSTree2016 root files.
class CMSAnalysis {
  public:
      CMSAnalysis();
      virtual ~CMSAnalysis();
	
      /// Initialize and add data sample to the analysis (always first sample to add)
      //void AddDataSample(const TString& id, const TString& file, double luminosity, int nfiles=1);
      void AddDataSampleFiles(const TString& id,  const TString& dir, const TString& file, double luminosity, int nfiles=1);
      //void AddDataSampleFiles(const string& id,  const string& dir, const string& file, double luminosity, int nfiles=1);
      /// Initialize and add the MC signal sample
     //void AddMCSignalSample(const TString& id, const TString& file, int maxevents, double xsection, int total_events_for_xsection=-1);
	//void AddMCSignalSample(const TString& id, const TString& file, int maxevents, double xsection, int total_events_for_xsection);
      /// Initialize and add another MC background sample. All samples with the same "id" are added to the same histogram components at the end (stacked)
      //void AddMCSample(const TString& id, const TString& file, int maxevents, double xsection, int total_events_for_xsection=-1);
      void AddMCSampleFiles(const TString& id, const TString& dir, const TString& file, int maxevents, double xsec, int total_events_for_xsection=-1, int nfiles=1);
	void AddMCSignalSampleFiles(const TString& id, const TString& dir, const TString& file, int maxevents, double xsec, int total_events_for_xsection, int nfiles) ;
    //  void ScalePlots();
	void AddFiles(TString tipo, TString name, double luminosity, double xsection, TString path, int ntrees, double maxevents, double numerofevents);
      /// Set tree to the sample with index i
      void SetTree(int i);
      bool SetTreeFile(int i, int fileJ); 

	TH1D* ReadingFileAndGettingTH1Histogram(TString path, TString histname);

	TH2D* ReadingFileAndGettingTH2Histogram(TString path, TString histname);
      /// Access and set the entry entryNumber from the current tree. It returns 0 if OK, and -1 if there is an error or we are beyond limits. Note that it does not read in the contents of the event branches. This behavior is useful if branches are going to be read individually later. To also read all active branches, use GetEntry(entryNumber) instead
      Int_t SetEntry(int entryNumber);
	double GettingSF_bTag(const TString& DeepCSV, const TString& WP , const TString& SysType, int Flavour, float pt);
	double bTagEventWeight(int nBtaggedJets, float bjetpt_1, int bjetflavour_1, float bjetpt_2, int bjetflavour_2, const TString& WP,const TString& SysType, int nBTags, const TString& DeepCSV);
      /// Access, set and read in the entry entryNumber from the current tree. It returns 0 if OK, and -1 if there is an error or we are beyond limits.
   //   Int_t GetEntry(int entryNumber);

      /// Get entry number of the tree currently in use
      Int_t GetCurrentEntryNumber(){return _NANOTREE->GetReadEntry();};

      /// True if the current file was declared to be real data
      bool isDataFile(){return _currentIndex==_dataIndex;};
      /// True if the current file was declared to be MC
      bool isMCFile(){return _currentIndex!=_dataIndex;};
      /// True if the current file was declared to be the signal MC
      bool isSignalFile(){return _signalFlag[_currentIndex]==true;};
      /// Get luminosity
      double GetLumi(){return _lumi;};
      /// Get the index to the data sample
      int GetDataIndex() {return _dataIndex;};
      /// Get pointer to current tree
      TTree* GetTree() {return _NANOTREE;};
      /// Get pointer to current file
      TFile* GetFile() {return _file;};
      /// Get the index of the sample currently in use
      int GetIndex() {return _currentIndex;};

      /// Get number of samples being processed
      int GetNumberOfSamples() {return _sampleId.size();};
 	/// Get number of files per sample being processed       
	int GetNumberOfFilesInSample(int i) {return _sampleNFiles[i];};
      /// Get number of events to be read in sample with index i
      int GetNumberOfEvents(int i) {return _sampleNevents[i];};
      /// Get sampleId of sample with index i
      TString GetSampleId(int i) {return _sampleId[i];};
      /// Get file name of sample with index i
      TString GetSampleFile(int i) {return _sampleFile[i];};
      /// Get cross section for sample with index i (in data this gives number_of_events/luminosity)
      double GetSampleXsection(int i) {return _sampleXsection[i];};
      /// Get effective luminosity provided by sample with index i (in data this gives just the declared luminosity)
      double GetSampleEquivalentLuminosity(int i) {return _lumi/_sampleWeight[i];};
      /// Get sample weight for sample with index i (to normalize it to the total integrated luminosity)
      double GetSampleWeight(int i) {return _sampleWeight[i];};

      /// Book 1D Plot for data and all bkgd components. It follows TH1 conventions
      void AddPlot1D(const TString& name, const TString& title, int nbins, double xmin, double xmax);
      /// Book 1D Plot for just one component
      //void AddPlot1D_bare(const TString& name, const TString& title, int nbins, double xmin, double xmax);
      /// Fill histogram for the 1D Plot. It follows TH1 conventions.
      void FillPlot1D(const TString& name, const SAMPLES &sample, double value, double weight=1.);
	double PileupReweighting(const TH1D* Ratio, const float Pileup_nTrueInt);
	double ScaleFactors(const TH2D* SFHistogram, float lepton_variable1, float lepton_variable2);
	///
	void SavingHistograms(const SAMPLES &sample, const TString& name, const TString& option);
      /// Fill histogram for the 1D Plot with only one component
      //void FillPlot1D_bare(const TString& name, double value, double weight=1.);
      /// Draw 1D plot with data and all bckg components. It follows TH1 conventions. It also produces .pdf, .png and .root versions of the histogram, with an optional suffix (to avoid overwriting other plots)
      void DrawPlot1D(const TString& name, const TString& suffix="", const TString& dir="png_00_generatorweight");

      /// Draw 1D plot for only one component
      //void DrawPlot1D_bare(const TString& name, const TString& suffix="", const TString& dir="");
	
      /// Normalization utility for histograms
    //  void Normalize_Histogram(const TString& name);
      /// Normalization utility for histogram (one component)
    //  void Normalize_Histogram_bare(const TString& name);
      /// Divide bare histograms
     // void Divide_Histogram_bare(const TString& name1, const TString& name2);
      /// Get an individual MC histogram
    //  TH1D* GetMCHistogram(const TString& name);
      /// Get an individual Data histogram
    //  TH1D* GetDataHistogram(const TString& name);
	long long int NumberOfEntries(const TString& path, int nfiles);
      /// True if the current event is a real data event
      bool isDataEvent(){return (!isMCFile());};
      /// True if the current event is a MC event
      bool isMCEvent(){return (isMCFile());};

      static TTree* _NANOTREE;
      static std::vector<TBranch*> _BUFFERBRANCHES;
      static std::vector<Char_t> _BUFFER;

      /// Get a NANOAOD variable of type T. For instance:
      /**   \verbatim 
                  auto nmu = Get<UInt_t>("nMuon");
                  printf("There are %ul muons in this event\n", nmu);
            \endverbatim 
       *    You can also access a given object in an array:
       *    \verbatim 
                  auto nmu = Get<UInt_t>("nMuon");
                  for (auto im=0; im<nmu; im++) {
                        auto ptmu = Get<Float_t>("Muon_pt",im);
                        printf("Muon %ul, pt[%ul]=%f\n", im, im, ptmu);
                  }
            \endverbatim 
       *  The method is optimized to read branch contents only the first time that is called in the event.
       */
      template <typename T> static T& Get(const std::string& branchName, unsigned int index=0) {

            TBranch* b = NULL;
            unsigned int bufferBranchesSize = _BUFFERBRANCHES.size();
            for (unsigned int i=0 ; i<bufferBranchesSize; i++) {
                  if (strcmp(_BUFFERBRANCHES[i]->GetName(),branchName.data())==0) {
                        b = _BUFFERBRANCHES[i];
                        break;
                  }
            }

            if (!b) {
                  b = _NANOTREE->GetBranch(branchName.data());
                  try {
                        if (!b) throw std::exception();
                  } catch (const std::exception&) {
                        printf("\n>>>>\n>>>> ERROR: Branch \"%s\" does not exist, throwing exception !!!\n>>>>\n\n", branchName.data());
                  }

                  _BUFFERBRANCHES.push_back(b);

                  // Strict test of data validity
                  TClass* clptr; EDataType this_type; 
                  b->GetExpectedType(clptr,this_type);
                  if (this_type!=TDataType::GetType(typeid(T))) {
                        TString this_chtype = TDataType::GetTypeName(this_type);
                        TString requested_chtype = TDataType::GetTypeName(TDataType::GetType(typeid(T)));
                        bool still_OK = (this_chtype=="UInt_t" && requested_chtype=="Int_t");
                        if (!still_OK) {
                              printf("\n>>>>\n>>>> WARNING IN BRANCH \"%s\": TypeId in root tree is \"%s\", which is not the requested one \"%s\". This can lead to problems!!!\n>>>>\n\n", branchName.data(), this_chtype.Data(), requested_chtype.Data());
                        }
                  }
            }

            Int_t new_entry = _NANOTREE->GetReadEntry();
            if (b->GetReadEntry()!=new_entry) b->GetEntry(new_entry);

            return reinterpret_cast<T &>(*(b->GetAddress()+index));

      };

      /// Get a NANOAOD array of type T. For instance:
      /**   \verbatim 
                  auto nmu = Get<UInt_t>("nMuon");
                  auto ptmu_array = GetP<Float_t>("Muon_pt");
                  for (auto im=0; im<nmu; im++) {
                        auto ptmu = ptmu_array[im];
                        printf("Muon %ul, pt[%ul]=%f\n", im, im, ptmu);
                  }
            \endverbatim 
       *    You can also access a given pointer to an object of the array:
       *    \verbatim 
                  auto nmu = Get<UInt_t>("nMuon");
                  for (auto im=0; im<nmu; im++) {
                        auto ptmu_pointer = GetP<Float_t>("Muon_pt",im);
                        printf("Muon %ul, pt[%ul]=%f\n", im, im, *ptmu_pointer);
                  }
            \endverbatim 
       *  The method is optimized to read branch contents only the first time that is called in the event.
       */
      template <typename T> static T* GetP(const std::string& branchName, unsigned int index=0) {

            TBranch* b = NULL;
            unsigned int bufferBranchesSize = _BUFFERBRANCHES.size();
            for (unsigned int i=0 ; i<bufferBranchesSize; i++) {
                  if (strcmp(_BUFFERBRANCHES[i]->GetName(),branchName.data())==0) {
                        b = _BUFFERBRANCHES[i];
                        break;
                  }
            }

            if (!b) {
                  b = _NANOTREE->GetBranch(branchName.data());
                  try {
                        if (!b) throw std::exception();
                  } catch (const std::exception&) {
                        printf("\n>>>>\n>>>> ERROR: Branch \"%s\" does not exist, throwing exception !!!\n>>>>\n\n", branchName.data());
                  }

                  _BUFFERBRANCHES.push_back(b);

                  // Strict test of data validity
                  TClass* clptr; EDataType this_type; 
                  b->GetExpectedType(clptr,this_type);
                  if (this_type!=TDataType::GetType(typeid(T))) {
                        TString this_chtype = TDataType::GetTypeName(this_type);
                        TString requested_chtype = TDataType::GetTypeName(TDataType::GetType(typeid(T)));
                        bool still_OK = (this_chtype=="UInt_t" && requested_chtype=="Int_t");
                        if (!still_OK) {
                              printf("\n>>>>\n>>>> WARNING IN BRANCH \"%s\": TypeId in root tree is \"%s\", which is not the requested one \"%s\". This can lead to problems!!!\n>>>>\n\n", branchName.data(), this_chtype.Data(), requested_chtype.Data());
                        }
                  }
            }

            Int_t new_entry = _NANOTREE->GetReadEntry();
            if (b->GetReadEntry()!=new_entry) b->GetEntry(new_entry);

            return reinterpret_cast<T *>(b->GetAddress()+index);

      };

      /// Set a NANOAOD variable of type T. It returns true if the branch exists, and false otherwise. Examples:
      /**   \verbatim 
                  UInt_t nmu; bool branch_exists = Set("nMuon", nmu);
                  if (branch_exists) printf("There are %ul muons in this event\n", nmu);
            \endverbatim 
       *    You can also access a given object in an array:
       *    \verbatim 
                  UInt_t nmu; Set("nMuon", nmu);
                  for (auto im=0; im<nmu; im++) {
                        Float_t ptmu; Set("Muon_pt", ptmu, im);
                        printf("Muon %ul, pt[%ul]=%f\n", im, im, ptmu);
                  }
            \endverbatim 
       *  The method is optimized to read branch contents only the first time that is called in the event.
       */
      template <typename T> static bool Set(const std::string& branchName, T& var, unsigned int index=0) {

            if (!_NANOTREE) return false;

            TBranch* b = NULL;
            unsigned int bufferBranchesSize = _BUFFERBRANCHES.size();
            for (unsigned int i=0 ; i<bufferBranchesSize; i++) {
                  if (strcmp(_BUFFERBRANCHES[i]->GetName(),branchName.data())==0) {
                        b = _BUFFERBRANCHES[i];
                        break;
                  }
            }

            if (!b) {
                  b = _NANOTREE->GetBranch(branchName.data());
                  if (!b) return false;

                  _BUFFERBRANCHES.push_back(b);

                  // Strict test of data validity
                  TClass* clptr; EDataType this_type; 
                  b->GetExpectedType(clptr,this_type);
                  if (this_type!=TDataType::GetType(typeid(T))) {
                        TString this_chtype = TDataType::GetTypeName(this_type);
                        TString requested_chtype = TDataType::GetTypeName(TDataType::GetType(typeid(T)));
                        bool still_OK = (this_chtype=="UInt_t" && requested_chtype=="Int_t");
                        if (!still_OK) {
                              printf("\n>>>>\n>>>> WARNING IN BRANCH \"%s\": TypeId in root tree is \"%s\", which is not the requested one \"%s\". This can lead to problems!!!\n>>>>\n\n", branchName.data(), this_chtype.Data(), requested_chtype.Data());
                        }
                  }
            }

            Int_t new_entry = _NANOTREE->GetReadEntry();
            if (b->GetReadEntry()!=new_entry) b->GetEntry(new_entry);

            var = reinterpret_cast<T &>(*(b->GetAddress()+index));

            return true;

      };

      /// Set a NANOAOD array of type T. It returns true if the branch exists, and false otherwise. Examples:
      /**   \verbatim 
                  UInt_t nmu; Set("nMuon", nmu);
                  Float_t* ptmu_array; bool branch_exists = SetP("Muon_pt", ptmu_array);
                  if (branch_exists) {
                        for (auto im=0; im<nmu; im++) {
                              auto ptmu = ptmu_array[im];
                              printf("Muon %ul, pt[%ul]=%f\n", im, im, ptmu);
                        }
                  }
            \endverbatim 
       *    You can also access a given pointer to an object of the array:
       *    \verbatim 
                  UInt_t nmu; Set("nMuon", nmu);
                  for (auto im=0; im<nmu; im++) {
                        Float_t* ptmu_pointer; SetP("Muon_pt", ptmu_pointer, im);
                        printf("Muon %ul, pt[%ul]=%f\n", im, im, *ptmu_pointer);
                  }
            \endverbatim 
       *  The method is optimized to read branch contents only the first time that is called in the event.
       */
      template <typename T> static bool SetP(const std::string& branchName, T*& pvar, unsigned int index=0) {

            if (!_NANOTREE) return false;

            TBranch* b = NULL;
            unsigned int bufferBranchesSize = _BUFFERBRANCHES.size();
            for (unsigned int i=0 ; i<bufferBranchesSize; i++) {
                  if (strcmp(_BUFFERBRANCHES[i]->GetName(),branchName.data())==0) {
                        b = _BUFFERBRANCHES[i];
                        break;
                  }
            }

            if (!b) {
                  b = _NANOTREE->GetBranch(branchName.data());
                  if (!b) return false;

                  _BUFFERBRANCHES.push_back(b);

                  // Strict test of data validity
                  TClass* clptr; EDataType this_type; 
                  b->GetExpectedType(clptr,this_type);
                  if (this_type!=TDataType::GetType(typeid(T))) {
                        TString this_chtype = TDataType::GetTypeName(this_type);
                        TString requested_chtype = TDataType::GetTypeName(TDataType::GetType(typeid(T)));
                        bool still_OK = (this_chtype=="UInt_t" && requested_chtype=="Int_t");
                        if (!still_OK) {
                              printf("\n>>>>\n>>>> WARNING IN BRANCH \"%s\": TypeId in root tree is \"%s*\", which is not the requested one \"%s*\". This can lead to problems!!!\n>>>>\n\n", branchName.data(), this_chtype.Data(), requested_chtype.Data());
                        }
                  }
            }

            Int_t new_entry = _NANOTREE->GetReadEntry();
            if (b->GetReadEntry()!=new_entry) b->GetEntry(new_entry);

            pvar = reinterpret_cast<T *>(b->GetAddress()+index);

            return true;

      };

 private:
      double _lumi; ///<
      std::vector<TString> _sampleFile; ///<
      std::vector<TString> _sampleId; ///<
      std::vector<int> _sampleNevents; ///<
      //std::vector<int> _sampleNeventsTree; ///<
      std::vector<int> _nDump10; ///<
      std::vector<double> _sampleXsection; ///<
      std::vector<double> _sampleWeight; ///<
      std::vector<bool> _signalFlag; ///<
      std::vector<int> _sampleNFiles; ///<
      int _dataIndex; ///<
      std::vector<int> _sampleNeventsSumFiles;
      double _totalLuminosity;
      int _currentIndex; ///<
    
	TFile* _file; ///<

      	std::vector<TH1D*> hists_1D; ///<
	
public:
	std::vector<SAMPLES> _SampleInfo;
};

#define Int_b(a) Int_t a; if (!CMSAnalysis::Set(#a,a)) {printf("Branch %s does not exists!!\n", #a); throw std::exception();};
#define UInt_b(a) UInt_t a; if (!CMSAnalysis::Set(#a,a)) {printf("Branch %s does not exists!!\n", #a); throw std::exception();};
#define Float_b(a) Float_t a; if (!CMSAnalysis::Set(#a,a)) {printf("Branch %s does not exists!!\n", #a); throw std::exception();};
#define Bool_b(a) Bool_t a; if (!CMSAnalysis::Set(#a,a)) {printf("Branch %s does not exists!!\n", #a); throw std::exception();};
#define Ulong64_b(a) ULong64_t a; if (!CMSAnalysis::Set(#a,a)) {printf("Branch %s does not exists!!\n", #a); throw std::exception();};

#define VInt_b(a) Int_t* a; if (!CMSAnalysis::SetP(#a,a)) {printf("Branch %s does not exists!!\n", #a); throw std::exception();};
#define VFloat_b(a) Float_t* a; if (!CMSAnalysis::SetP(#a,a)) {printf("Branch %s does not exists!!\n", #a); throw std::exception();};
#define VBool_b(a) Bool_t* a; if (!CMSAnalysis::SetP(#a,a)) {printf("Branch %s does not exists!!\n", #a); throw std::exception();};
#define String_b(a) UChar_t* a; if (!CMSAnalysis::SetP(#a,a)) {printf("Branch %s does not exists!!\n", #a); throw std::exception();};

