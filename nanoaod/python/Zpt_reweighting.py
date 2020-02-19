#
#	First python macro to plot the Missing ET using pyROOT
#When it says luminisity it means 1/luminosity

import ROOT
from ROOT import TTree, TH1D, TH2D, TFile, TCanvas
import struct
import numpy as np 
ROOT.objs = []

### PATHS ###
directory = '/eos/user/l/lurda/CMS/SamplesToMergeJanuary2020_madgraphMLM/'

DY_path  = directory+'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoA0-v1_NANOAODSIM/skimmed-nano_17.root'
DY1_path = directory+'DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM/skimmed-nano_63.root'
DY2_path = directory+'DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM/skimmed-nano_7.root'
DY3_path = directory+'DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM/skimmed-nano_6.root'
DY4_path = directory+'DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18Nano0-v1_NANOAODSIM/skimmed-nano_5.root'

### LISTS ####

List = ['DY', 'DY1', 'DY2', 'DY3', 'DY4']

Jets = [0, 1, 2, 3, 4]

Paths = [DY_path, DY1_path, DY2_path, DY3_path, DY4_path]

cross_sections_LO = [5343., 877.8, 304., 111.5, 44.03]

Colors = [ROOT.kBlue+1, ROOT.kYellow+1, ROOT.kPink+1, ROOT.kGreen+1, ROOT.kOrange+1]

### DICTIONARIES ###
Paths_dict = dict(zip(List, Paths))

Cross_sections_LO_dict = dict(zip(List, cross_sections_LO))

Colors_dict = dict(zip(List, Colors))

merge_weight = {}

original_generator_weight = {}

inclusive_weight_dict = dict(zip(Jets, List))
### CONSTANTS ###

cross_section_NNLO = 6077.22 #pb

LOToNNLOfactor = (cross_section_NNLO/Cross_sections_LO_dict.get("DY"))

Analyzed_luminosity = 59.266*1000.#1/pb

Total_merge_weight = 0.0

###
def Zpt_reweighting():
	"""Routine to calculate the Zpt_reweighting"""


	eventcount = {}
	original_numberofevents = {}
	Histo_pt = {}
	Histo_mass = {}

	for sample in List:
		Fin = ROOT.TFile.Open(Paths_dict.get(sample), 'read')
		ROOT.gDirectory.cd('PyROOT:/')

		#First of all I count the number of events per root file to know how many interesting event I do have in them
		eventcount[sample] = Fin.Get('eventcount').Clone()
		original_numberofevents[sample] = eventcount[sample].GetBinContent(0)

		#Now I define the histogrmas I wanna plot
		Histo_pt[sample] = ROOT.TH1D("Zpt_"+sample, "Distribution of the Zpt for the sample "+sample, 101, -0.5, 100.5)
		Histo_mass[sample] = ROOT.TH1D("Zmass_"+sample, "Distribution of the Z mass for the sample "+ sample, 201, -0.5, 200.5)

		#Now I get the branches of the tree that I need 
		mytree = Fin.Get('Events')
		mytree.SetBranchStatus('*', 0)
		mytree.SetBranchStatus('GenPart_pdgId', 1)
		mytree.SetBranchStatus('GenPart_mass', 1)
		mytree.SetBranchStatus('GenPart_pt', 1)
		mytree.SetBranchStatus('GenPart_status', 1)
		mytree.SetBranchStatus('GenPart_statusFlags', 1)
		mytree.SetBranchStatus('GenPart_genPartIdxMother', 1)
		mytree.SetBranchStatus('nGenPart', 1)
		entries = mytree.GetEntriesFast()
		print('Number of entries in sample '+sample, entries)

		eventcounter = -1

		ZFound = 0
		ZNotFound = 0
		ZFound_With_Status = 0 
		ZFound_Without_Status = 0 
		ZFound_With_Status_but_isNotLastCopy = 0
		ZFound_With_Status_and_isLastCopy = 0
		ZFoundEvent=0

		MuonPairFound = 0
		ElectronPairFound = 0
		TauonPairFound = 0

		MuonFound = 0
		ElectronFound = 0
		TauonFound = 0
		AntiMuonFound =0
		AntiElectronFound = 0
		AntiTauonFound = 0

		for event in mytree:
			eventcounter = eventcounter + 1
			#print('Number of nGenPart in the current event ', eventcounter, " is " , event.nGenPart)
			foundZ=False
			for particle in range(0, event.nGenPart):

				if(event.GenPart_pdgId[particle] == 23): 

					ZFound = ZFound + 1
					foundZ=True

					if(event.GenPart_status[particle] == 62): 

						ZFound_With_Status = ZFound_With_Status + 1
						binaryflag = np.binary_repr(event.GenPart_statusFlags[particle], width=15)

						if(binaryflag[13] == "1" or binaryflag[14] == "1"): #isLastCopy() == true

							ZFound_With_Status_and_isLastCopy = ZFound_With_Status_and_isLastCopy + 1
							Histo_pt[sample].Fill(event.GenPart_pt[particle])
							Histo_mass[sample].Fill(event.GenPart_mass[particle])

						else: ZFound_With_Status_but_isNotLastCopy = ZFound_With_Status_but_isNotLastCopy + 1

					else: ZFound_Without_Status = ZFound_Without_Status + 1

				if(abs(event.GenPart_pdgId[particle]) == 13 or abs(event.GenPart_pdgId[particle]) == 15 or abs(event.GenPart_pdgId[particle]) == 11):
		
					if(event.GenPart_genPartIdxMother[particle] == -1): continue

					if(event.GenPart_pdgId[event.GenPart_genPartIdxMother[particle]] == 23 and event.GenPart_status[event.GenPart_genPartIdxMother[particle]] == 62):

							if(event.GenPart_pdgId[particle]==13): MuonPairFound = MuonPairFound+1
							if(event.GenPart_pdgId[particle]==11): ElectronPairFound = ElectronPairFound+1
							if(event.GenPart_pdgId[particle]==15): TauonPairFound = TauonPairFound+1
							#print('Found a', event.GenPart_pdgId[particle], " in event ", eventcounter)

		   	if not foundZ: 
				ZNotFound = ZNotFound + 1
				if(abs(event.GenPart_pdgId[particle]) == 13 or abs(event.GenPart_pdgId[particle]) == 15 or abs(event.GenPart_pdgId[particle]) == 11):
		
					if(event.GenPart_genPartIdxMother[particle] == -1): continue
					#event.GenPart_pdgId[event.GenPart_genPartIdxMother[particle]] != event.GenPart_pdgId[particle] and
					if(event.GenPart_genPartIdxMother[particle] == 0):

						if(event.GenPart_pdgId[particle]==13): MuonFound = MuonFound+1
						if(event.GenPart_pdgId[particle]==11): ElectronFound = ElectronFound+1
						if(event.GenPart_pdgId[particle]==15): TauonFound = TauonFound+1
						if(event.GenPart_pdgId[particle]==-13): AntiMuonFound = AntiMuonFound+1
						if(event.GenPart_pdgId[particle]==-11): AntiElectronFound = AntiElectronFound+1
						if(event.GenPart_pdgId[particle]==-15): AntiTauonFound = AntiTauonFound+1

					#print('In NOT_ZFOUND: Found a', event.GenPart_pdgId[particle], " in event ", eventcounter)
			else: ZFoundEvent = ZFoundEvent +1
			
		print('_________eventcounter', eventcounter)
		#print('_________ZFound', ZFound)
		print('_________ZFoundEvent', ZFoundEvent, ZFoundEvent*100/eventcounter,'%')
		print('_________ZNotFoundEvent', ZNotFound, ZNotFound*100/eventcounter,'%')
		print('_________ZFound_With_Status', ZFound_With_Status, ZFound_With_Status*100/ZFoundEvent,'%')
		#print('_________ZFound_Without_Status', ZFound_Without_Status, ZFound_Without_Status*100/ZFoundEvent,'%')
		print('_________ZFound_With_Status_but_isNotLastCopy', ZFound_With_Status_but_isNotLastCopy, ZFound_With_Status_but_isNotLastCopy*100/ZFoundEvent,'%')
		print('_________ZFound_With_Status_and_isLastCopy', ZFound_With_Status_and_isLastCopy, ZFound_With_Status_and_isLastCopy*100/ZFoundEvent,'%')
		print('_________PARTICLES COUNT (%)')
		print('_________ZFound_MuonPairFound', MuonPairFound, MuonPairFound*100/ZFoundEvent,'%')
		print('_________ZFound_ElectronPairFound', ElectronPairFound, ElectronPairFound*100/ZFoundEvent,'%')
		print('_________ZFound_TauonPairFound', TauonPairFound, TauonPairFound*100/ZFoundEvent,'%')
		print('_________ZNotFoundEvent_MuonFound', MuonFound, MuonFound*100/ZNotFound,'%')
		print('_________ZNotFoundEvent_AntiMuonFound', AntiMuonFound, AntiMuonFound*100/ZNotFound,'%')
		print('_________ZNotFoundEvent_ElectronFound', ElectronFound, ElectronFound*100/ZNotFound,'%')
		print('_________ZNotFoundEvent_AntiElectronFound', AntiElectronFound, AntiElectronFound*100/ZNotFound,'%')
		print('_________ZNotFoundEvent_TauonFound', TauonFound, TauonFound*100/ZNotFound,'%')
		print('_________ZNotFoundEvent_AntiTauonFound', AntiTauonFound, AntiTauonFound*100/ZNotFound,'%')


		

		'''canvas1 = ROOT.TCanvas("cv", "Zpt", 900, 900)	
		canvas1.Draw()
		Histo_pt['DY'].Draw("hist")

		canvas2 = ROOT.TCanvas("cv", "Number of jets", 900, 900)	
		canvas2.Draw()
		Histo_mass['DY'].Draw("hist")
		ROOT.objs.append([canvas1,canvas2, Histo_mass,Histo_pt])'''
		



#Code to be run in order to produce the plots
if __name__ == '__main__':

	import sys

    	Zpt_reweighting()
	
	raw_input("Press ENTER... ")

