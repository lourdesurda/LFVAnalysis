#
#	First python macro get QCD background
#

import ROOT
import numpy as np 
from ROOT import TFile, TTree

ROOT.objs = []

###
def Getting_QCD_Histograms ():
	"""Routine to get QCD background"""

	directory = '../../TOPLOTQCD/TOPLOT_SS_IM_0jets_0bjets/'

	DataList = ['DataA', 'DataB', 'DataC', 'DataD']
	MCSignalList = ['GG', 'VBF']
	MCBckgList = ['ZZ', 'WZ', 'WW', 'WJets', 'TW', 'TbarW', 'TTbar', 'DY']
	MCBckgMergingWJetsList = ['WJets', 'WJets1', 'WJets2', 'WJets3', 'WJets4']

	TotalList = DataList + MCSignalList + MCBckgList + MCBckgMergingWJetsList

	#f = open('../../NameOfHistograms.txt', "r")
	
	#Openning the files to get the histogram

	Hin = dict()   # Dictionario: {}

	for sample in TotalList:
		try:
			Fin = ROOT.TFile.Open(directory+sample+'_ToPlot_QCD_SS_IM_0jets_0bjets.root')
			#Fin = ROOT.TFile.Open(directory+sample+'_ToPlot_QCD_SS_IM_1jet_0bjets.root')
			#Fin = ROOT.TFile.Open(directory+sample+'_ToPlot_QCD_SS_IM_2jets_0bjets.root')
			
		except:
			print('Could not open file')
			raise

		try:
			print(sample)
			ROOT.gDirectory.cd('PyROOT:/')
			Hin[sample] = Fin.Get(sample+'_hInvariantMassMuonElectron').Clone()

		except:
			print('Could not get the histogram '+sample)
			raise

		Fin.Close()

	##### MERGING WJETS SAMPLES ####

	hMerge = []
	msamp = 1 
	for sample in MCBckgMergingWJetsList:
		if(len(hMerge)==0):
			hMerge.append(Hin[sample])
		else:
			hMerge.append(Hin[sample]+ hMerge[msamp-1])
			msamp = msamp + 1

	print("Merging", len(hMerge))
	
	hMerge[-1].SetName("WJets_hInvariantMassMuonElectron")

	print(hMerge[-1])
	print(Hin["WJets"])

	Samples = {"GG":Hin["GG"],
		  "VBF":Hin["VBF"],
		  "WW": Hin["WW"], 
		  "WJets":hMerge[-1], 
		  "TW":Hin["TW"], 
		  "TbarW":Hin["TbarW"], 
		  "TTbar":Hin["TTbar"], 
		  "DY":Hin["DY"], 
		  "ZZ":Hin["ZZ"], 
		  "WZ":Hin["WZ"]}
	
	##### MC BACKGROUNDS #####

	hSumaBck = []
	xsmp = 1
	for sample in MCBckgList:
		if(len(hSumaBck)==0):
			hSumaBck.append(Samples.get(sample))
		else: 
			hSumaBck.append(Samples.get(sample) + hSumaBck[xsmp-1])
			xsmp = xsmp + 1

	print(len(hSumaBck))
	hSumaBck[-1].SetName("MCBck_hInvariantMassMuonElectron")		


	###### DATA #######
	
	hSumData = []
	xdata = 1
	for sample in DataList:
		if(len(hSumData)==0):
			hSumData.append(Hin[sample])
		else:
			hSumData.append(Hin[sample]+hSumData[xdata-1])
			xdata = xdata + 1

	print(len(hSumData))
	hSumData[-1].SetName("Data_hInvariantMassMuonElectron")
	
	###### QCD #######

	hQCD = hSumData[-1].Clone()
	hQCD.Add(hSumaBck[-1],-1)

	#hQCD.SetMinimum(0)
	hQCD.SetName("QCD_hInvariantMassMuonElectron")

	print(hQCD.GetName())
	f = TFile("QCD_ToPlot_V6_OS_IM_0jets_0bjets.root","UPDATE");
	hQCD.Write("QCD_hInvariantMassMuonElectron")

#Code to be run in order to produce the plots
if __name__ == '__main__':

    	Getting_QCD_Histograms()
	
	raw_input("Press ENTER... ")

