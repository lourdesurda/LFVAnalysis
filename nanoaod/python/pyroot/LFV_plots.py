#
#	First python macro to plot the Missing ET using pyROOT
#

import ROOT
from ROOT import TFile, TTree
ROOT.objs = []

Hin = dict()

#directory = '../../Plotting_14Feb2020/'
#directory = '../../Plotting_29Feb2020_26GeVmuonptcut/'
#directory = '../../Plotting_3March2020_FullSelection/'
directory = '../../Plotting_5March2020_vetos/'

DataList = ['DataA', 'DataB', 'DataC', 'DataD']
MCSignalList = ['ggHMuTau', 'VBFHMuTau']

MCBckgMergingWJetsList = ['WJetsInclusiveMerge', 'WJets1Merge', 'WJets2Merge', 'WJets3Merge', 'WJets4Merge']

MCBckgMergingDYList = ['DYInclusiveMerge', 'DY1Merge', 'DY2Merge', 'DY3Merge', 'DY4Merge']

ToTauTauList = ['ggHTauTau', 'VBFHTauTau', 'ZHTauTau', 'WminusHTauTau', 'WplusHTauTau']
ToWWList = ['ggHWW', 'VBFHWW']
SMHList = ToWWList+ToTauTauList

DiBosonsList = ['ZZ', 'WZ', 'WW2L2Nu']

NonpromtleptonsList = ['WWQQ', 'WJetsInclusiveMerge', 'QCD']#METER QCD AQUI y Wjets

ttbarList = ['TTbar', 'TTbarSemileptonic']
tWList = ['TW', 'TbarW']

MCList = SMHList + DiBosonsList + ['WWQQ'] + ttbarList + tWList + MCBckgMergingWJetsList + MCBckgMergingDYList

MCListQCD = MCList + ['QCD']

ListForBckgQCD = ['WJetsInclusiveMerge', 'DYInclusiveMerge', 'ZZ', 'TTbar', 'TW', 'ggHWW', 'WWQQ']

ALL = DataList + MCSignalList + MCListQCD

ShortList = ['ZZ', 'WWQQ', 'TTbar', 'TW', 'ggHWW', 'DYInclusiveMerge']

SMBR = {'ggHTauTau':0.06,
	'VBFHTauTau':0.06,
	'ZHTauTau':0.06,
	'WminusHTauTau':0.06,
	'WplusHTauTau':0.06,
	'ggHWW':0.2,
	'VBFHWW':0.2}

#nameofhist = '_hInvariantMassMuonElectron'
nameofhist = '_hCollinearMass'
#nameofhist = '_hnVertex'
#nameofhist = '_hMuonPt'
#nameofhist = '_hElectronPt'
#nameofhist = '_hCollinearMass_1jet'
#nameofhist = '_hCollinearMass_0jets'
#nameofhist = '_hCollinearMass_2jetsVBF'
#nameofhist = '_hCollinearMass_BoostingGGH'
#nameofhist = '_hnJetSel'
#nameofhist = '_hnbMediumSel'
#nameofhist = '_hJetPt1'

length_nameofhist = len(nameofhist)

###
def Getting_histograms(sign, muon, jets, bjets, lista):
	"""Routine to get the histogram to plot from the .root files"""

	for sample in lista:
		try:
                        Fin = ROOT.TFile.Open(directory+'rootfiles/'+sample+'_ToPlot_V6_'+sign+'_'+muon+'_'+jets+'_'+bjets+'.root')			
		except:
			print('Could not open file')
			raise
		try:
			print('Taking histogram from sample ',sample)
			ROOT.gDirectory.cd('PyROOT:/')
			Hin[sample] = Fin.Get(sample+nameofhist).Clone()

		except:
			print('Could not get the histogram '+sample)
			raise
		Fin.Close()

	return Hin

###
def LFV_plots (sign, muon, jets, bjets):
	"""Routine to plot the LFV analysis plots"""

	#Setting up the plot
	canvas = ROOT.TCanvas("plot"+nameofhist, "Plot of the "+nameofhist, 1000, 1000)	
	canvas.Draw()	
	pad1 = ROOT.TPad("pad1", "The pad 70% of the height", 0.0,0.3,1.0,1.0,31)
	pad1.Draw("same")
	pad1.SetFillColor(0)
	pad2 = ROOT.TPad("pad2", "The pad 30% of the height", 0.0,0.0,1.0,0.3,32)
	pad2.Draw("same")
	pad2.SetFillColor(0)

	#Here I have to define what I wanna plot

	colors = {"ggHMuTau":ROOT.kRed+1, 
		  "VBFHMuTau":ROOT.kBlue+1,
		  "dibosons":ROOT.kPink+1, 
		  "nonpromptleptons":ROOT.kGreen+1, 
		  "tW":ROOT.kOrange+1, 
		  "WJetsInclusive":ROOT.kAzure+1, 
		  "ttbar":ROOT.kMagenta+1, 
		  "DYMerge":ROOT.kCyan+1, 
		  "toTauTau":ROOT.kSpring+1, 
		  "smhiggs":ROOT.kViolet+1, 
		  "QCD":ROOT.kTeal+1, 
		  "toWW":ROOT.kRed+1}

	legend = {"ggHMuTau":"ggF h#rightarrow#mu#tau (BR = 20%)",
		  "VBFHMuTau":"VBF h#rightarrow#mu#tau (BR = 20%)",
		  "dibosons": "dibosons: ZZ+WW(2L2Nu)+WZ", 
		  "nonpromptleptons":"(W+NJets)+WWQQ+QCD", 
		  "WJetsInclusive":"W+NJets",
		  "tW":"tW+#bar{t}W", 
		  "ttbar":"t#bar{t}", 
		  "DYMerge":"Drell-Yan", 
		  "toTauTau":"hToTauTau: ggF,VBF,ZH,W-,W+",
		  "toWW":"hToWW: ggF, VBF",
		  "smhiggs":'SM Higgs'}
			
	pad1.cd()
	pad1.SetBottomMargin(0.02)
	pad1.SetTopMargin(0.1)

	
	#title = ROOT.TPaveText(.11, .95, .35, .99, "brndc")
	#title.AddText( "CMS preliminary - L = 59.266 fb^{-1} (13 TeV)")
	#title.Draw("same")

	#preliminary = ROOT.TLatex(-0.4304896,35834.14, "CMS preliminary")
	#luminosity = ROOT.TLatex(117.3359,35834.14, "L = 59.266 fb^{-1} (13 TeV)")

	#title = ROOT.TLatex(0.5,0.5, "CMS preliminary - L = 59.266 fb^{-1} (13 TeV)")

	Samples = {"ggHMuTau":Hin["ggHMuTau"],
		   "VBFHMuTau":Hin["VBFHMuTau"],
		   "ZZ": hDiBosons, 
		   "WWQQ":hNonPromptLeptons,
		   "TW":htW,  
		   "TTbar":httbar, 
		   "ggHWW":hSMHiggs,
		   "DYInclusiveMerge":hDYMerge}
	
	
	##### MC BACKGROUNDS #####

	hAux = []
	xsmp = 1
	for sample in ShortList:
		if(len(hAux)==0):
			print(sample)
			hAux.append(Samples.get(sample))
		else: 
			print(sample)
			hAux.append(Samples.get(sample) + hAux[xsmp-1])
			xsmp = xsmp + 1

	print(len(hAux))

	rebinvaluestring = "10"
	rebinvalue = int(rebinvaluestring)
	
	binmax = hAux[-1].GetMaximumBin()
	y = hAux[-1].GetBinContent(binmax)
	print('binmax, ', binmax)
	print('y, ', y)

	for xsample in reversed(hAux): #[::-1]:
		print('Painting : ' + xsample.GetName()[0:-length_nameofhist])
		xsample.SetFillColor(colors.get(xsample.GetName()[0:-length_nameofhist]))
		xsample.SetLineColor(ROOT.TColor.GetColorDark(colors.get(xsample.GetName()[0:-length_nameofhist])))
		xsample.GetXaxis().SetLabelSize(0)
		xsample.SetYTitle("Events / "+rebinvaluestring+" GeV")
		xsample.GetYaxis().SetTitleOffset(1.3)
		xsample.GetYaxis().SetTitleSize(0.04)
   		xsample.GetYaxis().SetLabelFont(42);
   		xsample.GetYaxis().SetLabelOffset(0.007);
   		xsample.GetYaxis().SetTitleFont(42);
		xsample.SetTitleOffset(1.4)
		#xsample.SetMaximum(50000)
		xsample.SetMaximum(y+y*15)
		xsample.RebinX(rebinvalue)
		xsample.Draw("hist same")
		ROOT.objs.append(xsample)			
	
	##### SIGNALS #### FOR THE BLINDING?
	
	hSumaSignals = None 
	for sample in MCSignalList:
		Hin[sample]=Hin[sample]*0.2
		Hin[sample].SetLineColor(colors.get(sample))
		Hin[sample].SetLineWidth(3)
		Hin[sample].RebinX(rebinvalue)
		Hin[sample].Draw("hist same")
		if (hSumaSignals is None): hSumaSignals = Hin[sample]
		else: hSumaSignals += Hin[sample]
		print('Painting sample ' + sample)

	###### DATA #######
	
	hSumData = None
	for sample in DataList:
		if (hSumData is None): hSumData = Hin[sample]
		else: hSumData = hSumData + Hin[sample]

	hSumData.RebinX(rebinvalue)
	hBlinding = hSumData.Clone()
	hBlinding.Divide(hSumaSignals)
	
	bine=hSumData.GetNbinsX()+1
	while (bine>=0):
		if(hSumaSignals.GetBinContent(bine)>hAux[-1].GetBinContent(bine)*0.1):
			hSumData.SetBinContent(bine, -10000)
		bine-=1

	hSumData.SetMinimum(0)
	hSumData.SetMarkerStyle(20)
	hSumData.SetMarkerSize(1.1)
	hSumData.Draw("e same F")

	##### LEGEND #####
	
	leg = ROOT.TLegend(0.42,0.32,0.88,0.88)
	leg.SetNColumns(2);
   	leg.SetFillColor(4001);
 	leg.SetBorderSize(0);
	leg.SetColumnSeparation(0.1)
	leg.AddEntry(hSumData, "Data", "P")
	leg.SetTextSize(0.025)
	leg.SetTextFont(62)
	leg.SetHeader("CMS preliminary - L = 59.266 fb^{-1} (13 TeV) "+jets+" "+bjets+" ("+sign+" and " + muon+ ")", "C") 

	for sample in MCSignalList:
		leg.AddEntry(Hin[sample], legend.get(sample), "F")

	for sample in reversed(hAux):
		leg.AddEntry(sample, legend.get(sample.GetName()[0:-length_nameofhist]), "F")
        leg.SetTextFont(42);
	leg.Draw("same")
	
	ROOT.objs.append(leg)




	##### ERRORS

	'''errorBand = None

	errorBand = AddingHistograms(ShortMCList, "errorband")
	errorBand.SetMarkerSize(0) 
	errorBand.SetFillColor(ROOT.kGray+2)
	errorBand.SetFillStyle(3001)
	errorBand.SetLineWidth(1)
	errorBand.Draw("e2same")'''

	##### RATIO DATA/MC #####
	pad2.cd()
	pad2.SetTopMargin(0.00)
	pad2.SetBottomMargin(0.3)
	hRatio = None
	hRatio = hSumData.Clone()
	hRatio.SetName("RatioDataMC")
	hRatio.Divide(hAux[-1])
	hRatio.Draw("")
	hRatio.SetYTitle("Ratio Data/MC")
	hRatio.GetYaxis().CenterTitle()
	#hRatio.SetXTitle("Electron pt [GeV]")
	#hRatio.SetXTitle("Muon pt [GeV]")
	#hRatio.SetXTitle("Number of vertexes [#]")
	#hRatio.SetXTitle("Collinear Mass Muon-Electron [GeV]")
	#hRatio.SetXTitle("Invariant Mass Muon-Electron [GeV]")
	hRatio.SetXTitle(nameofhist)
	hRatio.GetYaxis().SetTitleOffset(0.6)
	hRatio.SetMinimum(0.4)
	hRatio.SetMaximum(1.6)
	hRatio.SetMarkerStyle(20)
	hRatio.SetMarkerSize(1)
	hRatio.SetMarkerColor(ROOT.kRed)
	hRatio.SetLineColor(ROOT.kRed+1)
	hRatio.SetLabelSize(0.09 , "xy")
	hRatio.SetTitleSize(0.08, "xy")
	hLine = ROOT.TLine(hRatio.GetXaxis().GetXmin(), 1, hRatio.GetXaxis().GetXmax(), 1)
	hLine.SetLineColor(ROOT.kRed+1)
	hLine.SetLineStyle(6)
	hLine.Draw("P same")
	pad2.SetGridy(1)
	hRatio.GetYaxis().SetNdivisions(5)
	ROOT.objs.append(hRatio)
	ROOT.objs.append(hLine)
		
	ROOT.objs.append([canvas,hSumData,Hin,hRatio, hLine, leg])
	canvas.Update()
###
def AddingHistograms(mylist, nameforhistogram):
	"""Routine to calculate the sum of the merged samples"""

	hAux = []
	msamp = 1 
	for sample in mylist:

		if(SMBR.get(sample)!=None): Hin[sample] = Hin[sample]*SMBR.get(sample)
		else: Hin[sample] = Hin[sample]

		if(len(hAux)==0):
			hAux.append(Hin[sample])
		else:
			hAux.append(Hin[sample]+hAux[msamp-1])
			msamp = msamp + 1

	hAux[-1].SetName(nameforhistogram+nameofhist)

	print(hAux[-1])

	return hAux[-1]

###
def GettingQCD(sign, muon, jets, bjets):
	"""Routine to get QCD Histograms using OSSS method"""

	'''Filling Hin dictionary with SS IM features to calculate QCD background'''
	Getting_histograms("SS", "IM", jets, bjets, DataList+MCList)

	hSumData = AddingHistograms(DataList, "sumdata")
	
	Hin["WJetsInclusiveMerge"] = AddingHistograms(MCBckgMergingWJetsList, "WJetsMerge")
	Hin["DYInclusiveMerge"] = AddingHistograms(MCBckgMergingDYList, "DYMerge")
	Hin["ZZ"] = AddingHistograms(DiBosonsList, "dibosons")
	Hin["TTbar"] = AddingHistograms(ttbarList, "ttbar")
	Hin["TW"] = AddingHistograms(tWList, "tW")
	Hin["ggHWW"] = AddingHistograms(SMHList, "smhiggs")
	
	hSumBckg = AddingHistograms(ListForBckgQCD, "sumbckg")

	if(hSumData.GetXaxis().GetNbins()!=hSumBckg.GetXaxis().GetNbins()): 
		print('ERROR: DIFFERENT NUMBER OF BINS FOR DATA AND MC')
	else:
		hQCD = ROOT.TH1D("QCD", "QCD",hSumData.GetXaxis().GetNbins(),hSumData.GetXaxis().GetXmin(), hSumData.GetXaxis().GetXmax())

	for entryDAT in range(0,hSumData.GetXaxis().GetNbins()):
		if(hSumData.GetBinContent(entryDAT)<hSumBckg.GetBinContent(entryDAT)): 
                     	hQCD.SetBinContent(entryDAT, 0.0)
		else: 
			hQCD.SetBinContent(entryDAT,hSumData.GetBinContent(entryDAT)-hSumBckg.GetBinContent(entryDAT) )

	hQCD.SetName("QCD"+nameofhist)

	print(hQCD.GetName())
	f = TFile(directory+"rootfiles/QCD_ToPlot_V6_"+sign+"_"+muon+"_"+jets+"_"+bjets+".root","RECREATE");
	hQCD.Write("QCD"+nameofhist)
	
###
def setStyle ():
	"""Routine to setup the configuration of the style for the plots"""

	#Using the TDR Style from CMS

	import tdrstyle

	tdrstyle.setTDRStyle()

	#And some additional "fine-tuning" of the default plotting options

	ROOT.gStyle.SetHistMinimumZero(ROOT.kTRUE)

###
#Code to be run in order to produce the plots
if __name__ == '__main__':

   	import sys
    	setStyle()

	jets = "2jets"
	bjets = "0bjets"
	sign = "OS"
	muon = "IM"

	'''First of all I have to calculate the missing background (QCD) using OSSS Method'''
	GettingQCD(sign, muon, jets, bjets)

	'''I clear the dictionary because it has the histograms to calculate qcd using the OSS method'''
	Hin.clear()

	'''Once I have all the Data, Signal and MC Bacground I will start to get the histograms that I want to plot'''
	Getting_histograms(sign, muon,jets, bjets, ALL)

	'''Once I have the Hin dictionary filled, I will put together some of the background'''
	hDYMerge = AddingHistograms(MCBckgMergingDYList, "DYMerge")
	hSMHiggs = AddingHistograms(SMHList, "smhiggs")
	hDiBosons = AddingHistograms(DiBosonsList, "dibosons")
	'''I need to put together WJets first because it belongs to the nonpromptleptons list'''
	Hin['WJetsInclusiveMerge'] = AddingHistograms(MCBckgMergingWJetsList, "WJetsMerge")
	hNonPromptLeptons = AddingHistograms(NonpromtleptonsList, "nonpromptleptons")
	httbar = AddingHistograms(ttbarList, "ttbar")
	htW = AddingHistograms(tWList, "tW")

    	LFV_plots(sign, muon, jets, bjets)
	
	raw_input("Press ENTER... ")

