#
# Skimming macro for the NANOAOD within the pyCiemat framework to filter events
# requiring two muons or a muon and an electron intended for filtering events
# for the LFV Higgs search.
#
# Written by Maria Cepeda (2019_11_20) [+changes by Oscar Gonzalez]
#

from CMSSWCiemat.PyCiemat4NanoAOD.GeneralTools.SkimmingProcess_cff import SkimmingProcess
#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *


def skimmingMacro(filenames,goodlumis,maxevents=None, isData=False):
    """Routine to run the skimmer to select events containing leptons
    (electrons or muons)."""

    skimmer=SkimmingProcess('mariaMuELFVLeptonFilterSkimmer')

    # Configuration parameters are given as arguments, but they may be hardcoded for tests

    skimmer.inputFiles=filenames

    skimmer.goodLumis = goodlumis

    skimmer.isData=isData

    #skimmer.maxEvents=107 # Poner 107 para pruebas cortas en crab (just uncomment this line)
    if (maxevents is not None):
        skimmer.maxEvents=maxevents

    # -----------------------------------------------------------------------

    # Setting up the filters! (1) Basic lepton filter

    import CMSSWCiemat.PyCiemat4NanoAOD.GeneralTools.CiematNanoAODFilter_cff
    mufilter=CMSSWCiemat.PyCiemat4NanoAOD.GeneralTools.CiematNanoAODFilter_cff.CiematNanoAODFilter('muonfilter')
    dileptonfilter=CMSSWCiemat.PyCiemat4NanoAOD.GeneralTools.CiematNanoAODFilter_cff.CiematNanoAODFilter('dileptonfilter')

    #
    import CMSSWCiemat.PyCiemat4NanoAOD.Filters.CiematMuonFilter
    muonfilter = CMSSWCiemat.PyCiemat4NanoAOD.Filters.CiematMuonFilter.CiematMuonFilter('muonfilter')
    muonfilter.muonClassRequired=["looseId",'pTPFCut']
    muonfilter.minPtValue = 20
    muonfilter.nMinMuons = 1  # For checks!

    dimuonfilter = CMSSWCiemat.PyCiemat4NanoAOD.Filters.CiematMuonFilter.CiematMuonFilter('dimuonfilter')
    dimuonfilter.muonClassRequired=["looseId",'pTPFCut']
    dimuonfilter.minPtValue = 20
    dimuonfilter.nMinMuons = 2  # For checks!

    mufilter.subfilters.append(muonfilter)
    #dileptonfilter.subfilters.append(dimuonfilter)
    #dileptonfilter.subfilters.append(muonfilter)

    import CMSSWCiemat.PyCiemat4NanoAOD.Filters.CiematElectronFilter
    electronfilter=CMSSWCiemat.PyCiemat4NanoAOD.Filters.CiematElectronFilter.CiematElectronFilter('electronfilter')
    electronfilter.nMinElectrons = 1  # For checks!
    dileptonfilter.subfilters.append(electronfilter)

    skimmer.filterPath.append(mufilter)
    skimmer.filterPath.append(dileptonfilter)

    METCorrector = createJMECorrector(isMC=True, dataYear='2018', jesUncert="Total", jetType = "AK4PFchs", isFastSim=False)
    modules = [METCorrector()]

    p = PostProcessor('./', filenames, modules=modules, haddFileName="")
    p.run()
    skimmer.filterPath.append(METCorrector)
    # Setting up the skimming process

    skimmer.setup()

    # Running!
    skimmer.run()

    # We return the skimmer in case it is needed outside!
    return skimmer

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# For running as a shell command
if __name__ == '__main__':
    import sys

    nevents=-1
    if (len(sys.argv)>2):
        if (sys.argv[2]=='-'): 
		nevents=-1
        else: 
		print("HOLA",sys.argv[2])
		nevents=int(sys.argv[2])

    files=sys.argv[1].split(',')

    sys.argv=[sys.argv[0]]   # To avoid pyROOT to make strange things with the arguments.

    skimmingMacro(files,None,nevents)

# =======================================================================
