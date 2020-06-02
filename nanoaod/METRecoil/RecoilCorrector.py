import os
import re
import ROOT
import math 

class RecoilCorrector:

   def __init__(self, fileName):
      cmsswBase ="/afs/hep.wisc.edu/home/psiddire/CMSSW_10_2_16_UL"
      baseDir = cmsswBase + "/src/FinalStateAnalysis/TagAndProbe/data"
      _fileName = baseDir+"/"+fileName
      File = ROOT.TFile.Open(_fileName)
      if (File.IsZombie()):
        print "file not found"
      projH = File.Get("projH")
      if (projH == None):
        print "File should contain histogram with the name projH " 
      firstBinStr  = projH.GetXaxis().GetBinLabel(1)
      secondBinStr = projH.GetXaxis().GetBinLabel(2)
      paralZStr = firstBinStr
      perpZStr  = secondBinStr
      if ("Perp" in firstBinStr):
        paralZStr = secondBinStr
        perpZStr  = firstBinStr
      ZPtBinsH = File.Get("ZPtBinsH")
      if (ZPtBinsH == None):
        print "File should contain histogram with the name ZPtBinsH "
      nZPtBins = ZPtBinsH.GetNbinsX()
      ZPtBins=[]
      ZPtStr=[]
      for i in range(0, nZPtBins+1):
        ZPtBins.append(ZPtBinsH.GetXaxis().GetBinLowEdge(i+1))
        if (i < nZPtBins):
          ZPtStr.append(ZPtBinsH.GetXaxis().GetBinLabel(i+1))
      nJetBinsH = File.Get("nJetBinsH")
      if (nJetBinsH == None):
        print "File should contain histogram with the name nJetBinsH"
      nJetsBins = nJetBinsH.GetNbinsX()
      nJetsStr = []
      for i in range(0, nJetsBins):
        nJetsStr.append(nJetBinsH.GetXaxis().GetBinLabel(i+1))
      self.InitMEtWeights(File, perpZStr, paralZStr, nZPtBins, ZPtBins, ZPtStr, nJetsBins, nJetsStr)
      _epsrel = 5e-4
      _epsabs = 5e-4
      _range = 0.95


   def InitMEtWeights(self, _file, _perpZStr, _paralZStr, nZPtBins, ZPtBins, _ZPtStr, nJetsBins, _nJetsStr):
	newPerpZStr  = str(_perpZStr)
	newParalZStr = str(_paralZStr)
        newZPtBins = []
        newZPtStr = []
        newNJetsStr = []
	for idx in range(0, nZPtBins+1):
	  newZPtBins.append(ZPtBins[idx])
	for idx in range(0, nZPtBins):
	  newZPtStr.append(str(_ZPtStr[idx]))
	for idx in range(0, nJetsBins):
	  newNJetsStr.append(str(_nJetsStr[idx]))
	self.InitMEtWeightsC(_file, newZPtBins, newPerpZStr, newParalZStr, newZPtStr, newNJetsStr)


   def InitMEtWeightsC(self, _fileMet, ZPtBins, _perpZStr, _paralZStr, _ZPtStr, _nJetsStr):
      _nZPtBins = len(ZPtBins)-1
      self._nJetsBins = len(_nJetsStr)
      self._ZPtBins = ZPtBins
      _ZPtBins = ZPtBins
      _metZParalData=[[0 for x in range(3)] for y in range(5)]
      _xminMetZParalData=[[0 for x in range(3)] for y in range(5)]
      _xmaxMetZParalData=[[0 for x in range(3)] for y in range(5)]
      _metZPerpData=[[0 for x in range(3)] for y in range(5)]
      _xminMetZPerpData=[[0 for x in range(3)] for y in range(5)]
      _xmaxMetZPerpData=[[0 for x in range(3)] for y in range(5)]
      _metZParalMC=[[0 for x in range(3)] for y in range(5)]
      _xminMetZParalMC=[[0 for x in range(3)] for y in range(5)]
      _xmaxMetZParalMC=[[0 for x in range(3)] for y in range(5)]
      _metZPerpMC=[[0 for x in range(3)] for y in range(5)] 
      _xminMetZPerpMC=[[0 for x in range(3)] for y in range(5)] 
      _xmaxMetZPerpMC=[[0 for x in range(3)] for y in range(5)] 
      _xminMetZParal=[[0 for x in range(3)] for y in range(5)] 
      _xmaxMetZParal=[[0 for x in range(3)] for y in range(5)] 
      _xminMetZPerp=[[0 for x in range(3)] for y in range(5)] 
      _xmaxMetZPerp=[[0 for x in range(3)] for y in range(5)] 
      self._meanMetZParalData=[[0 for x in range(3)] for y in range(5)] 
      self._rmsMetZParalData=[[0 for x in range(3)] for y in range(5)] 
      self._rmsMetZPerpData=[[0 for x in range(3)] for y in range(5)] 
      self._meanMetZParalMC=[[0 for x in range(3)] for y in range(5)] 
      self._rmsMetZParalMC=[[0 for x in range(3)] for y in range(5)] 
      self._rmsMetZPerpMC=[[0 for x in range(3)] for y in range(5)] 
      for ZPtBin in range(0,_nZPtBins):
        for jetBin in range(0,self._nJetsBins):
          binStrPerpData  = _perpZStr  + "_" + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_data"
          binStrParalData = _paralZStr + "_" + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_data"
          binStrPerpMC    = _perpZStr  + "_" + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_mc"
          binStrParalMC   = _paralZStr + "_" + _nJetsStr[jetBin] + _ZPtStr[ZPtBin] + "_mc"
          _metZParalData[ZPtBin][jetBin] = _fileMet.Get(binStrParalData)
          _metZPerpData[ZPtBin][jetBin]  = _fileMet.Get(binStrPerpData)
          _metZParalMC[ZPtBin][jetBin]   = _fileMet.Get(binStrParalMC)
          _metZPerpMC[ZPtBin][jetBin]    = _fileMet.Get(binStrPerpMC)
          if (_metZParalData[ZPtBin][jetBin]==None):
    	      print  "Function with names not found in file "
          if (_metZPerpData[ZPtBin][jetBin]==None):
    	      print  "Function with names not found in file "
          if (_metZParalMC[ZPtBin][jetBin]==None):
    	      print  "Function with names not found in file "
          if (_metZPerpMC[ZPtBin][jetBin]==None):
    	      print  "Function with names not found in file "
          xminD=ROOT.Double(0.0000000000000)
          xmaxD=ROOT.Double(0.0000000000000)
          _metZParalData[ZPtBin][jetBin].GetRange(xminD,xmaxD)
          _xminMetZParalData[ZPtBin][jetBin] = float(xminD)
          _xmaxMetZParalData[ZPtBin][jetBin] = float(xmaxD)
          _metZPerpData[ZPtBin][jetBin].GetRange(xminD,xmaxD);
          _xminMetZPerpData[ZPtBin][jetBin] = float(xminD);
          _xmaxMetZPerpData[ZPtBin][jetBin] = float(xmaxD);
          _metZParalMC[ZPtBin][jetBin].GetRange(xminD,xmaxD);
          _xminMetZParalMC[ZPtBin][jetBin] = float(xminD);
          _xmaxMetZParalMC[ZPtBin][jetBin] = float(xmaxD);
          _metZPerpMC[ZPtBin][jetBin].GetRange(xminD,xmaxD);
          _xminMetZPerpMC[ZPtBin][jetBin] = float(xminD);
          _xmaxMetZPerpMC[ZPtBin][jetBin] = float(xmaxD);
          _xminMetZParal[ZPtBin][jetBin] = max(_xminMetZParalData[ZPtBin][jetBin],_xminMetZParalMC[ZPtBin][jetBin]);
          _xmaxMetZParal[ZPtBin][jetBin] = min(_xmaxMetZParalData[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]);
          _xminMetZPerp[ZPtBin][jetBin] = max(_xminMetZPerpData[ZPtBin][jetBin],_xminMetZPerpMC[ZPtBin][jetBin]);
          _xmaxMetZPerp[ZPtBin][jetBin] = min(_xmaxMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]);
          self._meanMetZParalData[ZPtBin][jetBin] = _metZParalData[ZPtBin][jetBin].Mean(_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin])
          self._rmsMetZParalData[ZPtBin][jetBin] = (_metZParalData[ZPtBin][jetBin].CentralMoment(2,_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin]))**(0.5)
          self._rmsMetZPerpData[ZPtBin][jetBin] =(_metZPerpData[ZPtBin][jetBin].CentralMoment(2,_xminMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpData[ZPtBin][jetBin]))**(0.5)
          self._meanMetZParalMC[ZPtBin][jetBin] = _metZParalMC[ZPtBin][jetBin].Mean(_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin])
          self._rmsMetZParalMC[ZPtBin][jetBin] = (_metZParalMC[ZPtBin][jetBin].CentralMoment(2,_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]))**(0.5)
          self._rmsMetZPerpMC[ZPtBin][jetBin] = (_metZPerpMC[ZPtBin][jetBin].CentralMoment(2,_xminMetZPerpMC[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]))**(0.5)


   def CorrectByMeanResolution(self, MetPx, MetPy, genVPx, genVPy, visVPx, visVPy, njets):
       Zpt = math.sqrt(genVPx*genVPx + genVPy*genVPy)
       self.U1 = 0
       self.U2 = 0
       self.metU1 = 0
       self.metU2 = 0
       self.CalculateU1U2FromMet(MetPx, MetPy, genVPx, genVPy, visVPx, visVPy)
       if (Zpt > 1000):
           Zpt = 999
       if (njets >= self._nJetsBins):
           njets = self._nJetsBins - 1
       MetCorrPx = 0.00000000
       MetCorrPy = 0.00000000
       for iB in range(0, len(self._ZPtBins)):
          if (Zpt >= self._ZPtBins[iB] and Zpt < self._ZPtBins[iB+1]):
            ZptBin = iB 
          else:
            ZptBin = 0
       self.U1U2CorrectionsByWidth(ZptBin, njets) 
       tmpMet = self.CalculateMetFromU1U2(genVPx, genVPy, visVPx, visVPy, MetCorrPx, MetCorrPy)
       return tmpMet


   def U1U2CorrectionsByWidth(self, ZptBin, njets):
       if (njets >= self._nJetsBins):
          njets = self._nJetsBins - 1
       width = self.U1 - self._meanMetZParalMC[ZptBin][njets]
       width *= self._rmsMetZParalData[ZptBin][njets]/self._rmsMetZParalMC[ZptBin][njets]
       self.U1 = self._meanMetZParalData[ZptBin][njets] + width
       width = self.U2
       width *= self._rmsMetZPerpData[ZptBin][njets]/self._rmsMetZPerpMC[ZptBin][njets]
       self.U2 = width


   def CalculateU1U2FromMet(self, metPx, metPy, genZPx, genZPy, diLepPx, diLepPy):  
       hadRecX = metPx + diLepPx - genZPx
       hadRecY = metPy + diLepPy - genZPy
       hadRecPt = (hadRecX*hadRecX + hadRecY*hadRecY)**(0.5)
       phiHadRec = math.atan2(hadRecY, hadRecX)
       phiDiLep = math.atan2(diLepPy, diLepPx)
       phiMEt =math.atan2(metPy, metPx)
       metPt = (metPx*metPx + metPy*metPy)**(0.5)
       phiZ = math.atan2(genZPy, genZPx)
       deltaPhiZHadRec  = phiHadRec - phiZ
       deltaPhiDiLepMEt = phiMEt - phiDiLep
       self.U1 = hadRecPt * math.cos(deltaPhiZHadRec)
       self.U2 = hadRecPt * math.sin(deltaPhiZHadRec)
       self.metU1 = metPt * math.cos(deltaPhiDiLepMEt)      
       self.metU2 = metPt * math.sin(deltaPhiDiLepMEt)


   def CalculateMetFromU1U2(self, genZPx, genZPy, diLepPx, diLepPy, metPx, metPy):  
       hadRecPt = math.sqrt(self.U1*self.U1 + self.U2*self.U2)
       deltaPhiZHadRec = math.atan2(self.U2, self.U1)
       phiZ = math.atan2(genZPy, genZPx)
       phiHadRec = phiZ + deltaPhiZHadRec
       hadRecX = hadRecPt*math.cos(phiHadRec)
       hadRecY = hadRecPt*math.sin(phiHadRec)
       metPx = hadRecX + genZPx - diLepPx
       metPy = hadRecY + genZPy - diLepPy
       return [metPx, metPy]


def Metcorrected(fileName):
    return RecoilCorrector(fileName)
