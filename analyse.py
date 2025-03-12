import ROOT
import yaml


ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro('inc/RooGausDExp.cxx+')

with open("configurations.yaml", "r") as stream:
  conf = yaml.safe_load(stream)

dataDfMB = ROOT.RDataFrame("tree", "snapshots/data_mb.root")
dataDf = ROOT.RDataFrame("tree", "snapshots/data.root")
primaryDf = ROOT.RDataFrame("tree", "snapshots/primary.root")
charmDf = ROOT.RDataFrame("tree", "snapshots/charm.root")
beautyDf = ROOT.RDataFrame("tree", "snapshots/beauty.root")

outFile = ROOT.TFile("output.root", "recreate")
purityDir = outFile.mkdir("purity")
mcDir = purityDir.mkdir("mc")
dataDir = purityDir.mkdir("data")

massMB = dataDfMB.Histo2D(("hMassMB", ";#it{p}_{T} (GeV/#it{c});#it{m} (GeV/#it{c}^{2})", conf['nPtBins'], conf['minPt'], conf['maxPt'], 70, 1.64, 1.71), "fCascPt", "fMassOmega")
massMC = primaryDf.Histo2D(("hMassMC", ";#it{p}_{T} (GeV/#it{c});#it{m} (GeV/#it{c}^{2})", conf['nPtBins'], conf['minPt'], conf['maxPt'], 70, 1.64, 1.71), "fCascPt", "fMassOmega")

wsPurity = ROOT.RooWorkspace("wsPurity")
wsPurity.factory("RooCrystalBall::signal(mass[1.64, 1.71], mu[1.68, 1.64, 1.71], sigma[0.01, 0.001, 0.1], alphaL[1.5, 0.1, 10], nL[1.5, 0.1, 10], alphaR[1.5, 0.1, 10], nR[1.5, 0.1, 10])")
wsPurity.factory("RooExponential::background(mass, lambda[-1, -40, 0])")
wsPurity.factory("SUM::model(fraction[0.5, 0, 1]*signal, background)")

hPurity = ROOT.TH1D("hPurity", ";#it{p}_{T} (GeV/#it{c});Fraction of signal", conf['nPtBins'], conf['minPt'], conf['maxPt'])
hMean = ROOT.TH1D("hMean", ";#it{p}_{T} (GeV/#it{c});#mu of the signal (GeV/#it{c}^{2})", conf['nPtBins'], conf['minPt'], conf['maxPt'])
hSigma = ROOT.TH1D("hSigma", ";#it{p}_{T} (GeV/#it{c});#sigma of the signal (GeV/#it{c}^{2})", conf['nPtBins'], conf['minPt'], conf['maxPt'])
for i in range(conf['nPtBins']):
  datasetMC = ROOT.RooDataHist(f"datasetMC_{i}", f"datasetMC_{i}", ROOT.RooArgList(wsPurity.var("mass")), massMC.ProjectionY(f"hMassMC_{i}", i+1, i+1))
  wsPurity.pdf("signal").fitTo(datasetMC)
  plotMC = wsPurity.var("mass").frame()
  datasetMC.plotOn(plotMC)
  wsPurity.pdf("signal").plotOn(plotMC)
  mcDir.cd()
  plotMC.Write(f"frameMC_{i}")

  wsPurity.var("alphaR").setConstant(True)
  wsPurity.var("alphaL").setConstant(True)
  wsPurity.var("nR").setConstant(True)
  wsPurity.var("nL").setConstant(True)

  datasetMB = ROOT.RooDataHist(f"datasetMB_{i}", f"datasetMB_{i}", ROOT.RooArgList(wsPurity.var("mass")), massMB.ProjectionY(f"hMassMB_{i}", i+1, i+1))
  wsPurity.pdf("model").fitTo(datasetMB)
  plotMB = wsPurity.var("mass").frame()
  datasetMB.plotOn(plotMB)
  wsPurity.pdf("model").plotOn(plotMB)
  wsPurity.pdf("model").paramOn(plotMB)

  wsPurity.var("alphaR").setConstant(False)
  wsPurity.var("alphaL").setConstant(False)
  wsPurity.var("nR").setConstant(False)
  wsPurity.var("nL").setConstant(False)

  mean = wsPurity.var("mu").getVal()
  sigma = wsPurity.var("sigma").getVal()
  wsPurity.var("mass").setRange(f"sigRange{i}", mean - 2 * sigma, mean + 2 * sigma)
  integralSignal = wsPurity.pdf("signal").createIntegral(ROOT.RooArgSet(wsPurity.var("mass")), ROOT.RooFit.NormSet(ROOT.RooArgSet(wsPurity.var("mass"))), ROOT.RooFit.Range(f"sigRange{i}"))
  integralModel = wsPurity.pdf("model").createIntegral(ROOT.RooArgSet(wsPurity.var("mass")), ROOT.RooFit.NormSet(ROOT.RooArgSet(wsPurity.var("mass"))), ROOT.RooFit.Range(f"sigRange{i}"))
  fractionSignal = integralSignal.getVal() * wsPurity.var("fraction").getVal() / integralModel.getVal()
  hMean.SetBinContent(i+1, mean)
  hMean.SetBinError(i+1, wsPurity.var("mu").getError())
  hSigma.SetBinContent(i+1, sigma)
  hSigma.SetBinError(i+1, wsPurity.var("sigma").getError())
  hPurity.SetBinContent(i+1, fractionSignal)

  dataDir.cd()
  plotMB.Write(f"frameMB_{i}")

purityDir.cd()
hMean.Fit("pol0")
averageMass = hMean.GetFunction("pol0").GetParameter(0)
hSigma.Fit("pol4")
sigPars = [hSigma.GetFunction("pol4").GetParameter(i) for i in range(5)]
sigmaFilter = f"({sigPars[0]} + {sigPars[1]} * fCascPt + {sigPars[2]} * fCascPt * fCascPt + {sigPars[3]} * fCascPt * fCascPt * fCascPt + {sigPars[4]} * fCascPt * fCascPt * fCascPt * fCascPt)"
hMean.Write()
hSigma.Write()
hPurity.Write()

signalDataDf = dataDf.Filter(f"fMassOmega > {averageMass} - 2 * {sigmaFilter} && fMassOmega < {averageMass} + 2 * {sigmaFilter}")
backgroundDataDf = dataDf.Filter(f"fMassOmega < {averageMass} - 5 * {sigmaFilter} || fMassOmega > {averageMass} + 5 * {sigmaFilter}")
hMassSignal = signalDataDf.Histo2D(("hMassSignal", ";#it{p}_{T} (GeV/#it{c});#it{m} (GeV/#it{c}^{2})", conf['nPtBins'], conf['minPt'], conf['maxPt'], 70, 1.64, 1.71), "fCascPt", "fMassOmega")
hMassBackground = backgroundDataDf.Histo2D(("hMassBackground", ";#it{p}_{T} (GeV/#it{c});#it{m} (GeV/#it{c}^{2})", conf['nPtBins'], conf['minPt'], conf['maxPt'], 70, 1.64, 1.71), "fCascPt", "fMassOmega")
purityDir.cd()
hMassSignal.Write()
hMassBackground.Write()

hDCAxyData = signalDataDf.Histo2D(("hDCAxyData", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")
hDCAxyBackground = backgroundDataDf.Histo2D(("hDCAxyBackground", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")
hDCAxyPrimary = primaryDf.Histo2D(("hDCAxyPrimary", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")
hDCAxyCharm = charmDf.Histo2D(("hDCAxyCharm", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")
hDCAxyBeauty = beautyDf.Histo2D(("hDCAxyBeauty", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")

dcaxy = ROOT.RooRealVar("dcaxy", "DCA_{xy}", 0., -conf['maxAbsDCAxy'], conf['maxAbsDCAxy'], "cm")
muDCAxy = ROOT.RooRealVar("muDCAxy", "#mu_{DCA_{xy}}", 0., -conf['maxAbsDCAxy'], conf['maxAbsDCAxy'], "cm")
sigmaDCAxy = ROOT.RooRealVar("sigmaDCAxy", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-4, 1.e-2, "cm")
resoDCAxy = ROOT.RooRealVar("resoDCAxy", "#sigma_{resolution xy}", 1.e-3, 1.e-4, 1.e-2, "cm")
sigmaDCAxyFit = ROOT.RooFormulaVar("sigmaDCAxyFit", "#sigma_{DCA_{xy}}", "TMath::Sqrt(@0 * @0 + @1 * @1)", ROOT.RooArgList(sigmaDCAxy, resoDCAxy))
alphaR = ROOT.RooRealVar("alphaR", "#alpha_{R}", 1.5, 0.1, 10.)
alphaL = ROOT.RooFormulaVar("alphaL", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaR))
gausExpMCxy = ROOT.RooGausDExp("gausExpMCxy", "GausDExp", dcaxy, muDCAxy, sigmaDCAxy, alphaL, alphaR)

sigmaPrimary = ROOT.TH1D("sigmaPrimary", ";#it{p}_{T} (GeV/#it{c});#sigma_{DCA_{xy}} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'])
sigmaCharm = ROOT.TH1D("sigmaCharm", ";#it{p}_{T} (GeV/#it{c});#sigma_{DCA_{xy}} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'])
sigmaBeauty = ROOT.TH1D("sigmaBeauty", ";#it{p}_{T} (GeV/#it{c});#sigma_{DCA_{xy}} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'])
primaryTplDir = outFile.mkdir("primaryTpl")
for i in range(conf['nPtBins']):
  datasetDCAxyData = ROOT.RooDataHist(f"datasetDCAxyData_{i}", f"datasetDCAxyData_{i}", ROOT.RooArgList(dcaxy), hDCAxyData.ProjectionY(f"hDCAxyData_{i}", i+1, i+1))
  datasetDCAxyBackground = ROOT.RooDataHist(f"datasetDCAxyBackground_{i}", f"datasetDCAxyBackground_{i}", ROOT.RooArgList(dcaxy), hDCAxyBackground.ProjectionY(f"hDCAxyBackground_{i}", i+1, i+1))
  datasetDCAxyPrimary = ROOT.RooDataHist(f"datasetDCAxyPrimary_{i}", f"datasetDCAxyPrimary_{i}", ROOT.RooArgList(dcaxy), hDCAxyPrimary.ProjectionY(f"hDCAxyPrimary_{i}", i+1, i+1))
  datasetDCAxyCharm = ROOT.RooDataHist(f"datasetDCAxyCharm_{i}", f"datasetDCAxyCharm_{i}", ROOT.RooArgList(dcaxy), hDCAxyCharm.ProjectionY(f"hDCAxyCharm_{i}", i+1, i+1))
  datasetDCAxyBeauty = ROOT.RooDataHist(f"datasetDCAxyBeauty_{i}", f"datasetDCAxyBeauty_{i}", ROOT.RooArgList(dcaxy), hDCAxyBeauty.ProjectionY(f"hDCAxyBeauty_{i}", i+1, i+1))

  gausExpMCxy.fitTo(datasetDCAxyPrimary)
  plotDCAxyPrimary = dcaxy.frame()
  datasetDCAxyPrimary.plotOn(plotDCAxyPrimary)
  gausExpMCxy.plotOn(plotDCAxyPrimary)
  gausExpMCxy.paramOn(plotDCAxyPrimary)
  sigmaPrimary.SetBinContent(i+1, sigmaDCAxy.getVal())
  sigmaPrimary.SetBinError(i+1, sigmaDCAxy.getError())
  primaryTplDir.cd()
  plotDCAxyPrimary.Write(f"frameDCAxyPrimary_{i}")
