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
sigmaDCAxy = ROOT.RooRealVar("sigmaDCAxy", "#sigma_{DCA_{xy}}", 1.e-3, 7.e-4, 1.e-2, "cm")
alphaR = ROOT.RooRealVar("alphaR", "#alpha_{R}", 1.5, 0.01, 10.)
alphaL = ROOT.RooFormulaVar("alphaL", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaR))
gausExpMCxy = ROOT.RooGausDExp("gausExpMCxy", "GausDExp", dcaxy, muDCAxy, sigmaDCAxy, alphaL, alphaR)

sigmaPrimary = ROOT.TH1D("sigmaPrimary", ";#it{p}_{T} (GeV/#it{c});#sigma_{DCA_{xy}} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'])
sigmaCharm = ROOT.TH1D("sigmaCharm", ";#it{p}_{T} (GeV/#it{c});#sigma_{DCA_{xy}} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'])
sigmaBeauty = ROOT.TH1D("sigmaBeauty", ";#it{p}_{T} (GeV/#it{c});#sigma_{DCA_{xy}} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'])
sigmaBkg = ROOT.TH1D("sigmaBkg", ";#it{p}_{T} (GeV/#it{c});#sigma_{DCA_{xy}} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'])

alphaRPrimary = ROOT.TH1D("alphaRPrimary", ";#it{p}_{T} (GeV/#it{c});#alpha_{R}", conf['nPtBins'], conf['minPt'], conf['maxPt'])
alphaRCharm = ROOT.TH1D("alphaRCharm", ";#it{p}_{T} (GeV/#it{c});#alpha_{R}", conf['nPtBins'], conf['minPt'], conf['maxPt'])
alphaRBeauty = ROOT.TH1D("alphaRBeauty", ";#it{p}_{T} (GeV/#it{c});#alpha_{R}", conf['nPtBins'], conf['minPt'], conf['maxPt'])
alphaRBkg = ROOT.TH1D("alphaRBkg", ";#it{p}_{T} (GeV/#it{c});#alpha_{R}", conf['nPtBins'], conf['minPt'], conf['maxPt'])

primaryTplDir = outFile.mkdir("primaryTpl")
primaryTplDir.cd()
hDCAxyPrimary.Write()
charmTplDir = outFile.mkdir("charmTpl")
charmTplDir.cd()
hDCAxyCharm.Write()
beautyTplDir = outFile.mkdir("beautyTpl")
beautyTplDir.cd()
hDCAxyBeauty.Write()
backgroundTplDir = outFile.mkdir("backgroundTpl")
backgroundTplDir.cd()
hDCAxyBackground.Write()

def get_and_store_templates(dataset, gausExpMCxy, dcaxy, sigmaHist, alphaRHist, tplDir, i, name):
    gausExpMCxy.fitTo(dataset)
    plot = dcaxy.frame()
    plot.SetName(f"frameDCAxy{name}_{i}")
    dataset.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
    gausExpMCxy.plotOn(plot, ROOT.RooFit.Name("model"))
    gausExpMCxy.paramOn(plot, ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plot.chiSquare("model", "data")}'))
    sigmaHist.SetBinContent(i+1, sigmaDCAxy.getVal())
    sigmaHist.SetBinError(i+1, sigmaDCAxy.getError())
    alphaRHist.SetBinContent(i+1, alphaR.getVal())
    alphaRHist.SetBinError(i+1, alphaR.getError())
    tplDir.cd()
    plot.Write()

hKSprob = ROOT.TH1D("hKSprob", ";#it{p}_{T} (GeV/#it{c});KS probability", conf['nPtBins'], conf['minPt'], conf['maxPt'])
for i in range(conf['nPtBins']):
  projPrimary = hDCAxyPrimary.ProjectionY(f"hDCAxyPrimary_{i}", i+1, i+1)
  projCharm = hDCAxyCharm.ProjectionY(f"hDCAxyCharm_{i}", i+1, i+1)
  hKSprob.SetBinContent(i+1, projPrimary.KolmogorovTest(projCharm))

  datasetDCAxyBackground = ROOT.RooDataHist(f"datasetDCAxyBackground_{i}", f"datasetDCAxyBackground_{i}", ROOT.RooArgList(dcaxy), hDCAxyBackground.ProjectionY(f"hDCAxyBackground_{i}", i+1, i+1))
  datasetDCAxyPrimary = ROOT.RooDataHist(f"datasetDCAxyPrimary_{i}", f"datasetDCAxyPrimary_{i}", ROOT.RooArgList(dcaxy), projPrimary)
  datasetDCAxyCharm = ROOT.RooDataHist(f"datasetDCAxyCharm_{i}", f"datasetDCAxyCharm_{i}", ROOT.RooArgList(dcaxy), projCharm)
  datasetDCAxyBeauty = ROOT.RooDataHist(f"datasetDCAxyBeauty_{i}", f"datasetDCAxyBeauty_{i}", ROOT.RooArgList(dcaxy), hDCAxyBeauty.ProjectionY(f"hDCAxyBeauty_{i}", i+1, i+1))

  alphaR.setRange(0.8, 3)
  alphaR.setVal(1.)
  get_and_store_templates(datasetDCAxyPrimary, gausExpMCxy, dcaxy, sigmaPrimary, alphaRPrimary, primaryTplDir, i, "Primary")
  get_and_store_templates(datasetDCAxyCharm, gausExpMCxy, dcaxy, sigmaCharm, alphaRCharm, charmTplDir, i, "Charm")
  get_and_store_templates(datasetDCAxyBackground, gausExpMCxy, dcaxy, sigmaBkg, alphaRBkg, backgroundTplDir, i, "Background")
  alphaR.setRange(0.01, 2)
  get_and_store_templates(datasetDCAxyBeauty, gausExpMCxy, dcaxy, sigmaBeauty, alphaRBeauty, beautyTplDir, i, "Beauty")
outFile.cd()
hKSprob.Write()

## resolution
resoDCAxy = ROOT.RooRealVar("resoDCAxy", "#sigma_{resolution xy}", 1.e-3, 4.e-4, 12.e-4, "cm")
## primary MC template
sigmaDCAxyPrim = ROOT.RooRealVar("sigmaDCAxyPrim", "#sigma_{DCA_{xy}}", 1.e-3, 7.e-4, 1.e-2, "cm")
sigmaDCAxyPrimConv = ROOT.RooFormulaVar("sigmaDCAxyPrimConv", "#sigma_{DCA_{xy}}", "sqrt(@0 * @0 + @1 * @1)", ROOT.RooArgList(sigmaDCAxyPrim, resoDCAxy))
alphaRDCAxyPrim = ROOT.RooRealVar("alphaRDCAPrim", "#alpha_{R}", 1.5, 0.1, 10.)
alphaLDCAxyPrim = ROOT.RooFormulaVar("alphaLDCAPrim", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaRDCAxyPrim))
pdfPrompt = ROOT.RooGausDExp("pdfPrompt", "GausDExp", dcaxy, muDCAxy, sigmaDCAxyPrimConv, alphaLDCAxyPrim, alphaRDCAxyPrim)
## charm and beauty MC templates
sigmaDCAxyCharm = ROOT.RooRealVar("sigmaDCAxyCharm", "#sigma_{DCA_{xy}}", 1.e-3, 7.e-4, 1.e-2, "cm")
sigmaDCAxyCharmConv = ROOT.RooFormulaVar("sigmaDCAxyCharmConv", "#sigma_{DCA_{xy}}", "sqrt(@0 * @0 + @1 * @1)", ROOT.RooArgList(sigmaDCAxyCharm, resoDCAxy))
alphaRDCAxyCharm = ROOT.RooRealVar("alphaRDCACharm", "#alpha_{R}", 1.5, 0.1, 10.)
alphaLDCAxyCharm = ROOT.RooFormulaVar("alphaLDCACharm", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaRDCAxyCharm))
pdfCharm = ROOT.RooGausDExp("pdfCharm", "GausDExp", dcaxy, muDCAxy, sigmaDCAxyCharmConv, alphaLDCAxyCharm, alphaRDCAxyCharm)
sigmaDCAxyBeauty = ROOT.RooRealVar("sigmaDCAxyBeauty", "#sigma_{DCA_{xy}}", 1.e-3, 7.e-4, 1.e-2, "cm")
sigmaDCAxyBeautyConv = ROOT.RooFormulaVar("sigmaDCAxyBeautyConv", "#sigma_{DCA_{xy}}", "sqrt(@0 * @0 + @1 * @1)", ROOT.RooArgList(sigmaDCAxyBeauty, resoDCAxy))
alphaRDCAxyBeauty = ROOT.RooRealVar("alphaRDCABeauty", "#alpha_{R}", 1.5, 0.1, 10.)
alphaLDCAxyBeauty = ROOT.RooFormulaVar("alphaLDCABeauty", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaRDCAxyBeauty))
pdfBeauty = ROOT.RooGausDExp("pdfBeauty", "GausDExp", dcaxy, muDCAxy, sigmaDCAxyBeautyConv, alphaLDCAxyBeauty, alphaRDCAxyBeauty)
## background template
sigmaDCAxyBkg = ROOT.RooRealVar("sigmaDCAxyBkg", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-4, 1.e-2, "cm")
alphaRDCAxyBkg = ROOT.RooRealVar("alphaRDCABkg", "#alpha_{R}", 1.5, 0.1, 10.)
alphaLDCAxyBkg = ROOT.RooFormulaVar("alphaLDCABkg", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaRDCAxyBkg))
pdfBkg = ROOT.RooGausDExp("pdfBkg", "GausDExp", dcaxy, muDCAxy, sigmaDCAxyBkg, alphaLDCAxyBkg, alphaRDCAxyBkg)
## total pdf
nBkg = ROOT.RooRealVar("nBkg", "nBkg", 10., 0., 1.e8)
nPrompt = ROOT.RooRealVar("nPrompt", "nPrompt", 1.e4, 0., 1.e8)
nCharm = ROOT.RooRealVar("nCharm", "nCharm", 1.e4, 100, 1.e8)
nBeauty = ROOT.RooRealVar("nBeauty", "nBeauty", 1.e4, 0., 1.e8)
totalPdf = ROOT.RooAddPdf("totalPdf", "totalPdf", ROOT.RooArgList(pdfPrompt, pdfCharm, pdfBeauty, pdfBkg), ROOT.RooArgList(nPrompt, nCharm, nBeauty, nBkg))

## now fit the data
outFile.mkdir("fit_data")
outFile.cd("fit_data")

hPromptFrac = ROOT.TH1D("hPromptFrac", ";#it{p}_{T} (GeV/#it{c});Prompt fraction", conf['nPtBins'], conf['minPt'], conf['maxPt'])
hCharmFrac = ROOT.TH1D("hCharmFrac", ";#it{p}_{T} (GeV/#it{c});Charm fraction", conf['nPtBins'], conf['minPt'], conf['maxPt'])

for i in range(conf['nPtBins']):

  sigmaDCAxyPrim.setVal(sigmaPrimary.GetBinContent(i+1))
  sigmaDCAxyPrim.setConstant(True)
  alphaRDCAxyPrim.setVal(alphaRPrimary.GetBinContent(i+1))
  alphaRDCAxyPrim.setConstant(True)
  sigmaDCAxyCharm.setVal(sigmaCharm.GetBinContent(i+1))
  sigmaDCAxyCharm.setConstant(True)
  alphaRDCAxyCharm.setVal(alphaRCharm.GetBinContent(i+1))
  alphaRDCAxyCharm.setConstant(True)
  sigmaDCAxyBeauty.setVal(sigmaBeauty.GetBinContent(i+1))
  sigmaDCAxyBeauty.setRange(max(0, sigmaBeauty.GetBinContent(i+1) - 1 * sigmaBeauty.GetBinError(i+1)), sigmaBeauty.GetBinContent(i+1) + 1 * sigmaBeauty.GetBinError(i+1))
  sigmaDCAxyBeauty.setVal(sigmaBeauty.GetBinContent(i+1))
  alphaRDCAxyBeauty.setRange(max(0, alphaRBeauty.GetBinContent(i+1) - 1 * alphaRBeauty.GetBinError(i+1)), alphaRBeauty.GetBinContent(i+1) + 1 * alphaRBeauty.GetBinError(i+1))
  sigmaDCAxyBkg.setVal(sigmaBkg.GetBinContent(i+1))
  sigmaDCAxyBkg.setConstant(True)
  alphaRDCAxyBkg.setVal(alphaRBkg.GetBinContent(i+1))
  alphaRDCAxyBkg.setConstant(True)

  datasetDCAxyData = ROOT.RooDataHist(f"datasetDCAxyData_{i}", f"datasetDCAxyData_{i}", ROOT.RooArgList(dcaxy), hDCAxyData.ProjectionY(f"hDCAxyData_{i}", i+1, i+1))
  nBkg.setVal((1 - hPurity.GetBinContent(i+1)) * datasetDCAxyData.sumEntries())
  nBkg.setConstant(True)
  norm = datasetDCAxyData.sumEntries() * hPurity.GetBinContent(i+1)

  totalPdf.fitTo(datasetDCAxyData)
  plot = dcaxy.frame()
  datasetDCAxyData.plotOn(plot, ROOT.RooFit.MarkerSize(0.5), ROOT.RooFit.Name("data"))
  totalPdf.plotOn(plot, ROOT.RooFit.Name("model"))
  totalPdf.paramOn(plot,ROOT.RooFit.Label(f'#chi^{{2}}/NDF = {plot.chiSquare("model", "data")}'))
  totalPdf.plotOn(plot, ROOT.RooFit.Components("pdfCharm"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kViolet))
  totalPdf.plotOn(plot, ROOT.RooFit.Components("pdfBeauty"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))
  totalPdf.plotOn(plot, ROOT.RooFit.Components("pdfBkg"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange + 3))
  plot.SetMinimum(0.1)
  plot.Write(f"frameDCAxyData_{i}")
  hPromptFrac.SetBinContent(i+1, nPrompt.getVal() / norm)
  hPromptFrac.SetBinError(i+1, nPrompt.getError() / norm)
  hCharmFrac.SetBinContent(i+1, nCharm.getVal() / norm)
  hCharmFrac.SetBinError(i+1, nCharm.getError() / norm)


outFile.cd()
sigmaPrimary.Write()
sigmaCharm.Write()
sigmaBeauty.Write()
alphaRPrimary.Write()
alphaRCharm.Write()
alphaRBeauty.Write()

hPromptFrac.Write()
hCharmFrac.Write()

outFile.Close()
