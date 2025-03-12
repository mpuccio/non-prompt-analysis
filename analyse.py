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
alphaR = ROOT.RooRealVar("alphaR", "#alpha_{R}", 1.5, 0.1, 10.)
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
    dataset.plotOn(plot)
    gausExpMCxy.plotOn(plot)
    gausExpMCxy.paramOn(plot)
    sigmaHist.SetBinContent(i+1, sigmaDCAxy.getVal())
    sigmaHist.SetBinError(i+1, sigmaDCAxy.getError())
    alphaRHist.SetBinContent(i+1, alphaR.getVal())
    alphaRHist.SetBinError(i+1, alphaR.getError())
    tplDir.cd()

for i in range(conf['nPtBins']):
  datasetDCAxyData = ROOT.RooDataHist(f"datasetDCAxyData_{i}", f"datasetDCAxyData_{i}", ROOT.RooArgList(dcaxy), hDCAxyData.ProjectionY(f"hDCAxyData_{i}", i+1, i+1))
  datasetDCAxyBackground = ROOT.RooDataHist(f"datasetDCAxyBackground_{i}", f"datasetDCAxyBackground_{i}", ROOT.RooArgList(dcaxy), hDCAxyBackground.ProjectionY(f"hDCAxyBackground_{i}", i+1, i+1))
  datasetDCAxyPrimary = ROOT.RooDataHist(f"datasetDCAxyPrimary_{i}", f"datasetDCAxyPrimary_{i}", ROOT.RooArgList(dcaxy), hDCAxyPrimary.GetPtr().ProjectionY(f"hDCAxyPrimary_{i}", i+1, i+1))
  datasetDCAxyCharm = ROOT.RooDataHist(f"datasetDCAxyCharm_{i}", f"datasetDCAxyCharm_{i}", ROOT.RooArgList(dcaxy), hDCAxyCharm.ProjectionY(f"hDCAxyCharm_{i}", i+1, i+1))
  datasetDCAxyBeauty = ROOT.RooDataHist(f"datasetDCAxyBeauty_{i}", f"datasetDCAxyBeauty_{i}", ROOT.RooArgList(dcaxy), hDCAxyBeauty.ProjectionY(f"hDCAxyBeauty_{i}", i+1, i+1))

  get_and_store_templates(datasetDCAxyPrimary, gausExpMCxy, dcaxy, sigmaPrimary, alphaRPrimary, primaryTplDir, i, "Primary")
  get_and_store_templates(datasetDCAxyCharm, gausExpMCxy, dcaxy, sigmaCharm, alphaRCharm, charmTplDir, i, "Charm")
  get_and_store_templates(datasetDCAxyBeauty, gausExpMCxy, dcaxy, sigmaBeauty, alphaRBeauty, beautyTplDir, i, "Beauty")
  get_and_store_templates(datasetDCAxyBackground, gausExpMCxy, dcaxy, sigmaBkg, alphaRBkg, backgroundTplDir, i, "Background")


## resolution
resoDCAxy = ROOT.RooRealVar("resoDCAxy", "#sigma_{resolution xy}", 1.e-3, 1.e-4, 1.e-2, "cm")
sigmaDCAxyFit = ROOT.RooFormulaVar("sigmaDCAxyFit", "#sigma_{DCA_{xy}}", "TMath::Sqrt(@0 * @0 + @1 * @1)", ROOT.RooArgList(sigmaDCAxy, resoDCAxy))
## primary MC template
sigmaDCAxyPrim = ROOT.RooRealVar("sigmaDCAxyPrim", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-4, 1.e-2, "cm")
sigmaDCAxyPrimConv = ROOT.RooFormulaVar("sigmaDCAxyPrimConv", "#sigma_{DCA_{xy}}", "TMath::Sqrt(@0 * @0 + @1 * @1)", ROOT.RooArgList(sigmaDCAxyPrim, resoDCAxy))
alphaRDCAxyPrim = ROOT.RooRealVar("alphaDCAPrim", "#alpha_{R}", 1.5, 0.1, 10.)
alphaLDCAxyPrim = ROOT.RooFormulaVar("alphaDCAPrim", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaRDCAxyPrim))
pdfPrompt = ROOT.RooGausDExp("pdfPrompt", "GausDExp", dcaxy, muDCAxy, sigmaDCAxyPrimConv, alphaLDCAxyPrim, alphaRDCAxyPrim)
promptFrac = ROOT.RooRealVar("promptFrac", "Prompt fraction", 0.5, 0., 1.)
## charm and beauty MC templates
sigmaDCAxyCharm = ROOT.RooRealVar("sigmaDCAxyCharm", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-4, 1.e-2, "cm")
sigmaDCAxyCharmConv = ROOT.RooFormulaVar("sigmaDCAxyCharmConv", "#sigma_{DCA_{xy}}", "TMath::Sqrt(@0 * @0 + @1 * @1)", ROOT.RooArgList(sigmaDCAxyCharm, resoDCAxy))
alphaRDCAxyCharm = ROOT.RooRealVar("alphaDCACharm", "#alpha_{R}", 1.5, 0.1, 10.)
alphaLDCAxyCharm = ROOT.RooFormulaVar("alphaDCACharm", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaRDCAxyCharm))
pdfCharm = ROOT.RooGausDExp("pdfCharm", "GausDExp", dcaxy, muDCAxy, sigmaDCAxyCharmConv, alphaLDCAxyCharm, alphaRDCAxyCharm)
charmFrac = ROOT.RooRealVar("charmFrac", "Charm fraction", 0.5, 0., 1.)
sigmaDCAxyBeauty = ROOT.RooRealVar("sigmaDCAxyBeauty", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-4, 1.e-2, "cm")
sigmaDCAxyBeautyConv = ROOT.RooFormulaVar("sigmaDCAxyBeautyConv", "#sigma_{DCA_{xy}}", "TMath::Sqrt(@0 * @0 + @1 * @1)", ROOT.RooArgList(sigmaDCAxyBeauty, resoDCAxy))
alphaRDCAxyBeauty = ROOT.RooRealVar("alphaDCABeauty", "#alpha_{R}", 1.5, 0.1, 10.)
alphaLDCAxyBeauty = ROOT.RooFormulaVar("alphaDCABeauty", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaRDCAxyBeauty))
pdfBeauty = ROOT.RooGausDExp("pdfBeauty", "GausDExp", dcaxy, muDCAxy, sigmaDCAxyBeautyConv, alphaLDCAxyBeauty, alphaRDCAxyBeauty)
## background template
sigmaDCAxyBkg = ROOT.RooRealVar("sigmaDCAxyBkg", "#sigma_{DCA_{xy}}", 1.e-3, 1.e-4, 1.e-2, "cm")
alphaRDCAxyBkg = ROOT.RooRealVar("alphaDCABkg", "#alpha_{R}", 1.5, 0.1, 10.)
alphaLDCAxyBkg = ROOT.RooFormulaVar("alphaDCABkg", "#alpha_{L}", "-@0", ROOT.RooArgList(alphaRDCAxyBkg))
pdfBkg = ROOT.RooGausDExp("pdfBkg", "GausDExp", dcaxy, muDCAxy, sigmaDCAxyBkg, alphaLDCAxyBkg, alphaRDCAxyBkg)
## total pdf
pdfNonPrompt = ROOT.RooAddPdf("pdfNonPrompt", "pdfNonPrompt", ROOT.RooArgList(pdfCharm, pdfBeauty), ROOT.RooArgList(charmFrac))
signalPdf = ROOT.RooAddPdf("signalPdf", "signalPdf", ROOT.RooArgList(pdfPrompt, pdfNonPrompt), ROOT.RooArgList(promptFrac))
purityFrac = ROOT.RooRealVar("purityFrac", "Purity fraction", 0.5, 0., 1.)
totalPdf = ROOT.RooAddPdf("totalPdf", "totalPdf", ROOT.RooArgList(signalPdf, pdfBkg), ROOT.RooArgList(purityFrac))


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
  sigmaDCAxyBeauty.setConstant(True)
  alphaRDCAxyBeauty.setVal(alphaRBeauty.GetBinContent(i+1))
  alphaRDCAxyBeauty.setConstant(True)
  sigmaDCAxyBkg.setVal(sigmaBkg.GetBinContent(i+1))
  alphaRDCAxyBkg.setVal(alphaRBkg.GetBinContent(i+1))
  alphaRDCAxyBkg.setConstant(True)

  purityFrac.setVal(hPurity.GetBinContent(i+1))
  purityFrac.setConstant(True)

  datasetDCAxyData = ROOT.RooDataHist(f"datasetDCAxyData_{i}", f"datasetDCAxyData_{i}", ROOT.RooArgList(dcaxy), hDCAxyData.ProjectionY(f"hDCAxyData_{i}", i+1, i+1))
  datasetDCAxyData.plotOn(dcaxy.frame())
  totalPdf.fitTo(datasetDCAxyData)
  plot = dcaxy.frame()
  datasetDCAxyData.plotOn(plot)
  totalPdf.plotOn(plot)
  totalPdf.paramOn(plot)
  plot.Write(f"frameDCAxyData_{i}")
  hPromptFrac.SetBinContent(i+1, promptFrac.getVal())
  hPromptFrac.SetBinError(i+1, promptFrac.getError())
  hCharmFrac.SetBinContent(i+1, charmFrac.getVal())
  hCharmFrac.SetBinError(i+1, charmFrac.getError())


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



