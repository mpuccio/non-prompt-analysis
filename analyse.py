import ROOT
import yaml


ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro('inc/RooGausDExp.cxx+')

with open("configurations.yaml", "r") as stream:
  confFile = yaml.safe_load(stream)
  conf = confFile["analysis"]

dataDfMB = ROOT.RDataFrame("tree", "snapshots/data_mb.root")
dataDf = ROOT.RDataFrame("tree", "snapshots/data.root")
primaryDf = ROOT.RDataFrame("tree", "snapshots/primary.root")
charmDf = ROOT.RDataFrame("tree", "snapshots/charm.root")
beautyDf = ROOT.RDataFrame("tree", "snapshots/beauty.root")

templateFile = ROOT.TFile(confFile["template_preparation"]["outputFilename"])
hPurity = templateFile.Get("purity/hPurity")
hMean = templateFile.Get("purity/hMean")
hSigma = templateFile.Get("purity/hSigma")

averageMass = hMean.GetFunction("pol0").GetParameter(0)
sigPars = [hSigma.GetFunction("pol4").GetParameter(i) for i in range(5)]
sigmaFilter = f"({sigPars[0]} + {sigPars[1]} * fCascPt + {sigPars[2]} * fCascPt * fCascPt + {sigPars[3]} * fCascPt * fCascPt * fCascPt + {sigPars[4]} * fCascPt * fCascPt * fCascPt * fCascPt)"

fAlphaRPrimary = templateFile.Get("template_shapes/hAlphaRPrimary").GetFunction("expo")
fAlphaRCharm = templateFile.Get("template_shapes/hAlphaRCharm").GetFunction("expo")
fAlphaRBeauty = templateFile.Get("template_shapes/hAlphaRBeauty").GetFunction("expo")
fAlphaRBkg = templateFile.Get("template_shapes/hAlphaRBkg").GetFunction("expo")
fSigmaPrimary = templateFile.Get("template_shapes/hSigmaPrimary").GetFunction("expo")
fSigmaCharm = templateFile.Get("template_shapes/hSigmaCharm").GetFunction("expo")
fSigmaBeauty = templateFile.Get("template_shapes/hSigmaBeauty").GetFunction("expo")
fSigmaBkg = templateFile.Get("template_shapes/hSigmaBkg").GetFunction("expo")

signalDataDf = dataDf.Filter(f"fMassOmega > {averageMass} - 2 * {sigmaFilter} && fMassOmega < {averageMass} + 2 * {sigmaFilter}")
backgroundDataDf = dataDf.Filter(f"fMassOmega < {averageMass} - 5 * {sigmaFilter} || fMassOmega > {averageMass} + 5 * {sigmaFilter}")

outFile = ROOT.TFile(conf["outputFilename"], "recreate")
hDCAxyData = signalDataDf.Histo2D(("hDCAxyData", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")
hDCAxyBackground = backgroundDataDf.Histo2D(("hDCAxyBackground", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")
hDCAxyPrimary = primaryDf.Histo2D(("hDCAxyPrimary", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")
hDCAxyCharm = charmDf.Histo2D(("hDCAxyCharm", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")
hDCAxyBeauty = beautyDf.Histo2D(("hDCAxyBeauty", ";#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)", conf['nPtBins'], conf['minPt'], conf['maxPt'], conf['nBinsDCAxy'], -conf['maxAbsDCAxy'], conf['maxAbsDCAxy']), "fCascPt", "fCascDCAxy")

dcaxy = ROOT.RooRealVar("dcaxy", "DCA_{xy}", 0., -conf['maxAbsDCAxy'], conf['maxAbsDCAxy'], "cm")
muDCAxy = ROOT.RooRealVar("muDCAxy", "#mu_{DCA_{xy}}", 0., -conf['maxAbsDCAxy'], conf['maxAbsDCAxy'], "cm")

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
  ptMax = hDCAxyData.GetXaxis().GetBinUpEdge(i+1)
  ptMin = hDCAxyData.GetXaxis().GetBinLowEdge(i+1)
  deltaPt = ptMax - ptMin

  sigmaDCAxyPrim.setVal(fSigmaPrimary.Integral(ptMin, ptMax) / deltaPt)
  sigmaDCAxyPrim.setConstant(True)
  alphaRDCAxyPrim.setVal(fAlphaRPrimary.Integral(ptMin, ptMax) / deltaPt)
  alphaRDCAxyPrim.setConstant(True)
  sigmaDCAxyCharm.setVal(fSigmaCharm.Integral(ptMin, ptMax) / deltaPt)
  sigmaDCAxyCharm.setConstant(True)
  alphaRDCAxyCharm.setVal(fAlphaRCharm.Integral(ptMin, ptMax) / deltaPt)
  alphaRDCAxyCharm.setConstant(True)
  sigmaDCAxyBeauty.setVal(fSigmaBeauty.Integral(ptMin, ptMax) / deltaPt)
  sigmaDCAxyBeauty.setRange(fSigmaBeauty.Integral(ptMin, ptMax) / deltaPt * (1 - conf['beautyTplUnc']), fSigmaBeauty.Integral(ptMin, ptMax) / deltaPt * (1 + conf['beautyTplUnc']))
  alphaRDCAxyBeauty.setRange(fAlphaRBeauty.Integral(ptMin, ptMax) / deltaPt * (1 - conf['beautyTplUnc']), fAlphaRBeauty.Integral(ptMin, ptMax) / deltaPt * (1 + conf['beautyTplUnc']))
  alphaRDCAxyBeauty.setConstant(False)
  sigmaDCAxyBkg.setVal(fSigmaBkg.Integral(ptMin, ptMax) / deltaPt)
  sigmaDCAxyBkg.setConstant(True)
  alphaRDCAxyBkg.setVal(fAlphaRBkg.Integral(ptMin, ptMax) / deltaPt)
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

hPromptFrac.Write()
hCharmFrac.Write()

outFile.Close()
