#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TChain.h>
#include <TKey.h>

using namespace ROOT::RDF;

void createChain(TChain& chainD, const std::string &filename, const std::string &treename) {
  TFile fileD(filename.c_str());
  TIter next(fileD.GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    std::string keyName = key->GetName();
    if (keyName.find("DF_") != std::string::npos) {
      chainD.Add((filename + "/" + keyName + "/" + chainD.GetName()).c_str());
    }
  }
}

auto adaptAndFilter(ROOT::RDataFrame &df, std::string filters) {
  auto columnNames = df.GetColumnNames();
  if (std::find(columnNames.begin(), columnNames.end(), "fBachKaonNClusTPC") != columnNames.end()) {
    return df.Define("fBachNClusTPC", "fBachKaonNClusTPC").Define("fBachHasTOF", "fBachKaonHasTOF").Filter(filters);
  } else {
    return df.Filter(filters);
  }
}

void createSnapshots() {
  std::string filters = "std::abs(fPvZ)<10 && std::abs(fProtonEta)<0.9 && std::abs(fPionEta)<0.9 && std::abs(fBachEta)<0.9 && fProtonNClusTPC>100 && fPionNClusTPC>100 && fBachNClusTPC>100 && std::abs(fProtonTPCNSigma)<4 && std::abs(fPionTPCNSigma)<4 && std::abs(fBachKaonTPCNSigma)<4";
  std::string mcOnlyFilters = "fMCcollisionMatch && std::abs(fPDGcode) == 3334";
  std::string dataOnlyFilters = "(fProtonHasTOF || fPionHasTOF || fBachHasTOF) && fCollisionTimeRes < 20 &&  std::abs(fMassXi-1.32171) > 0.008 && std::abs(fMassOmega-1.67245) < 0.3";

  ROOT::EnableImplicitMT();

  TChain chainDataMB("O2npcasctable");
  createChain(chainDataMB, "data/AO2D_mb.root", "O2npcasctable");
  auto dfDataMB_orig = ROOT::RDataFrame(chainDataMB);
  auto dfDataMB = adaptAndFilter(dfDataMB_orig, filters + " && " + dataOnlyFilters);

  TChain chainData("O2npcasctable");
  createChain(chainData, "data/AO2D.root", "O2npcasctable");
  auto dfData_orig = ROOT::RDataFrame(chainData);
  auto dfData = adaptAndFilter(dfData_orig, filters + " && " + dataOnlyFilters);

  TChain chainPrimary("O2npcasctablemc");
  createChain(chainPrimary, "mc/primary/AO2D.root", "O2npcastablemc");
  auto dfPrimary_orig = ROOT::RDataFrame(chainPrimary);
  auto dfPrimary = adaptAndFilter(dfPrimary_orig, filters + " && " + mcOnlyFilters).Define("fIsPrimary", "!fIsFromCharm && !fIsFromBeauty");

  TChain chainCharm("O2npcasctablemc");
  createChain(chainCharm, "mc/charm/AO2D.root", "O2npcastablemc");
  auto dfCharm_orig = ROOT::RDataFrame(chainCharm);
  auto dfCharm = adaptAndFilter(dfCharm_orig, filters + " && " + mcOnlyFilters).Define("fIsPrimary", "!fIsFromCharm && !fIsFromBeauty");

  TChain chainBeauty("O2npcasctablemc");
  createChain(chainBeauty, "mc/beauty/AO2D.root", "O2npcastablemc");
  auto dfBeauty_orig = ROOT::RDataFrame(chainBeauty);
  auto dfBeauty = adaptAndFilter(dfBeauty_orig, filters + " && " + mcOnlyFilters).Define("fIsPrimary", "!fIsFromCharm && !fIsFromBeauty");

  dfDataMB.Snapshot("tree", "snapshots/data_mb.root", {"fMassOmega", "fCascPt", "fCascDCAxy", "fCascDCAz"});
  dfData.Snapshot("tree", "snapshots/data.root", {"fMassOmega", "fCascPt", "fCascDCAxy", "fCascDCAz"});
  dfPrimary.Filter("fIsPrimary").Snapshot("tree", "snapshots/primary.root",{"fMassOmega", "fCascPt", "fCascDCAxy", "fCascDCAz", "fIsPrimary", "fIsFromCharm", "fIsFromBeauty"});
  dfCharm.Filter("fIsFromCharm && !fIsFromBeauty").Snapshot("tree", "snapshots/charm.root",{"fMassOmega", "fCascPt", "fCascDCAxy", "fCascDCAz", "fIsPrimary", "fIsFromCharm", "fIsFromBeauty"});
  dfBeauty.Filter("fIsFromBeauty").Snapshot("tree", "snapshots/beauty.root",{"fMassOmega", "fCascPt", "fCascDCAxy", "fCascDCAz", "fIsPrimary", "fIsFromCharm", "fIsFromBeauty"});
}
