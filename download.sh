mkdir -p data
mkdir -p mc/beauty
mkdir -p mc/charm
mkdir -p mc/primary
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0036/361688/84941/AO2D.root file:mc/beauty/LHC24g5.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0036/361503/84829/AO2D.root file:mc/beauty/LHC24e3.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0036/361532/84889/AO2D.root file:mc/beauty/LHC24d3.root
# alien_cp alien:/alice/cern.ch/user/m/mpuccio/non-prompt-omega-mc/AO2D.root file:mc/charm/AO2D_full.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0035/357838/82806/AO2D.root file:mc/charm/AO2D.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0035/356381/82317/AO2D.root file:mc/primary/AO2D.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0035/356560/82408/AO2D.root file:data/AO2D_mb.root
alien_cp alien:/alice/cern.ch/user/a/alihyperloop/outputs/0034/348595/79079/AO2D.root file:data/AO2D.root
