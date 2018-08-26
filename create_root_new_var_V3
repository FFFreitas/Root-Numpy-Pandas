#!/usr/bin/python
import sys
import ROOT 
import numpy as np
from ROOT import TLorentzVector
import csv
import pandas as pd
from ROOT import TFile, TTree
from rootpy.io import root_open
from rootpy.tree import Tree, TreeChain
from rootpy.plotting import Hist
from rootpy.extern.six.moves import range
from root_numpy import hist2array, root2array

if len(sys.argv) < 2:
  print " Usage: Example1.py input_file"
  sys.exit(1)

ROOT.gSystem.Load("/home/felipe/madanalysis5_1_5/tools/delphes/libDelphes")
inputFile = sys.argv[1]

# Create chain of root trees
chain1 = ROOT.TChain("Delphes")
chain1.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain1)
numberOfEntries = treeReader.GetEntries()

# create new root file
root_name = raw_input("name of new root: ")
csv_name = raw_input("name of new csv: ")
f = root_open(root_name, "recreate")
tree = Tree("test")
tree.create_branches({'DETA_jj': 'F',
	'Dphi_jj': 'F',
	'DPHI_aa': 'F',
	'DPHI_HJ1': 'F',
	'M_jj': 'F',
	'PT_a1': 'F',
	'PT_a2': 'F',
	'PT_H': 'F',
	'M_H': 'F',
	'Eta_H': 'F',
	'Phi_H': 'F',
	'Cos_Ha1': 'F',
	'PT_j1': 'F',
	'PT_j2': 'F'
	})

# Get pointers to branches used in this analysis
branchJet = treeReader.UseBranch("Jet")
branchElectron = treeReader.UseBranch("Electron")
branchMuon = treeReader.UseBranch("Muon")
branchPhoton = treeReader.UseBranch("Photon")
branchMET = treeReader.UseBranch("MissingET")
# Loop over all events
for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
	treeReader.ReadEntry(entry)
##################################################
	photons = []
	for an in xrange(branchPhoton.GetEntries()):
		photons.append(branchPhoton.At(an))

	photons = sorted(photons, key=lambda pa: pa.P4().Pt(), reverse=True)

	a1 = ROOT.TLorentzVector()
	a2 = ROOT.TLorentzVector()
	a1.SetPtEtaPhiE(photons[0].P4().Pt(),photons[0].P4().Eta(),photons[0].P4().Phi(),photons[0].P4().E())
	a2.SetPtEtaPhiE(photons[1].P4().Pt(),photons[1].P4().Eta(),photons[1].P4().Phi(),photons[1].P4().E())
	jato1 = ROOT.TLorentzVector()
	jato2 = ROOT.TLorentzVector()
####################################################################################
	bjets, ljets = [], []
	for n in xrange(branchJet.GetEntries()):
		if branchJet.At(n).BTag == 1:
			bjets.append(branchJet.At(n))
		else:
			ljets.append(branchJet.At(n))

	if len(ljets) >= 2:
		ljets = sorted(ljets, key=lambda Jet:  Jet.P4().Pt(), reverse=True)
	else:
		continue
####################################################################################
	jato1.SetPtEtaPhiE(ljets[0].P4().Pt(),ljets[0].P4().Eta(),ljets[0].P4().Phi(),ljets[0].P4().E())
	jato2.SetPtEtaPhiE(ljets[1].P4().Pt(),ljets[1].P4().Eta(),ljets[1].P4().Phi(),ljets[1].P4().E())
####################################################################################
	if((jato1 + jato2).M() > 450 and 115 < (a1 + a2).M() < 135) and np.abs(jato1.Eta() - jato2.Eta()) >= 2.8:
		tree.DETA_jj = np.abs(jato1.Eta() - jato2.Eta())
		tree.Dphi_jj = np.abs((jato1).DeltaPhi(jato2))
		tree.PT_j1 = (jato1).Pt()
		tree.PT_j2 = (jato2).Pt()
		tree.M_jj = (jato1 + jato2).M()
		H = ROOT.TLorentzVector()
		H = (a1 + a2)
		tree.DPHI_aa = np.abs((a1).DeltaPhi(a2))
		tree.PT_a1 = a1.Pt()
		tree.PT_a2 = a2.Pt()
		tree.PT_H = H.Pt()
		tree.M_H = H.M()
		tree.Eta_H = H.Eta()
		tree.Phi_H = H.Phi()
		tree.DPHI_HJ1 = np.abs((jato1).DeltaPhi(H))
################# boosted objects ################################################
                Htob = ROOT.TLorentzVector()
                Htob.SetPxPyPzE(H.Px(),H.Py(),H.Pz(),H.E())
                Hboost = ROOT.TVector3()
                Hboost = Htob.BoostVector()
                a1.Boost(-Hboost)
                v = Hboost.Unit()
                tree.Cos_Ha1 = np.cos(a1.Angle(v))
		tree.Fill()

# Show resulting histograms
#hist_PT_l1.Draw()
#raw_input("Press Enter to continue...")
tree.write()
f.close()
#create the csv output

to_convert = root2array(root_name,'test')

df_conv = pd.DataFrame(to_convert)

df_conv.to_csv( csv_name + '.csv', index=False, header= df_conv.keys(), mode='w', sep=' ')

