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
from rootpy.plotting import Hist2D
from rootpy.extern.six.moves import range
from root_numpy import hist2array, root2array
from itertools import combinations, permutations

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
tree.create_branches({'PT_l': 'F',
	'MT_VH': 'F',
	'PT_VH': 'F',
	'PT_W': 'F',
        'Cos_lw': 'F',
	'DPHI_lmet': 'F',
	'met': 'F',
	'PT_b1': 'F',
	'PT_b2': 'F',
	'PT_lj1': 'F',
	'PT_lj2': 'F',
        'Eta_H': 'F',
	'Phi_H': 'F',
        'M_H': 'F',
        'MT_W': 'F',
        'Cos_Hb1': 'F',
	'PT_H': 'F',
	})


# Get pointers to branches used in this analysis
branchJet = treeReader.UseBranch("Jet")
branchElectron = treeReader.UseBranch("Electron")
branchMuon = treeReader.UseBranch("Muon")
branchPhoton = treeReader.UseBranch("Photon")
branchMET = treeReader.UseBranch("MissingET")
####################################################################
# Loop over all events
for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
	treeReader.ReadEntry(entry)
##########################################################################################################
	eletrons = sorted(branchElectron, key=lambda Electron: Electron.PT, reverse=True)
        missing = sorted(branchMET, key=lambda MisingET: MisingET.MET, reverse=True)
	elec1 = eletrons[0]
        eletron1 = ROOT.TLorentzVector()
	eletron1.SetPtEtaPhiE(elec1.PT,elec1.Eta,elec1.Phi,elec1.P4().E())
        met = ROOT.TLorentzVector()
	met.SetPtEtaPhiE(missing[0].P4().Pt(),missing[0].P4().Eta(),missing[0].P4().Phi(),missing[0].P4().E())
        bjato1 = ROOT.TLorentzVector()
        bjato2 = ROOT.TLorentzVector()
        jato1 = ROOT.TLorentzVector()
        jato2 = ROOT.TLorentzVector()
####################################################################################
        bjets, ljets = [], []
        for n in xrange(branchJet.GetEntries()):
                if branchJet.At(n).BTag == 1:
                        bjets.append(branchJet.At(n))
                else:
                        ljets.append(branchJet.At(n))

        if len(bjets) >= 2:
                bjets = sorted(bjets, key=lambda BJet:  BJet.P4().Pt(), reverse=True)
        else:
                continue

        if len(ljets) >= 2:
                ljets = sorted(ljets, key=lambda Jet:  Jet.P4().Pt(), reverse=True)
        else:
                continue
####################################################################################
        jato1.SetPtEtaPhiE(ljets[0].P4().Pt(),ljets[0].P4().Eta(),ljets[0].P4().Phi(),ljets[0].P4().E())
        jato2.SetPtEtaPhiE(ljets[1].P4().Pt(),ljets[1].P4().Eta(),ljets[1].P4().Phi(),ljets[1].P4().E())
####################################################################################
        bjato1.SetPtEtaPhiE(bjets[0].P4().Pt(),bjets[0].P4().Eta(),bjets[0].P4().Phi(),bjets[0].P4().E())
        bjato2.SetPtEtaPhiE(bjets[1].P4().Pt(),bjets[1].P4().Eta(),bjets[1].P4().Phi(),bjets[1].P4().E())
####################################################################################
	if 115 < (bjato1 + bjato2).M() < 135:
		tree.PT_l = (eletron1).Pt()
		tree.met = np.abs(met.Mt())
		tree.PT_b1 = (bjato1).Pt()
		tree.PT_b2 = (bjato2).Pt()
		tree.PT_lj1 = jato1.Pt()
		tree.PT_lj2 = jato2.Pt()
		tree.PT_H = (bjato1 + bjato2).Pt()
		tree.Eta_H = (bjato1 + bjato2).Eta()
		W = ROOT.TLorentzVector()
		W = (eletron1 + met)
		tree.DPHI_lmet = np.abs(eletron1.DeltaPhi(met))
		tree.MT_W = np.sqrt(2*np.abs(met.Et())*np.abs(eletron1.Pt())*(1-np.cos(eletron1.DeltaPhi(met))))
		tree.PT_W = W.Pt()
                H = ROOT.TLorentzVector()
       	        H = (bjato1 + bjato2)
		tree.MT_VH = (W + H).Mt() #H.Mt() + np.sqrt(2*np.abs(met.Et())*np.abs(eletron1.Pt())*(1-np.cos(eletron1.DeltaPhi(met))))
		tree.PT_VH = ((bjato1 + bjato2) + (eletron1 + met)).Pt()
		tree.Phi_H = H.Phi()
		tree.M_H = H.M()
#########################boosted objects#########################################################
                Wtob = ROOT.TLorentzVector()
                Wtob.SetPxPyPzE(W.Px(),W.Py(),W.Pz(),W.E())
                Wboost = ROOT.TVector3()
                Wboost = Wtob.BoostVector()
                v = Wboost.Unit()
                Htob = ROOT.TLorentzVector()
                Htob.SetPxPyPzE(H.Px(),H.Py(),H.Pz(),H.E())
                Hboost = ROOT.TVector3()
                Hboost = Htob.BoostVector()
                ang = Hboost.Unit()
                bjato1.Boost(-Hboost)
                tree.Cos_Hb1 = np.cos(bjato1.Angle(ang))
                eletron1.Boost(-Wboost)
                tree.Cos_lw = np.cos(eletron1.Angle(v))
		tree.Fill()

###############################################
tree.write()
f.close()


#create the csv output

to_convert = root2array(root_name,'test')

df_conv = pd.DataFrame(to_convert)

df_conv.to_csv( csv_name + '.csv', index=False, header= df_conv.keys(), mode='w', sep=' ')
