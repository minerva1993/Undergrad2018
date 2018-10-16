import ROOT, os
from ROOT import *

if not os.path.exists('hist_pdf'):
  os.makedirs('hist_pdf')

gROOT.SetBatch();

f1 = TFile.Open('./ana_sig.root')
f2 = TFile.Open('./ana_ttbar.root')

#Get list of histogram's name
hist_names = [x.GetName() for x in f1.GetListOfKeys()]
hist_names.remove("tree");
c = TCanvas("a", "a", 400, 400)
l = TLegend(0.7,0.75,0.9,0.9)

for names in hist_names:
  h1 = f1.Get(names)
  h2 = f2.Get(names)
  h1.SetLineColor(ROOT.kRed)
  h2.SetLineColor(ROOT.kBlue)
  h1.SetLineWidth(2)
  h2.SetLineWidth(2)
  l.AddEntry(h1,"Signal","l")
  l.AddEntry(h2,"TTbar","l")

  h1_max = 0
  h2_max = 0
  if h1.Integral() > 0: h1_max = h1.GetMaximum()
  if h2.Integral() > 0: h2_max = h2.GetMaximum()
  h_max = max(h1_max, h2_max)
  h1.SetMaximum(h_max)
  h1.SetStats(0)

  if h1.Integral() > 0: h1.DrawNormalized("hist")
  else: h1.Draw("hist")
  if h2.Integral() > 0: h2.DrawNormalized("hist same")
  else: h2.Draw("hist same")
  l.Draw("same")

  c.Print('hist_pdf/' + names + '.pdf')
  l.Clear()
