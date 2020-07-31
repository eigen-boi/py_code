# Aidan Patton, June 2020
# Based on code by S. Langrock (plot_fit.py)
# Also heavily borrowing from
# https://indico.cern.ch/event/704163/contributions/2936719/attachments/1693833/2726445/Tutorial-PyROOT.pdf
# and Energy.C by Yang Zhang

import ROOT
from rat import dsreader
import math
import sys

def histMaker(input_files, nHitsRange, QRange, hist_display_title):

    """

    Description:
    Extract data from root file(s), apply retrigger filtering, and return histograms "h_nHits_name" and "h_totcharge_name"
    # contains commented out code to save histogram "h_nHits_name" as a root file 

    Parameters:
    input_files - root files to extract data from
    nHitsRange - range to plot on x axis (units of hits)
    QRange - range to plot on x axis (units of Coulombs)
    hist_display_title - list containing the titles of the histograms, unused

    Returns:
    h_nHits - the TH1D nHits histogram
    h_totcharge - TH1D charge histogram
    Counts - list  [CountTotal, CountTrigger, CountValid, CountRadius, LowECountFiltered]

    """

    fitName = "partialFitter"

    #ERange = [0.0, 1.5]             # ERange and ECut as outputs so can use in next script too
    #ERange = [0, 9]
    #ECut = 0.2                      # cutoff below which events are considered "low energy" [MeV]

    nbins = 100
    Counts = [0, 0, 0, 0, 0]        # CountTotal, CountTrigger, CountRadius, CountValid, LowECountFiltered

    h_nHits = ROOT.TH1D("h_nHits_name", "nHits Fit", nbins, nHitsRange[0], nHitsRange[1])
    h_nHits.SetDirectory(0)
    h_totcharge = ROOT.TH1D("h_totcharge_name", "Total Charge Fit", nbins, QRange[0], QRange[1])
    h_totcharge.SetDirectory(0)

    for fname in input_files :

        for ds, run in dsreader(fname): # ds is equivalent to dsreader.GetEntry(i) in c++

            #mc = ds.GetMC()
            #mcparticle = mc.GetMCParticle(0)
            #mcenergy = mcparticle.GetKineticEnergy()

            for iev in range(0, ds.GetEVCount()):

                Counts[0] += 1                                           # total event count

                if iev > 0 :
                    continue                                             # re-trigger filter
                Counts[1] += 1                                           # events that aren't re-triggers count

                ev = ds.GetEV(iev)

                if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                    continue

                fVertex = ev.GetFitResult(fitName).GetVertex(0)
                PosEV = fVertex.GetPosition()
                REV = PosEV.Mag()
                RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                # see if inside the AV
                if PosEV.Z() >= 6000 and RhoEV >= 730 :
                    continue
                if PosEV.Z() < 6000 and REV > 6000 :
                    continue
                Counts[2] += 1                                           # inside AV count

                if not fVertex.ValidEnergy():
                    continue
                Counts[3] += 1                                           # how many actually pass the gauntlet

                h_nHits.Fill(ev.GetNhits())
                h_totcharge.Fill(ev.GetTotalCharge())

                if fVertex.GetEnergy() < 0.2: #ECut:
                    Counts[4] += 1                                       # keep track of LowE which pass the filter

    return h_nHits, h_totcharge, Counts

def plot_nHits_totcharge(h_nHits, h_totcharge, Counts, nHitsRange, QRange):
    
    """ 

    Description:
    Plot the two input histograms on different pads, export single pdf
    
    Parameters:
    h_nHits     - pre-processed data in a histogram labelled "h_nHits_name"
    h_totcharge - pre-processed data in a histogram labelled "h_totcharge_name"
    Counts      - counts from the filtered data, as they make it through each of the filtering gauntlets 
                - CountTotal, CountTrigger, CountRadius, CountValid, LowECountFiltered
    nHitsRange - range to plot on x axis (units of hits)
    QRange - range to plot on x axis (units of Coulombs)

    Returns:
    return - The histogram plot TCanvas

    """

    nbins = 100                                                   # for histogram formatting
    nHitsbinsize = (nHitsRange[1]-nHitsRange[0])/nbins            # [number of hits]
    Qbinsize = (QRange[1]-QRange[0])/nbins                        # [C]

    h_nHits.SetDirectory(0)
    h_nHits.SetStats(0)
    h_totcharge.SetDirectory(0)
    h_totcharge.SetStats(0)

    if not h_nHits and not h_totcharge:
        print "Failed to get both data histograms"
        sys.exit(1)
    elif not h_nHits and h_totcharge:
        print "Failed to get nHits"
        sys.exit(1)
    elif h_nHits and not h_totcharge:
        print "Failed to get totcharge"

    canvas = ROOT.TCanvas("canvas")                                                    # make canvas
    canvas.Divide(1,2)

    # =================================================================
    canvas.cd(1)

    max_bin_value = h_nHits.GetMaximum()

    h_nHits.SetLineColor(ROOT.kRed)
    h_nHits.GetYaxis().SetTitleOffset(1.1)
    h_nHits.GetYaxis().SetTitle("Count per %.1f hits bin"%(nHitsbinsize))
    h_nHits.GetXaxis().SetTitle("nHits")
    h_nHits.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)

    gaussFit00 = ROOT.TF1("gaussfit", "gaus", nHitsRange[0], nHitsRange[1])            # auto-determine range for gauss fitting
    h_nHits.Fit(gaussFit00, "ERL")
    fit00 = []
    fit00.append(gaussFit00.GetParameter(1))                                           # mean
    fit00.append(gaussFit00.GetParameter(2))                                           # stdev

    gaussFit0 = ROOT.TF1("gaussfit", "gaus", fit00[0]-5*fit00[1], fit00[0]+5*fit00[1]) # actual gaussian fit
    h_nHits.Fit(gaussFit0, "ERL")

    fit0 = []
    fit0.append(gaussFit0.GetChisquare())  # chi20       0
    fit0.append(gaussFit0.GetNDF())        # ndf0        1
    fit0.append(gaussFit0.GetParameter(1)) # mean0       2
    fit0.append(gaussFit0.GetParError(1))  # mean0error  3
    fit0.append(gaussFit0.GetParameter(2)) # sigma0      4
    fit0.append(gaussFit0.GetParError(2))  # sigma0error 5

    latex_nHits = ROOT.TLatex()
    latex_nHits.SetTextFont(62)
    latex_nHits.SetNDC()
    latex_nHits.SetTextSize(0.04)
    xloc = 0.74
    yloc = 0.6
    latex_nHits.DrawText(xloc-0.01, yloc+0.04, "nHits")
    latex_nHits.SetTextSize(0.03)
    latex_nHits.DrawText(xloc, yloc, "%i Entries"%(h_nHits.GetEntries()))
    latex_nHits.DrawText(xloc, yloc-0.04, "Mean = %.3f +/- %.4f MeV"%(fit0[2], fit0[3]))
    latex_nHits.DrawText(xloc, yloc-0.08, "stdev = %.3f +/- %.4f"%(fit0[4], fit0[5]))
    latex_nHits.DrawText(xloc, yloc-0.12, "chi^2/ndf = %.1f/%i"%(fit0[0], fit0[1]))
    # =================================================================
    canvas.cd(2)

    max_bin_value = h_totcharge.GetMaximum()

    h_totcharge.SetLineColor(ROOT.kBlue)
    h_totcharge.GetYaxis().SetTitleOffset(1.1)
    h_totcharge.GetYaxis().SetTitle("Count per %.1f C bin"%(Qbinsize))
    h_totcharge.GetXaxis().SetTitle("Fitted Total Charge [C]")
    h_totcharge.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)

    gaussFit01 = ROOT.TF1("gaussfit", "gaus", QRange[0], QRange[1])                    # auto-determine range for gauss fitting
    h_totcharge.Fit(gaussFit01, "ERL")
    fit01 = []
    fit01.append(gaussFit01.GetParameter(1))                                           # mean
    fit01.append(gaussFit01.GetParameter(2))                                           # stdev

    gaussFit1 = ROOT.TF1("gaussfit", "gaus", fit01[0]-5*fit01[1], fit01[0]+5*fit01[1]) # actual gaussian fit
    h_totcharge.Fit(gaussFit1, "ERL")

    fit1 = []
    fit1.append(gaussFit1.GetChisquare())  # chi20       0
    fit1.append(gaussFit1.GetNDF())        # ndf0        1
    fit1.append(gaussFit1.GetParameter(1)) # mean0       2
    fit1.append(gaussFit1.GetParError(1))  # mean0error  3
    fit1.append(gaussFit1.GetParameter(2)) # sigma0      4
    fit1.append(gaussFit1.GetParError(2))  # sigma0error 5

    h_totcharge.Draw()

    latex_Q = ROOT.TLatex()
    latex_Q.SetTextFont(62)
    latex_Q.SetNDC()
    latex_Q.SetTextSize(0.04)
    xloc = 0.74
    yloc = 0.6
    latex_Q.DrawText(xloc-0.01, yloc+0.04, "totcharge")
    latex_Q.SetTextSize(0.03)
    latex_Q.DrawText(xloc, yloc, "%i Entries"%(h_totcharge.GetEntries()))
    latex_Q.DrawText(xloc, yloc-0.04, "Mean = %.3f +/- %.4f C"%(fit1[2], fit1[3]))
    latex_Q.DrawText(xloc, yloc-0.08, "stdev = %.3f +/- %.4f"%(fit1[4], fit1[5]))
    latex_Q.DrawText(xloc, yloc-0.12, "chi^2/ndf = %.1f/%i"%(fit1[0], fit1[1]))

    # =================================================================

    canvas.Print("/home/eigenboi/pdfs/nHits_totcharge_point_10MeV_AV_00.pdf")

    print "CountTotal = ", Counts[0]
    print "CountTrigger", Counts[1]
    print "CountRadius = ", Counts[2]
    print "CountValid =", Counts[3]
    print "LowECountFiltered = ", Counts[4]
    
    return canvas
    
    # =================================================================

if __name__ == '__main__':

    input_files_p = [['/home/eigenboi/root_files/sim_point_1MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_05.root'], 
                     ['/home/eigenboi/root_files/sim_point_2p5MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_05.root'],    
                     ['/home/eigenboi/root_files/sim_point_5MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_05.root'],    
                     ['/home/eigenboi/root_files/sim_point_7p5MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_05.root'],
                     ['/home/eigenboi/root_files/sim_point_10MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_05.root']]
    
    input_files = [['/home/eigenboi/root_files/sim_point_10MeV_AV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_05_tablecheck.root']]

    #nHitsRange = [0.0, 600] # 1 MeV
    nHitsRange = [1800.0, 2800] # 10 MeV
    QRange  = [60000, 160000] # 10 MeV

    hist_display_title = ["nHits and totcharge Fit", "nHits and totcharge Fit", "nHits and totcharge Fit", "nHits and totcharge Fit", "nHits and totcharge Fit"]

    h_nHits, h_totcharge, Counts = histMaker(input_files[0], nHitsRange, QRange, hist_display_title[4])
    plot_nHits_totcharge(h_nHits, h_totcharge, Counts, nHitsRange, QRange)

