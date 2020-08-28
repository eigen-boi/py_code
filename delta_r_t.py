# Aidan Patton, August 2020

import ROOT
from rat import dsreader
import math
import sys


# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def histMaker(input_files, delta_r_range, delta_t_range, retriggerfilter):

    """

    Description:
    Extract data from root file(s), and return histograms "h_delta_r_name" and "h_delta_t_name"

    Parameters:
    input_files : list of strs
        root files to extract data from
    delta_r_range list of floats
        Range to plot on x axis (units of mm)
    delta_t_range : list of floats
        Range to plot on x axis (units of micro s)
    retriggerfilter : bool 
        Whether to only use first event in an entry

    Returns:
    h_delta_r - the TH1D delta_r histogram
    h_delta_t - TH1D delta_t histogram

    """

    fitName = "partialFitter"

    nbins = 100
    Counts = [0, 0, 0, 0, 0]        # CountTotal, CountTrigger, CountRadius, CountValid, LowECountFiltered

    h_delta_r = ROOT.TH1D("h_delta_r_name", "#Delta r", nbins, delta_r_range[0], delta_r_range[1])
    h_delta_r.SetDirectory(0)
    h_delta_t = ROOT.TH1D("h_delta_t_name", "#Delta t", nbins, delta_t_range[0], delta_t_range[1])
    h_delta_t.SetDirectory(0)

    skip_pair = False

    for fname in input_files :

        dsRead = dsreader(fname) # really shouldn't be looping, could do ds, run = dsRead.next() to do next individual one

        for ds0, run in dsRead: # ds is equivalent to dsreader.GetEntry(i) in c++

            if skip_pair :
                print "skipp"
                skip_pair = False
                continue

            for iev0 in range(0, ds0.GetEVCount()):

                Counts[0] += 1                                           # total event count

                if is_mc and retriggerfilter and iev0 > 0 :
                    skip_pair = True
                    break

                bi_ev = ds0.GetEV(iev0)

                fVertex = bi_ev.GetFitResult(fitName).GetVertex(0)
                bi_fit_pos = fVertex.GetPosition()
                bi_ev_time = bi_ev.GetClockCount50()*20/1000

                # can delete if they're all zero, since they should have been checked already (just checking to make sure)
                if not fVertex.ValidEnergy():
                    Counts[1] += 1                                       # doesn't have valid energy
                    skip_pair = True
                    break

                if bi_fit_pos.Mag() > 6000 or bi_fit_pos.Z() < 747.5 :                  
                    Counts[2] += 1
                    skip_pair = True
                    break

                for ds1, run1 in dsRead : 
                    for iev1 in range(0, ds1.GetEVCount()) :

                        po_ev = ds1.GetEV(iev1)

                        po_fit_pos = po_ev.GetFitResult(fitName).GetVertex(0).GetPosition()
                        po_ev_time = po_ev.GetClockCount50()*20/1000

                        delta_r = (po_fit_pos - bi_fit_pos).Mag()
                        delta_t = po_ev_time - bi_ev_time

                        if delta_t > delta_t_range[1] or delta_t < delta_t_range[0] :
                            skip_pair = True
                            break

                        Counts[3] += 1 
                        h_delta_r.Fill(delta_r)
                        h_delta_t.Fill(delta_t)

                        break # break iev1 loop, should do anyways but to be safe
                    break # break ds1 loop, back to iev0 loop 
                break # break iev0 loop and start the cycle again at ds0

    return h_delta_r, h_delta_t
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def plot_delta_r_t(h_delta_r, h_delta_t, Counts, delta_r_range, delta_t_range, outFile):
    
    """ 

    Description:
    Plot the two input histograms on different pads, export single pdf
    
    Parameters:
    h_delta_r     - pre-processed data in a histogram labelled "h_delta_r_name"
    h_delta_t     - pre-processed data in a histogram labelled "h_delta_t_name"
    Counts        - counts from the filtered data, as they make it through each of the filtering gauntlets 
                  - CountTotal, bi invalid energy, bi outside FV
    delta_r_range - range to plot on x axis (units of mm)
    delta_t_range - range to plot on x axis (units of ns)
    outFile       - str, name of output pdf

    Returns:
    return - The histogram plot TCanvas

    """

    if not h_delta_r or not h_delta_t:
        print "Failed to get data histogram(s)"
        sys.exit(1)

    nbins = 100                                                       # for histogram formatting
    delta_r_binsize = (delta_r_range[1]-delta_r_range[0])/nbins       # [mm]
    delta_t_binsize = (delta_t_range[1]-delta_t_range[0])/(nbins)     # [mus]

    h_delta_r.SetDirectory(0)
    h_delta_r.SetStats(0)
    h_delta_t.SetDirectory(0)
    h_delta_t.SetStats(0)


    canvas = ROOT.TCanvas("canvas") # make canvas
    canvas.Divide(1,2)

    # =================================================================
    canvas.cd(1)

    max_bin_value = h_delta_r.GetMaximum()

    h_delta_r.SetLineColor(ROOT.kBlue)
    h_delta_r.GetYaxis().SetTitleOffset(1.1)
    h_delta_r.GetYaxis().SetTitle("Count per %.1f mm bin"%(delta_r_binsize))
    h_delta_r.GetXaxis().SetTitle("#Delta R")
    h_delta_r.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)

    '''
    expoFit00 = ROOT.TF1("expofit", "gaus", delta_r_range[0], delta_r_range[1])            # auto-determine range for expo fitting
    h_delta_r.Fit(expoFit00, "ERL")
    fit00 = []
    fit00.append(h_delta_r.GetMean())                                                  # mean
    fit00.append(expoFit00.GetParameter(2))                                           # stdev

    expoFit0 = ROOT.TF1("expofit", "gaus", fit00[0]-5*fit00[1], fit00[0]+5*fit00[1]) # actual exponential fit
    h_delta_r.Fit(expoFit0, "ERL")

    fit0 = []
    fit0.append(expoFit0.GetChisquare())  # chi20       0
    fit0.append(expoFit0.GetNDF())        # ndf0        1
    fit0.append(expoFit0.GetParameter(1)) # mean0       2
    fit0.append(expoFit0.GetParError(1))  # mean0error  3
    fit0.append(expoFit0.GetParameter(2)) # sigma0      4
    fit0.append(expoFit0.GetParError(2))  # sigma0error 5
    '''

    h_delta_r.Draw()

    latex_delta_r = ROOT.TLatex()
    latex_delta_r.SetTextFont(62)
    latex_delta_r.SetNDC()
    latex_delta_r.SetTextSize(0.04)
    xloc = 0.15
    yloc = 0.75
    latex_delta_r.DrawLatex(xloc-0.01, yloc+0.04, "#Delta r")
    latex_delta_r.SetTextSize(0.03)
    latex_delta_r.DrawText(xloc, yloc, "%i Entries"%(h_delta_r.GetEntries()))
    latex_delta_r.DrawLatex(xloc, yloc-0.04, "Mean = %.1f mm"%(h_delta_r.GetMean()))
    # =================================================================
    canvas.cd(2)

    max_bin_value = h_delta_t.GetMaximum()

    h_delta_t.SetLineColor(ROOT.kGreen) # hi
    h_delta_t.GetYaxis().SetTitleOffset(1.1)
    h_delta_t.GetYaxis().SetTitle("Count per %.1f #mus bin"%(delta_t_binsize))
    h_delta_t.GetXaxis().SetTitle("#Delta t [#mus]")
    h_delta_t.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)

    expoFit1 = ROOT.TF1("expofit", "expo", delta_t_range[0], delta_t_range[1]) # actual exponential fit
    h_delta_t.Fit(expoFit1, "ERL")

    fit1 = []
    fit1.append(expoFit1.GetParameter(0))
    fit1.append(expoFit1.GetParError(0))
    fit1.append(expoFit1.GetParameter(1)) # mean0       2
    fit1.append(expoFit1.GetParError(1))  # mean0error  3

    h_delta_t.Draw()

    latex_delta_t = ROOT.TLatex()
    latex_delta_t.SetTextFont(62)
    latex_delta_t.SetNDC()
    latex_delta_t.SetTextSize(0.04)
    xloc = 0.25
    yloc = 0.7
    latex_delta_t.DrawLatex(xloc-0.01, yloc+0.04, "#Delta t")
    latex_delta_t.SetTextSize(0.03)
    latex_delta_t.DrawText(xloc, yloc, "%i Entries"%(h_delta_t.GetEntries()))
    latex_delta_t.DrawLatex(xloc, yloc-0.04, "Mean = %.1f #mus"%(h_delta_t.GetMean()))
    latex_delta_t.DrawLatex(xloc, yloc-0.08, "#tau = %.1f #pm %.1f #mus"%(abs(1/fit1[2]), abs(fit1[3]/fit1[2]**2)))

    # =================================================================

    canvas.Print(outFile)

    print "\n"
    
    return canvas
    
    # =================================================================
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////