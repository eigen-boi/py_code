# Aidan Patton, June 2020

import ROOT
from rat import dsreader
import math
import sys


# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def histMaker(input_files, NhitsRange, QRange, hist_display_title, retriggerfilter):

    """

    Description:
    Extract data from root file(s), apply retrigger filtering, and return histograms "h_Nhits_name" and "h_totcharge_name"

    Parameters:
    input_files - root files to extract data from
    NhitsRange - range to plot on x axis (units of hits)
    QRange - range to plot on x axis (units of Coulombs)
    hist_display_title - list containing the titles of the histograms, unused
    retriggerfilter - bool 

    Returns:
    h_Nhits - the TH1D Nhits histogram
    h_totcharge - TH1D charge histogram
    Counts - list  [CountTotal, CountTrigger, CountValid, CountRadius, LowECountFiltered]

    """

    fitName = "partialFitter"

    nbins = 100
    Counts = [0, 0, 0, 0, 0]        # CountTotal, CountTrigger, CountRadius, CountValid, LowECountFiltered

    h_Nhits = ROOT.TH1D("h_Nhits_name", "Nhits Fit", nbins, NhitsRange[0], NhitsRange[1])
    h_Nhits.SetDirectory(0)
    h_totcharge = ROOT.TH1D("h_totcharge_name", "Total Charge Fit", nbins, QRange[0], QRange[1])
    h_totcharge.SetDirectory(0)

    for fname in input_files :

        for ds, run in dsreader(fname): # ds is equivalent to dsreader.GetEntry(i) in c++

            for iev in range(0, ds.GetEVCount()):

                Counts[0] += 1                                           # total event count

                if retriggerfilter and iev > 0 :
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

                should_be_charge = fVertex.GetPositiveEnergyError()
                if Counts[0] % 100 == 0 : print should_be_charge

                h_Nhits.Fill(ev.GetNhits())
                h_totcharge.Fill(should_be_charge)

                if fVertex.GetEnergy() < 0.2: #ECut:
                    Counts[4] += 1                                       # keep track of LowE which pass the filter

    return h_Nhits, h_totcharge, Counts
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def plot_Nhits_totcharge(h_Nhits, h_totcharge, Counts, NhitsRange, QRange):
    
    """ 

    Description:
    Plot the two input histograms on different pads, export single pdf
    
    Parameters:
    h_Nhits     - pre-processed data in a histogram labelled "h_Nhits_name"
    h_totcharge - pre-processed data in a histogram labelled "h_totcharge_name"
    Counts      - counts from the filtered data, as they make it through each of the filtering gauntlets 
                - CountTotal, CountTrigger, CountRadius, CountValid, LowECountFiltered
    NhitsRange - range to plot on x axis (units of hits)
    QRange - range to plot on x axis (units of Coulombs)

    Returns:
    return - The histogram plot TCanvas

    """

    nbins = 100                                                   # for histogram formatting
    Nhitsbinsize = (NhitsRange[1]-NhitsRange[0])/nbins            # [number of hits]
    Qbinsize = (QRange[1]-QRange[0])/nbins                        # [C]

    h_Nhits.SetDirectory(0)
    h_Nhits.SetStats(0)
    h_totcharge.SetDirectory(0)
    h_totcharge.SetStats(0)

    if not h_Nhits and not h_totcharge:
        print "Failed to get both data histograms"
        sys.exit(1)
    elif not h_Nhits and h_totcharge:
        print "Failed to get Nhits"
        sys.exit(1)
    elif h_Nhits and not h_totcharge:
        print "Failed to get totcharge"

    canvas = ROOT.TCanvas("canvas")                                                    # make canvas
    canvas.Divide(1,2)

    # =================================================================
    canvas.cd(1)

    max_bin_value = h_Nhits.GetMaximum()

    h_Nhits.SetLineColor(ROOT.kRed)
    h_Nhits.GetYaxis().SetTitleOffset(1.1)
    h_Nhits.GetYaxis().SetTitle("Count per %.1f hits bin"%(Nhitsbinsize))
    h_Nhits.GetXaxis().SetTitle("Nhits")
    h_Nhits.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)
    '''
    gaussFit00 = ROOT.TF1("gaussfit", "gaus", NhitsRange[0], NhitsRange[1])            # auto-determine range for gauss fitting
    h_Nhits.Fit(gaussFit00, "ERL")
    fit00 = []
    fit00.append(gaussFit00.GetParameter(1))                                           # mean
    fit00.append(gaussFit00.GetParameter(2))                                           # stdev

    gaussFit0 = ROOT.TF1("gaussfit", "gaus", fit00[0]-5*fit00[1], fit00[0]+5*fit00[1]) # actual gaussian fit
    h_Nhits.Fit(gaussFit0, "ERL")

    fit0 = []
    fit0.append(gaussFit0.GetChisquare())  # chi20       0
    fit0.append(gaussFit0.GetNDF())        # ndf0        1
    fit0.append(gaussFit0.GetParameter(1)) # mean0       2
    fit0.append(gaussFit0.GetParError(1))  # mean0error  3
    fit0.append(gaussFit0.GetParameter(2)) # sigma0      4
    fit0.append(gaussFit0.GetParError(2))  # sigma0error 5
    '''
    h_Nhits.Draw()

    latex_Nhits = ROOT.TLatex()
    latex_Nhits.SetTextFont(62)
    latex_Nhits.SetNDC()
    latex_Nhits.SetTextSize(0.04)
    xloc = 0.74
    yloc = 0.6
    latex_Nhits.DrawText(xloc-0.01, yloc+0.04, "Nhits")
    latex_Nhits.SetTextSize(0.03)
    latex_Nhits.DrawText(xloc, yloc, "%i Entries"%(h_Nhits.GetEntries()))
    #latex_Nhits.DrawLatex(xloc, yloc-0.04, "Mean = %.3f #pm %.4f MeV"%(fit0[2], fit0[3]))
    #latex_Nhits.DrawLatex(xloc, yloc-0.08, "#sigma = %.3f #pm %.4f"%(fit0[4], fit0[5]))
    #latex_Nhits.DrawLatex(xloc, yloc-0.12, "#Chi^{2}/ndf = %.1f/%i"%(fit0[0], fit0[1]))
    # =================================================================
    canvas.cd(2)

    max_bin_value = h_totcharge.GetMaximum()

    h_totcharge.SetLineColor(ROOT.kBlue)
    h_totcharge.GetYaxis().SetTitleOffset(1.1)
    h_totcharge.GetYaxis().SetTitle("Count per %.1f C bin"%(Qbinsize))
    h_totcharge.GetXaxis().SetTitle("Fitted Total Charge [C]")
    h_totcharge.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)

    '''
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
    '''
    h_totcharge.Draw()

    latex_Q = ROOT.TLatex()
    latex_Q.SetTextFont(62)
    latex_Q.SetNDC()
    latex_Q.SetTextSize(0.04)
    xloc = 0.74
    yloc = 0.6
    latex_Q.DrawText(xloc-0.01, yloc+0.04, "Total Charge")
    latex_Q.SetTextSize(0.03)
    latex_Q.DrawText(xloc, yloc, "%i Entries"%(h_totcharge.GetEntries()))
    #latex_Q.DrawLatex(xloc, yloc-0.04, "Mean = %.3f #pm %.4f C"%(fit1[2], fit1[3]))
    #latex_Q.DrawLatex(xloc, yloc-0.08, "#sigma = %.3f #pm %.4f"%(fit1[4], fit1[5]))
    #latex_Q.DrawLatex(xloc, yloc-0.12, "#Chi^{2}/ndf = %.1f/%i"%(fit1[0], fit1[1]))

    # =================================================================

    canvas.Print(outDir + outFile)

    print "CountTotal = ", Counts[0]
    print "CountTrigger", Counts[1]
    print "CountRadius = ", Counts[2]
    print "CountValid =", Counts[3]
    print "LowECountFiltered = ", Counts[4]
    
    return canvas
    
    # =================================================================
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def histMaker_combo(input_files_MC, input_files_data, NhitsRange, QRange, retriggerfilter):

    """

    Description:
    Extract data from root file(s), apply retrigger filtering, and return count normalized MC and data Nhits and charge histograms 

    Parameters:
    input_files_MC/data : list of str 
        root files to extract data from
    NhitsRange : list or floats
        Range to plot on x axis (units of hits)
    QRange : list of floats
        Range to plot on x axis (units of Coulombs)
    retriggerfilter : bool 
        Whether to only use first event in an entry

    Returns:
    h_Nhits - the TH1D Nhits histogram
    h_totcharge - TH1D charge histogram

    """

    fitName = "partialFitter"

    nbins = 100
    Counts = [0, 0, 0, 0, 0]        # CountTotal, CountTrigger, CountRadius, CountValid, LowECountFiltered

    h_Nhits_MC = ROOT.TH1D("h_Nhits_name", "Normalized Nhits Fit", nbins, NhitsRange[0], NhitsRange[1])
    h_Nhits_MC.SetDirectory(0)
    h_totcharge_MC = ROOT.TH1D("h_totcharge_name", "Normalized Charge Fit", nbins, QRange[0], QRange[1])
    h_totcharge_MC.SetDirectory(0)

    h_Nhits_data = ROOT.TH1D("h_Nhits_name", "Normalized Nhits Fit", nbins, NhitsRange[0], NhitsRange[1])
    h_Nhits_data.SetDirectory(0)
    h_totcharge_data = ROOT.TH1D("h_totcharge_name", "Normalized Charge Fit", nbins, QRange[0], QRange[1])
    h_totcharge_data.SetDirectory(0)

    for fname in input_files_MC :

        for ds, run in dsreader(fname): # ds is equivalent to dsreader.GetEntry(i) in c++

            for iev in range(0, ds.GetEVCount()):

                if retriggerfilter and iev > 0 :
                    continue                                             # re-trigger filter

                ev = ds.GetEV(iev)

                if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                    continue

                fVertex = ev.GetFitResult(fitName).GetVertex(0)
                PosEV = fVertex.GetPosition()
                REV = PosEV.Mag()
                RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                # see if inside the AV
                if REV > 6000 or PosEV.Z() < 747.5 :
                    continue

                if not fVertex.ValidEnergy(): continue
                finalE = fVertex.GetEnergy()
                if finalE < 0.5: continue

                h_Nhits_MC.Fill(ev.GetNhits())
                h_totcharge_MC.Fill(ev.GetTotalCharge())

    for fname in input_files_data :

        for ds, run in dsreader(fname): # ds is equivalent to dsreader.GetEntry(i) in c++

            for iev in range(0, ds.GetEVCount()):

                if retriggerfilter and iev > 0 :
                    continue                                             # re-trigger filter

                ev = ds.GetEV(iev)

                if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                    continue

                fVertex = ev.GetFitResult(fitName).GetVertex(0)
                PosEV = fVertex.GetPosition()
                REV = PosEV.Mag()
                RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                # see if inside the AV
                if REV > 6000 or PosEV.Z() < 747.5 :
                    continue

                if not fVertex.ValidEnergy():
                    continue

                should_be_charge = fVertex.GetPositiveEnergyError()

                h_Nhits_data.Fill(ev.GetNhits())
                h_totcharge_data.Fill(should_be_charge)

    h_Nhits_MC.Scale(1/h_Nhits_MC.GetEntries())
    h_Nhits_data.Scale(1/h_Nhits_data.GetEntries())
    h_totcharge_MC.Scale(1/h_totcharge_MC.GetEntries())
    h_totcharge_data.Scale(1/h_totcharge_data.GetEntries())

    return h_Nhits_MC, h_Nhits_data, h_totcharge_MC, h_totcharge_data
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def plot_Nhits_totcharge_combo(h_Nhits_MC, h_Nhits_data, h_totcharge_MC, h_totcharge_data, NhitsRange, QRange):
    
    """ 

    Description:
    Plot the two input histograms on different pads, save single pdf
    
    Parameters:
    -----------
    h_Nhits_MC/data : TH1D     
        Pre-processed data in a histogram labelled "h_Nhits_name"
    h_totcharge_MC/data : TH1D 
        Pre-processed data in a histogram labelled "h_totcharge_name"
    NhitsRange : list 
        Range to plot on x axis (units of hits)
    QRange : list 
        Range to plot on x axis (units of Coulombs)

    Returns:
    canvas : The histogram plot TCanvas

    """

    nbins = 100                                                   # for histogram formatting
    Nhitsbinsize = (NhitsRange[1]-NhitsRange[0])/nbins            # [number of hits]
    Qbinsize = (QRange[1]-QRange[0])/nbins                        # [C]

    h_Nhits_MC.SetDirectory(0)
    h_Nhits_MC.SetStats(0)
    h_totcharge_MC.SetDirectory(0)
    h_totcharge_MC.SetStats(0)

    h_Nhits_data.SetDirectory(0)
    h_Nhits_data.SetStats(0)
    h_totcharge_data.SetDirectory(0)
    h_totcharge_data.SetStats(0)


    canvas = ROOT.TCanvas("canvas")                                                    # make canvas
    canvas.Divide(1,2)

    # =================================================================
    canvas.cd(1)

    max_bin_value = max( h_Nhits_MC.GetMaximum(), h_Nhits_data.GetMaximum() )

    h_Nhits_MC.SetLineColor(ROOT.kRed)
    h_Nhits_data.SetLineColor(7)

    h_Nhits_data.GetYaxis().SetTitleOffset(1.1)
    h_Nhits_data.GetYaxis().SetTitle("Normalized Counts")
    h_Nhits_data.GetXaxis().SetTitle("Nhits")
    h_Nhits_data.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)
    h_Nhits_MC.GetYaxis().SetTitleOffset(1.1)
    h_Nhits_MC.GetYaxis().SetTitle("Normalized Counts")
    h_Nhits_MC.GetXaxis().SetTitle("Nhits")
    h_Nhits_MC.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)

    h_Nhits_data.Draw()
    h_Nhits_MC.Draw("same")

    legendhasitNhits = ROOT.TLegend(0.75, 0.79, 0.9, 0.9) 
    legendhasitNhits.AddEntry(h_Nhits_MC, "Nhits MC") 
    legendhasitNhits.AddEntry(h_Nhits_data, "Nhits Data")
    legendhasitNhits.SetLineWidth(0)
    legendhasitNhits.Draw("same")

    latex_Nhits = ROOT.TLatex()
    latex_Nhits.SetTextFont(62)
    latex_Nhits.SetNDC()
    xloc = 0.74
    yloc = 0.6

    latex_Nhits.SetTextSize(0.04)
    latex_Nhits.DrawText(xloc-0.01, yloc+0.04, "Nhits MC")
    latex_Nhits.SetTextSize(0.03)
    latex_Nhits.DrawText(xloc, yloc, "%i Entries"%(h_Nhits_MC.GetEntries()))

    latex_Nhits.SetTextSize(0.04)
    latex_Nhits.DrawText(xloc-0.01, yloc-0.04, "Nhits Data")
    latex_Nhits.SetTextSize(0.03)
    latex_Nhits.DrawText(xloc, yloc-0.08, "%i Entries"%(h_Nhits_data.GetEntries()))
    # =================================================================
    canvas.cd(2)

    max_bin_value = max( h_totcharge_MC.GetMaximum(), h_totcharge_data.GetMaximum() )

    h_totcharge_MC.SetLineColor(ROOT.kBlue)
    h_totcharge_data.SetLineColor(ROOT.kGreen)

    h_totcharge_data.GetYaxis().SetTitleOffset(1.1)
    h_totcharge_data.GetYaxis().SetTitle("Normalized Counts")
    h_totcharge_data.GetXaxis().SetTitle("Fitted Total Charge [C]")
    h_totcharge_data.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)
    h_totcharge_MC.GetYaxis().SetTitleOffset(1.1)
    h_totcharge_MC.GetYaxis().SetTitle("Normalized Counts")
    h_totcharge_MC.GetXaxis().SetTitle("Fitted Total Charge [C]")
    h_totcharge_MC.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)
    

    h_totcharge_data.Draw()
    h_totcharge_MC.Draw("same")

    legendhasitQ = ROOT.TLegend(0.75, 0.79, 0.9, 0.9) 
    legendhasitQ.AddEntry(h_totcharge_MC, "Charge MC")
    legendhasitQ.AddEntry(h_totcharge_data, "Charge Data")
    legendhasitQ.SetLineWidth(0)
    legendhasitQ.Draw("same")

    latex_Q = ROOT.TLatex()
    latex_Q.SetTextFont(62)
    latex_Q.SetNDC()
    xloc = 0.74
    yloc = 0.6

    latex_Q.SetTextSize(0.04)
    latex_Q.DrawText(xloc-0.01, yloc+0.04, "Total Charge MC")
    latex_Q.SetTextSize(0.03)
    latex_Q.DrawText(xloc, yloc, "%i Entries"%(h_totcharge_MC.GetEntries()))

    latex_Q.SetTextSize(0.04)
    latex_Q.DrawText(xloc-0.01, yloc-0.04, "Total Charge Data")
    latex_Q.SetTextSize(0.03)
    latex_Q.DrawText(xloc, yloc-0.08, "%i Entries"%(h_totcharge_data.GetEntries()))

    # =================================================================

    canvas.Print(outDir + outFile)
    
    return canvas
    
    # =================================================================
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////