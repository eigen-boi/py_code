# Aidan Patton, June 2020
# Based on code by S. Langrock (plot_fit.py)
# Also heavily influenced by 
# https://indico.cern.ch/event/704163/contributions/2936719/attachments/1693833/2726445/Tutorial-PyROOT.pdf
# and Energy.C by Yang Zhang

print "\n"

import ROOT
from rat import dsreader
import sys
import math

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def EHistMaker(input_files, ERange, ECut, hist_display_title, retriggerfilter): 

    """

    Description:
    Extract data from root file(s) fitted with the PartialFitter, apply cuts, and return energy histogram with name "h_energy_name" 

    Parameters:
    input_files : list of str 
        Name of rat root file to process
    ERange : list of floats
        The energy range for the histogram
    ECut : float
        Energy below which events are cut
    retriggerfilter : bool 
        Whether to only use first event in an entry


    Returns:
    h_energy : TH1D histogram
        The energy histogram
    Counts : list 
        CountTotal, CountTrigger, CountValid, CountRadius, LowECountFiltered

    """

    fitName = "partialFitter"

    nbins = 100
    Counts = [0, 0, 0, 0, 0]        # CountTotal, CountTrigger, CountRadius, CountValid, LowECountFiltered

    h_energy = ROOT.TH1D("h_energy_name", hist_display_title, nbins, ERange[0], ERange[1])
    h_energy.SetDirectory(0)

    for fname in input_files :

        for ds, run in dsreader(fname): # ds is equivalent to dsreader.GetEntry(i) in c++

            for iev in range(0, ds.GetEVCount()):

                Counts[0] += 1                                           # total event count

                if retriggerfilter and iev > 0 :
                    Counts[1] += 1 
                    continue                                             # re-trigger filter
                                                            

                ev = ds.GetEV(iev)

                if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                    continue # valid position cut

                fVertex = ev.GetFitResult(fitName).GetVertex(0)
                PosEV = fVertex.GetPosition()
                REV = PosEV.Mag() 
                RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                
                # FV cut
                if REV > 6000 or PosEV.Z() < 747.5 :                  
                    continue                                             # if outside AV, don't use
                Counts[2] += 1                                           # inside AV count
                
                '''
                # see if inside the AV (including the neck)
                if PosEV.Z() >= 6000 and RhoEV >= 730 :             
                    continue
                if PosEV.Z() < 6000 and REV > 6000 :
                    continue
                Counts[2] += 1                                           
                '''

                '''  
                if REV > 5500 : # FV cut instead will give better fit, boundary events are tricky 
                    continue
                Counts[2] += 1   
                '''  

                if not fVertex.ValidEnergy(): continue # valid energy cut
                finalE = fVertex.GetEnergy()
                if finalE < ECut: continue # low energy cut
                Counts[3] += 1                                           # how many actually pass the gauntlet

                h_energy.Fill(finalE) 

                if fVertex.GetEnergy() < ECut: 
                    Counts[4] += 1                                       # keep track of LowE which pass the filter 

    #h_energy.Scale(1/h_energy.GetEntries()) # scale by entries
    
    return h_energy, Counts
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def plot_fit_energy(input_files, h_energy, compare_uncleaned, Counts, ERange, ECut, plotTitle):
    
    """ 

    Description:
    Compare pre and post retrigger event filtering/screening, automatically determines a range for gauss fit of energy, save as pdf
    
    Parameters:
    h_energy    - pre-processed data in a histogram labelled "h_energy_name"
    input_files - list of name(s) of unfiltered root file(s) that h_energy was made from
    Counts      - counts from the filtered data, as they make it through each of the filtering gauntlets 
                - CountTotal, CountTrigger, CountRadius, CountValid, LowECountFiltered
    ERange      - energy range of the histogram h_energy_range [MeV]
    ECut        - arbitrary cutoff below which events are considered "low energy" [MeV]

    Returns:
    return - The histogram plot TCanvas

    """

    nbins = 100                                                   # for histogram formatting
    binsize = 1000*(ERange[1]-ERange[0])/nbins                    # [keV]

    h_fit_energy_0 = h_energy
    h_fit_energy_0.SetDirectory(0)
    h_fit_energy_0.SetStats(0)

    if not h_fit_energy_0 :
        print "Failed to get data histogram"
        sys.exit(1)

    # create low E histogram from original unfiltered data
    if compare_uncleaned :
        OriginalCount = 0
        OriginalCountValidFit = 0
        LowECountUnfiltered = 0
        LowEnbins = int(nbins*(ECut-ERange[0])/(ERange[1]-ERange[0]))
        h_fit_energy_3 = ROOT.TH1D("QpartialFitter", "Removing Retrigger Events Test", LowEnbins, ERange[0], ECut)
        h_fit_energy_3.SetDirectory(0)
        h_fit_energy_3.SetStats(0)
        for fname in input_files :
            for ds, run in dsreader(fname):
                for iev in range(0, ds.GetEVCount()):

                    ev = ds.GetEV(iev)
                    OriginalCount += 1

                    fVertex = ev.GetFitResult("partialFitter").GetVertex(0)

                    if not fVertex.ValidEnergy():
                        continue

                    OriginalCountValidFit += 1

                    if ev.GetDefaultFitVertex().GetEnergy() < ECut : 
                        LowECountUnfiltered += 1

                    h_fit_energy_3.Fill(ev.GetDefaultFitVertex().GetEnergy())

    canvas = ROOT.TCanvas("canvas")                                                    # make canvas       
    canvas.cd()

    if compare_uncleaned:
        max_bin_value = max([h_fit_energy_0.GetMaximum(), h_fit_energy_3.GetMaximum()])
    else:
        max_bin_value = h_fit_energy_0.GetMaximum()

    max_bin_value = h_fit_energy_0.GetMaximum()    

    if compare_uncleaned :
        h_fit_energy_3.SetLineColor(ROOT.kBlue)                                        # format histograms
    h_fit_energy_0.SetLineColor(ROOT.kRed)
    h_fit_energy_0.GetYaxis().SetTitleOffset(1.1)
    h_fit_energy_0.GetYaxis().SetTitle("Count per %.1f keV bin"%(binsize)) 
    h_fit_energy_0.GetXaxis().SetTitle("Fitted energy [MeV]")
    h_fit_energy_0.GetYaxis().SetRangeUser(0, max_bin_value + 0.1*max_bin_value)

    gaussFit00 = ROOT.TF1("gaussfit", "gaus", ECut, ERange[1])                         # auto-determine range for gauss fitting
    h_fit_energy_0.Fit(gaussFit00, "ERL")                    
    fit00 = []
    fit00.append(h_fit_energy_0.GetMean())                                             # mean 
    fit00.append(gaussFit00.GetParameter(2))                                           # stdev

    gaussFit0 = ROOT.TF1("gaussfit", "gaus", fit00[0]-5*fit00[1], fit00[0]+5*fit00[1]) # actual gaussian fit
    h_fit_energy_0.Fit(gaussFit0, "ERL")    
    
    fit0 = [] 
    fit0.append(gaussFit0.GetChisquare())  # chi20       0
    fit0.append(gaussFit0.GetNDF())        # ndf0        1
    fit0.append(gaussFit0.GetParameter(1)) # mean0       2
    fit0.append(gaussFit0.GetParError(1))  # mean0error  3
    fit0.append(gaussFit0.GetParameter(2)) # sigma0      4
    fit0.append(gaussFit0.GetParError(2))  # sigma0error 5

    h_fit_energy_0.Draw()
    if compare_uncleaned :
        h_fit_energy_3.Draw("same")

    legend = ROOT.TLegend(0.73, 0.73, 0.9, 0.9) 
    legend.AddEntry(h_fit_energy_0, "Re-Trigger Filtered Data") # hist made by retriggerFilter
    if compare_uncleaned :
        legend.AddEntry(h_fit_energy_3, "Low-E Events")         # just events below ECut MeV in original unprocessed data 
    legend.AddEntry(gaussFit0, "Gaussian Fit")
    legend.SetLineWidth(0)
    legend.Draw("same")

    latex = ROOT.TLatex()
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.02)
    xloc = 0.74
    yloc = 0.6
    latex.DrawText(xloc-0.01, yloc+0.04, "Re-Trigger Filtered Data") 
    latex.SetTextSize(0.015)
    latex.DrawText(xloc, yloc, "%i Entries"%(h_fit_energy_0.GetEntries()))
    latex.DrawLatex(xloc, yloc-0.04, "Mean = %.3f #pm %.4f MeV"%(fit0[2], fit0[3]))
    latex.DrawLatex(xloc, yloc-0.08, "#sigma = %.3f #pm %.4f"%(fit0[4], fit0[5]))
    latex.DrawLatex(xloc, yloc-0.12, "#Chi^{2}/ndf = %.1f/%i"%(fit0[0], fit0[1]))
    #if compare_uncleaned :
    #    latex.DrawText(xloc-0.5, yloc+0.04, "%i-%i=%i fewer entries below %.1f MeV"%(LowECountUnfiltered,Counts[4], LowECountUnfiltered-Counts[4], ECut))
    #    latex.DrawText(xloc-0.5, yloc-0.04, "%.2f percent reduction of entries below %.1f MeV"%((100*(LowECountUnfiltered-Counts[4])/LowECountUnfiltered, ECut)))
    #    latex.DrawText(xloc-0.5, yloc, "%i-%i=%i fewer events in total."%(OriginalCountValidFit, Counts[3], OriginalCountValidFit-Counts[3]))

    canvas.Print(plotTitle)

    if compare_uncleaned :
        print "h_fit_energy_0.GetEntries() = ", h_fit_energy_0.GetEntries()
        print "Original count = ", OriginalCount
        print "Original count valid fit = ", OriginalCountValidFit
        print "(original count below ECut): LowECountUnfiltered = ", LowECountUnfiltered
    print "CountTotal = ", Counts[0]
    print "CountTrigger", Counts[1]
    print "CountRadius = ", Counts[2]
    print "CountValid =", Counts[3]
    print "LowECountFiltered = ", Counts[4]
    
    return canvas
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def plot_dual_multiple_fit_E(h_energy_MC, h_energy_data, ERange, plotTitle):
    
    ''' 
    Description:
    2 pads, 2 peaks in each, fit both peaks in each canvas with Gaussians
    Save as pdf
    
    Parameters:
    h_energy_MC : ROOT TH1D
        Histogram containing MC data
    h_energy_data : ROOT TH1D
        Histogram containing data
    ERange : list of lists
        structured like [[total energy range], [range of 1st MC peak], [range of 2nd MC peak], [range of 1st data peak], [range of 2nd data peak]]
    plotTitle : str
        What to title the output file

    Returns:
    canvas : ROOT TCanvas

    '''

    nbins = 100                                                   # for histogram formatting
    binsize = 1000*(ERange[0][1]-ERange[0][0])/nbins              # [keV]

    h_energy_MC.SetDirectory(0)
    h_energy_MC.SetStats(0)
    h_energy_data.SetDirectory(0)
    h_energy_data.SetStats(0)
    
    if not h_energy_MC or not h_energy_data :
        print "Failed to get data histogram"
        sys.exit(1)

    canvas = ROOT.TCanvas("canvas")     
    canvas.Divide(1,2)  

    # ===============================================================================================                                                  
    canvas.cd(1)

    max_bin_value_MC = h_energy_MC.GetMaximum()

    h_energy_MC.SetLineColor(ROOT.kBlue)
    h_energy_MC.GetYaxis().SetTitleOffset(1.1)
    h_energy_MC.GetYaxis().SetTitle("Count per %.1f keV bin"%(binsize)) 
    h_energy_MC.GetXaxis().SetTitle("Fitted energy [MeV]")
    h_energy_MC.GetYaxis().SetRangeUser(0, max_bin_value_MC + 0.1*max_bin_value_MC)

    gaussFit0MC = ROOT.TF1("gaussFit0MC", "gaus", ERange[1][0], ERange[1][1]) # temp fit with built in range
    h_energy_MC.Fit(gaussFit0MC, "ERL") # fit Po peak                
    fitMCtemppeakPo = []
    fitMCtemppeakPo.append(gaussFit0MC.GetParameter(1))                                    # mean 
    fitMCtemppeakPo.append(gaussFit0MC.GetParameter(2))                                    # stdev

    gaussFit1MC = ROOT.TF1("gaussFit1MC", "gaus", ERange[2][0], ERange[2][1]) # temp fit with built in range
    h_energy_MC.Fit(gaussFit1MC, "ERL") # fit Bi peak                
    fitMCtemppeakBi = []
    fitMCtemppeakBi.append(gaussFit1MC.GetParameter(1))                                    # mean 
    fitMCtemppeakBi.append(gaussFit1MC.GetParameter(2))                                    # stdev


    gaussFit0MC = ROOT.TF1("gaussFit0MC", "gaus", fitMCtemppeakPo[0]-2.0*fitMCtemppeakPo[1], fitMCtemppeakPo[0]+2.0*fitMCtemppeakPo[1]) # actual gaussian fit
    h_energy_MC.Fit(gaussFit0MC, "ERL")   
    fitMCtemppeakPo = []
    fitMCtemppeakPo.append(gaussFit0MC.GetParameter(1))                                    # mean 
    fitMCtemppeakPo.append(gaussFit0MC.GetParameter(2))                                    # stdev
 
    gaussFit1MC = ROOT.TF1("gaussFit1MC", "gaus", fitMCtemppeakBi[0]-2.0*fitMCtemppeakBi[1], fitMCtemppeakBi[0]+2.0*fitMCtemppeakBi[1]) # actual gaussian fit
    h_energy_MC.Fit(gaussFit1MC, "ERL")                
    fitMCtemppeakBi = []
    fitMCtemppeakBi.append(gaussFit1MC.GetParameter(1))                                    # mean 
    fitMCtemppeakBi.append(gaussFit1MC.GetParameter(2))                                    # stdev


    # last time, symmetrical 
    gaussFit0MC = ROOT.TF1("gaussFit0MC", "gaus", fitMCtemppeakPo[0]-2.0*fitMCtemppeakPo[1], fitMCtemppeakPo[0]+2.0*fitMCtemppeakPo[1]) # actual gaussian fit
    h_energy_MC.Fit(gaussFit0MC, "ERL")    
 
    gaussFit1MC = ROOT.TF1("gaussFit1MC", "gaus", fitMCtemppeakBi[0]-2.0*fitMCtemppeakBi[1], fitMCtemppeakBi[0]+2.0*fitMCtemppeakBi[1]) # actual gaussian fit
    h_energy_MC.Fit(gaussFit1MC, "ERL+")  
    
    fit0 = [] 
    fit0.append(gaussFit0MC.GetChisquare())  # chi20       0
    fit0.append(gaussFit0MC.GetNDF())        # ndf0        1
    fit0.append(gaussFit0MC.GetParameter(1)) # mean0       2
    fit0.append(gaussFit0MC.GetParError(1))  # mean0error  3
    fit0.append(gaussFit0MC.GetParameter(2)) # sigma0      4
    fit0.append(gaussFit0MC.GetParError(2))  # sigma0error 5
    fit1 = [] 
    fit1.append(gaussFit1MC.GetChisquare())  # chi20       0
    fit1.append(gaussFit1MC.GetNDF())        # ndf0        1
    fit1.append(gaussFit1MC.GetParameter(1)) # mean0       2
    fit1.append(gaussFit1MC.GetParError(1))  # mean0error  3
    fit1.append(gaussFit1MC.GetParameter(2)) # sigma0      4
    fit1.append(gaussFit1MC.GetParError(2))  # sigma0error 5
    

    h_energy_MC.Draw()

    legendMC = ROOT.TLegend(0.73, 0.73, 0.9, 0.9) 
    legendMC.AddEntry(h_energy_MC, "Fitted BiPo MC  ")        
    legendMC.AddEntry(gaussFit0MC, "Gaussian Fit")
    legendMC.SetLineWidth(0)
    legendMC.Draw("same")

    latex = ROOT.TLatex()
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.04)
    xloc0 = 0.74
    yloc0 = 0.6
    yloc1 = 0.38

    
    latex.DrawText(xloc0-0.01, yloc0+0.04, "Peak with Mean %.2f MeV"%(fit0[2])) 
    latex.SetTextSize(0.03)
    latex.DrawText(xloc0, yloc0, "n = %i (Entries)"%(h_energy_MC.GetEntries()))
    latex.DrawLatex(xloc0, yloc0-0.04, "Mean = %.3f #pm %.4f MeV"%(fit0[2], fit0[3]))
    latex.DrawLatex(xloc0, yloc0-0.08, "#sigma = %.3f #pm %.4f"%(fit0[4], fit0[5]))
    latex.DrawLatex(xloc0, yloc0-0.12, "#Chi^{2}n/ndf = %.1f/%i"%(fit0[0]*h_energy_MC.GetEntries(), fit0[1]))

    latex.SetTextSize(0.04)
    latex.DrawText(xloc0-0.01, yloc1, "Peak with Mean %.2f MeV"%(fit1[2])) 
    latex.SetTextSize(0.03)
    latex.DrawLatex(xloc0, yloc1-0.04, "Mean = %.3f #pm %.4f MeV"%(fit1[2], fit1[3]))
    latex.DrawLatex(xloc0, yloc1-0.08, "#sigma = %.3f #pm %.4f"%(fit1[4], fit1[5]))
    latex.DrawLatex(xloc0, yloc1-0.12, "#Chi^{2}n/ndf = %.1f/%i"%(fit1[0]*h_energy_MC.GetEntries(), fit1[1]))
    
    # ===============================================================================================     
    canvas.cd(2)

    max_bin_value_data = h_energy_data.GetMaximum()

    h_energy_data.SetLineColor(ROOT.kGreen)
    h_energy_data.GetYaxis().SetTitleOffset(1.1)
    h_energy_data.GetYaxis().SetTitle("Count per %.1f keV bin"%(binsize)) 
    h_energy_data.GetXaxis().SetTitle("Fitted energy [MeV]")
    h_energy_data.GetYaxis().SetRangeUser(0, max_bin_value_data + 0.1*max_bin_value_data) 

    
    gaussFit0data = ROOT.TF1("gaussFit0data", "gaus", ERange[3][0], ERange[3][1]) # temp fit with hardcoded range
    h_energy_data.Fit(gaussFit0data, "ERL") # fit Po peak                
    fitDatatemppeakPo = []
    fitDatatemppeakPo.append(gaussFit0data.GetParameter(1))                                    # mean 
    fitDatatemppeakPo.append(gaussFit0data.GetParameter(2))                                    # stdev

    gaussFit1data = ROOT.TF1("gaussFit1data", "gaus", ERange[4][0], ERange[4][1]) # temp fit with hard coded range
    h_energy_data.Fit(gaussFit1data, "ERL") # fit Po peak                
    fitDatatemppeakBi = []
    fitDatatemppeakBi.append(gaussFit1data.GetParameter(1))                                    # mean 
    fitDatatemppeakBi.append(gaussFit1data.GetParameter(2))                                    # stdev

    gaussFit0data = ROOT.TF1("gaussFit0data", "gaus", fitDatatemppeakPo[0]-2.0*fitDatatemppeakPo[1], fitDatatemppeakPo[0]+2.0*fitDatatemppeakPo[1]) # actual gaussian fit
    h_energy_data.Fit(gaussFit0data, "ERL")  
    fitDatatemppeakPo = []
    fitDatatemppeakPo.append(gaussFit0data.GetParameter(1))                                    # mean 
    fitDatatemppeakPo.append(gaussFit0data.GetParameter(2))                                    # stdev  
 
    gaussFit1data = ROOT.TF1("gaussFit1data", "gaus", fitDatatemppeakBi[0]-2.0*fitDatatemppeakBi[1], fitDatatemppeakBi[0]+2.0*fitDatatemppeakBi[1]) # actual gaussian fit
    h_energy_data.Fit(gaussFit1data, "ERL")     
    fitDatatemppeakBi = []
    fitDatatemppeakBi.append(gaussFit1data.GetParameter(1))                                    # mean 
    fitDatatemppeakBi.append(gaussFit1data.GetParameter(2))                                    # stdev


    gaussFit0data = ROOT.TF1("gaussFit0data", "gaus", fitDatatemppeakPo[0]-2.0*fitDatatemppeakPo[1], fitDatatemppeakPo[0]+2.0*fitDatatemppeakPo[1]) # actual gaussian fit
    h_energy_data.Fit(gaussFit0data, "ERL")    
 
    gaussFit1data = ROOT.TF1("gaussFit1data", "gaus", fitDatatemppeakBi[0]-2.0*fitDatatemppeakBi[1], fitDatatemppeakBi[0]+2.0*fitDatatemppeakBi[1]) # actual gaussian fit
    h_energy_data.Fit(gaussFit1data, "ERL+")     

    
    fit0D = [] 
    fit0D.append(gaussFit0data.GetChisquare())  # chi20       0
    fit0D.append(gaussFit0data.GetNDF())        # ndf0        1
    fit0D.append(gaussFit0data.GetParameter(1)) # mean0       2
    fit0D.append(gaussFit0data.GetParError(1))  # mean0error  3
    fit0D.append(gaussFit0data.GetParameter(2)) # sigma0      4
    fit0D.append(gaussFit0data.GetParError(2))  # sigma0error 5
    fit1D = [] 
    fit1D.append(gaussFit1data.GetChisquare())  # chi20       0
    fit1D.append(gaussFit1data.GetNDF())        # ndf0        1
    fit1D.append(gaussFit1data.GetParameter(1)) # mean0       2
    fit1D.append(gaussFit1data.GetParError(1))  # mean0error  3
    fit1D.append(gaussFit1data.GetParameter(2)) # sigma0      4
    fit1D.append(gaussFit1data.GetParError(2))  # sigma0error 5
    

    h_energy_data.Draw()

    legend = ROOT.TLegend(0.73, 0.73, 0.9, 0.9) 
    legend.AddEntry(h_energy_data, "Fitted BiPo Data")         
    legend.AddEntry(gaussFit0data, "Gaussian Fit")
    legend.SetLineWidth(0)
    legend.Draw("same")

    latex = ROOT.TLatex()
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.04)
    xloc0 = 0.74
    yloc0 = 0.6
    yloc1 = 0.38
    
    latex.DrawText(xloc0-0.01, yloc0+0.04, "Peak with Mean %.2f MeV"%(fit0D[2])) 
    latex.SetTextSize(0.03)
    latex.DrawText(xloc0, yloc0, "n = %i (Entries)"%(h_energy_data.GetEntries()))
    latex.DrawLatex(xloc0, yloc0-0.04, "Mean = %.3f #pm %.4f MeV"%(fit0D[2], fit0D[3]))
    latex.DrawLatex(xloc0, yloc0-0.08, "#sigma = %.3f #pm %.4f"%(fit0D[4], fit0D[5]))
    latex.DrawLatex(xloc0, yloc0-0.12, "#Chi^{2}n/ndf = %.1f/%i"%(fit0D[0]*h_energy_data.GetEntries(), fit0D[1]))

    latex.SetTextSize(0.04)
    latex.DrawText(xloc0-0.01, yloc1, "Peak with Mean %.2f MeV"%(fit1D[2])) 
    latex.SetTextSize(0.03)
    latex.DrawLatex(xloc0, yloc1-0.04, "Mean = %.3f #pm %.4f MeV"%(fit1D[2], fit1D[3]))
    latex.DrawLatex(xloc0, yloc1-0.08, "#sigma = %.3f #pm %.4f"%(fit1D[4], fit1D[5]))
    latex.DrawLatex(xloc0, yloc1-0.12, "#Chi^{2}n/ndf = %.1f/%i"%(fit1D[0]*h_energy_data.GetEntries(), fit1D[1]))
    
    # =============================================================================================== 


    canvas.Print(plotTitle)


    return canvas
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def plot_dual_energy(h_energy_MC, h_energy_data, ERange, plotTitle):
    
    """ 

    plot MC and scaled BiPo data overlapping, BiPo with points

    """

    if not h_energy_MC or not h_energy_data :
        print "Failed to get data histogram"
        sys.exit(1)

    nbins = 100                                                   # for histogram formatting
    binsize = 1000*(ERange[0][1]-ERange[0][0])/nbins              # [keV]

    h_energy_MC.SetDirectory(0)
    h_energy_MC.SetStats(0)
    h_energy_data.SetDirectory(0)
    h_energy_data.SetStats(0)

    canvas = ROOT.TCanvas("canvas")
    canvas.cd()

    h_energy_MC.SetLineColor(ROOT.kBlue)

    h_energy_data.SetMarkerStyle(2)
    h_energy_data.GetYaxis().SetTitleOffset(1.1)
    h_energy_data.GetYaxis().SetTitle("Normalized count per %.1f keV bin"%(binsize)) 
    h_energy_data.GetXaxis().SetTitle("Fitted energy [MeV]")

    h_energy_MC.Draw()
    h_energy_data.Draw("SAME P")

    legend = ROOT.TLegend(0.75, 0.79, 0.9, 0.9) 
    legend.AddEntry(h_energy_MC, "MC Energy") 
    legend.AddEntry(h_energy_data, "Data Energy, Ratio Scaled") 
    legend.SetLineWidth(0)
    legend.Draw("same")

    xloc0 = 0.44
    latex = ROOT.TLatex()
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.02)
    latex.DrawText(xloc0-0.01, 0.68, "Data Scaled By 0.7369")#Ratio Function (2nd Order Polynomial)") 
    

    canvas.Print(outDir + plotTitle)

    return canvas
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def plot_dual_fit_energy(input_files, plotTitles):
    
    """ 

    Plot the fitted energy of two root files on same plot.
    Automatically determines a range to apply gaussian fit to.

    Parameters 
    input_files : list
        Path to the RAT DS file(s) to plot.
    plotTitles : list
        What to label the histograms.

    Returns - The histogram plot as a pdf

    """

    file_name_0 = input_files[0]
    file_name_1 = input_files[1]
    plot0Title = plotTitle[0]
    plot1Title = plotTitle[1]

    Erange = [0.0, 1.5] # range to plot
    nbins = 100 
    binsize = 1000*(Erange[1]-Erange[0])/nbins # keV

    event_count_0 = 0
    zeroptwo = 0

    h_fit_energy_0 = ROOT.TH1D(plot0Title, plot0Title, nbins, Erange[0], Erange[1])
    for ds, run in dsreader(file_name_0):
        for iev in range(0, ds.GetEVCount()):
            ev = ds.GetEV(iev)
            if not ev.DefaultFitVertexExists() or not ev.GetDefaultFitVertex().ContainsEnergy() or not ev.GetDefaultFitVertex().ValidEnergy():
                continue  # if didn't fit succesfully
            h_fit_energy_0.Fill(ev.GetDefaultFitVertex().GetEnergy())
            event_count_0 += 1
            if ev.GetDefaultFitVertex().GetEnergy() < 0.2 :
                zeroptwo += 1
    
    event_count_1 = 0
    h_fit_energy_1 = ROOT.TH1D(plot1Title, plot1Title, nbins, Erange[0], Erange[1])
    for ds, run in dsreader(file_name_1):
        for iev in range(0, ds.GetEVCount()):
            ev = ds.GetEV(iev)
            if not ev.DefaultFitVertexExists() or not ev.GetDefaultFitVertex().ContainsEnergy() or not ev.GetDefaultFitVertex().ValidEnergy():
                continue  # if didn't fit succesfully
            h_fit_energy_1.Fill(ev.GetDefaultFitVertex().GetEnergy())
            event_count_1 += 1
            
    h_fit_energy_0.SetDirectory(0)
    h_fit_energy_1.SetDirectory(0)
    h_fit_energy_0.SetStats(0)
    h_fit_energy_1.SetStats(0)

    canvas = ROOT.TCanvas("canvas")
    canvas.cd()

    h_fit_energy_0.SetLineColor(ROOT.kBlue)

    h_fit_energy_1.SetLineColor(ROOT.kRed)
    h_fit_energy_1.GetYaxis().SetTitleOffset(1.1)
    h_fit_energy_1.GetYaxis().SetTitle("Count per %.1f keV bin"%(binsize)) 
    h_fit_energy_1.GetXaxis().SetTitle("Fitted energy [MeV]")

    # do this real quick first, find mean and sigma, then choose new range for fit from ~ 5 sigma on either side of mean
    gaussFit00 = ROOT.TF1("gaussfit", "gaus", Erange[0], Erange[1]) 
    h_fit_energy_0.Fit(gaussFit0, "ERL")                    
    gaussFit11 = ROOT.TF1("gaussfit", "gaus", Erange[0], Erange[1])
    h_fit_energy_1.Fit("gaussfit", "ERL")

    fit00 = [] # just for computing range to apply real fit to
    fit00.append(h_fit_energy_0.GetMean()) # mean 
    fit00.append(gaussFit00.GetParameter(2)) # stdev
    fit01 = []
    fit01.append(h_fit_energy_0.GetMean()) # mean 
    fit01.append(gaussFit01.GetParameter(2)) # stdev

    gaussFit0 = ROOT.TF1("gaussfit", "gaus", fit00[0]-5*fit00[1], fit00[0]+5*fit00[1]) 
    h_fit_energy_0.Fit(gaussFit0, "ERL")                    
    gaussFit1 = ROOT.TF1("gaussfit", "gaus", fit01[0]-5*fit01[1], fit01[0]+5*fit01[1])
    h_fit_energy_1.Fit("gaussfit", "ERL")
    
    fit0 = [] # real fit parameters now
    fit0.append(gaussFit0.GetChisquare())  # chi20       0
    fit0.append(gaussFit0.GetNDF())        # ndf0        1
    fit0.append(gaussFit0.GetParameter(1)) # mean0       2
    fit0.append(gaussFit0.GetParError(1))  # mean0error  3
    fit0.append(gaussFit0.GetParameter(2)) # sigma0      4
    fit0.append(gaussFit0.GetParError(2))  # sigma0error 5
    fit1 = [] 
    fit1.append(gaussFit1.GetChisquare())  # chi21       0
    fit1.append(gaussFit1.GetNDF())        # ndf1        1
    fit1.append(gaussFit1.GetParameter(1)) # mean1       2
    fit1.append(gaussFit1.GetParError(1))  # mean1error  3
    fit1.append(gaussFit1.GetParameter(2)) # sigma1      4
    fit1.append(gaussFit1.GetParError(2))  # sigma1error 5

    h_fit_energy_1.Draw()
    h_fit_energy_0.Draw("same")

    legend = ROOT.TLegend(0.75, 0.79, 0.9, 0.9) 
    legend.AddEntry(h_fit_energy_0, plot0Title) 
    legend.AddEntry(gaussFit0, "Gaussian Fit")
    legend.AddEntry(h_fit_energy_1, plot1Title) 
    legend.AddEntry(gaussFit1, "Gaussian Fit")
    legend.SetLineWidth(0)
    legend.Draw("same")

    latex = ROOT.TLatex()
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.02)
    latex.DrawText(0.69, 0.48, plot0Title) 
    latex.DrawText(0.44, 0.68, plot1Title) 
    latex.SetTextSize(0.015)
    latex.DrawText(0.7, 0.45, "%i Entries"%(event_count_0))
    latex.DrawText(0.45, 0.65, "%i Entries"%(event_count_1))
    latex.DrawText(0.7, 0.41, "Mean = %.3f +/- %.4f MeV"%(fit0[2], fit0[3]))
    latex.DrawText(0.7, 0.37, "stdev = %.3f +/- %.4f"%(fit0[4], fit0[5]))
    latex.DrawText(0.7, 0.33, "chi^2/ndf = %.1f/%i"%(fit0[0], fit0[1]))
    latex.DrawText(0.45, 0.61, "Mean = %.3f +/- %.4f MeV"%(fit1[2], fit1[3]))
    latex.DrawText(0.45, 0.57, "stdev = %.3f +/- %.4f"%(fit1[4], fit1[5]))
    latex.DrawText(0.45, 0.53, "chi^2/ndf = %.1f/%i"%(fit1[0], fit1[1]))

    canvas.Print(outDir + outFilename)

    return canvas
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////