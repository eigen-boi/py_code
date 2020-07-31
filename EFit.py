# Aidan Patton, June 2020
# Based on code by S. Langrock (plot_fit.py)
# Also heavily influenced by 
# https://indico.cern.ch/event/704163/contributions/2936719/attachments/1693833/2726445/Tutorial-PyROOT.pdf
# and Energy.C by Yang Zhang

import ROOT
from rat import dsreader
import sys

#outDir = "/SNO+/pdfs/" # for local VM
outDir = "/home/eigenboi/pdfs/" # for cedar
outFileName = "example.pdf"

index = 1



def BiPoEnergy():

    """ For finding energy of BiPo data """

    input_files = ["uhhh.root"]
    ERange = [0, 10]
    ECut = 0.015 # already baked in to QPartialEnergy method
    hist_display_title = "BiPo214 Energy"
    plotTitle = outDir + outFileName
    compare_uncleaned = True

    h_energy, Counts = retriggerFilter(input_files, ERange, ECut, hist_display_title) 
    plot_fit_energy(input_files, h_energy, compare_uncleaned, Counts, ERange, ECut, plotTitle)



def retriggerFilter(input_files, ERange, ECut, hist_display_title): 

    """

    Description:
    Extract data from root file(s), apply retrigger filtering, and return histogram "h_energy_name" as root file "hist_energy.root"
    # contains commented out code to save histogram "h_energy_name" as root file "hist_energy.root"

    Parameters:
    fname - name of rat root file to process

    Returns:
    h_energy - the TH1D energy histogram
    Counts - CountTotal, CountTrigger, CountValid, CountRadius, LowECountFiltered
    ERange - energy range of h_energy_name [MeV]
    ECut - arbitrary cutoff below which events are considered "low energy" [MeV]

    """

    fitName = "partialFitter"

    #ERange = [0.0, 1.5]             # ERange and ECut as outputs so can use in next script too 
    #ERange = [0, 9]
    #ECut = 0.2                      # cutoff below which events are considered "low energy" [MeV]

    nbins = 100
    Counts = [0, 0, 0, 0, 0]        # CountTotal, CountTrigger, CountRadius, CountValid, LowECountFiltered

    #h_energy = ROOT.TH1D("h_energy_name","Retrigger Event Filter Test", nbins, ERange[0], ERange[1])
    h_energy = ROOT.TH1D("h_energy_name", hist_display_title, nbins, ERange[0], ERange[1])
    h_energy.SetDirectory(0)
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

                
                # radius cut, replaced below with inside AV cut
                if REV > 6000 :                  # Energy.C used 5500 mm
                    continue                                             # if outside AV, don't use
                Counts[2] += 1                                           # inside AV count
                
                '''
                # see if inside the AV
                if PosEV.Z() >= 6000 and RhoEV >= 730 :             
                    continue
                if PosEV.Z() < 6000 and REV > 6000 :
                    continue
                Counts[2] += 1                                           # inside AV count
                '''
                '''  
                if REV > 5500 : # FV cut instead will give better fit, boundary events are tricky 
                    continue
                Counts[2] += 1   
                '''  

                if not fVertex.ValidEnergy():
                    continue
                Counts[3] += 1                                           # how many actually pass the gauntlet

                h_energy.Fill(fVertex.GetEnergy()) 

                if fVertex.GetEnergy() < ECut: 
                    Counts[4] += 1                                       # keep track of LowE which pass the filter 

                #nHits = ev.GetNhits()
                #totcharge = ev.GetTotalCharge()

    # save file
    '''
    fileout =  ROOT.TFile("/SNO+/root_files/hist_energy_point_00.root", "RECREATE") 
    fileout.cd()
    h_energy.Write() 
    fileout.Close()
    '''
    return h_energy, Counts



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

    #histFile = ROOT.TFile.Open(file_name_0, "READ")               # read in filtered histogram
    #h_fit_energy_0 = histFile.Get("h_energy_name")
    h_fit_energy_0 = h_energy
    h_fit_energy_0.SetDirectory(0)
    h_fit_energy_0.SetStats(0)
    #histFile.Close()
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
    latex.DrawText(xloc, yloc-0.04, "Mean = %.3f +/- %.4f MeV"%(fit0[2], fit0[3]))
    latex.DrawText(xloc, yloc-0.08, "stdev = %.3f +/- %.4f"%(fit0[4], fit0[5]))
    latex.DrawText(xloc, yloc-0.12, "chi^2/ndf = %.1f/%i"%(fit0[0], fit0[1]))
    #if compare_uncleaned :
    #    latex.DrawText(xloc-0.5, yloc+0.04, "%i-%i=%i fewer entries below %.1f MeV"%(LowECountUnfiltered,Counts[4], LowECountUnfiltered-Counts[4], ECut))
    #    latex.DrawText(xloc-0.5, yloc-0.04, "%.2f percent reduction of entries below %.1f MeV"%((100*(LowECountUnfiltered-Counts[4])/LowECountUnfiltered, ECut)))
    #    latex.DrawText(xloc-0.5, yloc, "%i-%i=%i fewer events in total."%(OriginalCountValidFit, Counts[3], OriginalCountValidFit-Counts[3]))

    canvas.Print(plotTitle)
    #canvas.Print("/home/eigenboi/pdfs/tablecheck_fill_10MeV_00.pdf")
    #canvas.Print("/home/eigenboi/pdfs/tablecheck_point_2p5MeV_02.pdf")
    #canvas.Print("/home/eigenboi/pdfs/tablecheck_point_10MeV_AV_03.pdf")

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



def plot_dual_fit_energy(input_files, plotTitles):
    
    """ 

    Plot the fitted energy of two root files on same plot.
    Automatically determines a range to apply gaussian fit to.

    Parameters 
    input_files : list
        Path to the RAT DS file(s) to plot.
    plotTitles : list
        What to label the top and bottom histograms (respectively).

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



def fileArr(pickedfile):

    #cedar
    
    input_files = [['/home/eigenboi/root_files/sim_point_1MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_05.root'],
    ['/home/eigenboi/root_files/sim_point_2p5MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_05.root'],
    ['/home/eigenboi/root_files/sim_point_5MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_05.root'],
    ['/home/eigenboi/root_files/sim_point_7p5MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_05.root'],
    ['/home/eigenboi/root_files/sim_point_10MeV_AV_01.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_02.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_03.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_04.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_05.root']]
    

    #input_files = ['/home/eigenboi/root_files/sim_fill_1MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_fill_1MeV_02_tablecheck.root', '/home/eigenboi/root_files/sim_fill_1MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_fill_1MeV_04_tablecheck.root', '/home/eigenboi/root_files/sim_fill_1MeV_05_tablecheck.root']
    #input_files = ['/home/eigenboi/root_files/sim_fill_1MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_fill_1MeV_02_tablecheck.root', '/home/eigenboi/root_files/sim_fill_1MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_fill_1MeV_05_tablecheck.root']   
    
    input_files_C = [['/home/eigenboi/root_files/sim_point_1MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_1MeV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_1MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_1MeV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_1MeV_05_tablecheck.root'],
    ['/home/eigenboi/root_files/sim_point_2p5MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_2p5MeV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_2p5MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_2p5MeV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_2p5MeV_05_tablecheck.root'],
    ['/home/eigenboi/root_files/sim_point_5MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_5MeV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_5MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_5MeV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_5MeV_05_tablecheck.root'],
    ['/home/eigenboi/root_files/sim_point_7p5MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_7p5MeV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_7p5MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_7p5MeV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_7p5MeV_05_tablecheck.root'],
    ['/home/eigenboi/root_files/sim_point_10MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_05_tablecheck.root']]
    
    input_files_AV = [['/home/eigenboi/root_files/sim_point_1MeV_AV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_1MeV_AV_05_tablecheck.root'],
    ['/home/eigenboi/root_files/sim_point_2p5MeV_AV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_2p5MeV_AV_05_tablecheck.root'],
    ['/home/eigenboi/root_files/sim_point_5MeV_AV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_5MeV_AV_05_tablecheck.root'],
    ['/home/eigenboi/root_files/sim_point_7p5MeV_AV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_7p5MeV_AV_05_tablecheck.root'],
    ['/home/eigenboi/root_files/sim_point_10MeV_AV_01_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_02_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_03_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_04_tablecheck.root', '/home/eigenboi/root_files/sim_point_10MeV_AV_05_tablecheck.root']]

    input_files_fill = [["empty"],
                        ['/home/eigenboi/root_files/sim_fill_2p5MeV_01.root', '/home/eigenboi/root_files/sim_fill_2p5MeV_02.root', '/home/eigenboi/root_files/sim_fill_2p5MeV_03.root', '/home/eigenboi/root_files/sim_fill_2p5MeV_04.root', '/home/eigenboi/root_files/sim_fill_2p5MeV_05.root'],
                        ['/home/eigenboi/root_files/sim_fill_5MeV_01.root', '/home/eigenboi/root_files/sim_fill_5MeV_02.root', '/home/eigenboi/root_files/sim_fill_5MeV_03.root', '/home/eigenboi/root_files/sim_fill_5MeV_04.root', '/home/eigenboi/root_files/sim_fill_5MeV_05.root'],
                        ['/home/eigenboi/root_files/sim_fill_7p5MeV_01.root', '/home/eigenboi/root_files/sim_fill_7p5MeV_02.root', '/home/eigenboi/root_files/sim_fill_7p5MeV_03.root', '/home/eigenboi/root_files/sim_fill_7p5MeV_04.root'],
                        ['/home/eigenboi/root_files/sim_fill_10MeV_03.root', '/home/eigenboi/root_files/sim_fill_10MeV_04.root', '/home/eigenboi/root_files/sim_fill_10MeV_05.root']]

    input_files_fill2 = [["empty"],
                        ['/home/eigenboi/root_files/sim_fill_2p5MeV_21.root', '/home/eigenboi/root_files/sim_fill_2p5MeV_22.root', '/home/eigenboi/root_files/sim_fill_2p5MeV_23.root', '/home/eigenboi/root_files/sim_fill_2p5MeV_24.root', '/home/eigenboi/root_files/sim_fill_2p5MeV_25.root'],
                        ['/home/eigenboi/root_files/sim_fill_5MeV_24.root', '/home/eigenboi/root_files/sim_fill_5MeV_25.root'],
                        ['/home/eigenboi/root_files/sim_fill_7p5MeV_22.root', '/home/eigenboi/root_files/sim_fill_7p5MeV_23.root', '/home/eigenboi/root_files/sim_fill_7p5MeV_24.root', '/home/eigenboi/root_files/sim_fill_7p5MeV_25.root'],
                        ['/home/eigenboi/root_files/sim_fill_10MeV_21.root', '/home/eigenboi/root_files/sim_fill_10MeV_22.root', '/home/eigenboi/root_files/sim_fill_10MeV_23.root', '/home/eigenboi/root_files/sim_fill_10MeV_24.root', '/home/eigenboi/root_files/sim_fill_10MeV_25.root']]


    input_files_fill_T = [["empty"],
                        ['/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_01_tablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_02_tablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_03_tablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_04_tablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_05_tablecheck.root'],
                        ['/home/eigenboi/root_files/sim_fill_5MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_fill_5MeV_02_tablecheck.root', '/home/eigenboi/root_files/sim_fill_5MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_fill_5MeV_04_tablecheck.root', '/home/eigenboi/root_files/sim_fill_5MeV_05_tablecheck.root'],
                        ['/home/eigenboi/root_files/sim_fill_7p5MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_fill_7p5MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_fill_7p5MeV_04_tablecheck.root', '/home/eigenboi/root_files/sim_fill_7p5MeV_05_tablecheck.root'],
                        ['/home/eigenboi/root_files/sim_fill_10MeV_01_tablecheck.root', '/home/eigenboi/root_files/sim_fill_10MeV_02_tablecheck.root', '/home/eigenboi/root_files/sim_fill_10MeV_03_tablecheck.root', '/home/eigenboi/root_files/sim_fill_10MeV_04_tablecheck.root', '/home/eigenboi/root_files/sim_fill_10MeV_05_tablecheck.root']]

    input_files_fill_T2 = [["N/A"],
                        ['/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_01_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_02_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_03_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_04_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_2p5MeV_05_newtablecheck.root'],
                        ['/home/eigenboi/scratch/root_files/sim_fill_5MeV_01_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_5MeV_02_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_5MeV_03_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_5MeV_04_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_5MeV_05_newtablecheck.root'],
                        ['/home/eigenboi/scratch/root_files/sim_fill_7p5MeV_01_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_7p5MeV_02_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_7p5MeV_03_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_7p5MeV_04_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_7p5MeV_05_newtablecheck.root'],
                        ['/home/eigenboi/scratch/root_files/sim_fill_10MeV_01_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_10MeV_02_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_10MeV_03_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_10MeV_04_newtablecheck.root', '/home/eigenboi/scratch/root_files/sim_fill_10MeV_05_newtablecheck.root']]

    filedict = {0:input_files, 1:input_files_C, 2:input_files_AV, 3:input_files_fill, 4:input_files_fill2, 5:input_files_fill_T, 6:input_files_fill_T2}

    return filedict[pickedfile]



def CompERanges():

    """ Useful for comparing the energy of a single configuaration of RAT over an energy range 0f 1MeV, 2.5MeV, 5MeV, 7.5MeV, 10MeV """

    input_files = fileArr(6)[index] # just change index to get each energy of a series of runs
    ERange = [[0.0, 1.5], [0, 4], [0, 8], [0, 11], [0, 14]]
    ERange = ERange[index]
    ECut = 0.2

    hist_display_title = ["1 MeV Fill Dist. Energy Fit", "2.5 MeV Fill Dist. Energy Fit", "5 MeV Fill Dist. Energy Fit", "7.5 MeV Fill Dist. Energy Fit", "10 MeV Fill Dist. Energy Fit"]
    #hist_display_title = ["1 MeV Point Dist. Energy Fit", "2.5 MeV Point Dist. Energy Fit, (rho=0m)", "5 MeV Point Dist. Energy Fit", "7.5 MeV Point Dist. Energy Fit", "10 MeV Point Dist. Energy Fit"]

    #plotTitle = ["/home/eigenboi/pdfs/fill_1MeV_01.pdf",  "/home/eigenboi/pdfs/fill_2p5MeV_03.pdf", "/home/eigenboi/pdfs/fill_5MeV_03.pdf", "/home/eigenboi/pdfs/fill_7p5MeV_03.pdf",  "/home/eigenboi/pdfs/fill_10MeV_03.pdf"]
    plotTitle = [outDir+"fill_1MeV_04.pdf",  outDir+"fill_2p5MeV_05.pdf", outDir+"fill_5MeV_04.pdf", outDir+"fill_7p5MeV_04.pdf", outDir+"fill_10MeV_04.pdf"]

    compare_uncleaned = True
    h_energy, Counts = retriggerFilter(input_files, ERange, ECut, hist_display_title[index]) 
    plot_fit_energy(input_files, h_energy, compare_uncleaned, Counts, ERange, ECut, plotTitle[index])



if __name__ == '__main__':

    pass