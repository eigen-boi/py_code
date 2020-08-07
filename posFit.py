
# Aidan Patton, June 2020
# Based on code by S. Langrock (plot_fit.py)
# Also heavily influenced by 
# https://indico.cern.ch/event/704163/contributions/2936719/attachments/1693833/2726445/Tutorial-PyROOT.pdf

import ROOT
from rat import dsreader
import sys
import math
from array import array

#outDir = "/SNO+/pdfs/" # for local VM
outDir = "/home/eigenboi/pdfs/" # for cedar
outFileName = "zrho_data_TH2D.pdf"
ECut = 0


def ZRHO_in_AV(file_name):
    """ 
    Description:
    Plot the Z vs Rho for MC and Fitted events THAT ARE INSIDE THE ACRYLIC ITSELF
    See where MC events are fitted, and destinguish between MC events that are fitted inside and outside the AV
    Compare valid energy and valid position fits 

    Parameters: 
    file_name - Path to the RAT DS file to play around with
    Return: nothing
    """

    fitName = "partialFitter"
    CountEV = 0
    just_compare_neck = False

    # ======================================================================= for z vs rho

    ArrayRho_EV = array("d")
    ArrayZ_EV = array("d")

    ArrayRho_EVO = array("d")
    ArrayZ_EVO = array("d")

    ArrayRho_MC = array("d")
    ArrayZ_MC = array("d")

    ArrayRho_MCO = array("d")
    ArrayZ_MCO = array("d")

    # ======================================================================= 

    AVCountMC = 0
    AVCountEV = 0
    OCountEV = 0
    OCountMC = 0
    filtercuts = [0, 0, 0, 0] # total fit count, retrigger filtered, posvalid, Evalid

    for ds, run in dsreader(file_name):

        for iev in range(0, ds.GetEVCount()): # ds is an entry, multiple EV events in one MC event possible so loop through them
            
            filtercuts[0] += 1
            if iev > 0: # retrigger filter
                continue

            ev = ds.GetEV(iev) # single EV event

            filtercuts[1] += 1

            if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                continue # valid position filter

            filtercuts[2] += 1

            fVertex = ev.GetFitResult(fitName).GetVertex(0)     
            PosEV = fVertex.GetPosition()                  # get position vector 
            REV = PosEV.Mag()                              # and radius
            RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

            mc = ds.GetMC() # individual MC event

            PosMC = mc.GetMCParticle(0).GetPosition() # get position vector and radius 
            RMC = PosMC.Mag()
            RhoMC = math.sqrt(PosMC.X()**2 + PosMC.Y()**2)

            if PosEV.Z() > 6000 and RhoEV > 730 and RhoEV < 785 :
                in_AV = True
            elif PosEV.Z() > 0 and RhoEV > 785 and REV > 6000 and REV < 6055 :
                in_AV = True  
            elif PosEV.Z() < 0 and REV > 6000 and REV < 6055 :
                in_AV = True  
            else :
                in_AV = False

            if PosEV.Z() >= 6055 and RhoEV >= 755 :
                out_AV = True
            elif PosEV.Z() < 6055 and REV > 6055 :
                out_AV = True
            else:
                out_AV = False

            if in_AV :

                # store values in arrays for plotting
                ArrayZ_EV.append(PosEV.Z()) 
                ArrayRho_EV.append(RhoEV)
                AVCountEV += 1  

                # store values in arrays for plotting
                ArrayZ_MC.append(PosMC.Z())
                ArrayRho_MC.append(RhoMC)
                AVCountMC += 1 # keep count of mc events

                if not fVertex.ValidEnergy():
                    continue          
                filtercuts[3] += 1  

            if out_AV :

                # store values in arrays for plotting
                ArrayZ_EVO.append(PosEV.Z()) 
                ArrayRho_EVO.append(RhoEV)
                OCountEV += 1  

                # store values in arrays for plotting
                ArrayZ_MCO.append(PosMC.Z())
                ArrayRho_MCO.append(RhoMC)
                OCountMC += 1

                if not fVertex.ValidEnergy():
                    continue           
                filtercuts[3] += 1  

    # =======================================================================
    
    canvas1 = ROOT.TCanvas("canvas1") # make canvas
    canvas1.cd()

    ZRhoPlotMC = ROOT.TGraph(AVCountMC, ArrayRho_MC, ArrayZ_MC) # make and format TGraph for Fitted events
    ZRhoPlotMC.SetTitle("MC Fitted Inside Acrylic")
    ZRhoPlotMC.SetMarkerStyle(6)
    ZRhoPlotMC.SetMarkerColorAlpha(ROOT.kRed, 1)
    ZRhoPlotMC.SetLineColorAlpha(ROOT.kRed, 1)
    ZRhoPlotMC.SetFillStyle(0)
    ZRhoPlotMC.Draw()

    ZRhoPlotEV = ROOT.TGraph(AVCountEV, ArrayRho_EV, ArrayZ_EV) # make and format TGraph for Fitted events
    ZRhoPlotEV.SetTitle("Fitted Events in Acrylic")
    ZRhoPlotEV.SetMarkerStyle(6)
    ZRhoPlotEV.SetMarkerColorAlpha(ROOT.kBlue, 1)
    ZRhoPlotEV.SetLineColorAlpha(ROOT.kBlue, 1)
    ZRhoPlotEV.SetFillStyle(0)
    ZRhoPlotEV.Draw()

    ZRhoPlotEVO = ROOT.TGraph(OCountEV, ArrayRho_EVO, ArrayZ_EVO) # make and format TGraph for validE stuff 
    ZRhoPlotEVO.SetTitle("Outside AV")
    ZRhoPlotEVO.SetMarkerStyle(6)
    ZRhoPlotEVO.SetMarkerColorAlpha(3, 1)
    ZRhoPlotEVO.SetLineColorAlpha(3, 1)
    ZRhoPlotEVO.SetFillStyle(0)
    ZRhoPlotEVO.Draw()

    ZRhoPlotMCO = ROOT.TGraph(OCountMC, ArrayRho_MCO, ArrayZ_MCO) # make and format TGraph for validE stuff 
    ZRhoPlotMCO.SetTitle("MC Fitted Outside AV")
    ZRhoPlotMCO.SetMarkerStyle(6)
    ZRhoPlotMCO.SetMarkerColorAlpha(6, 1)
    ZRhoPlotMCO.SetLineColorAlpha(6, 1)
    ZRhoPlotMCO.SetFillStyle(0)
    ZRhoPlotMCO.Draw()

    r = 6000 # save to file eventually
    theta = -math.pi/2
    x, y = array("d"), array("d")
    for i in range(100):
        x.append(r*math.cos(theta)), y.append(r*math.sin(theta))
        theta += 0.0302
    for i in range(50):
        x.append(730), y.append(6000+i*140)
    boundaryplot = ROOT.TGraph(150, x, y)
    boundaryplot.SetTitle("Inner AV")
    boundaryplot.SetLineStyle(3)
    boundaryplot.SetMarkerStyle(0)
    boundaryplot.SetFillStyle(0)
    boundaryplot.Draw()
    
    lineplot = ROOT.TGraph(2, array("d", [0, 8000]), array("d", [747.5, 747.5])) # a line at Z = 747.5 mm to tell where scint ends
    lineplot.SetTitle("747.5 mm")
    lineplot.SetMarkerStyle(0)
    lineplot.SetLineStyle(7)
    lineplot.SetFillStyle(0)
    lineplot.Draw()

    ZRho = ROOT.TMultiGraph() # draw the graphs on a multigraph with a legend
    ZRho.SetTitle("Z vs. Rho for Fitted Events; Rho [mm]; Z [mm]")
    ZRho.Add(ZRhoPlotMC, "AP")
    ZRho.Add(ZRhoPlotMCO, "AP")
    ZRho.Add(ZRhoPlotEV, "AP")
    ZRho.Add(ZRhoPlotEVO, "AP")
    ZRho.Add(lineplot, "L")
    ZRho.Add(boundaryplot, "C")
    ZRho.Draw("A")
    leg = canvas1.BuildLegend(0.67, 0.76, 0.9, 0.9)
    canvas1.Update()
    ZRho.GetYaxis().SetTitleOffset(1.36)

    latex = ROOT.TLatex() # write messages
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.02)
    xloc = 0.27
    yloc = 0.81
    latex.DrawText(xloc, yloc, "%i Events Outside the Inner Boundary Have Valid Energy Fits"%(filtercuts[3]))
    latex.DrawText(xloc, yloc-0.05, "%i Red and Blue Events"%(AVCountMC))
    latex.DrawText(xloc, yloc-0.1, "%i Pink and Green Events"%(OCountMC))
    
    canvas1.Print(outDir + outFileName) # save pdf
    
    print AVCountMC 
    print AVCountEV 
    print filtercuts 
    # =======================================================================



def CompareMCValidInvalidFitsECut(input_files, ECut):
    """ 

    Description:
    Plot Z vs Rho for MC and Fitted for energies greater than ECut
    Destinguish between MC events that are fitted inside and outside the AV
    
    Parameters: 
    file_name - Path to the RAT DS file to play around with
    Return: nothing

    """

    fitName = "partialFitter"
    CountEV = 0
    just_compare_neck = False

    # ======================================================================= for z vs rho

    ArrayRho_EV = array("d")
    ArrayZ_EV = array("d")

    ArrayRho_EVO = array("d")
    ArrayZ_EVO = array("d")

    ArrayRho_MC = array("d")
    ArrayZ_MC = array("d")

    ArrayRho_MCO = array("d")
    ArrayZ_MCO = array("d")

    # ======================================================================= 

    AVCountMC = 0
    AVCountEV = 0
    OCountEV = 0
    OCountMC = 0
    filtercuts = [0, 0, 0, 0] # total fit count, retrigger filtered, posvalid, Evalid

    for file_name in input_files :

        for ds, run in dsreader(file_name):

            for iev in range(0, ds.GetEVCount()): # ds is an entry, multiple EV events in one MC event possible so loop through them
                
                filtercuts[0] += 1
                if iev > 0: # retrigger filter
                    continue

                ev = ds.GetEV(iev) # single EV event
                fVertex = ev.GetFitResult(fitName).GetVertex(0)     

                filtercuts[1] += 1

                if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                    continue # valid position filter
                else:
                    filtercuts[2] += 1
                if not fVertex.ValidEnergy():
                    continue          
                filtercuts[3] += 1  
   
                PosEV = fVertex.GetPosition()                  # get position vector 
                REV = PosEV.Mag()                              # and radius
                RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                mc = ds.GetMC() # individual MC event

                PosMC = mc.GetMCParticle(0).GetPosition() # get position vector and radius 
                RMC = PosMC.Mag()
                RhoMC = math.sqrt(PosMC.X()**2 + PosMC.Y()**2)

                if PosEV.Z() >= 6000 and RhoEV >= 730 :    
                    in_AV = False         
                    out_AV = True
                elif PosEV.Z() < 6000 and REV > 6000 :
                    in_AV = False
                    out_AV = True
                else :
                    in_AV = True
                    out_AV = False

                if in_AV and fVertex.GetEnergy() > ECut:

                    # store values in arrays for plotting
                    ArrayZ_EV.append(PosEV.Z()) 
                    ArrayRho_EV.append(RhoEV)
                    AVCountEV += 1  

                    # store values in arrays for plotting
                    ArrayZ_MC.append(PosMC.Z())
                    ArrayRho_MC.append(RhoMC)
                    AVCountMC += 1 # keep count of mc events

                elif out_AV and fVertex.GetEnergy() > ECut:

                    # store values in arrays for plotting
                    ArrayZ_EVO.append(PosEV.Z()) 
                    ArrayRho_EVO.append(RhoEV)
                    OCountEV += 1  

                    # store values in arrays for plotting
                    ArrayZ_MCO.append(PosMC.Z())
                    ArrayRho_MCO.append(RhoMC)
                    OCountMC += 1  

    # =======================================================================
    
    canvas1 = ROOT.TCanvas("canvas1") # make canvas
    canvas1.cd()

    ZRhoPlotMC = ROOT.TGraph(AVCountMC, ArrayRho_MC, ArrayZ_MC) # make and format TGraph for Fitted events
    ZRhoPlotMC.SetTitle("MC Fitted Inside AV")
    ZRhoPlotMC.SetMarkerStyle(6)
    ZRhoPlotMC.SetMarkerColorAlpha(ROOT.kRed, 1)
    ZRhoPlotMC.SetLineColorAlpha(ROOT.kRed, 1)
    ZRhoPlotMC.SetFillStyle(0)
    ZRhoPlotMC.Draw()

    ZRhoPlotEV = ROOT.TGraph(AVCountEV, ArrayRho_EV, ArrayZ_EV) # make and format TGraph for Fitted events
    ZRhoPlotEV.SetTitle("Fitted Events in AV")
    ZRhoPlotEV.SetMarkerStyle(6)
    ZRhoPlotEV.SetMarkerColorAlpha(ROOT.kBlue, 1)
    ZRhoPlotEV.SetLineColorAlpha(ROOT.kBlue, 1)
    ZRhoPlotEV.SetFillStyle(0)
    ZRhoPlotEV.Draw()

    ZRhoPlotEVO = ROOT.TGraph(OCountEV, ArrayRho_EVO, ArrayZ_EVO) # make and format TGraph for validE stuff 
    ZRhoPlotEVO.SetTitle("Fitted Outside AV")
    ZRhoPlotEVO.SetMarkerStyle(6)
    ZRhoPlotEVO.SetMarkerColorAlpha(3, 1)
    ZRhoPlotEVO.SetLineColorAlpha(3, 1)
    ZRhoPlotEVO.SetFillStyle(0)
    ZRhoPlotEVO.Draw()

    ZRhoPlotMCO = ROOT.TGraph(OCountMC, ArrayRho_MCO, ArrayZ_MCO) # make and format TGraph for validE stuff 
    ZRhoPlotMCO.SetTitle("MC Fitted Outside AV")
    ZRhoPlotMCO.SetMarkerStyle(6)
    ZRhoPlotMCO.SetMarkerColorAlpha(6, 1)
    ZRhoPlotMCO.SetLineColorAlpha(6, 1)
    ZRhoPlotMCO.SetFillStyle(0)
    ZRhoPlotMCO.Draw()

    r = 6000 # save to file eventually
    theta = -math.pi/2
    x, y = array("d"), array("d")
    for i in range(100):
        x.append(r*math.cos(theta)), y.append(r*math.sin(theta))
        theta += 0.0302
    for i in range(50):
        x.append(730), y.append(6000+i*140)
    boundaryplot = ROOT.TGraph(150, x, y)
    boundaryplot.SetTitle("Inner AV")
    boundaryplot.SetLineStyle(3)
    boundaryplot.SetMarkerStyle(0)
    boundaryplot.SetFillStyle(0)
    boundaryplot.Draw()
    
    lineplot = ROOT.TGraph(2, array("d", [0, 8000]), array("d", [747.5, 747.5])) # a line at Z = 747.5 mm to tell where scint ends
    lineplot.SetTitle("747.5 mm")
    lineplot.SetMarkerStyle(0)
    lineplot.SetLineStyle(7)
    lineplot.SetFillStyle(0)
    lineplot.Draw()

    ZRho = ROOT.TMultiGraph() # draw the graphs on a multigraph with a legend
    ZRho.SetTitle("Z vs. Rho for Fitted Events with E > %.1f MeV; Rho [mm]; Z [mm]"%(ECut))
    ZRho.Add(ZRhoPlotMC, "AP")
    ZRho.Add(ZRhoPlotMCO, "AP")
    ZRho.Add(ZRhoPlotEV, "AP")
    ZRho.Add(ZRhoPlotEVO, "AP")
    ZRho.Add(lineplot, "L")
    ZRho.Add(boundaryplot, "C")
    ZRho.Draw("A")
    leg = canvas1.BuildLegend(0.67, 0.76, 0.9, 0.9)
    canvas1.Update()
    ZRho.GetYaxis().SetTitleOffset(1.36)

    latex = ROOT.TLatex() # write messages
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.02)
    xloc = 0.27
    yloc = 0.81
    latex.DrawText(xloc, yloc, "%i Events Outside the Inner Boundary Have Valid Energy Fits"%(filtercuts[3]))
    latex.DrawText(xloc, yloc-0.05, "%i Red and Blue Events"%(AVCountMC))
    latex.DrawText(xloc, yloc-0.1, "%i Pink and Green Events"%(OCountMC))
    
    canvas1.Print(outDir + outFileName) # save pdf
    
    print AVCountMC 
    print AVCountEV 
    print filtercuts 
    # =======================================================================



def CompareMCValidInvalidFits(input_files):

    """ 

    Description:
    Plot the Z vs Rho for MC and Fitted events
    Compare valid energy and invalid energy events

    Parameters: 
    file_name - Path to the RAT DS file to play around with
    Return: nothing

    """

    fitName = "partialFitter"
    CountEV = 0
    CountMC = 0
    just_compare_neck = False

    # ======================================================================= for z vs rho

    ArrayRho_EV = array("d")
    ArrayZ_EV = array("d")

    ArrayRho_F = array("d")
    ArrayZ_F = array("d")

    ArrayRho_MC = array("d")
    ArrayZ_MC = array("d")

    # ======================================================================= 

    FitBelowCount = 0 # keep track of fitted events below 747.5 mm
    OutsideEV = 0
    NeckCountMC = 0
    NeckCountEV = 0
    outsideradius = 0
    filtercuts = [0, 0, 0, 0] # total fit count, retrigger filtered, posvalid, Evalid

    for file_name in input_files :

        for ds, run in dsreader(file_name):

            mc = ds.GetMC() # individual MC event
            CountMC += 1 # keep count of mc events

            PosMC = mc.GetMCParticle(0).GetPosition() # get position vector and radius 
            RMC = PosMC.Mag()
            RhoMC = math.sqrt(PosMC.X()**2 + PosMC.Y()**2)

            if RMC > 6000 and RhoMC < 730 :
                NeckCountMC += 1

            # store values in arrays for plotting
            ArrayZ_MC.append(PosMC.Z())
            ArrayRho_MC.append(RhoMC)

            for iev in range(0, ds.GetEVCount()): # ds is an entry, multiple EV events in one MC event possible so loop through them
                filtercuts[0] += 1

                #if iev > 0: # retrigger filter
                #    continue

                ev = ds.GetEV(iev) # single EV event

                filtercuts[1] += 1

                if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                    continue # valid position filter

                filtercuts[2] += 1 # keep count of valid pos fitted events

                fVertex = ev.GetFitResult(fitName).GetVertex(0)     
                    
                PosEV = fVertex.GetPosition()                  # get position vector 
                REV = PosEV.Mag()                              # and radius
                RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                if PosEV.Z() < 747.5 :
                    FitBelowCount += 1

                if PosEV.Z() > 6000 and RhoEV < 730 :
                    NeckCountEV += 1
                
                if PosEV.Z() >= 6000 and RhoEV >= 730 :
                    OutsideEV += 1
                if PosEV.Z() < 6000 and REV > 6000 :
                    OutsideEV += 1
                if REV > 6000 :
                    outsideradius += 1

                # store values in arrays for plotting
                ArrayZ_EV.append(PosEV.Z()) 
                ArrayRho_EV.append(RhoEV)

                if not fVertex.ValidEnergy():
                    continue      
                ArrayZ_F.append(PosEV.Z()) 
                ArrayRho_F.append(RhoEV)      
                filtercuts[3] += 1   # valid pos and energy

    # =======================================================================
    
    canvas1 = ROOT.TCanvas("canvas1") # make canvas
    canvas1.cd()

    ZRhoPlotMC = ROOT.TGraph(NeckCountMC, ArrayRho_MC, ArrayZ_MC) # make and format TGraph for Fitted events
    ZRhoPlotMC.SetTitle("MC Events")
    ZRhoPlotMC.SetMarkerStyle(6)
    ZRhoPlotMC.SetMarkerColorAlpha(ROOT.kRed, 1)
    ZRhoPlotMC.SetLineColorAlpha(ROOT.kRed, 1)
    ZRhoPlotMC.SetFillStyle(0)
    ZRhoPlotMC.Draw()

    ZRhoPlotEV = ROOT.TGraph(filtercuts[2], ArrayRho_EV, ArrayZ_EV) # make and format TGraph for Fitted events
    ZRhoPlotEV.SetTitle("Valid Position Events")
    ZRhoPlotEV.SetMarkerStyle(6)
    ZRhoPlotEV.SetMarkerColorAlpha(ROOT.kBlue, 1)
    ZRhoPlotEV.SetLineColorAlpha(ROOT.kBlue, 1)
    ZRhoPlotEV.SetFillStyle(0)
    ZRhoPlotEV.Draw()

    ZRhoPlotF = ROOT.TGraph(filtercuts[3], ArrayRho_F, ArrayZ_F) # make and format TGraph for validE stuff 
    ZRhoPlotF.SetTitle("Valid Pos & Energy Events")
    ZRhoPlotF.SetMarkerStyle(6)
    ZRhoPlotF.SetMarkerColorAlpha(3, 1)
    ZRhoPlotF.SetLineColorAlpha(3, 1)
    ZRhoPlotF.SetFillStyle(0)
    ZRhoPlotF.Draw()

    r = 6000 # save to file eventually
    theta = -math.pi/2
    x, y = array("d"), array("d")
    for i in range(100):
        x.append(r*math.cos(theta)), y.append(r*math.sin(theta))
        theta += 0.0302
    for i in range(50):
        x.append(730), y.append(6000+i*140)
    boundaryplot = ROOT.TGraph(150, x, y)
    boundaryplot.SetTitle("Inner AV")
    boundaryplot.SetLineStyle(3)
    boundaryplot.SetMarkerStyle(0)
    boundaryplot.SetFillStyle(0)
    boundaryplot.Draw()
    
    lineplot = ROOT.TGraph(2, array("d", [0, 8000]), array("d", [747.5, 747.5])) # a line at Z = 747.5 mm to tell where scint ends
    lineplot.SetTitle("747.5 mm")
    lineplot.SetMarkerStyle(0)
    lineplot.SetLineStyle(7)
    lineplot.SetFillStyle(0)
    lineplot.Draw()

    ZRho = ROOT.TMultiGraph() # draw the graphs on a multigraph with a legend
    ZRho.SetTitle("Z vs. Rho for Fitted Events; Rho [mm]; Z [mm]")
    ZRho.Add(ZRhoPlotMC, "AP")
    ZRho.Add(ZRhoPlotEV, "AP")
    ZRho.Add(ZRhoPlotF, "AP")
    ZRho.Add(lineplot, "L")
    ZRho.Add(boundaryplot, "C")
    ZRho.Draw("A")
    leg = canvas1.BuildLegend(0.68, 0.79, 0.9, 0.9)
    canvas1.Update()
    ZRho.GetYaxis().SetTitleOffset(1.36)

    latex = ROOT.TLatex() # write messages
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.02)
    xloc = 0.3
    yloc = 0.81
    latex.DrawText(xloc, yloc, "%i/%i MC Events In Neck "%(NeckCountMC, CountMC)) 
    latex.DrawText(xloc, yloc-0.05, "%i/%i Events With Valid Pos Fit Have Valid Energy Fit"%(filtercuts[3], filtercuts[2]))
    #latex.DrawText(xloc, yloc-0.1, "%i Fitted Into Neck / %i MC Neck Events"%(NeckCountEV, NeckCountMC))
    #latex.DrawText(xloc, yloc-0.1, "%i/%i Fitted With Z < 747.5 mm"%(FitBelowCount, CountEV))
    #latex.DrawText(xloc, yloc-0.1, "%i Fitted Outside AV, %i Radius > 6 m"%(OutsideEV, outsideradius))
    
    canvas1.Print(outDir + outFileName) # save pdf
    
    print FitBelowCount # keep track of fitted events below 747.5 mm
    print OutsideEV 
    print NeckCountMC 
    print NeckCountEV 
    print filtercuts 
    # =======================================================================



def NeckEventPositions(input_files):
    """ 
    Description:
    Plot the Z vs Rho for MC and Fitted 
    See where MC events in neck are fitted
    Compare valid energy and valid position fits

    Parameters: 
    file_name - Path to the RAT DS file to play around with
    Return: nothing
    """

    fitName = "partialFitter"
    CountEV = 0
    CountMC = 0
    just_compare_neck = True

    # ======================================================================= for z vs rho

    ArrayRho_EV = array("d")
    ArrayZ_EV = array("d")

    ArrayRho_F = array("d")
    ArrayZ_F = array("d")

    ArrayRho_MC = array("d")
    ArrayZ_MC = array("d")

    # ======================================================================= 

    FitBelowCount = 0 # keep track of fitted events below 747.5 mm
    OutsideEV = 0
    NeckCountMC = 0
    NeckCountEV = 0
    outsideradius = 0
    filtercuts = [0, 0, 0, 0] # total fit count, retrigger filtered, posvalid, Evalid

    for file_name in input_files:
        for ds, run in dsreader(file_name):

            mc = ds.GetMC() # individual MC event
            CountMC += 1 # keep count of mc events

            PosMC = mc.GetMCParticle(0).GetPosition() # get position vector and radius 
            RMC = PosMC.Mag()
            RhoMC = math.sqrt(PosMC.X()**2 + PosMC.Y()**2)

            if RMC > 6000 and RhoMC < 730 and PosMC.Z() < 8000: # "and PosMC.Z() < 8000" added just to see bottom of neck
                NeckCountMC += 1

                # store values in arrays for plotting
                ArrayZ_MC.append(PosMC.Z())
                ArrayRho_MC.append(RhoMC)

                if just_compare_neck :

                    for iev in range(0, ds.GetEVCount()): # ds is an entry, multiple EV events in one MC event possible so loop through them
                        filtercuts[0] += 1
                        if iev > 0: # retrigger filter
                            continue

                        ev = ds.GetEV(iev) # single EV event

                        filtercuts[1] += 1

                        if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                            continue # valid position filter

                        filtercuts[2] += 1

                        fVertex = ev.GetFitResult(fitName).GetVertex(0)     
                            
                        PosEV = fVertex.GetPosition()                  # get position vector 
                        REV = PosEV.Mag()                              # and radius
                        RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                        if PosEV.Z() < 747.5 :
                            FitBelowCount += 1

                        if PosEV.Z() > 6000 and RhoEV < 730 :
                            NeckCountEV += 1
                        
                        if PosEV.Z() >= 6000 and RhoEV >= 730 :
                            OutsideEV += 1
                        elif PosEV.Z() < 6000 and REV > 6000 :
                            OutsideEV += 1
                        if REV > 6000 :
                            outsideradius += 1

                        # store values in arrays for plotting
                        ArrayZ_EV.append(PosEV.Z()) 
                        ArrayRho_EV.append(RhoEV)
                        CountEV += 1 # keep count of fitted events

                        if not fVertex.ValidEnergy():
                            continue      
                        ArrayZ_F.append(PosEV.Z()) 
                        ArrayRho_F.append(RhoEV)      
                        filtercuts[3] += 1  

            if not just_compare_neck :

                for iev in range(0, ds.GetEVCount()): # ds is an entry, multiple EV events in one MC event possible so loop through them
                    filtercuts[0] += 1
                    if iev > 0: # retrigger filter
                        continue

                    ev = ds.GetEV(iev) # single EV event

                    filtercuts[1] += 1

                    if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                        continue # valid position filter

                    filtercuts[2] += 1 # keep count of valid pos fitted events

                    fVertex = ev.GetFitResult(fitName).GetVertex(0)     
                        
                    PosEV = fVertex.GetPosition()                  # get position vector 
                    REV = PosEV.Mag()                              # and radius
                    RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                    if PosEV.Z() < 747.5 :
                        FitBelowCount += 1

                    if PosEV.Z() > 6000 and RhoEV < 730 :
                        NeckCountEV += 1
                    
                    if PosEV.Z() >= 6000 and RhoEV >= 730 :
                        OutsideEV += 1
                    if PosEV.Z() < 6000 and REV > 6000 :
                        OutsideEV += 1
                    if REV > 6000 :
                        outsideradius += 1

                    # store values in arrays for plotting
                    ArrayZ_EV.append(PosEV.Z()) 
                    ArrayRho_EV.append(RhoEV)

                    if not fVertex.ValidEnergy():
                        continue      
                    ArrayZ_F.append(PosEV.Z()) 
                    ArrayRho_F.append(RhoEV)      
                    filtercuts[3] += 1   # valid pos and energy

    # =======================================================================
    
    canvas1 = ROOT.TCanvas("canvas1") # make canvas
    canvas1.cd()

    ZRhoPlotMC = ROOT.TGraph(NeckCountMC, ArrayRho_MC, ArrayZ_MC) # make and format TGraph for Fitted events
    ZRhoPlotMC.SetTitle("MC Neck Events")
    ZRhoPlotMC.SetMarkerStyle(6)
    ZRhoPlotMC.SetMarkerColorAlpha(ROOT.kRed, 1)
    ZRhoPlotMC.SetLineColorAlpha(ROOT.kRed, 1)
    ZRhoPlotMC.SetFillStyle(0)
    ZRhoPlotMC.Draw()

    ZRhoPlotEV = ROOT.TGraph(filtercuts[2], ArrayRho_EV, ArrayZ_EV) # make and format TGraph for Fitted events
    ZRhoPlotEV.SetTitle("Valid Position Events")
    ZRhoPlotEV.SetMarkerStyle(6)
    ZRhoPlotEV.SetMarkerColorAlpha(ROOT.kBlue, 1)
    ZRhoPlotEV.SetLineColorAlpha(ROOT.kBlue, 1)
    ZRhoPlotEV.SetFillStyle(0)
    ZRhoPlotEV.Draw()

    ZRhoPlotF = ROOT.TGraph(filtercuts[3], ArrayRho_F, ArrayZ_F) # make and format TGraph for validE stuff 
    ZRhoPlotF.SetTitle("Valid Pos & Energy Events")
    ZRhoPlotF.SetMarkerStyle(6)
    ZRhoPlotF.SetMarkerColorAlpha(3, 1)
    ZRhoPlotF.SetLineColorAlpha(3, 1)
    ZRhoPlotF.SetFillStyle(0)
    ZRhoPlotF.Draw()

    r = 6000 # save to file eventually
    theta = -math.pi/2
    x, y = array("d"), array("d")
    for i in range(100):
        x.append(r*math.cos(theta)), y.append(r*math.sin(theta))
        theta += 0.0302
    for i in range(50):
        x.append(730), y.append(6000+i*140)
    boundaryplot = ROOT.TGraph(150, x, y)
    boundaryplot.SetTitle("Inner AV")
    boundaryplot.SetLineStyle(3)
    boundaryplot.SetMarkerStyle(0)
    boundaryplot.SetFillStyle(0)
    boundaryplot.Draw()
    
    lineplot = ROOT.TGraph(2, array("d", [0, 8000]), array("d", [747.5, 747.5])) # a line at Z = 747.5 mm to tell where scint ends
    lineplot.SetTitle("747.5 mm")
    lineplot.SetMarkerStyle(0)
    lineplot.SetLineStyle(7)
    lineplot.SetFillStyle(0)
    lineplot.Draw()

    ZRho = ROOT.TMultiGraph() # draw the graphs on a multigraph with a legend
    ZRho.SetTitle("Z vs. Rho for Fitted Events; Rho [mm]; Z [mm]")
    ZRho.Add(ZRhoPlotMC, "AP")
    ZRho.Add(ZRhoPlotEV, "AP")
    ZRho.Add(ZRhoPlotF, "AP")
    ZRho.Add(lineplot, "L")
    ZRho.Add(boundaryplot, "C")
    ZRho.Draw("A")
    leg = canvas1.BuildLegend(0.68, 0.79, 0.9, 0.9)
    canvas1.Update()
    ZRho.GetYaxis().SetTitleOffset(1.36)

    latex = ROOT.TLatex() # write messages
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.02)
    xloc = 0.3
    yloc = 0.81
    #latex.DrawText(xloc, yloc, "%i/%i MC Events In Neck "%(NeckCountMC, CountMC)) 
    latex.DrawText(xloc, yloc, "%i/%i MC Neck Events Betwenn 6-8 m"%(NeckCountMC, CountMC)) 
    latex.DrawText(xloc, yloc-0.05, "%i/%i Events With Valid Pos Fit Have Valid Energy Fit"%(filtercuts[3], filtercuts[2]))
    latex.DrawText(xloc, yloc-0.1, "%i Fitted Into Neck / %i MC Neck Events"%(NeckCountEV, NeckCountMC))
    #latex.DrawText(xloc, yloc-0.1, "%i/%i Fitted With Z < 747.5 mm"%(FitBelowCount, CountEV))
    #latex.DrawText(xloc, yloc-0.1, "%i Fitted Outside AV, %i Radius > 6 m"%(OutsideEV, outsideradius))
    
    canvas1.Print(outDir + outFileName) # save pdf
    
    print FitBelowCount # keep track of fitted events below 747.5 mm
    print OutsideEV 
    print NeckCountMC 
    print NeckCountEV 
    print filtercuts 
    # =======================================================================



def PlotDataPosHist(input_files):

    """ 

    Description
    -----------
    Plot a 2D Z vs Rho histogram for Fitted evs in data

    Parameters
    ----------
    input_files : list of str 
        Path to the RAT DS file(s) to play around with
    

    """

    fitName = "partialFitter"
    CountEV = 0
    just_compare_neck = False

    # ======================================================================= for z vs rho

    ZRhoPlotEV = ROOT.TH2D("hist", "hist", 100, 0, 8000, 100, -8000, 8000)
    ZRhoPlotEV.SetDirectory(0)
    ZRhoPlotEV.SetStats(0)

    # ======================================================================= 

    FitBelowCount = 0 # keep track of fitted events below 747.5 mm
    OutsideEV = 0
    NeckCountEV = 0
    outsideradius = 0
    filtercuts = [0, 0, 0, 0] # total fit count, retrigger filtered, posvalid, Evalid

    for file_name in input_files :

        dsread = dsreader(file_name)
        print "\n\tSuccessfully ran dsreader(", file_name, ")\n"

        for ds, run in dsread:
        
            for iev in range(0, ds.GetEVCount()): # ds is an entry, multiple EV events in one MC event possible so loop through them

                filtercuts[0] += 1
                if filtercuts[0] < 10000: continue
                if filtercuts[0] % 1000 == 0: print "Event", filtercuts[0]

                if iev > 0: # retrigger filter
                    continue

                ev = ds.GetEV(iev) # single EV event

                filtercuts[1] += 1

                if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                    continue # valid position filter

                filtercuts[2] += 1 # keep count of valid pos fitted events

                fVertex = ev.GetFitResult(fitName).GetVertex(0)     
                    
                PosEV = fVertex.GetPosition()                  # get position vector 
                PosEV_Z = PosEV.Z()
                REV = PosEV.Mag()                              # and radius
                RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                if PosEV_Z < 747.5 :
                    FitBelowCount += 1

                if PosEV_Z > 6000 and RhoEV < 730 :
                    NeckCountEV += 1
                
                if PosEV_Z >= 6000 and RhoEV >= 730 :
                    OutsideEV += 1
                if PosEV_Z < 6000 and REV > 6000 :
                    OutsideEV += 1
                if REV > 6000 :
                    outsideradius += 1

                # plot in hist
                ZRhoPlotEV.Fill(RhoEV, PosEV_Z)
            
            if filtercuts[0] >= 20000: break



    # =======================================================================
    
    canvas1 = ROOT.TCanvas("canvas1") # make canvas
    canvas1.cd()

    #ZRhoPlotEV = ROOT.TGraph(filtercuts[2], ArrayRho_EV, ArrayZ_EV) # make and format TGraph for Fitted events
    ZRhoPlotEV.SetTitle("Z vs. Rho for Fitted Events; Rho [mm]; Z [mm]")
    #ZRhoPlotEV.SetMarkerStyle(6)
    #ZRhoPlotEV.SetMarkerColorAlpha(ROOT.kBlue, 1)
    #ZRhoPlotEV.SetLineColorAlpha(ROOT.kBlue, 1)
    #ZRhoPlotEV.SetFillStyle(0)
    ZRhoPlotEV.Draw("COL")

    '''
    r = 6000 # save to file eventually
    theta = -math.pi/2
    x, y = array("d"), array("d")
    for i in range(100):
        x.append(r*math.cos(theta)), y.append(r*math.sin(theta))
        theta += 0.0302
    for i in range(50):
        x.append(730), y.append(6000+i*140)
    boundaryplot = ROOT.TGraph(150, x, y)
    boundaryplot.SetTitle("Inner AV")
    boundaryplot.SetLineStyle(3)
    boundaryplot.SetMarkerStyle(0)
    boundaryplot.SetFillStyle(0)
    boundaryplot.Draw()
    
    lineplot = ROOT.TGraph(2, array("d", [0, 8000]), array("d", [747.5, 747.5])) # a line at Z = 747.5 mm to tell where scint ends
    lineplot.SetTitle("747.5 mm")
    lineplot.SetMarkerStyle(0)
    lineplot.SetLineStyle(7)
    lineplot.SetFillStyle(0)
    lineplot.Draw()

    ZRho = ROOT.TMultiGraph() # draw the graphs on a multigraph with a legend
    ZRho.SetTitle("Z vs. Rho for Fitted Events; Rho [mm]; Z [mm]")
    ZRho.Add(ZRhoPlotEV, "AP")
    ZRho.Add(lineplot, "L")
    ZRho.Add(boundaryplot, "C")
    ZRho.Draw("A")
    leg = canvas1.BuildLegend(0.68, 0.79, 0.9, 0.9)
    canvas1.Update()
    ZRho.GetYaxis().SetTitleOffset(1.36)
    '''

    latex = ROOT.TLatex() # write messages
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.02)
    xloc = 0.5
    yloc = 0.9
    #latex.DrawText(xloc, yloc-0.05, "%i/%i Events With Valid Pos Fit Have Valid Energy Fit"%(filtercuts[3], filtercuts[2]))
    latex.DrawText(xloc, yloc-0.1, "%i/%i Fitted With Z < 747.5 mm"%(FitBelowCount, CountEV))
    latex.DrawText(xloc, yloc-0.15, "%i Fitted Outside AV, %i Radius > 6 m"%(OutsideEV, outsideradius))
    
    canvas1.Print(outDir + outFileName) # save pdf
    
    print FitBelowCount, "fitted events below z = 747.5 mm" # keep track of fitted events below 747.5 mm
    print OutsideEV, "events outside the AV"
    print NeckCountEV, "events in the neck"
    print filtercuts, " total fit count, retrigger filtered, posvalid, Evalid"
    # =======================================================================



def PlotDataPosHist(input_files):

    """ 

    Description
    -----------
    Plot a 2D Z vs Rho histogram for Fitted evs in data

    Parameters
    ----------
    input_files : list of str 
        Path to the RAT ROOT file(s) to play around with
    

    """

    fitName = "partialFitter"
    CountEV = 0
    just_compare_neck = False

    # ======================================================================= for z vs rho

    ZRhoPlotEV = ROOT.TH2D("hist", "hist", 100, 0, 8000, 100, -8000, 8000)
    ZRhoPlotEV.SetDirectory(0)
    ZRhoPlotEV.SetStats(0)

    # ======================================================================= 

    FitBelowCount = 0 # keep track of fitted events below 747.5 mm
    OutsideEV = 0
    NeckCountEV = 0
    outsideradius = 0
    filtercuts = [0, 0, 0, 0] # total fit count, retrigger filtered, posvalid, Evalid

    for file_name in input_files :

        oldFile = ROOT.TFile(fname) 
        oldTree_T = oldFile.Get("T")
        #print "\n\tSuccessfully ran dsreader(", file_name, ")\n"

        for ds, run in dsread:
        
            for iev in range(0, ds.GetEVCount()): # ds is an entry, multiple EV events in one MC event possible so loop through them

                filtercuts[0] += 1
                if filtercuts[0] < 10000: continue
                if filtercuts[0] % 1000 == 0: print "Event", filtercuts[0]

                if iev > 0: # retrigger filter
                    continue

                ev = ds.GetEV(iev) # single EV event

                filtercuts[1] += 1

                if not ev.FitResultExists(fitName) or not ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() or not ev.GetFitResult(fitName).GetVertex(0).ValidPosition() :
                    continue # valid position filter

                filtercuts[2] += 1 # keep count of valid pos fitted events

                fVertex = ev.GetFitResult(fitName).GetVertex(0)     
                    
                PosEV = fVertex.GetPosition()                  # get position vector 
                PosEV_Z = PosEV.Z()
                REV = PosEV.Mag()                              # and radius
                RhoEV = math.sqrt(PosEV.X()**2 + PosEV.Y()**2)

                if PosEV_Z < 747.5 :
                    FitBelowCount += 1

                if PosEV_Z > 6000 and RhoEV < 730 :
                    NeckCountEV += 1
                
                if PosEV_Z >= 6000 and RhoEV >= 730 :
                    OutsideEV += 1
                if PosEV_Z < 6000 and REV > 6000 :
                    OutsideEV += 1
                if REV > 6000 :
                    outsideradius += 1

                # plot in hist
                ZRhoPlotEV.Fill(RhoEV, PosEV_Z)
            
            if filtercuts[0] >= 20000: break



    # =======================================================================
    
    canvas1 = ROOT.TCanvas("canvas1") # make canvas
    canvas1.cd()

    #ZRhoPlotEV = ROOT.TGraph(filtercuts[2], ArrayRho_EV, ArrayZ_EV) # make and format TGraph for Fitted events
    ZRhoPlotEV.SetTitle("Z vs. Rho for Fitted Events; Rho [mm]; Z [mm]")
    #ZRhoPlotEV.SetMarkerStyle(6)
    #ZRhoPlotEV.SetMarkerColorAlpha(ROOT.kBlue, 1)
    #ZRhoPlotEV.SetLineColorAlpha(ROOT.kBlue, 1)
    #ZRhoPlotEV.SetFillStyle(0)
    ZRhoPlotEV.Draw("COL")

    '''
    r = 6000 # save to file eventually
    theta = -math.pi/2
    x, y = array("d"), array("d")
    for i in range(100):
        x.append(r*math.cos(theta)), y.append(r*math.sin(theta))
        theta += 0.0302
    for i in range(50):
        x.append(730), y.append(6000+i*140)
    boundaryplot = ROOT.TGraph(150, x, y)
    boundaryplot.SetTitle("Inner AV")
    boundaryplot.SetLineStyle(3)
    boundaryplot.SetMarkerStyle(0)
    boundaryplot.SetFillStyle(0)
    boundaryplot.Draw()
    
    lineplot = ROOT.TGraph(2, array("d", [0, 8000]), array("d", [747.5, 747.5])) # a line at Z = 747.5 mm to tell where scint ends
    lineplot.SetTitle("747.5 mm")
    lineplot.SetMarkerStyle(0)
    lineplot.SetLineStyle(7)
    lineplot.SetFillStyle(0)
    lineplot.Draw()

    ZRho = ROOT.TMultiGraph() # draw the graphs on a multigraph with a legend
    ZRho.SetTitle("Z vs. Rho for Fitted Events; Rho [mm]; Z [mm]")
    ZRho.Add(ZRhoPlotEV, "AP")
    ZRho.Add(lineplot, "L")
    ZRho.Add(boundaryplot, "C")
    ZRho.Draw("A")
    leg = canvas1.BuildLegend(0.68, 0.79, 0.9, 0.9)
    canvas1.Update()
    ZRho.GetYaxis().SetTitleOffset(1.36)
    '''

    latex = ROOT.TLatex() # write messages
    latex.SetTextFont(62)
    latex.SetNDC()
    latex.SetTextSize(0.02)
    xloc = 0.5
    yloc = 0.9
    #latex.DrawText(xloc, yloc-0.05, "%i/%i Events With Valid Pos Fit Have Valid Energy Fit"%(filtercuts[3], filtercuts[2]))
    latex.DrawText(xloc, yloc-0.1, "%i/%i Fitted With Z < 747.5 mm"%(FitBelowCount, CountEV))
    latex.DrawText(xloc, yloc-0.15, "%i Fitted Outside AV, %i Radius > 6 m"%(OutsideEV, outsideradius))
    
    canvas1.Print(outDir + outFileName) # save pdf
    
    print FitBelowCount, "fitted events below z = 747.5 mm" # keep track of fitted events below 747.5 mm
    print OutsideEV, "events outside the AV"
    print NeckCountEV, "events in the neck"
    print filtercuts, " total fit count, retrigger filtered, posvalid, Evalid"
    # =======================================================================



if __name__ == '__main__':
    #globals()[sys.argv[1]](sys.argv[2])
    PlotDataPosHist(["/home/eigenboi/scratch/root_files/data_bipo214_results_test_of_concept.root"])