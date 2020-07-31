# attempt to translate BiPo tagger to python
# original code by Ana Sofia and Fady Shaker
# Aidan Patton, 07-2020

'''
Apply cuts to isolate BiPo214 coincidence events, save to a new .root file.
'''

import ROOT
import rat 
import tagger_settings

plot_dir = tagger_settings.plot_dir
root_dir = tagger_settings.root_dir
#in_data_file = tagger_settings.in_data_file
data_file_arr = tagger_settings.data_file_arr
mc_file_arr = tagger_settings.mc_file_arr
identifier = tagger_settings.identifier

r_cut = tagger_settings.r_cut                 # [mm]
tres_cut_min = tagger_settings.tres_cut_min   # [ns]
tres_cut_max = tagger_settings.tres_cut_max   # [ns]

bi_z_min = tagger_settings.bi_z_min # [mm]
bi_r_min = tagger_settings.bi_r_min # [mm]
bi_r_max = tagger_settings.bi_r_max # [mm] 
bi_nhit_cleaned_min = tagger_settings.bi_nhit_cleaned_min

po_z_min = tagger_settings.po_z_min # [mm]
po_r_max = tagger_settings.po_r_max # [mm] 
po_nhit_cleaned_min = tagger_settings.po_nhit_cleaned_min
po_nhit_cleaned_max = tagger_settings.po_nhit_cleaned_max

bipo_delta_r_max = tagger_settings.bipo_delta_r_max # [mm]
bipo_delta_t_min = tagger_settings.bipo_delta_t_min # [ns]
bipo_delta_t_max = tagger_settings.bipo_delta_t_max # [ns]

is_mc = tagger_settings.is_mc
time_res_cut = tagger_settings.time_res_cut
fitName = tagger_settings.fitName 



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def bipo214_comb(is_mc, time_res_cut, in_files): # main function to call

    ''' 
    Description
    -----------
    Combs through files and finds all BiPo214 coincidence events in all of them, 
    saving them to an output root file (which may be pre-existing).

    Parameters
    ----------
    is_mc : bool 
        Whether the input files are MC simulated data
    time_res_cut : bool
        Whether or not to apply a time residual cut
    in_files : list
        List of input root files to search through
    
    Returns
    -------
    No returns, but return is used to exit function.

    '''

    if is_mc :
        out_file_prefix = "mc_"
    else :
        out_file_prefix = "data_"
    
    if in_files in ("", None, [], [""]) or type(in_files) is not list :
        print "Invalid input files. Must be a list of strings."
        return 
    
    # First load all the utilities: PMT information, lightpath calculator. Without this call the calculated time residual will be BIASED at the first call of the corresponding lightpath and group velocity calculations 
    rat.utility().Get().LoadDBAndBeginRun() 
    
    fname_count = 0 # counter (starts at 1) for what file is being analysed
    nfiles = len(in_files) # total number of files 
    coincicount = [0] * nfiles # list of number of bipo events found in each file.

    for fname in in_files :

        fname_count += 1

        dsread = rat.dsreader(fname) # (case sensitive)
        print "\n\tSuccessfully ran rat.dsreader(", fname, ") !\n"
        #print "\t", dsread.GetEntryCount(), "entries in total.\n" # breaks

        # the output file to store the bipo events, will save in root_dir with name data_bipo214_results.root
        # Create a new file + a clone of old tree in new file
        oldFile = ROOT.TFile(fname) # READ mode by default
        oldTree = oldFile.Get("T") # main tree with EV info, https://snopl.us/docs/rat/user_manual/html/data_structure.html
        outFile = ROOT.TFile(str(root_dir + out_file_prefix + "bipo214_results_" + identifier + ".root"), "UPDATE")
        newTree = oldTree.CloneTree(0) # 0 parameter ensures no entries are coppied, all branches are copied though

        print "\tis_mc ==", is_mc, "\n"     
        coincicount[fname_count-1] = build_bipo214_results(dsread, time_res_cut, newTree) 
        
        #newTree.Print() # see the tree structure (number of entries, the branches, and the leaves)
        newTree.AutoSave()
        outFile.Write() # saves the file

        print "\n\tFinal coincidence event count for file", fname_count, "/", nfiles, ":", coincicount[fname_count-1]

    print "\n\tFinal coincidence event count for all", nfiles, "file(s):", sum(coincicount), "\n"

    outFile.Close() # worst case scenario, save iEntry number for each root file

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def build_bipo214_results(dsReader, time_res_cut, newTree):

    '''

    Description
    -----------
    Find all BiPo coincidence events in an input file and save them to newTree.

    Parameters
    ----------
    dsReader : rat.RAT::DU::DSReader generator
        RAT data accessed with the dsReader.
    time_res_cut : bool
        Whether or not to apply a time residual cut.
    newTree : ROOT.TTree
        TTree in an output file to save events to.

    '''

    iEntry = -1
    coincicount = 0 # count of bipo214 coincidence events

    # loop to identify the Bi candidates. This currently needs a RatDB connection every time it's run - any way to keep tables loaded?
    for ds0, run in dsReader : # ds0 is equivalent to rat.dsreader.GetEntry(i) in c++, dsReader = rat.dsreader(input_files)

        iEntry += 1
        if iEntry % 1000 == 0 : print "Entry", iEntry

        for iev in range(0, ds0.GetEVCount()) :

            bi_ev = ds0.GetEV(iev)

            beta_pibo214 = is_bipo214_beta(dsReader, bi_ev, iEntry, time_res_cut)
         
            if beta_pibo214 == False : continue

            coincicount += 1
            print coincicount, "coincidence events!    ( Entry =", iEntry, ")"

            newTree.Fill(bi_ev) # add event to new file to save if is a bipo event, don't know best way to do this

    return coincicount

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def is_bipo214_beta(dsReader, bi_ev, entry_indx, time_res_cut):

    '''

    Description
    -----------
    Given a single input Bi event, determine whether it is a bipo coincidence event by searching whole file for its respective Po.

    Parameters
    ----------
    dsReader : rat.RAT::DU::DSReader generator
        RAT data accessed with the dsReader.
    bi_ev : rat.RAT::DS::EV
        bipo event candidate.
    entry_indx : int
        Index of bi_ev in its file.
    time_res_cut : bool
        Whether or not to apply a time residual cut.
    
    Returns
    -------
    is_bipo214_beta_ev : bool
        Whether given event is a BiPo214 beta event

    '''
    
    if not bi_ev.FitResultExists(fitName) or not bi_ev.GetFitResult(fitName).GetVertex(0).ValidPosition() or not bi_ev.GetFitResult(fitName).GetVertex(0).ValidTime() :
        return False

    is_bipo214_beta_ev = False # False by default
    fill_bi_ev_OK = False
    fill_po_ev_OK = False

    bi_fit_pos = bi_ev.GetFitResult(fitName).GetVertex(0).GetPosition()
    bi_ev_time = bi_ev.GetClockCount50()*20 # to search between different events we have to use the clk not the fitted time (fit time is within the event window [0,400]ns)

    # Bi position
    bi_z = bi_fit_pos.Z()
    bi_r = bi_fit_pos.Mag()

    """ Apply Cuts """
    if not bi_z > bi_z_min :                        # FV cut (in scintillator cap)
        return False
    if not bi_r > bi_r_min or not bi_r < bi_r_max : # FV cut (Bi  2m < R < 6 m)
        return False
    # r and z passed FV cuts

    # Bi nhitsCleaned cut  
    if bi_ev.GetNhitsCleaned() < bi_nhit_cleaned_min :
        return False

    if is_mc == False : # the data cleaning (dc) cuts are only used for DATA, using it for MC will create a crash!
        latestPass = bi_ev.GetDataCleaningFlags().GetLatestPass()
        bi_dcApplied = bi_ev.GetDataCleaningFlags().GetApplied(latestPass).GetULong64_t(0)
        bi_dcFlagged = bi_ev.GetDataCleaningFlags().GetFlags(latestPass).GetULong64_t(0) 
        if ((bi_dcApplied & 0x210000000242) & bi_dcFlagged) != (bi_dcApplied & 0x210000000242) : # dc bitmask 
            return False
    """ All Bi Cuts Successful! """
    
    fill_bi_ev_OK = fill_ev_time_res(bi_ev, True, time_res_cut, bi_fit_pos) # checks pmt time residual cuts. To get here, every other cut must be passed, so this is the final test

    quit_search = False
    iEntryPo = -1
    # loop to identify the Po candidates
    for ds1, run in dsReader : #dsReader is rat.dsreader(input_files)

        iEntryPo += 1
        if iEntryPo < entry_indx : continue # want to start at entry_indx

        for iev in range(0, ds1.GetEVCount()) :

            po_ev = ds1.GetEV(iev)
            if not po_ev.FitResultExists(fitName) or not po_ev.GetFitResult(fitName).GetVertex(0).ValidPosition() or not po_ev.GetFitResult(fitName).GetVertex(0).ValidTime() :
                continue

            po_fit_pos = po_ev.GetFitResult(fitName).GetVertex(0).GetPosition()
            po_ev_time = po_ev.GetClockCount50()*20 # to search between different events we have to use the clk not the fitted time (fit time is within the event window [0,400]ns)

            # Po position
            po_z = po_fit_pos.Z()
            po_r = po_fit_pos.Mag()

            ''' Apply Po Cuts '''
            if po_z < po_z_min : continue # Po z > 0.85 m (in scintillator cap)
            if po_r > po_r_max : continue # po  R < 6 m
            # r and z cuts both passed if made it here

            if po_ev.GetNhitsCleaned() < po_nhit_cleaned_min or po_ev.GetNhitsCleaned() > po_nhit_cleaned_max : 
                continue # nhits cut

            if is_mc == False : # unnecessary to apply data cleaning bitmask to MC
                latestPass = po_ev.GetDataCleaningFlags().GetLatestPass()
                po_dcApplied = po_ev.GetDataCleaningFlags().GetApplied(latestPass).GetULong64_t(0)
                po_dcFlagged = po_ev.GetDataCleaningFlags().GetFlags(latestPass).GetULong64_t(0) 
                if ((po_dcApplied & 0x210000000242) & po_dcFlagged) != (po_dcApplied & 0x210000000242) : # dc bitmask
                    continue
            ''' All Cuts Passed! '''

            delta_r = (po_fit_pos - bi_fit_pos).Mag()
            delta_r_cut = delta_r < bipo_delta_r_max
            delta_t = po_ev_time - bi_ev_time
            delta_t_cut = delta_t > bipo_delta_t_min and delta_t < bipo_delta_t_max # coincidence decay time + some flat bkg

            if delta_r_cut and delta_t_cut : 
                is_bipo214_beta_ev = True # else stays False

            if is_bipo214_beta_ev == True :
                fill_po_ev_OK = fill_ev_time_res(po_ev, True, time_res_cut) # if doesn't open, stays false
            
            if is_bipo214_beta_ev == True or delta_t > bipo_delta_t_max: # quit search if found coincidence or if didn't but delta t is too big
                quit_search = True
                break # no need to check the rest of the file (break inner loop) the iev loop in a ds entry
               
        if quit_search == True : break
    
    is_bipo214_beta_ev = fill_bi_ev_OK and fill_po_ev_OK # need to pass all cuts to get to fill_OK, need to fill okay to be bipo coincidence

    return is_bipo214_beta_ev # final bool: answers question "is it a bipo214 beta event?"

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def fill_ev_time_res(ev, fv_cut, time_res_cut, fit_pos):

    '''

    Description
    -----------
    Apply FV cut and fitted time residual cut if their respective parameters are True. 
    Determine whether there are any valid time residualls.

    Parameters
    ----------
    ev : rat.RAT::DS::EV
        Input event to apply cuts to
    fv_cut : bool
        Whether or not to apply a fiducial volume cut 
    time_res_cut : bool
        Whether or not to apply a time residual cut
    fit_pos : ROOT.TVector
        Fitted position of the event ev.
    
    Returns
    -------
    fill_ok : bool
        Whether there were any valid fitted time residuals

    '''
    
    # rat.utility().GetLightPath() must be called *after* the rat.dsreader constructor.
    lightPath = rat.utility().GetLightPathCalculator()
    groupVelocity = rat.utility().GetGroupVelocity()
    pmtInfo = rat.utility().GetPMTInfo()
    
    timeres_list = [] # time residual list 
    fill_ok = False #stays False if no valid time residuals

    fit_time = ev.GetFitResult(fitName).GetVertex(0).GetTime()
    pos_in_fv = fit_pos.Mag() < r_cut

    if fv_cut == True and pos_in_fv == False : return False # do not process this event if fv cut was enabled and its position is not within the fv
    
    calibratedPMTs = ev.GetCalPMTs()
    for iPMT in range(calibratedPMTs.GetCount()) : # if breaks, check here 

        pmtCal = calibratedPMTs.GetPMT( iPMT )
        pmt_pos = pmtInfo.GetPosition(pmtCal.GetID())
        lightPath.CalcByPosition(fit_pos , pmt_pos)
        distInInnerAV = lightPath.GetDistInInnerAV()
        distInAV = lightPath.GetDistInAV()
        distInWater = lightPath.GetDistInWater()
        transitTime = groupVelocity.CalcByDistance( distInInnerAV, distInAV, distInWater )
        fit_time_res = pmtCal.GetTime() - transitTime - fit_time
        timeres_inwindow = (fit_time_res > tres_cut_min) and (fit_time_res < tres_cut_max)
        if time_res_cut == True and timeres_inwindow == False: continue # do not fill pmt time res outside of the provided window 
        timeres_list.append(fit_time_res) # append time residuals

    if timeres_list : fill_ok = True # if list is not empty

    return fill_ok 

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////