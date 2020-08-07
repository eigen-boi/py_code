# An attempt to translate a BiPo tagger to python.
# Original code by Ana Sofia and Fady Shaker.
# Aidan Patton, 07-2020

'''
Apply cuts to isolate BiPo214 coincidence events, save to a new .root file.
'''

import ROOT
import rat 
import tagger_settings
import gc

plot_dir = tagger_settings.plot_dir # str
root_dir = tagger_settings.root_dir # str

data_file_arr = tagger_settings.data_file_arr # list of strs
mc_file_arr = tagger_settings.mc_file_arr     # list of strs
identifier = tagger_settings.identifier       # str

bi_z_min = tagger_settings.bi_z_min # [mm] # float
bi_r_min = tagger_settings.bi_r_min # [mm] # float
bi_r_max = tagger_settings.bi_r_max # [mm] # float
bi_nhit_cleaned_min = tagger_settings.bi_nhit_cleaned_min # int

po_z_min = tagger_settings.po_z_min # [mm] # float
po_r_max = tagger_settings.po_r_max # [mm] # float
po_nhit_cleaned_min = tagger_settings.po_nhit_cleaned_min # int
po_nhit_cleaned_max = tagger_settings.po_nhit_cleaned_max # int

bipo_delta_r_max = tagger_settings.bipo_delta_r_max # [mm] # float
bipo_delta_t_min = tagger_settings.bipo_delta_t_min # [ns] # float
bipo_delta_t_max = tagger_settings.bipo_delta_t_max # [ns] # float

is_mc = tagger_settings.is_mc     # bool
fitName = tagger_settings.fitName # str
counts = tagger_settings.counts   # list of ints



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def bipo214_comb(is_mc, in_files): # main function to call

    ''' 
    Description
    -----------
    Combs through files and finds all BiPo214 coincidence events in all of them, 
    saving them to an output root file (which may be pre-existing).

    Parameters
    ----------
    is_mc : bool 
        Whether the input files are MC simulated data
    in_files : list
        List of input root files to search through.
        Requires PMT-level data, so ntuples have not been made acceptable 
    
    Returns
    -------
    No returns, but return is used to exit function in case of invalid input.

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

    # just do one for now: #FIXME
    in_files = [in_files.pop(0)]
    
    fname_count = 0 # counter (starts at 1) for what file is being analysed
    nfiles = len(in_files) # total number of files 
    coincicount = [0] * nfiles # to count of number of bipo events found in each file.
    
    for fname in in_files :

        fname_count += 1

        if ".ntuple.root" in fname :
            print "\n", fname, "is an invalid input file."
            print "ROOT ntuples are not acceptable input as they do not contain PMT-level data. Skipping."
            continue # skip to next file  

        # the output file to store the bipo events, will save in root_dir with name data_bipo214_results.root
        # Create a new file + a clone of old tree in new file https://root-forum.cern.ch/t/appending-branch-of-user-class-to-existing-tree/10432/5
        
        oldFile = ROOT.TFile(fname) # READ mode by default
        if oldFile.IsZombie(): 
            print "\n", fname, "is a zombie. Skipping."
            continue # skip to next file  

        oldTree_T = oldFile.Get("T") # main tree with EV info, https://snopl.us/docs/rat/user_manual/html/data_structure.html
        oldTree_runT = oldFile.Get("runT") # required to run dsreader on the output file
        outFile = ROOT.TFile(str(root_dir + out_file_prefix + "bipo214_results_" + identifier + ".root"), "UPDATE")

        if outFile.IsZombie(): # protection against corrupted files, don't want everything to break
            identifier = identifier + identifier # will technically always work
            print "\n", str(root_dir + out_file_prefix + "bipo214_results_" + identifier + ".root"), "is a zombie."
            outFile = ROOT.TFile(str(root_dir + out_file_prefix + "bipo214_results_" + identifier + ".root"), "UPDATE")
            print "Saving events to new file,", str(root_dir + out_file_prefix + "bipo214_results_" + identifier + ".root")
            continue # skip to next file 

        newTree_T = oldTree_T.CloneTree(0) # 0 parameter ensures no entries are coppied, all branches are copied though
        newTree_runT = oldTree_runT.CloneTree() # all entries copied

        counts[0] = 0 # reset individual file event counter
        print "\n\tis_mc ==", is_mc   
        coincicount[fname_count-1] = build_bipo214_results(fname, newTree_T) # main functionality
        
        newTree_T.AutoSave() # "If the acquisition fails to complete, you can recover the file and all the 
        newTree_runT.AutoSave() # contents since the last Autosave."
        del oldFile, oldTree_T, oldTree_runT # no trace of these in the outFile (might be uneeded FIXME)
        outFile.Write() # saves the file

        print "\n\t", counts[0], "events were analysed in this file."
        print "\n\tFinal coincidence event count for file", fname_count, "/", nfiles, ":", coincicount[fname_count-1]

    print "\n\n\tTotal Bi event count across all", nfiles, "file(s):", counts[1], "events."
    print "\n\tFinal coincidence event count for all", nfiles, "file(s):", sum(coincicount), "coincidence events.\n\n"

    outFile.Close() # hopefully, the file will be safely closed. All done!

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def build_bipo214_results(fname, newTree_T):

    '''

    Description
    -----------
    Find all BiPo coincidence events in an input file and save them to newTree_T.

    Parameters
    ----------
    fname : str
        full path to RAT data to access with the dsReader.
    newTree : ROOT.TTree
        TTree in an output file to save events to.

    Returns
    -------
    coincicount : int
        number of BiPo coincidences found in the file.

    '''

    print "\n\tAnalysing file", fname, "\n"

    try : # If this fails for any reason, we want to avoid corrupting the TFile, so exit the function.
        dsRead = rat.dsreader(fname) # unfortunately requires an additional opening of the file with dsreader, but avoids loosing all progress to a bad input file.
    except :
        print fname, "is unable to be read by the dsreader. Skipping this file."
        return 0

    iEntryBi = -1 # counting entries and events 
    coincicount = 0 # count of bipo214 coincidence events

    # loop to identify the Bi candidates. This needs a RatDB connection every time it's run for data. "You could store all [events] in a list, but that's a guarantee you'll run out of memory"
    for ds0, run0 in dsRead : # ds0 is equivalent to rat.dsreader.GetEntry(i) in c++, dsReader = rat.dsreader(input_files)

        iEntryBi += 1
        if iEntryBi % 1000 == 0 : print "Entry #", iEntryBi # print every thousand entries
        if coincicount >= 2 : break # just to test, don't take hours FIXME

        for iev0 in range(0, ds0.GetEVCount()) :

            bi_ev = ds0.GetEV(iev0)

            beta_bipo214, iEntryBi = is_bipo214_beta(dsRead, bi_ev, iEntryBi)
            if beta_bipo214 == False : continue

            coincicount += 1
            print "\t", coincicount, "coincidence event(s)!\n\n"

            newTree_T.SetBranchAddress("evs", ROOT.AddressOf(bi_ev))
            newTree_T.Fill() # add bi_ev to new file to save if is a bipo event, like https://root.cern/doc/master/copytree3_8C.html

    return coincicount

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def is_bipo214_beta(dsRead, bi_ev, iEntryBi):

    '''

    Description
    -----------
    Given a single input Bi event, determine whether it is a bipo coincidence event by seeking its respective Po in the events after itself.

    Parameters
    ----------
    dsRead : rat.RAT::DU::DSReader generator
        RAT data accessed with the dsReader. Every entry touched in this loop is inaccesible in the Bi loop afterwards, 
        but this is okay since the Po is only ever a few entries away (because of delta_t_max), making it unlikely 
        (though possible) that any true Bi will be missed; even if it is missed, reloading the DB tables every time 
        a potential Po candidate is found just to cover every event is not worth the extra time/resources. 
    bi_ev : rat.RAT::DS::EV
        Bi event candidate.
    iEntryBi : int
        Entry index of bi_ev in its file.
    
    Returns
    -------
    is_bipo214_beta_ev : bool
        Whether given input event bi_ev is a BiPo214 beta event.
    iEntryPo : int
        The highet entry index reached. The main Bi for loop will continue off after this point.

    '''

    counts[0] += 1 # counter for all events that are are tested in THIS file
    counts[1] += 1 # counter for all events across ALL files 

    if not bi_ev.FitResultExists(fitName) or not bi_ev.GetFitResult(fitName).GetVertex(0).ValidPosition() or not bi_ev.GetFitResult(fitName).GetVertex(0).ValidTime() :
        return (False, iEntryBi) # must have valid fit otherwise we cannot extract the necessary information to apply cuts

    quit_search = False        # redefine upon each opening of the function, for breaking loops
    is_bipo214_beta_ev = False # False by default

    bi_fit_pos = bi_ev.GetFitResult(fitName).GetVertex(0).GetPosition()
    bi_ev_time = bi_ev.GetClockCount50()*20 # to search between different events we have to use the clk not the fitted time (fit time is within the event window (0, 400) [ns])

    # Bi position
    bi_z = bi_fit_pos.Z()
    bi_r = bi_fit_pos.Mag()

    """ Apply Cuts """
    if bi_z < bi_z_min : return (False, iEntryBi)                               # FV cut (in scintillator cap)
    if bi_r < bi_r_min or bi_r > bi_r_max : return (False, iEntryBi)            # FV cut (Bi  2m < R < 6 m)
    if bi_ev.GetNhitsCleaned() < bi_nhit_cleaned_min : return (False, iEntryBi) # Bi nhitsCleaned cut 
    if is_mc == False : # the data cleaning (dc) cuts are only used for DATA, using it for MC will create a crash!
        latestPass = bi_ev.GetDataCleaningFlags().GetLatestPass()
        bi_dcApplied = bi_ev.GetDataCleaningFlags().GetApplied(latestPass).GetULong64_t(0)
        bi_dcFlagged = bi_ev.GetDataCleaningFlags().GetFlags(latestPass).GetULong64_t(0) 
        if ((bi_dcApplied & 0x210000000242) & bi_dcFlagged) != (bi_dcApplied & 0x210000000242) : # dc bitmask 
            return (False, iEntryBi)
    """ All Bi Cuts Successful! """

    print "Entry #", iEntryBi # to make it clear when a potential Bi candidate is found
    iEntryPo = iEntryBi # Po loop will automatically start with ds1 equal to the ds0 at iEntryBi+1, set Po index accordingly 
    
    for ds1, run1 in dsRead : 
        iEntryPo += 1
        for iev1 in range(0, ds1.GetEVCount()) :

            po_ev = ds1.GetEV(iev1)

            if not po_ev.FitResultExists(fitName) or not po_ev.GetFitResult(fitName).GetVertex(0).ValidPosition() or not po_ev.GetFitResult(fitName).GetVertex(0).ValidTime() :
                continue

            po_fit_pos = po_ev.GetFitResult(fitName).GetVertex(0).GetPosition()
            po_ev_time = po_ev.GetClockCount50()*20 

            # Po position
            po_z = po_fit_pos.Z()
            po_r = po_fit_pos.Mag()

            ''' Apply Po Cuts '''
            if po_z < po_z_min : continue # Po z > 0.85 m (in scintillator cap)
            if po_r > po_r_max : continue # po  R < 6 m
            if po_ev.GetNhitsCleaned() < po_nhit_cleaned_min or po_ev.GetNhitsCleaned() > po_nhit_cleaned_max : continue # nhits cut
            if is_mc == False : # unnecessary to apply data cleaning bitmask to MC
                latestPass = po_ev.GetDataCleaningFlags().GetLatestPass()
                po_dcApplied = po_ev.GetDataCleaningFlags().GetApplied(latestPass).GetULong64_t(0)
                po_dcFlagged = po_ev.GetDataCleaningFlags().GetFlags(latestPass).GetULong64_t(0) 
                if ((po_dcApplied & 0x210000000242) & po_dcFlagged) != (po_dcApplied & 0x210000000242) : # dc bitmask
                    continue
            ''' All Po Cuts Passed! '''

            # final cuts
            delta_r = (po_fit_pos - bi_fit_pos).Mag()
            delta_r_cut = delta_r < bipo_delta_r_max
            delta_t = po_ev_time - bi_ev_time
            delta_t_cut = delta_t > bipo_delta_t_min and delta_t < bipo_delta_t_max # coincidence decay time + some flat bkg

            if delta_r_cut and delta_t_cut : 
                print "\n\n\tIs BiPo214! iEntryBi =", iEntryBi, "iEntryPo =", iEntryPo
                is_bipo214_beta_ev = True # else stays False
            
            if is_bipo214_beta_ev == True or delta_t > bipo_delta_t_max : # quit search if found coincidence 
                quit_search = True # OR if didn't but delta t is too big for any later events to be bipo
                break 
               
        if quit_search == True : break # breaks for ds1, run in dsReader

    #del dsReadPo # forceful deletion of dsReadPo to ensure no chance of keeping it in memory
    #gc.collect() # forceful garbage collection

    return (is_bipo214_beta_ev, iEntryPo) # tuple - answers question "is it a bipo214 beta event?", and highet entry index reached in search for Po

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////