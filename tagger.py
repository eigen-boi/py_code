# just for running functions from tagger_functions.py
# run on an interactive job pls
print "\n" # easier on the eyes
import ROOT
import rat 
import tagger_functions
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

if __name__ == "__main__":
    tagger_functions.bipo214_comb(is_mc, time_res_cut, data_file_arr) # bool, bool, list
    # bipo214_comb() calls build_bipo214_results --> is_bipo214_beta() --> fill_ev_time_res()  