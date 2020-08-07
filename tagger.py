# tagger.py
# Aidan Patton, 07/2020

# just for running functions from tagger_functions.py
# run on an interactive job pls

print "\n" # easier on the eyes

import ROOT
import rat 
import tagger_functions
import tagger_settings

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

if __name__ == "__main__":
    tagger_functions.bipo214_comb(is_mc, data_file_arr) # bool, list of strs
    # bipo214_comb() --> build_bipo214_results --> is_bipo214_beta()   