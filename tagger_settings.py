# tagger_settings.py

# (re)define global variables
plot_dir = "/home/eigenboi/pdfs/"
root_dir = "/home/eigenboi/scratch/root_files/"
#in_data_file = "/home/eigenboi/scratch/data/testData.root" # one subrun, ~5mins of data
#data_file_arr =["/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257301_s000_p000.root", "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257532_s011_p000.root", "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257953_s011_p001.root",
#                "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257301_s001_p000.root",  "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257533_s000_p000.root",  "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257954_s000_p000.root",
#                "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257301_s002_p000.root",  "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257533_s001_p000.root",  "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257954_s001_p000.root",
#                "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257301_s003_p000.root",  "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257533_s002_p000.root",  "/scratch/djauty/data/full/Jan/6_18_0/Analysis40_r0000257954_s002_p000.root"]
# Physics, Bubblers ON, Cavity Recirculation ON, List(s) None
#data_file_arr = ["/home/eigenboi/scratch/data/258988.root", "/home/eigenboi/scratch/data/258988.root", "/home/eigenboi/scratch/data/testData.root"]
data_file_arr = ["/home/eigenboi/scratch/data/testData.root"]
mc_file_arr = []
#in_data_files = "/home/eigenboi/scratch/data/*.root"
identifier = "test_of_concept"
#identifier = "gold_line_00"

r_cut = 3500          # [mm]
tres_cut_min = -5.0   # [ns]
tres_cut_max = 15.0   # [ns]

bi_z_min = 1400 # [mm]
bi_r_min = 2000 # [mm]
bi_r_max = 6000 # [mm] # TO DO change this to 5000 and see 
bi_nhit_cleaned_min = 250

po_z_min = 850  # [mm]
po_r_max = 6000 # [mm] 
po_nhit_cleaned_min = 150
po_nhit_cleaned_max = 350

bipo_delta_r_max = 1000    # [mm]
bipo_delta_t_min = 3690    # [ns]
bipo_delta_t_max = 1800000 # [ns]

is_mc = False
time_res_cut = False
fitName = "partialFitter" # for both data and mc 