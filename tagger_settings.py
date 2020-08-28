# tagger_settings.py

# define global variables
root_dir = ""
ev_dir = ""
data_dir = ""

data_file_arr = []
new_data_file_arr = []
mc_file_arr = []               
identifier = ""

bi_z_min = 1400.0            # [mm] # float
bi_r_min = 2000.0            # [mm] # float
bi_r_max = 6000.0            # [mm] # float
bi_nhit_cleaned_min = 250    # int

po_z_min = 850.0             # [mm] # float
po_r_max = 6000.0            # [mm] # float
po_nhit_cleaned_min = 150    # int
po_nhit_cleaned_max = 350    # int

bipo_delta_r_max = 1000.0    # [mm]
bipo_delta_t_min = 3690.0    # [ns]
bipo_delta_t_max = 1800000.0 # [ns]

is_mc = True
fitName = "partialFitter" 
counts = [0, 0] # count in single file (most recent), count across all files