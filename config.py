import numpy as np
import get_data

wv = 700
unc_sw = 0.0224
DIST_THRES = 5 # km
TIME_THRES = 6/24 # days
DTEMP_THRES = 0.5 # deg C
DSAL_THRES = 0.1 # psu

## Preliminary stuff to control output of this code
PLOT_FIG = True
SAVE_FIG = True
OUTPUT_INFO = True
SAVE_INFO = True
SAVE_INFO_FN = 'comparisons.csv'
SAVE_ALIGNED_INFO_FN = 'comparisons_aligned.csv'

xx = [0,2e-3]

label = dict(
    FLBBRR = 'FLBB (RR)',
    FLBBSR = 'FLBB (SR)',
    BB9SR  = 'BB9 (SR)',
    HS6RR  = 'HS6 (RR)',
    BB9RR  = 'BB9 (RR)',
    FLBBLF = 'FLBB (LF)',
    BBFL2WW= 'BBFL2 (WW)',
    BB2FLSG= 'BB2FL (SG)',
    MCOMSBGC = 'MCOMS (BGC)'
)

sca = dict(
    FLBBRR = 142,
    FLBBSR = 142,
    BB9SR  = 124,
    HS6RR  = 140,
    BB9RR  = 124,
    FLBBLF = 142,
    BBFL2WW= 124,
    BB2FLSG= 124,
    MCOMSBGC = 150
)

chi = dict(
    FLBBRR = 1.14,
    FLBBSR = 1.14,
    BB9SR  = 1.04,
    HS6RR  = 1.14,
    BB9RR  = 1.04,
    FLBBLF = 1.14,
    BBFL2WW= 1.04,
    BB2FLSG= 1.04,
    MCOMSBGC = 1.16
)

unc_chi = dict(
    FLBBRR = 0.01,
    FLBBSR = 0.01,
    BB9SR  = 0.005,
    HS6RR  = 0.01,
    BB9RR  = 0.005,
    FLBBLF = 0.01,
    BBFL2WW= 0.005,
    BB2FLSG= 0.005,
    MCOMSBGC = 0.015
)

res_vsf = dict( # one count (except for HS6RR)
    FLBBRR = 5.98e-6,
    FLBBSR = 1.40e-6,
    BB9SR  = 3.99e-6, # at 717 nm
    HS6RR  = 2.28e-6, # from standard deviation on dark cast
    BB9RR  = 3.16e-6, # at 715 nm
    FLBBLF = 6.51e-6,
    BBFL2WW= 3.25e-6,
    BB2FLSG= 3.22e-6,
    MCOMSBGC = 3.27e-7
)

unc_sf = dict(
    FLBBRR = 0.020,
    FLBBSR = 0.020,
    BB9SR  = 0.020,
    HS6RR  = 0.020,
    BB9RR  = 0.022,
    FLBBLF = 0.021,
    BBFL2WW= 0.021,
    BB2FLSG= 0.020,
    MCOMSBGC = 0.020
)

unc_dk = dict(
    FLBBRR = 2.99e-6, # from standard deviation on dark cast
    FLBBSR = 1.11e-6, # from standard deviation on dark cast
    BB9SR  = 1.17e-5, # from standard deviation on dark cast (at 717 nm)
    HS6RR  = 2.28e-6, # standard deviation on dark
    BB9RR  = 2.53e-6, # from standard deviation on dark cast (at 715 nm)
    FLBBLF = 2.86e-5, # one count
    BBFL2WW= 2.93e-6, # from standard deviation on factory calibration
    BB2FLSG= 2.76e-6, # from standard deviation on factory calibration
    MCOMSBGC = 4.64e-7 # from standard deviation on factory calibration
)

max_dep = dict(
    FLBBRR = 1000,
    FLBBSR = 1000,
    BB9SR  = 200,
    HS6RR  = 200,
    BB9RR  = 200,
    FLBBLF = 200,
    BBFL2WW= 500,
    BB2FLSG= 500,
    MCOMSBGC = 1000
)

data_folder = '/Users/erickson/Documents/Data/SeaBASS_EXPORTSNP/'
data_files = dict(
    FLBBRR = '%s/UCSB/CRSEO/EXPORTS/EXPORTSNP/archive/EXPORTS_EXPORTSNP_CTDbinned_rosette_process_*_R2.sb'%data_folder,
    FLBBSR = '%s/UCSB/CRSEO/EXPORTS/EXPORTSNP/archive/EXPORTS_EXPORTSNP_CTDbinned_rosette_survey_*_R2.sb'%data_folder,
    BB9SR  = '%s/NASA_GSFC/EXPORTS/EXPORTSNP/archive/EXPORTSNP_BB9_*_R0.sb'%data_folder,
    HS6RR  = '%s/MAINE/boss/EXPORTS/exportsnp/archive/EXPORTS-EXPORTSNP_HS6_process_*_bbp_R1.sb'%data_folder,
    BB9RR  = '%s/MAINE/boss/EXPORTS/exportsnp/archive/EXPORTS-EXPORTSNP_ECO-BB9_process_*_bbp_R1.sb'%data_folder,
    FLBBLF = '%s/UWASH/dasaro/EXPORTS/EXPORTSNP/archive/EXPORTS-EXPORTSNP_bb_Seabird_float_20180814_R1.sb'%data_folder,
    BBFL2WW= '%s/URI/omand/EXPORTS/EXPORTSNP/archive/EXPORTS-EXPORTSNP_Wirewalker_*_R2.sb'%data_folder,
    BB2FLSG= '%s/UWASH/clee/EXPORTS/EXPORTSNP/associated/EXPORTS-EXPORTSNP_seaglider_L1toL3_netcdf/sg219_EXPORTS_Jul18_level2.nc'%data_folder,
    MCOMSBGC = '%s/5905988_Sprof.nc'%data_folder
)

data_function = dict(
    FLBBRR = get_data.get_FLBB_RR,
    FLBBSR = get_data.get_FLBB_SR,
    BB9SR  = get_data.get_BB9_SR,
    HS6RR  = get_data.get_HS,
    BB9RR  = get_data.get_BB9_RR,
    FLBBLF = get_data.get_FLBB_LF,
    BBFL2WW= get_data.get_WW,
    BB2FLSG= get_data.get_SG,
    MCOMSBGC = get_data.get_MCOMS_Argo
)
