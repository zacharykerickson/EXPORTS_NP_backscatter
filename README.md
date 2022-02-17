# EXPORTS_NP_backscatter
Code on EXPORTS North Pacific backscatter alignment, in support of publication Erickson et al., 2022 in Elementa.

## Steps to generate all of the analyses used in this paper:

1. Download data from SeaBASS. There is not, to my knowledge, a convenient way to download subsets of data, so the easiest thing to do is just download everything, which will be very large and require multiple .tar downloads, and then merge the complicated folder structures together. If you are on a Mac, the Finder application merges folders together well. The filenames used here are given in the `config.py` script (the `data_files` variable), and are also reproduced below.

2. Define locations for the data in the `config.py` file. This is the `data_folder` variable and the `data_files` dictionary. Make sure they point to the correct places on your computer.

3. Ensure you have all appropriate Python packages downloaded. Everything used here is easily downloaded from 'conda' or 'pip'. The exact environment I used is included here as `zke_elem22.yml`.

4. Run the `run_comparisons.sh` bash script. By default, this will print statistical results to the terminal window. You can port this into its own file. To duplicate exactly what I did, in Terminal type: `bash run_comparisons.sh >> out.log`, and you should be able to exactly recreate the 'out.log' file provided here.

5. Depending on what is in the `config.py file`, this will plot figures (if `PLOT_FIG = True`), save figures (if `SAVE_FIG = True`), output info (if `OUTPUT_INFO = True`), and save the alignment statistics  (if `SAVE_INFO = True`) in `comparisons.csv` (or whatever `SAVE_INFO_FN` is set to). These will give all of the information found in Table 2 and Figures 3, S1--S28.

6. Run the `run_aligned_comparisons.sh` bash script. This will go into the `comparisons.csv` file, make recommended alignments, and then re-do the comparison. By default, this will print statistical results to the terminal window. You can port this into its own file. To duplicate exactly what I did, in Terminal type: `bash run_aligned_comparisons.sh >> out_aligned.log`, and you should be able to exactly recreate the `out_aligned.log` file provided here.

7. This will similarly create figures and a .csv file, by default `comparisons_aligned.csv`, which will give all of the information found in Figure 4.

## Data files
- FLBB-RR: R2 dataset, filenames matching 'EXPORTS_EXPORTSNP_CTDbinned_rosette_process\_\*\_R2.sb', submitted by UCSB/CRSEO.  
- FLBB-SR: R2 dataset, filenames matching 'EXPORTS_EXPORTSNP_CTDbinned_rosette_survey\_\*\_R2.sb', submitted by UCSB/CRSEA.  
- BB9-SR: R0 dataset, filenames matching 'EXPORTSNP_BB9\_\*\_R0.sb', submitted by NASA_GSFC.  
- BB9-RR: R1 dataset, filenames matching 'EXPORTS-EXPORTSNP_ECO-BB9_process\_\*\_bbp_R1.sb', submitted by MAINE/boss.  
- HS6-RR: R1 dataset, filenames matching 'EXPORTS-EXPORTSNP_HS6_process\_\*\_bbp_R1.sb', submitted by MAINE/boss.
- FLBB-LF: R1 dataset, filename EXPORTS-EXPORTSNP_bb_Seabird_float_20180814_R1.sb, submitted by UWASH/dasaro.  
- BBFL2-WW: R2 dataset, filenames matching 'EXPORTS-EXPORTSNP_Wirewalker\_\*\_R2.sb', submitted by URI/omand.  
- BB2FL-SG: netCDF file (level 2) submitted by UWASH/clee as an 'associated' file, filename 'sg219_EXPORTS_Jul18_level2.nc'.  
- MCOMS-BGC: must be downloaded separately as a netCDF, float ID is 5905988, filename is '5905988_Sprof.nc'.

# Table of Contents

## Bash scripts

`run_comparisons.sh` : runs through all possible comparisons of instruments using `compare_inst.py`. Sample usage: `bash run_comparisons.sh` 

`run_aligned_comparisons.sh` : runs through all possible comparisons of instruments using `compare_aligned_inst.py`. Sample usage: `bash run_aligned_comparison.sh`

## Python scripts (main)

`compare_inst.py`: compare two sensors (i.e., `python compare_inst.py FLBBRR FLBBSR`).

`compare_aligned_isnt.py`: compre two sensors, using the values in `comparisons.csv` to decide how to align each of them to the reference sensor (i.e., `python compare_aligned_inst.py FLBBRR FLBBSR`)

## Python scripts (functions)

`get_bestfit.py`: includes functions to find nearest profiles between two platforms, interpolate data with respect to depth and density, and fit a line through data using orthogonal distance regression (ODR).

`get_data.py`: code to load data from each of the instruments (after first saving SeaBASS files onto computer)

## Python scripts (deep functions)

`seawater_scattering.py` was originally written by Xiaodong Zhang (xiaodong.zhang@ums.edu) and translated into Python by me while I was a postdoc at NASA/GSFC. It is used in various Python scripts here.

`SB_support.py` was written by Joel Scott (NASA/GSFC) and is used to load information from SeaBASS files ('.sb').

## Output files

`comparisons.csv`: Statistics for all comparisons, as given in Table 2 and Figure 3.

`comparisons_aligned.csv`: Statistics for all comparisons after the alignments as suggested in this paper, as given in Figure 4.

## Log files

`out.log`: Output of the `compare_inst.py` function (normally would be printed to the screen, but can be shunted to this file).

`out_aligned.log`: Output of the `compare_aligned_inst.py` function (normally would be printed to the screen, but can be shunted to this file).



