# EXPORTS_NP_backscatter
Code on EXPORTS North Pacific backscatter alignment, in support of publication Erickson et al., 2021 in Elementa.

# Table of Contents

`get_bestfit.py`: includes functions to find nearest profiles between two platforms, interpolate data with respect to depth and density, and fit a line through data using orthogonal distance regression (ODR).

`get_data.py`: code to load data from each of the instruments (after first saving SeaBASS files onto computer)

`prop_error.py`: propagates error to calculate best-fit lines between parameters after they are aligned.

`bbp_plot.ipynb`: construct Figure 2 of paper

`compare_XXXxx_YYYyy.ipynb`: compare sensor XXX on platform xx with sensor YYY on platform yy.
