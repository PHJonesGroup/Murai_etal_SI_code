# p53wt_MLE
Maximum likelihood estimate (MLE) calculation of the single-progenitor model parameters for WT epidermis.
Run MLE-calculator-WT-epidermis.m to calculate the MLE values based on the experimental data on clone sizes.

### Dependencies
- TotalCloneSizes-raw-data.mat : contains the experimental data on total clone sizes
- gillespie-EPC-total-paramest.m : runs the Gillespie algorithm
- logLike-calc : calculates the Log-likelihood
- pruning-outlierClones : excludes extremely large outlier clones from the analysis
- size2freq : converts a set of clone sizes into clone frequencies
- size2freqbinned : converts a set of clone sizes into clone frequencies (with sizes binned in powers of 2)

### Requirements
Matlab R2016b