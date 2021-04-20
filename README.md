# cmip6_plot_check

Usage example:
conda activate cmip6_plot_check
python cmip6_plot_check.py /Users/afridlin/data/CMIP6/CMIP/NASA-GISS/GISS-E3-G/amip /Users/afridlin/data/CMIP6/CMIP/NASA-GISS/GISS-E3-G2/amip

Output:
./GISS-E3-G_amip_vs_GISS-E3-G2_amip.pdf

Images are generated as follows:
- plot first sample from all data files in the first directory (recursively)
    - if matching data file exists in second directory, plot that alongside it
- plot first sample from any data files that exist in the second directory that do not appear in the first directory
