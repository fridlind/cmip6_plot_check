# cmip6_plot_check

Usage example on Discover:
> module load python/GEOSpyD/Ana2019.10_py3.7  
> ipython  
> run cmip6_plot_check.py 2014-09 /discover/nobackup/aackerma/cmip6_postprocess/CMIP6/CMIP/NASA-GISS/GISS-E2-G/amip/ /discover/nobackup/aackerma/cmip6_postprocess/CMIP6/CMIP/NASA-GISS/GISS-E3-G/amip/ --first  

Output example filename format:
> ls ./GISS-E2-G_amip_vs_GISS-E3-G2_amip.pdf  

Images are generated as follows:
- if option "--first" is added, operate on the earliest version directory (otherwise default to last)
- plot YYYY-MM field from each data file in the first source (recursively)
    - if matching YYYY-MM exists in second source, plot that alongside it
- then plot first sample from any data files that exist in the second source that do not appear in the first source
