import os
import sys
import argparse
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob2 as glob
from matplotlib.backends.backend_pdf import PdfPages as pdf

def main(inargs):
    """Run the program."""
    
    # specify source files from a model run, and year and month to plot

    dir_1 = inargs.dir_1
    tim_1 = inargs.mon_1
    if len(tim_1) != 7:
        print('INPUT ERROR: Date(s) must be in YYYY-MM format.')
        sys.exit() # abort    
    print('First source: ', dir_1)
    mod_1 = dir_1.split('/')[-4]
    run_1 = dir_1.split('/')[-3]
    sce_1 = dir_1.split('/')[-2]
    dir_1_var = glob.glob(dir_1+'*/*/*', recursive=True) # identify all reported variables
    if dir_1_var == []:
        print('INPUT ERROR: Directory specification error (no variables found).')
        print(dir_1)
        sys.exit()
    if mod_1[0:7] == 'GISS-E2': # skip *fx* and *fy* variables
        dir_1_var = [i for i in dir_1_var if 'fx' not in i]
        dir_1_var = [i for i in dir_1_var if 'fy' not in i]
        print('WARNING: Skipping fx class from E2 owing to dimensionality error.')
    if inargs.include != None:
        dir_1_var_incl = []
        for i_1, d_1 in enumerate(inargs.include.split(',')): # iterate over comma-separated list
            var_incl = [i for i in dir_1_var if d_1 in i] # include only specified variable(s)
            dir_1_var_incl.extend(var_incl)
        dir_1_var = list(dir_1_var_incl)
        if dir_1_var == []:
            print('DATA ERROR: First directory output is missing specified variable(s).')
            print(inargs.include)
            sys.exit()
    elif inargs.exclude != None:
        for i_1, d_1 in enumerate(inargs.exclude.split(',')): # iterate over comma-separated list
            dir_1_var = [i for i in dir_1_var if d_1 not in i] # exclude specified variable(s)
        if dir_1_var == []:
            print('DATA ERROR: First directory output is empty beyond excluded variable(s).')
            print(inargs.exclude)
            sys.exit()
    ncs_1 = []
    for i_1, d_1 in enumerate(dir_1_var):
        v_all = sorted(glob.glob(dir_1_var[i_1]+'/*', recursive=True)) # all versions
        if v_all != []:
            if inargs.first: # optionally use first version in the output
                f_all = sorted(glob.glob(v_all[0]+'/*.nc')) # all files in first version
            else: # otherwise use default last (most recent version)
                f_all = sorted(glob.glob(v_all[-1]+'/*.nc'))
            f_tim = [i for i in f_all if '2014' in i]
            if f_tim == []:
                f_tim = [i for i in f_all if tim_1.replace('-','') in i]
            if f_tim != []:
                ncs_1.append(f_tim[0])
    try: # make sure that year and month output exist
        dat_1 = xr.open_dataset(ncs_1[0]).sel(time=tim_1)
    except:
        print('DATA ERROR: First directory output is missing specified year and month.')
        print(ncs_1)
        sys.exit() # abort if date and month are not found

    # optionally specify source files from a second model run

    if inargs.compare:
        dir_2 = inargs.dir_2
        if inargs.mon_2 != None:
            tim_2 = inargs.mon_2 # optional different year and month from second run
        else:
            tim_2 = tim_1 # or same year and month from second run
        print('Second source: ', dir_2, tim_2)
        mod_2 = dir_2.split('/')[-4]
        run_2 = dir_2.split('/')[-3]
        sce_2 = dir_2.split('/')[-2]
        dir_2_var = glob.glob(dir_2+'*/*/*', recursive=True)
        if dir_2_var == []:
            print('INPUT ERROR: Directory specification error (no variables found).')
            print(dir_2)
            sys.exit()
        if mod_2[0:7] == 'GISS-E2': # skip *fx* and *fy* variables
            dir_2_var = [i for i in dir_2_var if 'fx' not in i]
            dir_2_var = [i for i in dir_2_var if 'fy' not in i]
            print('WARNING: Skipping fx class from E2 owing to dimensionality error.')
        dir_2_var = [i for i in dir_2_var if 'fx' not in i]
        if inargs.include != None:
            dir_2_var_incl = []
            for i_2, d_2 in enumerate(inargs.include.split(',')): # iterate over comma-separated list
                var_incl = [i for i in dir_2_var if d_2 in i] # include only specified variable(s)
                dir_2_var_incl.extend(var_incl)
            dir_2_var = list(dir_2_var_incl)
            if dir_2_var == []:
                print('DATA ERROR: Second directory output is missing specified variable(s).')
                print(inargs.include)
                sys.exit()
        elif inargs.exclude != None:
            for i_2, d_2 in enumerate(inargs.exclude.split(',')): # iterate over comma-separated list
                dir_2_var = [i for i in dir_2_var if d_2 not in i] # exclude specified variable(s)
            if dir_2_var == []:
                print('DATA ERROR: Second directory output is empty beyond excluded variable(s).')
                print(inargs.exclude)
                sys.exit()
        ncs_2 = []
        for i_2, d_2 in enumerate(dir_2_var):
            v_all = sorted(glob.glob(dir_2_var[i_2]+'/*', recursive=True))
            if v_all != []:
                if inargs.first:
                    f_all = sorted(glob.glob(v_all[0]+'/*.nc'))
                else:
                    f_all = sorted(glob.glob(v_all[-1]+'/*.nc'))
                f_tim = [i for i in f_all if '2014' in i]
                if f_tim == []:
                    f_tim = [i for i in f_all if tim_1.replace('-','') in i]
                if f_tim != []:
                    ncs_2.append(f_tim[0])
        try:
            dat_2 = xr.open_dataset(ncs_2[0]).sel(time=tim_2)
        except:
            print('DATA ERROR: Second directory output is missing specified year and month.')
            print(ncs_2[0])
            sys.exit()
                    
    # specify local destination for output comparison plots

    if inargs.compare:
        out_pdf = mod_1+'_'+run_1+'_'+sce_1+'_vs_'+mod_2+'_'+run_2+'_'+sce_2+'.pdf'
    else:
        out_pdf = mod_1+'_'+run_1+'_'+sce_1+'.pdf'

    # loop over source files in first model run

    pp = pdf('multipage.pdf') # initialize multipage package to receive sequential images

    print('Processing first source ...')
    for i_1, f_1 in enumerate(ncs_1): # loop over variables identified above
        print(f_1)
        dat_1 = xr.open_dataset(f_1)
        var_1 = list(dat_1.data_vars.keys())[-1] # identify the variable name
        ndims = len(dat_1[var_1].dims) # determine dimensionality
        if ndims!=2: dat_1 = xr.open_dataset(f_1).sel(time=tim_1) # most fields have time
        if ndims==1: # data is a scalar (dummy plot)
            fig = plt.figure(figsize=[8.5,11])
            ax = fig.add_subplot(211)
            fld_1 = dat_1[var_1].isel(time=0) # scalar value
            ax.annotate('SCALAR VALUE',xy=(0.4,0.5),xycoords='axes fraction')
            path, fname = os.path.split(f_1)
            parr = path.split(mod_1)
            title = parr[0]+mod_1+'\n'+parr[1]+'/\n'+fname+'\n'+fld_1.attrs['long_name']
            # parse directory name for title
            plt.title(title)
            val_str = ("value = "+"{:.5e}".format(fld_1.data))
            ax.annotate(tim_1+' '+val_str,xy=(0,-0.15),xycoords='axes fraction')

            if inargs.compare: # optionally search for matching variable in second directory
                search_str = '/'+fname.split('_')[0]+'_'+fname.split('_')[1]
                matching_file = [i for i in ncs_2 if search_str in i]
            else: matching_file = []
            if matching_file != []:
                f_2 = matching_file[0]
                print(f_2)
                dat_2 = xr.open_dataset(f_2)
                var_2 = list(dat_2.data_vars.keys())[-1]
                ax = fig.add_subplot(212)
                fld_2 = dat_2[var_2].isel(time=0)
                ax.annotate('SCALAR VALUE',xy=(0.4,0.5),xycoords='axes fraction')
                path, fname = os.path.split(f_2)
                parr = path.split(mod_2)
                title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']
                plt.title(title)
                val_str = ("value = "+"{:.5e}".format(fld_2.data))
                ax.annotate(tim_2+' '+val_str,xy=(0,-0.15),xycoords='axes fraction')
                fig.tight_layout(pad=6)
                
            pp.savefig() # completed page
        elif ndims==2:
            fig = plt.figure(figsize=[8.5,11])
            ax = fig.add_subplot(211,projection=ccrs.PlateCarree(central_longitude=180))
            fld_1 = dat_1[var_1]
            fld_1.plot(ax=ax,transform=ccrs.PlateCarree(),
                        cbar_kwargs={'label': fld_1.units},rasterized=True)
            ax.coastlines()
            path, fname = os.path.split(f_1)
            parr = path.split(mod_1)
            title = parr[0]+mod_1+'\n'+parr[1]+'/\n'+fname+'\n'+fld_1.attrs['long_name']
            plt.title(title)
            val_str = ("min, max, avg = "+"{:.5e}".format(fld_1.min().data)+", "
                "{:.5e}".format(fld_1.max().data)+", "+"{:.5e}".format(fld_1.mean().data))
            ax.annotate(tim_1+' '+val_str,xy=(0,-0.25),xycoords='axes fraction')

            if inargs.compare: # optionally search for matching variable in second directory
                search_str = '/'+fname.split('_')[0]+'_'+fname.split('_')[1]
                matching_file = [i for i in ncs_2 if search_str in i]
            else: matching_file = []
            if matching_file != []:
                f_2 = matching_file[0]
                print(f_2)
                dat_2 = xr.open_dataset(f_2)
                var_2 = list(dat_2.data_vars.keys())[-1]
                ax = fig.add_subplot(212,projection=ccrs.PlateCarree(central_longitude=180))
                fld_2 = dat_2[var_2]
                fld_2.plot(ax=ax,transform=ccrs.PlateCarree(),
                           cbar_kwargs={'label': fld_2.units},rasterized=True)
                ax.coastlines()
                path, fname = os.path.split(f_2)
                parr = path.split(mod_2)
                title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']
                plt.title(title)
                val_str = ("min, max, mean = "+"{:.5e}".format(fld_2.min().data)+", "
                    "{:.5e}".format(fld_2.max().data)+", "+"{:.5e}".format(fld_2.mean().data))
                ax.annotate(tim_2+' '+val_str,xy=(0,-0.25),xycoords='axes fraction')
                fig.tight_layout(pad=6)

            pp.savefig()
        elif ndims==3: # data is lat/lon (simplest case to plot)
            fig = plt.figure(figsize=[8.5,11]) # initialize letter-size page         
            if dat_1[var_1].dims[1]=='basin':
                # initialize top subplot with line plot
                ax = fig.add_subplot(211)
                fld_1 = dat_1[var_1].isel(basin=0,time=0) # data to plot
                subtit = ' (basin=0)'
                fld_1.plot(ax=ax)
            else:
                # initialize top subplot with a mapping projection
                ax = fig.add_subplot(211,projection=ccrs.PlateCarree(central_longitude=180))
                fld_1 = dat_1[var_1].isel(time=0) # data to plot
                subtit = ''
                # plot on specified projection with default color bar, rasterize to reduce file size
                fld_1.plot(ax=ax,transform=ccrs.PlateCarree(),
                           cbar_kwargs={'label': fld_1.units},rasterized=True)
                ax.coastlines()
            # parse directory name for title
            path, fname = os.path.split(f_1)
            parr = path.split(mod_1)
            title = parr[0]+mod_1+'\n'+parr[1]+'/\n'+fname+'\n'+fld_1.attrs['long_name']+subtit
            plt.title(title)
            # calculate statistics and report below figure
            val_str = ("min, max, mean = "+"{:.5e}".format(fld_1.min().data)+", "
                       "{:.5e}".format(fld_1.max().data)+", "+"{:.5e}".format(fld_1.mean().data))
            ax.annotate(tim_1+' '+val_str,xy=(0,-0.25),xycoords='axes fraction')

            if inargs.compare: # optionally search for matching variable in second directory
                search_str = '/'+fname.split('_')[0]+'_'+fname.split('_')[1]
                matching_file = [i for i in ncs_2 if search_str in i]
            else: matching_file = []
            if matching_file != []: # if it exists, execute same procedure for matching data
                f_2 = matching_file[0]
                print(f_2)                
                dat_2 = xr.open_dataset(f_2).sel(time=tim_2)
                var_2 = list(dat_2.data_vars.keys())[-1]
                if dat_2[var_2].dims[1]=='basin':
                    ax = fig.add_subplot(212)
                    fld_2 = dat_2[var_2].isel(basin=0,time=0)
                    subtit = ' (basin=0)'
                    fld_2.plot(ax=ax)
                else:
                    ax = fig.add_subplot(212,projection=ccrs.PlateCarree(central_longitude=180))
                    fld_2 = dat_2[var_2].isel(time=0)
                    subtit = ''
                    fld_2.plot(ax=ax,transform=ccrs.PlateCarree(),
                               cbar_kwargs={'label': fld_2.units},rasterized=True)
                    ax.coastlines()
                path, fname = os.path.split(f_2)
                parr = path.split(mod_2)
                title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']+subtit
                plt.title(title)
                val_str = ("min, max, mean = "+"{:.5e}".format(fld_2.min().data)+", "
                           "{:.5e}".format(fld_2.max().data)+", "+"{:.5e}".format(fld_2.mean().data))
                ax.annotate(tim_2+' '+val_str,xy=(0,-0.25),xycoords='axes fraction')
                fig.tight_layout(pad=6)

            pp.savefig() # completed page
        elif ndims==4: # narrow down to either one basin or longitude for plotting
            fig = plt.figure(figsize=[8.5,11])
            ax = fig.add_subplot(211)
            if dat_1[var_1].dims[1]=='basin':
                fld_1 = dat_1[var_1].isel(basin=0,time=0)
                subtit = ' (basin=0)'
            else:
                fld_1 = dat_1[var_1].isel(lon=0,time=0)
                subtit = ' (lon=0)'
            fld_1.plot(ax=ax,cbar_kwargs={'label': fld_1.units},rasterized=True)
            if dat_1[var_1].dims[1]==('lev') or dat_1[var_1].dims[1]==('plev'): ax.invert_yaxis()
            if dat_1[var_1].dims[2]==('lev'): ax.invert_yaxis() # ocean basin case
            path, fname = os.path.split(f_1)
            parr = path.split(mod_1)
            title = parr[0]+mod_1+'\n'+parr[1]+'/\n'+fname+'\n'+fld_1.attrs['long_name']+subtit
            plt.title(title)
            val_str = ("min, max, mean = "+"{:.5e}".format(fld_1.min().data)+", "
                       "{:.5e}".format(fld_1.max().data)+", "+"{:.5e}".format(fld_1.mean().data))
            ax.annotate(tim_1+' '+val_str,xy=(0,-0.25),xycoords='axes fraction')

            if inargs.compare:
                search_str = '/'+fname.split('_')[0]+'_'+fname.split('_')[1]
                matching_file = [i for i in ncs_2 if search_str in i]
            else: matching_file = []
            if matching_file != []:
                f_2 = matching_file[0]
                print(f_2)
                dat_2 = xr.open_dataset(f_2).sel(time=tim_2)
                var_2 = list(dat_2.data_vars.keys())[-1]
                ax = fig.add_subplot(212)
                if dat_2[var_2].dims[1]=='basin':
                    fld_2 = dat_2[var_2].isel(basin=0,time=0)
                    subtit = ' (basin=0)'
                else:
                    fld_2 = dat_2[var_2].isel(lon=0,time=0)
                    subtit = ' (lon=0)'
                fld_2.plot(ax=ax,cbar_kwargs={'label': fld_2.units},rasterized=True)
                if dat_2[var_2].dims[1]==('lev') or dat_2[var_2].dims[1]==('plev'): ax.invert_yaxis()
                if dat_2[var_2].dims[2]==('lev'): ax.invert_yaxis() # ocean basin case
                path, fname = os.path.split(f_2)
                parr = path.split(mod_2)
                title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']+subtit
                val_str = ("min, max, mean = "+"{:.5e}".format(fld_2.min().data)+", "
                           "{:.5e}".format(fld_2.max().data)+", "+"{:.5e}".format(fld_2.mean().data))
                ax.annotate(tim_2+' '+val_str,xy=(0,-0.25),xycoords='axes fraction')
                plt.title(title)
                fig.tight_layout(pad=6)
                
            pp.savefig()
        else: # more than 4 dimensions: also choose a latitude
            fig = plt.figure(figsize=[8.5,11])
            ax = fig.add_subplot(211)
            fld_1 = dat_1[var_1].isel(lat=0,lon=0,time=0)
            subtit = ' (Lat/Lon=0/0)'
            fld_1.plot(ax=ax,cbar_kwargs={'label': fld_1.units},rasterized=True)
            path, fname = os.path.split(f_1)
            parr = path.split(mod_1)
            title = parr[0]+mod_1+'\n'+parr[1]+'/\n'+fname+'\n'+fld_1.attrs['long_name']+subtit
            plt.title(title)
            val_str = ("min, max, mean = "+"{:.5e}".format(fld_1.min().data)+", "
                "{:.5e}".format(fld_1.max().data)+", "+"{:.5e}".format(fld_1.mean().data))
            ax.annotate(tim_1+' '+val_str,xy=(0,-0.2),xycoords='axes fraction')

            if inargs.compare:
                search_str = '/'+fname.split('_')[0]+'_'+fname.split('_')[1]
                matching_file = [i for i in ncs_2 if search_str in i]
            else: matching_file = []
            if matching_file != []:
                f_2 = matching_file[0]
                print(f_2)
                dat_2 = xr.open_dataset(f_2)
                var_2 = list(dat_2.data_vars.keys())[-1]
                ax = fig.add_subplot(212)
                fld_2 = dat_1[var_2].isel(lat=0,lon=0,time=0)
                subtit = ' (Lat/Lon=0/0)'
                fld_2.plot(ax=ax,cbar_kwargs={'label': fld_2.units},rasterized=True)
                path, fname = os.path.split(f_2)
                parr = path.split(mod_2)
                title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']+subtit
                val_str = ("min, max, mean = "+"{:.5e}".format(fld_2.min().data)+", "
                    "{:.5e}".format(fld_2.max().data)+", "+"{:.5e}".format(fld_2.mean().data))
                ax.annotate(tim_2+' '+val_str,xy=(0,-0.2),xycoords='axes fraction')
                plt.title(title)
                fig.tight_layout(pad=6)
                
            pp.savefig()
        
        plt.close() # clear matplotlib for next page (to avoid overflows)

    # loop over source files in second model run (plot only any missing from first run)

    if inargs.compare: 
        print('Processing second source ...')
        for i_2, f_2 in enumerate(ncs_2):
            path, fname = os.path.split(f_2)
            matching_file = [i for i in ncs_1 if fname.split('_')[0]+'_'+fname.split('_')[1] in i]
            if matching_file == []:
                print(f_2)
                dat_2 = xr.open_dataset(f_2)
                var_2 = list(dat_2.data_vars.keys())[-1]
                ndims = len(dat_2[var_2].dims)
                if ndims!=2: dat_2 = xr.open_dataset(f_2).sel(time=tim_2) # usually time is a dimension
                if ndims==1:
                    fig = plt.figure(figsize=[8.5,11])
                    ax = fig.add_subplot(212)
                    fld_2 = dat_2[var_2].isel(time=0)
                    ax.annotate('SCALAR VALUE',xy=(0.4,0.5),xycoords='axes fraction')
                    path, fname = os.path.split(f_2)
                    parr = path.split(mod_2)
                    title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']
                    plt.title(title)
                    val_str = ("value = "+"{:.5e}".format(fld_2.data))
                    ax.annotate(tim_2+' '+val_str,xy=(0,-0.15),xycoords='axes fraction')
                    pp.savefig()
                elif ndims==2:
                    fig = plt.figure(figsize=[8.5,11])
                    ax = fig.add_subplot(212,projection=ccrs.PlateCarree(central_longitude=180))
                    fld_2 = dat_2[var_2]
                    fld_2.plot(ax=ax,transform=ccrs.PlateCarree(),
                               cbar_kwargs={'label': fld_2.units},rasterized=True)
                    ax.coastlines()
                    path, fname = os.path.split(f_2)
                    parr = path.split(mod_2)
                    title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']
                    plt.title(title)
                    val_str = ("min, max, mean = "+"{:.5e}".format(fld_2.min().data)+", "
                        "{:.5e}".format(fld_2.max().data)+", "+"{:.5e}".format(fld_2.mean().data))
                    ax.annotate(tim_2+' '+val_str,xy=(0,-0.25),xycoords='axes fraction')
                    plt.title(title)
                    pp.savefig()
                elif ndims==3:
                    fig = plt.figure(figsize=[8.5,11])
                    if dat_2[var_2].dims[1]=='basin':
                        ax = fig.add_subplot(212)
                        fld_2 = dat_2[var_2].isel(basin=0,time=0)
                        subtit = ' (basin=0)'
                        fld_2.plot(ax=ax)
                    else:
                        ax = fig.add_subplot(211,projection=ccrs.PlateCarree(central_longitude=180))
                        fld_2 = dat_2[var_2].isel(time=0)
                        subtit = ''
                        fld_2.plot(ax=ax,transform=ccrs.PlateCarree(),
                                   cbar_kwargs={'label': fld_2.units},rasterized=True)
                        ax.coastlines()
                    path, fname = os.path.split(f_2)
                    parr = path.split(mod_2)
                    title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']+subtit
                    plt.title(title)
                    val_str = ("min, max, mean = "+"{:.5e}".format(fld_2.min().data)+", "
                            "{:.5e}".format(fld_2.max().data)+", "+"{:.5e}".format(fld_2.mean().data))
                    ax.annotate(tim_2+' '+val_str,xy=(0,-0.25),xycoords='axes fraction')
                    pp.savefig()
                elif ndims==4:
                    fig = plt.figure(figsize=[8.5,11])
                    ax = fig.add_subplot(212)
                    if dat_2[var_2].dims[1]=='basin':
                        fld_2 = dat_2[var_2].isel(basin=0,time=0)
                        subtit = ' (basin=0)'
                    else:
                        fld_2 = dat_2[var_2].isel(lon=0,time=0)
                        subtit = ' (lon=0)'
                    fld_2.plot(ax=ax,cbar_kwargs={'label': fld_2.units},rasterized=True)
                    if dat_2[var_2].dims[1]==('lev') or dat_2[var_2].dims[1]==('plev'): ax.invert_yaxis()
                    if dat_2[var_2].dims[2]==('lev'): ax.invert_yaxis() # ocean basin case
                    path, fname = os.path.split(f_2)
                    parr = path.split(mod_2)
                    title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']+subtit
                    plt.title(title)
                    val_str = ("min, max, mean = "+"{:.5e}".format(fld_2.min().data)+", "
                            "{:.5e}".format(fld_2.max().data)+", "+"{:.5e}".format(fld_2.mean().data))
                    ax.annotate(tim_2+' '+val_str,xy=(0,-0.25),xycoords='axes fraction')
                    pp.savefig()
                else:
                    fig = plt.figure(figsize=[8.5,11])
                    ax = fig.add_subplot(212)
                    fld_2 = dat_2[var_2].isel(lat=0,lon=0,time=0)
                    subtit = ' (Lat/Lon=0/0)'
                    fld_2.plot(ax=ax,cbar_kwargs={'label': fld_2.units},rasterized=True)
                    path, fname = os.path.split(f_2)
                    parr = path.split(mod_2)
                    title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']+subtit
                    plt.title(title)
                    val_str = ("min, max, mean = "+"{:.5e}".format(fld_2.min().data)+", "
                        "{:.5e}".format(fld_2.max().data)+", "+"{:.5e}".format(fld_2.mean().data))
                    ax.annotate(tim_2+' '+val_str,xy=(0,-0.2),xycoords='axes fraction')
                    pp.savefig()

    pp.close() # multipage document complete
    os.popen('mv multipage.pdf '+out_pdf) # save document to descriptive file name
    
    print('Output file: ', out_pdf)

if __name__ == '__main__':

    description='Plot specified month and year from all CMIP6 variables submitted.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("dir_1", type=str, help="Source directory")
    parser.add_argument("mon_1", type=str, help="YYYY-MM to plot")
    parser.add_argument("--compare", action="store_true", default=False, 
                        help="Compare with a second source directory?")
    parser.add_argument("--dir_2", type=str, help="Optional second source directory")
    parser.add_argument("--mon_2", type=str, help="Optional different YYYY-MM from second source")
    parser.add_argument("--include", type=str, help="Include only these variables or classes (comma separated)")
    parser.add_argument("--exclude", type=str, help="Exclude these variables or classes (comma separated)")
    parser.add_argument("--first", action="store_true", default=False, 
                        help="Plot first version (earliest) instead of last (latest=most recent)")

    args = parser.parse_args()            

    main(args)    
