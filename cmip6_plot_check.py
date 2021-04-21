import os
import argparse
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob2 as glob
from matplotlib.backends.backend_pdf import PdfPages as pdf

def main(inargs):
    """Run the program."""

    # specify source files from two model runs

    dir_1 = inargs.dirs[0]
    ncs_1 = glob.glob(dir_1+'/**/*.nc', recursive=True)
    mod_1 = dir_1.split('/')[-3]
    run_1 = dir_1.split('/')[-2]
    print('First source: ', dir_1)

    if len(inargs.dirs)==2:
        dir_2 = inargs.dirs[1]
        ncs_2 = glob.glob(dir_2+'/**/*.nc', recursive=True)
        mod_2 = dir_2.split('/')[-3]
        run_2 = dir_2.split('/')[-2]
        print('Second source: ', dir_2)

    # specify local destination for comparison plots

    if len(inargs.dirs)==2:
        out_pdf = mod_1+'_'+run_1+'_vs_'+mod_2+'_'+run_2+'.pdf'
    else:
        out_pdf = mod_1+'_'+run_1+'.pdf'

    # loop over source files in first model run

    print('Processing ...')
    pp = pdf('multipage.pdf')

    for i_1, f_1 in enumerate(ncs_1):
        print(f_1)
        dat_1 = xr.open_dataset(f_1)
        var_1 = list(dat_1.data_vars.keys())[-1]
        ndims = len(dat_1[var_1].dims)
        if ndims==3:
            fig = plt.figure(figsize=[8.5,11])
            ax = fig.add_subplot(211,projection=ccrs.PlateCarree(central_longitude=180))
            fld_1 = dat_1[var_1].isel(time=0)
            fld_1.plot(ax=ax,transform=ccrs.PlateCarree(),
                       cbar_kwargs={'label': fld_1.units},rasterized=True)
            path, fname = os.path.split(f_1)
            parr = path.split(mod_1)
            title = parr[0]+mod_1+'\n'+parr[1]+'/\n'+fname+'\n'+fld_1.attrs['long_name']
            plt.title(title)
            val_str = ("min, max, mean = "+"{:.6e}".format(fld_1.min().data)+", "
                       "{:.6e}".format(fld_1.max().data)+", "+"{:.6e}".format(fld_1.mean().data))
            ax.annotate(val_str,xy=(0,-0.1),xycoords='axes fraction')
            ax.coastlines()

            if len(inargs.dirs)==2:
                search_str = fname.split('_')[0]+'_'+fname.split('_')[1]
                matching_file = [i for i in ncs_2 if search_str in i]
            else: matching_file = []
            if matching_file != []:
                f_2 = matching_file[0]
                dat_2 = xr.open_dataset(f_2)
                var_2 = list(dat_2.data_vars.keys())[-1]
                ax = fig.add_subplot(212,projection=ccrs.PlateCarree(central_longitude=180))
                fld_2 = dat_2[var_2].isel(time=0)
                fld_2.plot(ax=ax,transform=ccrs.PlateCarree(),
                       cbar_kwargs={'label': fld_2.units},rasterized=True)
                path, fname = os.path.split(f_2)
                parr = path.split(mod_2)
                title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']
                plt.title(title)
                val_str = ("min, max, mean = "+"{:.6e}".format(fld_2.min().data)+", "
                           "{:.6e}".format(fld_2.max().data)+", "+"{:.6e}".format(fld_2.mean().data))
                ax.annotate(val_str,xy=(0,-0.1),xycoords='axes fraction')
                ax.coastlines()

            pp.savefig()
        elif ndims==4:
            fig = plt.figure(figsize=[8.5,11])
            ax = fig.add_subplot(211)
            fld_1 = dat_1[var_1].isel(time=0,lon=0)
            fld_1.plot(ax=ax,
                       cbar_kwargs={'label': fld_1.units},rasterized=True)
            if dat_1[var_1].dims[1]==('lev') or dat_1[var_1].dims[1]==('plev'): ax.invert_yaxis()
            title = fld_1.attrs['long_name'] + ' (Longitude=0)'
            path, fname = os.path.split(f_1)
            parr = path.split(mod_1)
            title = parr[0]+mod_1+'\n'+parr[1]+'/\n'+fname+'\n'+fld_1.attrs['long_name']+' (Longitude=0)'
            plt.title(title)
            val_str = ("min, max, mean = "+"{:.6e}".format(fld_1.min().data)+", "
                       "{:.6e}".format(fld_1.max().data)+", "+"{:.6e}".format(fld_1.mean().data))
            ax.annotate(val_str,xy=(0,-0.25),xycoords='axes fraction')

            if len(inargs.dirs)==2:
                search_str = fname.split('_')[0]+'_'+fname.split('_')[1]
                matching_file = [i for i in ncs_2 if search_str in i]
            else: matching_file = []
            if matching_file != []:
                f_2 = matching_file[0]
                dat_2 = xr.open_dataset(f_2)
                var_2 = list(dat_2.data_vars.keys())[-1]
                ax = fig.add_subplot(212)
                fld_2 = dat_2[var_2].isel(time=0,lon=0)
                fld_2.plot(ax=ax,
                       cbar_kwargs={'label': fld_2.units},rasterized=True)
                if dat_2[var_2].dims[1]==('lev') or dat_2[var_2].dims[1]==('plev'): ax.invert_yaxis()
                path, fname = os.path.split(f_2)
                parr = path.split(mod_2)
                title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']+' (Longitude=0)'
                val_str = ("min, max, mean = "+"{:.6e}".format(fld_2.min().data)+", "
                           "{:.6e}".format(fld_2.max().data)+", "+"{:.6e}".format(fld_2.mean().data))
                ax.annotate(val_str,xy=(0,-0.25),xycoords='axes fraction')
                plt.title(title)

            fig.tight_layout(pad=6)
            pp.savefig()
        plt.close()

    # loop over source files in second model run (plot only any missing from first run)

    if len(inargs.dirs)==2: 
        for i_2, f_2 in enumerate(ncs_2):
            path, fname = os.path.split(f_2)
            matching_file = [i for i in ncs_1 if fname.split('_')[0]+'_'+fname.split('_')[1] in i]
            if matching_file == []:
                print(f_2)
                var_2 = list(dat_2.data_vars.keys())[-1]
                ndims = len(dat_2[var_2].dims)
                if ndims==3:
                    fig = plt.figure(figsize=[8.5,11])
                    ax = fig.add_subplot(212,projection=ccrs.PlateCarree(central_longitude=180))
                    fld_2 = dat_2[var_2].isel(time=0)
                    fld_2.plot(ax=ax,transform=ccrs.PlateCarree(),
                    cbar_kwargs={'label': fld_2.units},rasterized=True)
                    path, fname = os.path.split(f_2)
                    parr = path.split(mod_2)
                    title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']
                    plt.title(title)
                    val_str = ("min, max, mean = "+"{:.6e}".format(fld_2.min().data)+", "
                            "{:.6e}".format(fld_2.max().data)+", "+"{:.6e}".format(fld_2.mean().data))
                    ax.annotate(val_str,xy=(0,-0.1),xycoords='axes fraction')
                    ax.coastlines()
                    pp.savefig()
                elif ndims==4:
                    fig = plt.figure(figsize=[8.5,11])
                    ax = fig.add_subplot(212)
                    fld_2 = dat_2[var_2].isel(time=0,lon=0)
                    fld_2.plot(ax=ax,cbar_kwargs={'label': fld_2.units},rasterized=True)
                    if dat_2[var_2].dims[1]==('lev') or dat_2[var_2].dims[1]==('plev'): ax.invert_yaxis()
                    title = fld_1.attrs['long_name'] + ' (Longitude=0)'
                    path, fname = os.path.split(f_2)
                    parr = path.split(mod_2)
                    title = parr[0]+mod_2+'\n'+parr[1]+'/\n'+fname+'\n'+fld_2.attrs['long_name']+' (Longitude=0)'
                    plt.title(title)
                    val_str = ("min, max, mean = "+"{:.6e}".format(fld_2.min().data)+", "
                            "{:.6e}".format(fld_2.max().data)+", "+"{:.6e}".format(fld_2.mean().data))
                    ax.annotate(val_str,xy=(0,-0.25),xycoords='axes fraction')
                    fig.tight_layout(pad=6)
                    pp.savefig()
                plt.close()

    pp.close()
    os.popen('mv multipage.pdf '+out_pdf)
    
    print('Output file: ', out_pdf)

if __name__ == '__main__':

    description='Print the input arguments to the screen.'
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument("dirs", type=str, nargs='*', help="One or two directory names")

    args = parser.parse_args()            

    main(args)    
