
import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import grid_set as gs
import CMIP_inputs as Ci
import CMIP_functions as Cf
from netCDF4 import Dataset
import xarray as xr
import glob

import sys



### input options are (proc_input) Out_file
### proc_input contains all information needed to find the input data, and save the post processed in a tidy way

inFile = sys.argv[1]
print_file = sys.argv[2]
print('inFile is '+inFile)
print('printfile is '+print_file)
pFile = open(print_file, 'w')
sys.stdout = pFile

### pre declare extra options (usually as false)
swflux_grid = False

with open(inFile,'r') as f:
    f.readline()
    CMIP6_path = f.readline().split(' ')[1].split('\n')[0]
    spath = f.readline().split(' ')[1].split('\n')[0]+'/'
    orient_opt = int(f.readline().split(' ')[1])
    print('OPTION: grid orientation: '+str(orient_opt),file=pFile,flush=True)
    #### try to read extra lines until exit
    more_lines = True
    while more_lines:
        options = f.readline()
        ### options to find here
        ### extra options if you want them
        ### for example some of my processing needs a seperate grid for surface flux
        ### say if you need to read from atmospheric or ocean variables for some models, but not for others.
        if 'swflux_grid' in options:
            swflux_grid = True
            print('OPTION: making a seperate swflux grid',file=pFile,flush=True)
        if len(options) == 0:
            more_lines = False
### orient_opt
### 1. lon/lat gets transposed, so does data
### 2. lon/lat needs to be mesh gridded, data is kept normal
### 3. lon/lat/data all unstructured, all left as list

### make this true to get extra information on output
load_verbos=False

### you can hack the CMIP6_path here if you need (say to switch scenario)
### this can alternately be done with a new proc_input file
### FORCE ssp126 scenario from 585
# CMIP6_path = CMIP6_path.replace('ssp585','ssp126')

### data object
### the data object can access multiple data variables. But they need to be consistent in shape, orientation and time period
### the buffer parameter here sets how many time slices to load with each call.
### this is a performance setting, loading full files can be a bad idea (memory limits)
### loading single slices can be a bad idea (slow), so we buffer a chunk. Don't worry, the object sorts this for you. You can push the parameter up to get more performance.
C = Ci.CMIP6_monthly(CMIP6_path,buffer=12,orient_opt = orient_opt)
### for the example option given above, I needed a new data object for surface fluxes (but not all the time)
# if swflux_grid: 
#     Cflux = Ci.CMIP6_monthly(CMIP6_path,buffer=12,orient_opt = orient_opt_flux)
# else:
#     Cflux = C

## dates, we will find the dates we want to analyse using start and end points
ts = dt.datetime(2015,1,1)
te = dt.datetime(2016,12,31)
### the get_dates method, uses a target variable to scan to see which time points we can use. This bit sorts out all the file structures for you.
C.get_dates(ts,te,'SImon/siconc')
### example dual data object option
# if swflux_grid: 
#     Cflux.get_dates(ts,te,'Amon/rsus')
# else: 
#     Cflux = C

#### SETTING UP GRIDS
## first a projection
# m = ccrs.RotatedPole(pole_latitude=00.0)
m = ccrs.LambertAzimuthalEqualArea(central_latitude=90)

### get grid information from the first file found
G = Cf.gs_from_CMIP(m,C.file_names[0],orient_opt)
### and for the example option
# if swflux_grid: 
#     Gflux = Cf.gs_from_CMIP(m,Cflux.file_names[0],orient_opt_flux)


### setting the saving grid
### this using my grid_set module. Formalises lots of gridding issues
GN = gs.grid_set(m)
GN.load_grid('./Ease_reduced_for_Eco.npz')

### regridder object using grid_set again
G2GN = gs.Gs2Gs(G,GN)
### example option
# if swflux_grid: 
#     Gflux2GN = gs.Gs2Gs(Gflux,GN)
# else:
#     Gflux2GN = G2GN


### MAKE THE PROCESSING LOOP
### to process data we will use a list of dates
### two examlples, first all the dates in the range
d_list = [d for d in C.dates if d>=ts and d<=te]

dstart = d_list[ 0]
dend   = d_list[-1]

### second was a case where I only wanted a few months in the year
#### find all the JFMA
# month_start = 1 ## Jan
# month_end = 4 ## Aprit
# ### finding start and end dates
# for d in C.dates:
#     if d.month == month_start:
#         dstart = d
#         break
# for d in C.dates[::-1]:
#     if d.month == month_end:
#         dend = d
#         break
# print(dstart.strftime('%Y-%m')+dend.strftime('_%Y-%m'),
#       file=pFile,flush=True)
# ### finding the number of years
# n_yrs = dend.year - dstart.year

# month_range = [dstart + relativedelta(months = mn) 
#                for mn in range(month_end-month_start+1)]
# ### make d_list 
# d_list = []
# for y in range(n_yrs):
#     duse = [d + relativedelta(years = y) for d in month_range]
#     d_list.extend(duse)

n_p = len(d_list)  
### set arrays
### outputs
all_ice_conc = np.empty([n_p,GN.n,GN.m])


### fill arrays and process
for n,d in enumerate(d_list):
    print(d.strftime('%Y-%m'),
      file=pFile,flush=True)
    ### load variables
    sic = G2GN.rg_array(C.get_var('SImon/siconc'    ,[d],verbos=load_verbos)[0])/100.0
    #### this uses a scalar array regridding method.
    #### a vector method also exists (rg_vecs), but this is much harder to use
    #### ask me for an example
    #### here is the example optional alternate grid too.
#     if swflux_grid:
#         ### sometimes the fluxes are on a different grid
#         swu = Gflux2GN.rg_array(Cflux.get_var('Amon/rsus',[d],verbos=load_verbos)[0])
    ### print summary
    print('Ice Concentration av: '+'{:.3}'.format(np.nanmean(sic)),file=pFile,flush=True)
    print('Ice Extent (10^6 Km^2): '+'{:.3}'.format(
        np.nansum(sic*GN.xdist*GN.ydist)*10e-12
                    ),file=pFile,flush=True)
    all_ice_conc[n] = sic
    
    ### you can obviously do loads more here in this loop, what do you want to calculate
    
    ### update log file
    sys.stdout.flush()    
    
### saving the output
d0 = dt.datetime(2000,1,1)
times = [(d-d0).days for d in d_list]
times = np.asarray(times,dtype="i4")
#### save it all
file = spath+'CMIP_proc_'+'_'.join([C.name,C.ensemble,C.run_type,str(dstart.year),str(dend.year)])+'.nc'
print('Saving in '+file,file=pFile,flush=True)
ds1=xr.Dataset(data_vars={
        'ice_concentration':(['time','x','y'],all_ice_conc,
                            {'units':'Fraction'}),
                            },
                   coords =  {
        'lon':(['x','y'],GN.lons,{'units':'Degrees E'}),
        'lat':(['x','y'],GN.lats,{'units':'Degrees N'}),
        'time':(['time'],times,{'units':'days since '+d0.strftime('%Y-%m-%d')})},
                  attrs = {
       'description':"Arctic EcoLight " ,
        'comment':"Processing code H. Heorton, and .....",
        'CMIP6_model':C.BIGNAME,
        'CMIP6_name':C.name,
        'CMIP6_ensemble':C.ensemble,
        'CMIP6_type':C.run_type,
        'Created_on':dt.date.today().strftime('%Y-%m-%d'),   
                  })
          
ds1.to_netcdf(file)
print('Success!',file=pFile,flush=True)
sys.stdout.flush()    