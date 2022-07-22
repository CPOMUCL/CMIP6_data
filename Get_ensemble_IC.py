
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

import gc
import sys
import os



### Find all the scenarios
home = os.getcwd()
os.chdir('/badc/cmip6/data/CMIP6/ScenarioMIP/')
base = os.getcwd()
models = {}
variable = 'siconc'
period   = 'SImon'
run_type   = 'ssp245' # 'ssp370'  'ssp585'

for dir1 in glob.glob('*'): #model centre
    os.chdir(dir1)
    for dir2 in glob.glob('*'): #modelID
        if os.path.exists(dir2+'/'+run_type):
            os.chdir(dir2+'/'+run_type)
            check = False
#             models[dir1+'/'+dir2]={}
            model={}
            model['path'] = base+'/'+dir1+'/'+dir2+'/'+run_type+'/'
            model['name'] = dir1+'/'+dir2
            for dir3 in glob.glob('*'): #ripfs
                path = glob.glob(dir3+'/'+period+'/'+str(variable)+'/gn/files/*/*.nc')
                ### example path
                ### '/badc/cmip6/data/CMIP6/ScenarioMIP/CCCma/CanESM5/ssp585/r1i1p2f1'
                if len(path)>0:
                    check=True
                    model.setdefault('ensembles', []).append(dir3)
            if check: models[dir1+'/'+dir2]=model
            os.chdir(base+'/'+dir1)
    os.chdir(base)

os.chdir(home)

### filtering perhaps???
n_ensembles = 1
n_models = 400


keep_models = [
#     'AS-RCEC/TaiESM1' ,
#     'BCC/BCC-CSM2-MR' ,
#     'CAMS/CAMS-CSM1-0',
#     'CAS/FGOALS-f3-L',
#     'CAS/FGOALS-g3',
#     'CCCma/CanESM5', #### just collapses
#     'CCCma/CanESM5-CanOE',
#     'CMCC/CMCC-CM2-SR5', ### no ice at all
#     'CMCC/CMCC-ESM2', ### no ice at all
#     'CNRM-CERFACS/CNRM-CM6-1', ### tripole stripe
#     'CNRM-CERFACS/CNRM-CM6-1-HR',
#     'CNRM-CERFACS/CNRM-ESM2-1', ### tripole stripe
#     'CSIRO/ACCESS-ESM1-5',
#     'CSIRO-ARCCSS/ACCESS-CM2',
#     'DKRZ/MPI-ESM1-2-HR',
#     'EC-Earth-Consortium/EC-Earth3',
#     'EC-Earth-Consortium/EC-Earth3-CC',
#     'EC-Earth-Consortium/EC-Earth3-Veg',
#     'EC-Earth-Consortium/EC-Earth3-Veg-LR',
#     'FIO-QLNM/FIO-ESM-2-0',
#     'IPSL/IPSL-CM6A-LR',### tripole stripe, low ice
#     'MIROC/MIROC-ES2L',
#     'MIROC/MIROC6',
#     'MOHC/HadGEM3-GC31-LL',
#     'MOHC/UKESM1-0-LL',
#     'MPI-M/MPI-ESM1-2-LR',
#     'MRI/MRI-ESM2-0',
#     'NCAR/CESM2',
#     'NCAR/CESM2-WACCM',
#     'NCC/NorESM2-LM',
#     'NCC/NorESM2-MM',
#     'NOAA-GFDL/GFDL-CM4',
#     'NOAA-GFDL/GFDL-ESM4',
#     'NUIST/NESM3', ## very little ice
#     'THU/CIESM', ### no ice at all
]
# models = {key:models[key] for key in models if key in keep_models}

del_models = [
    'CAS/CAS-ESM2-0', ### dodgy grid and no ice
    'CMCC/CMCC-CM2-SR5', ### no ice at all
    'CMCC/CMCC-ESM2', ### no ice at all
    'THU/CIESM', ### no ice at all
    'CCCma/CanESM5', #### just collapses
    'CNRM-CERFACS/CNRM-CM6-1', ### tripole stripe
    'CNRM-CERFACS/CNRM-CM6-1-HR',
    'CNRM-CERFACS/CNRM-ESM2-1', ### tripole stripe
    'IPSL/IPSL-CM6A-LR',### tripole stripe, low ice
        ]
models = {key:models[key] for key in models if key not in del_models}


print('MODELS FOUND')
for kn,key in enumerate(models):
    model = models[key]
    print(model['name'],model['ensembles'],flush=True)
    

### might need something here - can go through models list to do something perhaps
### orient_opt
### 1. lon/lat gets transposed, so does data
### 2. lon/lat needs to be mesh gridded, data is kept normal
### 3. lon/lat/data all unstructured, all left as list

### make this true to get extra information on output
load_verbos=False



### data object
### the data object can access multiple data variables. But they need to be consistent in shape, orientation and time period
### the buffer parameter here sets how many time slices to load with each call.
### this is a performance setting, loading full files can be a bad idea (memory limits)
### loading single slices can be a bad idea (slow), so we buffer a chunk. Don't worry, the object sorts this for you. You can push the parameter up to get more performance.
## make all the objects add to the dictionary
## dates, we will find the dates we want to analyse using start and end points
ts = dt.datetime(2015,1,1)
te = dt.datetime(2039,12,31)
# te = dt.datetime(2099,12,31)
for kn,key in enumerate(models):
    if kn>=n_models: break
    model = models[key]
    model['Cilist']=[]
    ### models default orientation option for later
    model['orient_opt'] = 1
    for ens in model['ensembles'][:n_ensembles]:
        en_path = model['path']+ens
        model['Cilist'].append(Ci.CMIP6_monthly(en_path,orient_opt = model['orient_opt']))
        model['Cilist'][-1].get_dates(ts,te,period+'/'+variable)
### the get_dates method, uses a target variable to scan to see which time points we can use. This bit sorts out all the file structures for you.


#### SETTING UP GRIDS
## first a projection NSIDC SIC
m = ccrs.NorthPolarStereo(central_longitude=-45)
G = gs.grid_set(m)
G.load_grid('./NSIDC_gs.npz')
# G.load_mask('./Ease_reduced_for Eco_mask.npz')
G.get_grid_mask()

#### hardwire some example models that need specific orient options
### those that need #1
OrOp2 = [
#     'CAS/CAS-ESM2-0',
        ]
### those that need #2
OrOp3 = [
        ]

### get grid information from the first file found, for one ensemble per model
for kn,key in enumerate(models):
    if kn>=n_models: break
    model = models[key]
    ### check the orientations
    for check in OrOp2:
        if check in model['name']: 
            model['orient_opt'] = 2
            print('Setting orientation to 2 for '+model['name'],flush=True)
    for check in OrOp3:
        if check in model['name']: 
            model['orient_opt'] = 3
            print('Setting orientation to 3 for '+model['name'],flush=True)
    model['grid']=Cf.gs_from_CMIP(m,model['Cilist'][0].file_names[0],model['orient_opt'])
    model['regridder']=gs.Gs2Gs(model['grid'],G,NaN_avoid=True)

### MAKE THE PROCESSING LOOP
### to process data we will use a list of dates
### use the last accessed model
d_list = [d for d in model['Cilist'][0].dates if d>=ts and d<=te]

#### SEPTEMBERS ONLY
# d_list = [d for d in d_list if d.month == 9]

dstart = d_list[ 0]
dend   = d_list[-1]


n_p = len(d_list)  
### set arrays
### outputs
###1 average the area
# all_ice_conc = np.empty([n_p,G.n,G.m])
ens_ice_conc = Cf.accumulator([n_p,G.m,G.n])
all_ice_conc = Cf.accumulator([n_p,G.m,G.n])

###2 average the extent
# all_ice_extn = np.empty([n_p,G.n,G.m])
ens_ice_extn = Cf.accumulator([n_p,G.m,G.n])
all_ice_extn = Cf.accumulator([n_p,G.m,G.n])

### but we also need
sic_temp = np.empty([n_p,G.m,G.n])

## count to accumulate the mean
# load_count = 0

### fill arrays and process
print('Data loading time',flush=True)

for kn,key in enumerate(models):
    ### accumulate each models ensembles
    ens_ice_conc.clean()
    ens_ice_extn.clean()
    if kn>=n_models: break
    model = models[key]
    print('Model: '+model['name'],flush=True)
    for ne, (ensCi,ens) in enumerate(zip(model['Cilist'],model['ensembles'])):
        if ne >= n_ensembles: break
        print('Ensemble: '+ens,flush=True)

        ### load variables
        sic_load = ensCi.get_var('SImon/siconc'    ,d_list,verbos=load_verbos)
        ### fill missing values (we will mask later)
#         sic_load[np.isnan(sic_load)] = 0.0
        
        print('regridding . . . . . . . .',flush=True)
        for n in range(n_p):
            sic_temp[n] = model['regridder'].rg_array(sic_load[n])*G.mask/100.
        
        ### accumulate
        ens_ice_conc.update(sic_temp,count_all = True)
        ens_ice_extn.update(sic_temp>0.15,count_all = True)
        
        ### print summary
        print('Ice Concentration av: '+'{:.3}'.format(np.nanmean(sic_temp)),flush=True)
        print('Ice Area av (10^6 Km^2): '+'{:.3}'.format(
            np.nansum(np.nanmean(sic_temp)*G.xdist*G.ydist)*10e-12
                        ),flush=True)
        ensCi.clean_var('SImon/siconc')
        gc.collect()
    ### accumulate all the models ensemble mean
    all_ice_conc.update(ens_ice_conc.mean(),count_all = True)
    all_ice_extn.update(ens_ice_extn.mean(),count_all = True)
    ### update log file
    sys.stdout.flush()    
    
    
### saving the output
d0 = dt.datetime(2000,1,1)
times = [(d-d0).days for d in d_list]
times = np.asarray(times,dtype="i4")
#### save it all

### info for the NC file
attrs = {
       'description':"Arctic CMIP6 Ensemble Mean ssp245 " ,
        'comment':"Processing code H. Heorton, and .....",
        'Created_on':dt.date.today().strftime('%Y-%m-%d')}
for kn,key in enumerate(models):
    if kn>=n_models: break
    model = models[key]
    attrs[model['name'].replace('/','--')+'_ensembles_used'] = ', '.join(model['ensembles'][:n_ensembles])
    
file = "IC_ensemble_save_+"+run_type+".nc"
print('Saving in '+file,flush=True)
ds1=xr.Dataset(data_vars={
        'ice_area_ensemble':(['time','x','y'],all_ice_conc.mean(),
                            {'units':'Fraction'}),
        'ice_extent_ensemble':(['time','x','y'],all_ice_extn.mean(),
                            {'units':'Fraction'}),
                            },
                   coords =  {
        'lon':(['x','y'],G.lons,{'units':'Degrees E'}),
        'lat':(['x','y'],G.lats,{'units':'Degrees N'}),
        'time':(['time'],times,{'units':'days since '+d0.strftime('%Y-%m-%d')})},
                  attrs = attrs)
          
ds1.to_netcdf(file)
print('Success!',flush=True)
sys.stdout.flush()    