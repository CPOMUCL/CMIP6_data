
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

# from guppy import hpy 
# h = hpy() 
import gc
import sys
import os



### Find all the scenarios
def get_scenlist(directory,runtype,variable,period):
    home = os.getcwd()
    os.chdir(directory)
    base = os.getcwd()
    models = {}
    var_per = '/'.join([period,variable,])

    for dir1 in glob.glob('*'): #model centre
        os.chdir(dir1)
        for dir2 in glob.glob('*'): #modelID
            if os.path.exists(dir2+'/'+run_type):
                os.chdir(dir2+'/'+run_type)
                check = False
    #             models[dir1+'/'+dir2]={}
                model={}
                model['path'] = base+'/'+dir1+'/'+dir2+'/'+runtype+'/'
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
    return models

def comp_mod_list(x1,y1,new_path = False,new_ens=False):
    def comp_mod(x2,y2):
        xy = set(x2['ensembles']) & set(y2['ensembles'])
        if xy == set():
            return  False
        else:
            return xy.copy()
    
    models={}
    for k in x1:
        if k in y1:
            xym = comp_mod(x1[k],y1[k])
            if xym is not False:
                models[k] = x1[k].copy()
                if new_path:
                    models[k]['path2'] = y1[k]['path']
                models[k]['ensembles'] = list(xym)
                if new_ens:
                    models[k]['ensembles2'] = y1[k]['ensembles']
    return models

if __name__=='__main__':
    
    period   = 'SImon'
    run_type   = 'ssp585' # 'ssp370'  'ssp585'
#     search_dir = '/badc/cmip6/data/CMIP6/ScenarioMIP/'
    search_dir = '/home/users/hheorton/CMIP_source/CMIP6/ScenarioMIP/'
    variable = 'sithick'
    models1_si = get_scenlist(search_dir,run_type,variable,period)
    search_dir = '/home/users/hheorton/CMIP_source/CMIP6/ScenarioMIP/'
    variable = 'sisnthick'
    models1_sn = get_scenlist(search_dir,run_type,variable,period)
    models1 = comp_mod_list(models1_si,models1_sn)
    
    print('MODELS FOUND '+run_type)
    for kn,key in enumerate(models1):
        model = models1[key]
        print(model['name'],model['ensembles'],flush=True)
    print('------------')


    run_type   = 'historical' # 'ssp370'  'ssp585'
#     search_dir = '/badc/cmip6/data/CMIP6/CMIP/'
    search_dir = '/home/users/hheorton/CMIP_source/CMIP6/CMIP/'
    variable = 'sithick'
    models2_si = get_scenlist(search_dir,run_type,variable,period)
    search_dir = '/home/users/hheorton/CMIP_source/CMIP6/CMIP/'
    variable = 'sisnthick'
    models2_sn = get_scenlist(search_dir,run_type,variable,period)
    models2 = comp_mod_list(models2_si,models2_sn)
    print('MODELS FOUND '+run_type)
    for kn,key in enumerate(models2):
        model = models2[key]
        print(model['name'],model['ensembles'],flush=True)
    print('------------')
    
    var_read = ['sithick','sisnthick']
    var_per = ['/'.join(['SImon',vr]) for vr in var_read]
    ### actually make sure the historical is done first
    models = comp_mod_list(models2,models1,new_path = True,new_ens=True)
    ### filtering perhaps???
    n_ensembles = 1
    n_models = 400
    ### load stride to keep data handling sensbile
    load_strd = 60
    
    save_model = True  ### file per model ensemble - CAREFUL NOW
    save_ens   = False ### accumulate the ensemble


    keep_models = [
    ]

    del_models = [
        'NASA-GISS/GISS-E2-1-G', ### good but unstructure and tricky right now.
        'AWI/AWI-CM-1-1-MR', ### good but unstructure and tricky right now.   
        'CAS/CAS-ESM2-0', ### dodgy grid and no ice
        'CMCC/CMCC-CM2-SR5', ### no ice at all
        'CMCC/CMCC-ESM2', ### no ice at all
        'THU/CIESM', ### no ice at all
        'CCCma/CanESM5', #### just collapses
        'DKRZ/MPI-ESM1-2-HR', #### only data from 2000
#         'CNRM-CERFACS/CNRM-CM6-1', ### tripole stripe
#         'CNRM-CERFACS/CNRM-CM6-1-HR',
#         'CNRM-CERFACS/CNRM-ESM2-1', ### tripole stripe
#         'IPSL/IPSL-CM6A-LR',### tripole stripe, low ice
            ]
    models = {key:models[key] for key in models if key not in del_models}
    
    alpha_check = 'NCB'
    
    models = {key:models[key] for key in models if sorted([alpha_check,key])[0]!=key}


    print('MODELS FOUND')
    for kn,key in enumerate(models):
        model = models[key]
        print(model['name'],model['ensembles'],flush=True)
    print('------------')


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
    ts = dt.datetime(1980,1,1)
#     ts = dt.datetime(2015,1,1)
    te = dt.datetime(2018,12,31)
#     te = dt.datetime(2049,12,31)
    # te = dt.datetime(2099,12,31)


    #### SETTING UP GRIDS
    ## first a projection NSIDC SIC
#     m = ccrs.NorthPolarStereo(central_longitude=-45)
    m = ccrs.LambertAzimuthalEqualArea(central_latitude=90)
    G = gs.grid_set(m)
#     G.load_grid('./NSIDC_gs.npz')
    G.load_grid('./Pathfinder_gs.npz')
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

    print('List size: ',sys.getsizeof(models))
    # print(h.heap())

    ### NEW ONLY loop models once - to save data overheads
    if save_ens:
        all_data1 = Cf.accumulator([n_p,G.n,G.m])
        all_data2 = Cf.accumulator([n_p,G.n,G.m])
    for kn,key in enumerate(models):
        if kn>=n_models: break
        model = models[key]
        model['Cilist']=[]
        ### models default orientation option for later
        model['orient_opt'] = 1
        for ens in model['ensembles'][:n_ensembles]:
            ### first should be historical
            en_path = model['path']+ens
            model['Cilist'].append(Ci.CMIP6_monthly(en_path,orient_opt = model['orient_opt']))
            model['Cilist'][-1].get_dates(ts,te,period+'/'+variable)
            ### and then second is the scenario, appending the Ci object dates
            ### see if ens is in the list, use if so, if not, then the first
            if ens in model['ensembles2']: ens2 = ens
            else: ens2 = model['ensembles2'][0]
            model['Cilist'][-1].path = model['path2']+ens2
            model['Cilist'][-1].get_dates(ts,te,period+'/'+variable,append=True)
        ### the get_dates method, uses a target variable to scan to see which time points we can use. This bit sorts out all the file structures for you.

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

        print('Model size: ',sys.getsizeof(model))

        if kn == 0: #### on the first time - check all the dates
            ### MAKE THE PROCESSING LOOP
            ### to process data we will use a list of dates
            ### use the last accessed model
            d_list = [d for d in model['Cilist'][0].dates if d>=ts and d<=te]

            #### SEPTEMBERS ONLY
            # d_list = [d for d in d_list if d.month == 9]

            ### loading stride list
            dstrd_list = [d_list[n:n+load_strd] for n in range(0,len(d_list),load_strd)]
            ### fill slices
            dstrd_slice = [range(n,n+load_strd) for n in range(0,len(d_list),load_strd)]
            ### fix last slice
            dstrd_slice[-1] = range(dstrd_slice[-1].start,
                                    np.minimum(dstrd_slice[-1].stop,len(d_list)))

            dstart = d_list[ 0]
            dend   = d_list[-1]


            n_p = len(d_list)  
            ### set arrays
            ### outputs
            ###1 average the area
            if save_ens:
                ens_data1 = Cf.accumulator([n_p,G.m,G.n])
                ens_data2 = Cf.accumulator([n_p,G.m,G.n])


        ### fill arrays and process
        print('Data loading time',flush=True)

        ### accumulate each models ensembles
        if save_ens:
            ens_data1.clean()
            ens_data2.clean()
        if kn>=n_models: break
        model = models[key]
        print('Model: '+model['name'],flush=True)
        for ne, (ensCi,ens) in enumerate(zip(model['Cilist'],model['ensembles'])):
            if ne >= n_ensembles: break
            print('Ensemble: '+ens,flush=True)

            ### but we also need
            data_temp1 = np.empty([n_p,G.m,G.n])
            data_temp2 = np.empty([n_p,G.m,G.n])

            ### load variables (STRIDING)
            for d_load,d_slice in zip(dstrd_list,dstrd_slice):
                load1 = ensCi.get_var(var_per[0],d_load,verbos=load_verbos)
                load2 = ensCi.get_var(var_per[1],d_load,verbos=load_verbos)
                ### fill missing values (we will mask later)
        #         sic_load[np.isnan(sic_load)] = 0.0

                print('regridding . . . . . . . .',flush=True)
                for n in d_slice:
                    data_temp1[n] = model['regridder'].rg_array(load1[n-d_slice.start])*G.mask
                    data_temp2[n] = model['regridder'].rg_array(load2[n-d_slice.start])*G.mask
            
            ### info for the NC file
            if save_model:
                d0 = dt.datetime(2000,1,1)
                times = [(d-d0).days for d in d_list]
                times = np.asarray(times,dtype="i4")
                file = '_'.join(var_read+[model['name'].replace('/','--'),ens,ts.strftime('%Y%m%d'),te.strftime('%Y%m%d')])+".nc"
                print('Saving in '+file,flush=True)
                attrs = {
                       'description':"Arctic CMIP6 standard grid" ,
                        'comment':"Processing code H. Heorton, and .....",
                        'Created_on':dt.date.today().strftime('%Y-%m-%d')}
                ## need model name in attributes
                model = models[key]
                attrs[model['name'].replace('/','--')+'_'+model['path'].replace('/','--')] = model['ensembles'][ne]
#                 attrs[model['name'].replace('/','--')+'_'+model['path2'].replace('/','--')] = model['ensembles2'][kn]
                attrs['data_1'] = var_read[0]
                attrs['data_2'] = var_read[1]

                ds1=xr.Dataset(data_vars={
                        var_read[0]:(['time','x','y'],data_temp1,
                                            {'units':'m'}),
                        var_read[1]:(['time','x','y'],data_temp2,
                                            {'units':'m'}),
                                            },
                                   coords =  {
                        'lon':(['x','y'],G.lons,{'units':'Degrees E'}),
                        'lat':(['x','y'],G.lats,{'units':'Degrees N'}),
                        'time':(['time'],times,{'units':'days since '+d0.strftime('%Y-%m-%d')})},
                                  attrs = attrs)

                ds1.to_netcdf(file)
                print('Success!',flush=True)

            ### accumulate
            if save_ens:
                ens_data1.update(data_temp1,count_all = True)
                ens_data2.update(data_temp1>0.15,count_all = True)

            ### print summary
#             print('Ice Concentration av: '+'{:.3}'.format(np.nanmean(sic_temp)),flush=True)
#             print('Ice Area av (10^6 Km^2): '+'{:.3}'.format(
#                 np.nansum(np.nanmean(sic_temp)*G.xdist*G.ydist)*10e-12
#                             ),flush=True)
            for vp in var_per: ensCi.clean_var(vp)
            gc.collect()
    #     print(h.heap())
        ### accumulate all the models ensemble mean
        if save_ens:
            all_data1.update(ens_data1.mean(),count_all = True)
            all_data2.update(ens_data2.mean(),count_all = True)
        ### clean all the grid info
        del model['regridder']
        del model['grid']
        gc.collect()
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

    file = variable+"_ensemble_save_"+".nc"
    file = '_'.join(var_read+['ensemble_save',ts.strftime('%Y%m%d'),te.strftime('%Y%m%d')])+".nc"
    print('Saving in '+file,flush=True)
    ds1=xr.Dataset(data_vars={
            var_read[0]:(['time','x','y'],all_data1.mean(),
                                {'units':'m'}),
            var_read[1]:(['time','x','y'],all_data2.mean(),
                                {'units':'m'}),
                                },
                       coords =  {
            'lon':(['x','y'],G.lons,{'units':'Degrees E'}),
            'lat':(['x','y'],G.lats,{'units':'Degrees N'}),
            'time':(['time'],times,{'units':'days since '+d0.strftime('%Y-%m-%d')})},
                      attrs = attrs)

    ds1.to_netcdf(file)
    print('Success!',flush=True)
    sys.stdout.flush()    