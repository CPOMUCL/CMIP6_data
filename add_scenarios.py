import numpy as np
import glob
import os
import shutil
import datetime as dt


#### set list of models and ensembles to get
models = [
#     'NCAR/CESM2-WACCM/r1i1p1f1',### FIRST batch - dowloaded ssp585!!!
#     'NCAR/CESM2/r4i1p1f1',
#     'THU/CIESM/r1i1p1f1',
#     'NUIST/NESM3/r1i1p1f1',
#     'CSIRO-ARCCSS/ACCESS-CM2/r1i1p1f1',
    'CMCC/CMCC-ESM2/r1i1p1f1',
    'CMCC/CMCC-CM2-SR5/r1i1p1f1',
#     'NUIST/NESM3/r2i1p1f1',
#     'MRI/MRI-ESM2-0/r1i1p1f1',
#     'AWI/AWI-CM-1-1-MR/r1i1p1f1', ### second batch more to find...
#     'CCCma/CanESM5/r10i1p1f1',
#     'BCC/BCC-CSM2-MR/r1i1p1f1',
#     'NASA-GISS/GISS-E2-1-G/r1i1p1f2',
#     'CAS/CAS-ESM2-0/r1i1p1f1',
#     'CAS/FGOALS-g3/r1i1p1f1',
#     'DKRZ/MPI-ESM1-2-HR/r1i1p1f1', 
#     'NCC/NorESM2-MM/r1i1p1f1',
#     'NCC/NorESM2-LM/r1i1p1f1',
#     'CAS/FGOALS-g3/r1i1p1f1',
#     'CNRM-CERFACS/CNRM-CM6-1/r1i1p1f2',
#     'CCCma/CanESM5/r10i1p1f1', ### another batch for rsds etc...
#     'MOHC/HadGEM3-GC31-MM/r1i1p1f3',
#     'NCC/NorESM2-MM/r1i1p1f1',
#     'NCC/NorESM2-LM/r1i1p1f1',
#     'CAS/FGOALS-g3/r1i1p1f1',
#     'CNRM-CERFACS/CNRM-CM6-1/r1i1p1f2',
#     'BCC/BCC-CSM2-MR/r1i1p1f1',
#     'NASA-GISS/GISS-E2-1-G/r1i1p1f2',
#     'CNRM-CERFACS/CNRM-CM6-1/r1i1p1f2',
    'IPSL/IPSL-CM6A-LR/r2i1p1f1',
    'NCC/NorESM2-LM/r1i1p1f1',
#     'EC-Earth-Consortium/EC-Earth3-CC/r1i1p1f1',
    'MPI-M/MPI-ESM1-2-LR/r1i1p1f1',
#     'NUIST/NESM3/r1i1p1f1',
    'MIROC/MIROC6/r1i1p1f1',
    'BCC/BCC-CSM2-MR/r1i1p1f1',
    'CCCma/CanESM5/r1i1p2f1',
    'NCAR/CESM2-WACCM/r2i1p1f1',
    'NCAR/CESM2/r11i1p1f1',
    'MRI/MRI-ESM2-0/r1i1p1f1',
    
]



n_check = 56
runs_found = 0
#### base dir
base_dir = '/badc/cmip6/data/CMIP6/'
extra_dir = '/home/users/hheorton/CMIP_source/CMIP6/'
group = 'ScenarioMIP' ### also ScenarioMIP
run_type = 'ssp126' ### others for SMIP


today_str = dt.datetime.today().strftime('d%Y%m%d')
### variables to seach for
# ncvarlist=['siconc','sithick','sisnthick','siflswutop','siflswdtop']
ncvarlist=['siconc','sithick','sisnthick','rsus','rsds','siflswutop','siflswdtop']

#### loop throught models
for md in models:
    ### /badc/cmip6/data/CMIP6/ScenarioMIP/AWI/AWI-CM-1-1-MR/ssp585/r1i1p1f1/SImon/
    md_s = md.split('/')
    md_dir = base_dir+'/'.join([group,md_s[0],md_s[1],run_type,md_s[2]])
    md_new = extra_dir+'/'.join([group,md_s[0],md_s[1],run_type,md_s[2]])
    save_dir = os.path.dirname(md_new)
    if not os.path.exists(save_dir):
        # make if it doesn't
        os.makedirs(save_dir)
        print('Creating directory: ',save_dir)
    vars_found = []
    vars_toget = []
    for var in ncvarlist:
        if 'si' in var:
            type_dir = '/SImon/'
        else:
            type_dir = '/Amon/'
        extra_time_dir = os.path.dirname(md_new+type_dir)
        if not os.path.exists(extra_time_dir):
            # make if it doesn't
            os.makedirs(extra_time_dir)
            print('Creating directory: ',extra_time_dir)
        path = glob.glob(md_dir+type_dir+var+'/gn/files/*/*.nc')
        new_path = glob.glob(md_new+type_dir+var+'/gn/files/*/*.nc')
        more_data = True
        ### also check the new path
        if len(path)>0:
            vars_found.append(var)
            more_data = False
            ### make symbolic link from md_dir+var to md_new+var 
            try:
                ## check dir before linking
                if not os.path.exists(md_new+type_dir):
                    os.makedirs(md_new+type_dir)
                os.symlink(### 'mytarget', 'mylink' )
                    md_dir+type_dir+var+'/',
                    md_new+type_dir+var)
                print(' ------ LINKING DATA ------ ')
            except FileExistsError:
                print(' ------ DATA ALREADY LINKED ------ ')
                pass
        elif len(new_path)>0:
            #### do we have them all?
            ### assume yes for now
            print(' ------ DOWNLOADED ALREADY ------ ')
            more_data = False
        if more_data:
            vars_toget.append(var)
            ### download to md_new+var+'/gn/files/today_str/ 
            save_dir = os.path.dirname(md_new+type_dir+var+'/gn/files/'+today_str+'/')
            print(' ------ DOWNLOADING DATA ------ ')
            if not os.path.exists(save_dir):
                # make if it doesn't
                os.makedirs(save_dir)
                print('Creating directory: ',save_dir)
            os.system('python cmip6_downloader_ALL_NODES.py --variable_id='+var+' --frequency=mon --experiment_id='+run_type+' --source_id='+md_s[1]+' --variant_label='+md_s[2])
            ### move data
            ### example data/s_CanESM5_e_ssp585_vl_r10i1p1f1_f_mon_v_sisnthick
            print(' ------ MOVING DATA ------ ')
            model_info = md.split('/')[1]+'*'
            for f in glob.glob('./data/*'+model_info+var+'/*.nc'):
                ### check file name - if it's in the date window, if not delete
                ### data/siconc_SImon_ACCESS-ESM1-5_ssp585_r10i1p1f1_gn_201501-210012.nc
                ### data/rsds_Amon_CNRM-CM6-1_ssp585_r1i1p1f2_gr_201501-210012.nc
                try:
                    fdates = f.split('_gn_')[1].split('.nc')[0]
                except IndexError:
                    fdates = f.split('_gr_')[1].split('.nc')[0]
                print(fdates)
                d0 = dt.datetime.strptime(fdates.split('-')[0],'%Y%m')
                d1 = dt.datetime.strptime(fdates.split('-')[1],'%Y%m')
#                 if d0 <= dt.datetime(2100,1,1):
                if d1 >= dt.datetime(2000,1,1):
                    dest =  md_new+type_dir+var+'/gn/files/'+today_str+'/'+f.split(var+'/')[1]
                    print('Moving to '+dest)
                    shutil.move(f, dest)
                else:
                    os.remove(f)

