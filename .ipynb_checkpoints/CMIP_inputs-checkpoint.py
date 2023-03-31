
import numpy as np
import datetime as dt
import struct
import glob
import gc
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
from os.path import exists


class CMIP6_daily():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath,orient_opt = 1,buffer = False,regrid_read=None,search_files =True):
        ### get name from path
        ### buffing assumes that all files are same size - iceconc, icethick etc
        psplit = ppath.split('CMIP6/')[1] ### gets CMIP/Scenario
        self.group = psplit.split('/')[0]
        self.BIGNAME = psplit.split('/')[1]
        self.name = psplit.split('/')[2]
        self.run_type = psplit.split('/')[3]
        self.ensemble = psplit.split('/')[4]
        self.path = ppath
        self.orient_opt = orient_opt
        self.buff_vars = False
        self.vars_loaded = False
        self.buff_length = 1
        self.search_files = search_files
        if buffer is not False:
            self.buff_vars = True
            self.buff_length = buffer
        if regrid_read is not None:
            self.regrid_read = True
            self.regrid_GS2GS = regrid_read
        else:
            self.regrid_read = False

    def get_dates(self,time_start,time_end,target,append=False):
        """
        returns the all encompassing date list for use with the forcing object
        we find all the start and end time for each file found in the directory
        Using these times we add the relevant daily dates also considering tiem start and end
        
        Parameters
        ---------
        time_start: datetime.datetime the start time point we want to find time points from
        time_end:   datetime.datetime the end time point to stop collecting data points
        append:   logical, option, default = False. 
                setting this to True  allows us to add a scenarioMIP to the end of a historical run.
                The list of files and dates generated from a first call, can then be append with a second directory. Do this carefully to make sure the runs are consecutive.
        """
        #### search the directories for all the files
        #### added 2021-08 recursive digging through flies/d837583 structure
        #### then split to sort
        self.target = target
        file_names = glob.glob(self.path+'/'+target+'/gn/files/d*/*.nc',recursive=True)
        try:
            file_names.sort(key = lambda x: x.split('_gn_')[1])
        except IndexError:
            file_names.sort(key = lambda x: x.split('_gr_')[1])
        if len(file_names) == 0:
            file_names = glob.glob(self.path+'/'+target+'/gr/files/d*/*.nc',recursive=True)
            try:
                file_names.sort(key = lambda x: x.split('_gr_')[1])
            except IndexError:
                file_names.sort(key = lambda x: x.split('_gn_')[1])
        ###parent_time_units: days since 1850-01-01-00-00-00
        if append:
            dates= self.dates
            ds = self.ds
            de = self.de
        else:
            dates= []
            ds = []
            de = []
        for file in file_names:
            both_dates = file.split('_'+self.ensemble)[1].split('_')[2].split('.')[0]
            d0 = dt.datetime.strptime(both_dates.split('-')[0],'%Y%m%d')
            d1 = dt.datetime.strptime(both_dates.split('-')[1],'%Y%m%d')
            ds.append(d0)
            de.append(d1)
            ### using the d0 and de get
            if time_start<=d1 and time_end>=d0:
                ### then this files contains relevant dates 
                d_start =np.maximum((time_start - d0).days,0)
                d_end   =np.minimum((time_end - d0).days,(d1 - d0).days)+1
                [dates.append(d0 + relativedelta(days = d)) 
                         for d in range(d_start, d_end)]
        self.dates = dates
        self.ds = ds
        self.de = de
        if append:
            self.file_names.extend(file_names)
        else:
            self.file_names = file_names
        print(self.name+' Found '+str(np.shape(dates)[0])+' dates')

    # daily points in yearly files

    def get_var(self,var,dates_u,verbos = False):
        ### use info about the var buffer and dates 
        ### either return the next points in the buffer, or 
        ### read the next chunk from the file
        set_var_now = False
        join_data = False
        if self.buff_vars and not hasattr(self,var):
            ### set it all up later
            set_var_now = True
            d_buff = dates_u[0]
            setattr(self,var+'_d_buff',d_buff)
            next_points = [(d - d_buff).days for d in dates_u] 
            print('Filling '+var+' buffer from '+d_buff.strftime('%Y-%m-%d')+' + ',self.buff_length)
        if not self.buff_vars:
            ### always load
            set_var_now = True
            next_dates = dates_u
        ### if we're buffering, get the data from the buffer
        if self.buff_vars and not set_var_now:
            d_buff = getattr(self,var+'_d_buff')
            d_buff_len = getattr(self,var+'_d_buff_len')
            data_points = [(d - d_buff).days for d in dates_u] 
            ### get all data_points within buffer length
            data_load = [p for p in data_points if p < d_buff_len]
            if len(data_load) == len(data_points):
                ### then all desired data is in the buffer
                x = getattr(self,var)
                if verbos: print(data_load)
                return x[data_load]
            elif len(data_load)>0:
                ### record all the useful stuff
                x = getattr(self,var)
                if verbos: print(data_load)
                x_use =  x[data_load]
                join_data = True
            next_points = [p - d_buff_len 
                           for p in data_points if p >= d_buff_len]
            ### then we need new data
            set_var_now = True
            ### record the new buffer point
            d_buff = dates_u[len(data_load)]
            setattr(self,var+'_d_buff',d_buff)
            print('Filling '+var+' buffer from '+d_buff.strftime('%Y-%m-%d')+' + ',self.buff_length)
        if set_var_now:
            #### need to set var_buff_length
            setattr(self,var+'_d_buff_len',self.buff_length)
            #### we want the dates to load
            if self.buff_vars:
                d_buff = getattr(self,var+'_d_buff')
                ## we want d_buff plus the next load from self.dates
                ### just use days for now
                next_dates = [d_buff+relativedelta(days = i) 
                                 for i in range(self.buff_length)]
                ### actaully we need to check if we've reached the end of the file list
                ### compare end point of next_dates with last entry in de
                if next_dates[-1] > self.de[-1]:
                    next_dates = [nd for nd in next_dates if nd <= self.de[-1]]
                    next_points = next_points[:len(next_dates)]
                    setattr(self,var+'_d_buff_len',len(next_dates))
                    print('Shrinking buffer to '+str(len(next_dates))+
                          ' to fit last date: '+self.de[-1].strftime('%Y-%m-%d'))
            data = self.load_var(var,next_dates,verbos = verbos)
            if not self.buff_vars:
                return data
            if self.buff_vars:
                setattr(self,var,data)
                ### set _d_buff_len from data shape
                setattr(self,var+'_d_buff_len',data.shape[0])
        if join_data:
            return np.vstack((x_use,data[next_points]))
        else:
            return data[next_points]
        
    def clean_var(self,var):
        ### to delete to save memory in big loops
        if hasattr(self,var):
            delattr(self,var)
            print('Cleaning: '+var+', from '+self.name)
        else:
            print('Cleaning: '+self.name+' does not have '+var)
        

    # next function will take a list of dates and return an appropriately orientated arrays
    def load_var(self,var,dates_u,verbos=False):
        ### uses the file_names as they are
        fvar = var.split('/')[1]
        targfvar = self.target.split('/')[1]
        d0  = dates_u[0]
        dn = np.shape(dates_u)[0]
        fn = 0
        looking = True
        load_again = False
        while looking:
            ### look through seached ds,de to get filename
            if d0>=self.ds[fn] and d0<=self.de[fn]:
                looking = False
                file = self.file_names[fn]
                if verbos: print('Target file: '+file)
                file = file.replace(self.target,var)
                file = file.replace(targfvar,fvar)
                if self.search_files:
                    #### replace the d345702347 word from the file with *
                    d_name = file.split('/files/')[1].split('/')[0]
                    if verbos: print('d_name ='+d_name)
                    #### glob.glob the wildcard to find the files
                    file = glob.glob(file.replace(d_name,'*'))[0]
                print('Accessing file: '+file)
                f_nc = Dataset(file)
                if verbos: print(file)
                t_ps = [(d - self.ds[fn]).days for d in dates_u]
                ### we overshoot the file
                if t_ps[-1] > f_nc[fvar].shape[0]:
                    load_again = True
                    ### the extra dates we want
                    dates_extra = [d for d in dates_u if d>self.de[fn]]
                    t_ps = [t for t in t_ps if t < f_nc[fvar].shape[0]]
                ### use orient_opt to load correctly 
                ### now load the next few days
                if self.orient_opt == 1:
                    data_now = f_nc[fvar][t_ps].transpose(0,2,1)
                if self.orient_opt == 2:
                    data_now = f_nc[fvar][t_ps]
                if self.orient_opt == 3:
                    #### unstructured ie AWI
                    data_now = f_nc[fvar][t_ps]
                data_now[data_now.mask]=np.nan
                f_nc.close()
                del f_nc
                gc.collect()
            else:
                fn+=1
        while load_again:
            file = self.file_names[fn+1]
            file = file.replace(self.target,var)
            file = file.replace(targfvar,fvar)
            if self.search_files:
                #### replace the d345702347 word from the file with *
                d_name = file.split('/files/')[1].split('/')[0]
                #### glob.glob the wildcard to find the files
                file = glob.glob(file.replace(d_name,'*'))[0]
            f_nc = Dataset(file)
            if verbos: print(file)
            t_ps = [(d - self.ds[fn+1]).days for d in dates_extra]
            ### we overshoot the file
            if t_ps[-1] > f_nc[fvar].shape[0]:
                load_again = True
                ### the extra dates we want
                dates_extra = [d for d in dates_u if d>self.de[fn+1]]
                t_ps = [t for t in t_ps if t < f_nc[fvar].shape[0]]
                fn+=1
            else: 
                load_again = False
            if self.orient_opt == 1:
                data_extra = f_nc[fvar][t_ps].transpose(0,2,1)
            if self.orient_opt == 2:
                data_extra = f_nc[fvar][t_ps]
            if self.orient_opt == 3:
                #### unstructured ie AWI
                data_extra = f_nc[fvar][t_ps]
            data_extra[data_extra.mask]=np.nan
            f_nc.close()
            del f_nc
            gc.collect()
            data_now = np.vstack((data_now,data_extra))
            
        if looking: 
            print('Something went wrong for '+self.name)
            print('we didn\'t find the array for '+var+' on '+d.strftime('%Y%m%d'))
            return np.array([False])
        else:
            return data_now

    def get_vector(self,vec_names,dates_u,verbos=False):
        data_x =  self.get_var(vec_names[0],dates_u,verbos=verbos)
        data_y =  self.get_var(vec_names[1],dates_u,verbos=verbos)
        if self.orient_opt == 2:
            ux = -data_y
            uy =  data_x
            data_x = ux
            data_y = uy
        if verbos: 
            data = np.hypot(data_x,data_y)
            print('vector: min/mean/max:',['{:.3}'.format(n) 
            for n in [np.nanmin(data),np.nanmean(data),np.nanmax(data),]])
        if self.regrid_read:
            #### loop through time vector
            data_rg_x = []
            data_rg_y = []
            for t in range(np.shape(data_x)[0]):
                x,y = self.regrid_GS2GS.rg_vecs(data_x[t],data_y[t])
                data_rg_x.append(x)
                data_rg_y.append(y)
            return np.array(data_rg_x),np.array(data_rg_y)
        else:
            return data_x,data_y

    def get_aice(self,dates_u,verbos=False):
        data = self.get_var('SIday/siconc',dates_u,verbos=verbos)/100
        if verbos: print('siconc: min/mean/max:',['{:.3}'.format(n) 
            for n in [np.nanmin(data),np.nanmean(data),np.nanmax(data),]])
        if self.regrid_read:
            #### loop through time vector
            data_rg = []
            for t in range(np.shape(data)[0]):
                data_rg.append(self.regrid_GS2GS.rg_array(data[t]))
            return np.array(data_rg)
        else:
            return data

    def get_hi(self,dates_u,verbos=False):
        data = self.get_var('SIday/sithick',dates_u,verbos=verbos)
        if verbos: print('sithick: min/mean/max:',['{:.3}'.format(n) 
            for n in [np.nanmin(data),np.nanmean(data),np.nanmax(data),]])
        if self.regrid_read:
            #### loop through time vector
            data_rg = []
            for t in range(np.shape(data)[0]):
                data_rg.append(self.regrid_GS2GS.rg_array(data[t]))
            return np.array(data_rg)
        else:
            return data

    def get_vels(self,dates_u,verbos=False):
        data_x =  self.get_var('SIday/siu',dates_u,verbos=verbos)
        data_y =  self.get_var('SIday/siv',dates_u,verbos=verbos)
        if self.orient_opt == 2:
            ux = -data_y
            uy =  data_x
            data_x = ux
            data_y = uy
        if verbos: 
            data = np.hypot(data_x,data_y)
            print('sispeed: min/mean/max:',['{:.3}'.format(n) 
            for n in [np.nanmin(data),np.nanmean(data),np.nanmax(data),]])
        if self.regrid_read:
            #### loop through time vector
            data_rg_x = []
            data_rg_y = []
            for t in range(np.shape(data_x)[0]):
                x,y = self.regrid_GS2GS.rg_vecs(data_x[t],data_y[t])
                data_rg_x.append(x)
                data_rg_y.append(y)
            return np.array(data_rg_x),np.array(data_rg_y)
        else:
            return data_x,data_y
        
        

class CMIP6_monthly():
    """
    forcing class for the budget
    lets the forcing load efficiently
    
    """
    def __init__(self,ppath,orient_opt = 1,buffer = False,regrid_read=None,search_files =True):
        ### get name from path
        ### buffing assumes that all files are same size - iceconc, icethick etc
        psplit = ppath.split('CMIP6/')[1] ### gets CMIP/Scenario
        self.group = psplit.split('/')[0]
        self.BIGNAME = psplit.split('/')[1]
        self.name = psplit.split('/')[2]
        self.run_type = psplit.split('/')[3]
        self.ensemble = psplit.split('/')[4]
        self.path = ppath
        self.orient_opt = orient_opt
        self.buff_vars = False
        self.vars_loaded = False
        self.buff_length = 1
        self.search_files = search_files
        if buffer is not False:
            self.buff_vars = True
            self.buff_length = buffer
        if regrid_read is not None:
            self.regrid_read = True
            self.regrid_GS2GS = regrid_read
        else:
            self.regrid_read = False

    def get_dates(self,time_start,time_end,target,append=False):
        """
        returns the all encompassing date list for use with the forcing object
        we find all the start and end time for each file found in the directory
        Using these times we add the relevant daily dates also considering tiem start and end
        
        Parameters
        ---------
        time_start: datetime.datetime the start time point we want to find time points from
        time_end:   datetime.datetime the end time point to stop collecting data points
        target: variable used to investigate the dates available, pick a variable that is likely to be indicative for all availabel data. For example SImon/siconc for monthly sea ice data.
        append:   logical, option, default = False. 
                setting this to True  allows us to add a scenarioMIP to the end of a historical run.
                The list of files and dates generated from a first call, can then be append with a second directory. Do this carefully to make sure the runs are consecutive.
        """
        #### search the directories for all the files
        #### added 2021-08 recursive digging through flies/d837583 structure
        #### then split to sort
        self.target = target
        file_names = glob.glob(self.path+'/'+target+'/gn/files/d*/*.nc',recursive=True)
        try:
            file_names.sort(key = lambda x: x.split('_gn_')[1])
        except IndexError:
            file_names.sort(key = lambda x: x.split('_gr_')[1])
        if len(file_names) == 0:
            file_names = glob.glob(self.path+'/'+target+'/gr/files/d*/*.nc',recursive=True)
            try:
                file_names.sort(key = lambda x: x.split('_gr_')[1])
            except IndexError:
                file_names.sort(key = lambda x: x.split('_gn_')[1])
        ###parent_time_units: days since 1850-01-01-00-00-00
        if append:
            dates= self.dates
            ds = self.ds
            de = self.de
        else:
            dates= []
            ds = []
            de = []
        for file in file_names:
            both_dates = file.split('_'+self.ensemble)[1].split('_')[2].split('.')[0]
            d0 = dt.datetime.strptime(both_dates.split('-')[0],'%Y%m')
            d1 = dt.datetime.strptime(both_dates.split('-')[1],'%Y%m')
            ds.append(d0)
            de.append(d1)
            ### using the d0 and de get
            if time_start<=d1 and time_end>=d0:
                ### then this files contains relevant dates 
                yn = time_start.year- d0.year
                mn = time_start.month-d0.month
                m_start =np.maximum(yn*12 + mn,0)
                yn = time_end.year- d0.year
                mn = time_end.month-d0.month+1
                yn0= d1.year- d0.year
                mn0= d1.month-d0.month+1
                m_end =np.minimum(yn*12 + mn,yn0*12 + mn0)
                [dates.append(d0 + relativedelta(months = m)) 
                         for m in range(m_start, m_end)]
#                 d_start =np.maximum((time_start - d0).days,0)
#                 d_end   =np.minimum((time_end - d0).days,(d1 - d0).days)+1
#                 [dates.append(d0 + relativedelta(days = d)) 
#                          for d in range(d_start, d_end)]
        self.dates = dates
        self.ds = ds
        self.de = de
        if append:
            self.file_names.extend(file_names)
        else:
            self.file_names = file_names
        print(self.name+' Found '+str(np.shape(dates)[0])+' dates')

    # daily points in yearly files

    def get_var(self,var,dates_u,verbos = False):
        ### use info about the var buffer and dates 
        ### either return the next points in the buffer, or 
        ### read the next chunk from the file
        set_var_now = False
        join_data = False
        if self.buff_vars and not hasattr(self,var):
            ### set it all up later
            set_var_now = True
            d_buff = dates_u[0]
            setattr(self,var+'_d_buff',d_buff)
            ### next_points = [(d - d_buff).days for d in dates_u] 
            next_points = []
            for d in dates_u:
                yn = d.year- d_buff.year
                mn = d.month-d_buff.month
                next_points.append(yn*12 + mn )
            print('Filling '+var+' buffer from '+d_buff.strftime('%Y-%m-%d')+' + ',self.buff_length)
        if not self.buff_vars:
            ### always load
            set_var_now = True
            next_dates = dates_u
        ### if we're buffering, get the data from the buffer
        if self.buff_vars and not set_var_now:
            d_buff = getattr(self,var+'_d_buff')
            d_buff_len = getattr(self,var+'_d_buff_len')
            ### data_points = [(d - d_buff).days for d in dates_u] 
            data_points = []
            for d in dates_u:
                yn = d.year- d_buff.year
                mn = d.month-d_buff.month
                data_points.append(yn*12 + mn )
            ### get all data_points within buffer length
            data_load = [p for p in data_points if p < d_buff_len]
            if len(data_load) == len(data_points):
                ### then all desired data is in the buffer
                x = getattr(self,var)
                if verbos: print(data_load)
                return x[data_load]
            elif len(data_load)>0:
                ### record all the useful stuff
                x = getattr(self,var)
                if verbos: print(data_load)
                x_use =  x[data_load]
                join_data = True
            next_points = [p - d_buff_len 
                           for p in data_points if p >= d_buff_len]
            ### then we need new data
            set_var_now = True
            ### record the new buffer point
            d_buff = dates_u[len(data_load)]
            setattr(self,var+'_d_buff',d_buff)
            print('Filling '+var+' buffer from '+d_buff.strftime('%Y-%m-%d')+' + ',self.buff_length)
        if set_var_now:
            #### need to set var_buff_length
            setattr(self,var+'_d_buff_len',self.buff_length)
            #### we want the dates to load
            if self.buff_vars:
                d_buff = getattr(self,var+'_d_buff')
                ## we want d_buff plus the next load from self.dates
                ### just use months for now
                next_dates = [d_buff+relativedelta(months = i) 
                                 for i in range(self.buff_length)]
                ### actaully we need to check if we've reached the end of the file list
                ### compare end point of next_dates with last entry in de
                if next_dates[-1] > self.de[-1]:
                    next_dates = [nd for nd in next_dates if nd <= self.de[-1]]
                    next_points = next_points[:len(next_dates)]
                    setattr(self,var+'_d_buff_len',len(next_dates))
                    print('Shrinking buffer to '+str(len(next_dates))+
                          ' to fit last date: '+self.de[-1].strftime('%Y-%m-%d'))
            data = self.load_var(var,next_dates,verbos = verbos)
            if not self.buff_vars:
                return data
            if self.buff_vars:
                setattr(self,var,data)
                ### set _d_buff_len from data shape
                setattr(self,var+'_d_buff_len',data.shape[0])
        if join_data:
            return np.vstack((x_use,data[next_points]))
        else:
            return data[next_points]

    def clean_var(self,var):
        ### to delete to save memory in big loops
        if hasattr(self,var):
            delattr(self,var)
            print('Cleaning: '+var+', from '+self.name)
        else:
            print('Cleaning: '+self.name+' does not have '+var)
        

    # next function will take a list of dates and return an appropriately orientated arrays
    def load_var(self,var,dates_u,verbos=False):
        ### uses the file_names as they are
        fvar = var.split('/')[1]
        targfvar = self.target.split('/')[1]
        d0  = dates_u[0]
        dn = np.shape(dates_u)[0]
        fn = 0
        looking = True
        load_again = False
        while looking:
            ### look through seached ds,de to get filename
            if d0>=self.ds[fn] and d0<=self.de[fn]:
                looking = False
                file = self.file_names[fn]
                if verbos: print('Target file: '+file)
                file = file.replace(self.target,var)
                file = file.replace(targfvar,fvar)
                if self.search_files:
                    #### replace the d345702347 word from the file with *
                    d_name = file.split('/files/')[1].split('/')[0]
                    if verbos: print('d_name ='+d_name)
                    #### glob.glob the wildcard to find the files
                    file = glob.glob(file.replace(d_name,'*'))[0]
                if verbos: print('Accessing file: '+file)
                f_nc = Dataset(file)
                #### convert to monthly point
                ## t_ps = [(d - self.ds[fn]).days for d in dates_u]
                t_ps = []
                for d in dates_u:
                    yn = d.year- self.ds[fn].year
                    mn = d.month-self.ds[fn].month
                    t_ps.append(yn*12 + mn )
                ### we overshoot the file
                if t_ps[-1] > f_nc[fvar].shape[0]:
                    load_again = True
                    ### the extra dates we want
                    dates_extra = [d for d in dates_u if d>self.de[fn]]
                    t_ps = [t for t in t_ps if t < f_nc[fvar].shape[0]]
                ### use orient_opt to load correctly 
                if verbos: print(t_ps)
                ### now load the next few days
                ## if len(t_ps) == 1:
                ##     t_ps = t_ps[0]
                if self.orient_opt == 1:
                    data_now = f_nc[fvar][t_ps,:,:].transpose(0,2,1)
                if self.orient_opt == 2:
                    data_now = f_nc[fvar][t_ps]
                if self.orient_opt == 3:
                    #### unstructured ie AWI
                    data_now = f_nc[fvar][t_ps,:]
                data_now[data_now.mask]=np.nan
                f_nc.close()
                del f_nc
                gc.collect()
            else:
                fn+=1
        while load_again:
            file = self.file_names[fn+1]
            file = file.replace(self.target,var)
            file = file.replace(targfvar,fvar)
            if self.search_files:
                #### replace the d345702347 word from the file with *
                d_name = file.split('/files/')[1].split('/')[0]
                #### glob.glob the wildcard to find the files
                file = glob.glob(file.replace(d_name,'*'))[0]
            f_nc = Dataset(file)
            if verbos: print(file)
            #### convert to monthly point
            ## t_ps = [(d - self.ds[fn+1]).days for d in dates_extra]
            t_ps = []
            for d in dates_extra:
                yn = d.year- self.ds[fn+1].year
                mn = d.month-self.ds[fn+1].month
                t_ps.append(yn*12 + mn )
            ### we overshoot the file
            if t_ps[-1] > f_nc[fvar].shape[0]:
                load_again = True
                ### the extra dates we want
                dates_extra = [d for d in dates_u if d>self.de[fn+1]]
                t_ps = [t for t in t_ps if t < f_nc[fvar].shape[0]]
                fn+=1
            else: 
                load_again = False
            if verbos: print(t_ps)
            ## if len(t_ps) == 1:
            ##     t_ps = t_ps[0]
            if self.orient_opt == 1:
                data_extra = f_nc[fvar][t_ps,:,:].transpose(0,2,1)
            if self.orient_opt == 2:
                data_extra = f_nc[fvar][t_ps]
            if self.orient_opt == 3:
                #### unstructured ie AWI
                data_extra = f_nc[fvar][t_ps,:]
            data_extra[data_extra.mask]=np.nan
            f_nc.close()
            del f_nc
            gc.collect()
            data_now = np.vstack((data_now,data_extra))
            
        if looking: 
            print('Something went wrong for '+self.name)
            print('we didn\'t find the array for '+var+' on '+d.strftime('%Y%m%d'))
            return np.array([False])
        else:
            return data_now



    def get_vector(self,vec_names,dates_u,verbos=False):
        data_x =  self.get_var(vec_names[0],dates_u,verbos=verbos)
        data_y =  self.get_var(vec_names[1],dates_u,verbos=verbos)
        if self.orient_opt == 2:
            ux = -data_y
            uy =  data_x
            data_x = ux
            data_y = uy
        if verbos: 
            data = np.hypot(data_x,data_y)
            print('vector: min/mean/max:',['{:.3}'.format(n) 
            for n in [np.nanmin(data),np.nanmean(data),np.nanmax(data),]])
        if self.regrid_read:
            #### loop through time vector
            data_rg_x = []
            data_rg_y = []
            for t in range(np.shape(data_x)[0]):
                x,y = self.regrid_GS2GS.rg_vecs(data_x[t],data_y[t])
                data_rg_x.append(x)
                data_rg_y.append(y)
            return np.array(data_rg_x),np.array(data_rg_y)
        else:
            return data_x,data_y
