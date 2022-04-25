# here are the classes for holding data
# the first is a class that is a single gridded array
# it's class so that we can access the time variables happily
# first data class is a data_month
# in the dimension [time,x,y] format
# why a class not just any old array?
# so we can have the methods that take averages according to the time scale we want
# yearlies... monthlies... runnning means... all nicely indexed

import numpy as np
import datetime as dt
import copy
from dateutil.relativedelta import relativedelta
# from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset


class data_year:
    # this will init with an input array
    # take the datetime time dims
    # construct a yrpd masked array to allow indexing
    def __init__(self,data,dates,periods = 12):
        """
        data map all tidy like
        initialise with how many data points per year
        default is periods = 12 (monthly)
        the data is the data - just links the object to methods 
        # the date are a list of datetimes, seperated by months
        it's up to the user to supply a list of datetimes, spaced by the correct periods
        ----------
        the datetimes and periods need to comply, if you want 3 monthly seasons
        make sure the datetime lists work for this
        likewise if you have daily data, give the correct periods,
        or data every 10 days, etc etc
        """
        # check if the data timeis the same shape as the date
        # get range of years from dates list
        # 
        [n_t,m,n] = data.shape
        self.n_t = n_t
        self.m = m
        self.n = n
        self.nyrs = dates[-1].year - dates[0].year + 1 
        self.yrpd = np.ma.empty([self.nyrs,periods],dtype = int)
        self.yrpd[:,:] = -1
        self.periods = periods
        self.data = data
        self.mask = False
        self.dates = dates
        self.saved = False
        
        if periods == 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
                mt = self.dates[tt].month - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        elif periods < 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
                mt = int((self.dates[tt].month-1)*periods/12)
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        elif periods > 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
#                 day_of_year = (self.dates[tt] - dt.datetime(self.dates[tt].year,1,1)).days + 1
                day_of_year = self.dates[tt].timetuple().tm_yday - 1
#                 print(day_of_year)
                mt = int(day_of_year*periods/366) - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        self.yrpd.mask = self.yrpd < 0
        self.files = True
        self.years = [self.dates[0].year + n for n in range(self.yrpd.shape[0])]
        
#     def resize(t=[0,0],m=[0,0],n=[0,0]):
#         """
#         reshapes an dy, cutting off t=[beggining,end] from either end of the time array
#         same for m and n dimension.
#         sort out the date strings, n_t, nyrs, and yrpd to be consistent
#         """
        
    
    def get_index(self,dt1,dt2=None):
        """
        if one datetime is given, the nearest index is returned
        if two are given, an array bracketed by the two nearest is returned
        output is time index of data array, for use in time_set for example
        """
        if dt2 is None:
            time_diff = np.abs([d - dt1 for d in self.dates])
            return time_diff.argmin(0)
        else:
            time_diff = np.abs([d - dt1 for d in self.dates])
            p1 = time_diff.argmin(0)
            time_diff = np.abs([d - dt2 for d in self.dates])
            p2 = time_diff.argmin(0)
            return p1,p2
        

    def get_periods(self,dt1,dt2=None,wrap_years = True):
        """
        if one datetime is given, the nearest periods is returned
        if two are given, an array bracketed by the two nearest is returned
        out put is periods, for use in periods for clim_map for example
        default is wrap years = True, so if the dates are on seperate years the periods that are between are selected, if more than a year seperates them, all years are selected.
        """
        if dt2 is None:
            time_diff = np.abs([d - dt1 for d in self.dates])
            x,y = np.where(self.yrpd == time_diff.argmin(0))
            return y[0]
        else:
            time_diff = np.abs([d - dt1 for d in self.dates])
            x1,y1 = np.where(self.yrpd == time_diff.argmin(0))
            time_diff = np.abs([d - dt2 for d in self.dates])
            x2,y2 = np.where(self.yrpd == time_diff.argmin(0))
            if x1==x2:
                return np.arange(y1,y2+1)
            elif x2==x1+1:
                if wrap_years:
                    return np.concatenate((np.arange(0,y2+1),
                                       np.arange(y1,self.periods)))
                else:
                    return np.arange(y1,y2+1)
            else:
                return np.arange(0,self.periods)

    
    def get_years_index(self,years):
        """
        Finds the year index to access the correct row of the yrpd attribute
        Checks to see if the years in question are available

        Parameters
        ---------

        years: list[int] the actual years we find the indexs for 
            typically [2010,2011,2014] for example

        Returns
        -------
        year_indexs: the access indexes for the years in question for the yrpd attribute
            typically [1,2,5] for example

        """
        year0 = self.dates[0].year
        years = [y for y in years if y>= year0]
        return [y-year0 for y in years if y-year0 < self.nyrs]

    def __getitem__(self,indx):
        if type(indx) == dt.datetime:
            t_p = self.get_index(indx)
            m = slice(None)
            n = slice(None)
        elif (type(indx) == int) or (type(indx) == np.int64) or (type(indx) == slice):
            t_p = indx
            m = slice(None)
            n = slice(None)
        else:
            t_p,m,n = indx
        if type(self.mask) == bool:
            # just return the data because it's not masked
            return self.data[t_p,m,n]
        else:
            temp = copy.copy(self.data[t_p,m,n])
            try: temp_mask = self.mask[t_p,m,n]
            except IndexError:
                return self.data[t_p,m,n]
            else: pass
        if type(temp) == np.ndarray:
            temp[temp_mask==False] = np.nan
            return temp
        elif type(temp) == np.ma.core.MaskedArray:
            temp[temp_mask==False] = np.nan
            return temp
        elif temp_mask:
            return temp
        else:
            return np.nan
        
    def print_date(self,t,string='auto',year_only=False):
        """
        Quickly return a date lable from a given data_year time point
        return format can be overidden by setting string to datetime string format
        otherwise it is 'auto'
        year_only = true overides everything and just gives the year
        """
        if type(t) == np.ma.core.MaskedConstant:
            return ''
        # simply get a date string for a time point
        if year_only: # year_only overides
            str_option = '%Y'
        elif string=='auto':
            # auto generate the strftime option from no. of periods
            # if periods = 4 then year + JFM etc...
            # Add this option later use yrpd to find quarter 
                # manually set JFM etc
            # if periods < 12 then months only
            if self.periods <= 12:
                str_option = '%Y-%m'
            elif self.periods <= 366:
                str_option = '%Y-%m-%d'
            # longer then days too
            else:
                str_option = '%Y-%m-%d-T%H'
        else:
            str_option = string
        return self.dates[t].strftime(str_option)
 
    def date_range(self,periods,unit = 'seconds'):
        """
        returns the time difference covered by consecutive periods given in 'periods'
        """
        p0 = periods[0]
        p1 = periods[-1]
        addp = False
        if p1 < self.periods-1: p1+=1
        else: addp = True
        ys = True
        y0 = 0
        while ys:
            if self.yrpd.mask[y0,p0]:
                y0+=1
            else: 
                ys = False
                d0 = self.dates[self.yrpd[y0,p0]]
            if y0>self.nyrs-1:
                print('dates available inconsistent with selected periods, check dy.yrpd')
                return False
        ys = True
        y1 = 0
        while ys:
            if addp:
                if self.yrpd.mask[y1+1,0]:
                    y1+=1
                else:
                    d1 = self.dates[self.yrpd[y1+1,0]]
                if y1>self.nyrs-2 and self.yrpd.mask[y1,-1]:
                    print('dates available inconsistent with selected periods, check dy.yrpd')
                    return False
                elif y1>self.nyrs-2:
                    d1 = (self.dates[self.yrpd[y1,-1]] +
                          (self.dates[self.yrpd[y1,-1]] - 
                           self.dates[self.yrpd[y1,-2]]))
            else:
                if self.yrpd.mask[y1,p1]:
                    y1+=1
                else: 
                    ys = False
                    d1 = self.dates[self.yrpd[y1,p1]]
                if y1>self.nyrs-1:
                    print('dates available inconsistent with selected periods, check dy.yrpd')
                    return False
        if addp and y0==y1-1: d1=d1
        elif addp:
            d1 = d1.replace(year = d0.year+1)
        elif y0==y1: d1=d1
        else:
            d1 = d1.replace(year = d0.year)
        if unit=='seconds':
            return (d1-d0).days*60*60*24
        elif unit=='days':
            return (d1-d0).days

    def build_mask(self):
        """
        fills the mask array with ones if it isn't there yet
        """
        if type(self.mask) == bool:
            self.mask = np.ones(self.data.shape,dtype=bool)
            print('Building new mask array')
        

    def build_static_mask(self,mask,points = False,overwrite=False):
        """
        Makes a satic 2d mask for all data
        or only the time points listed in points option
        mask is 1/True for good, nan/False bad
        overwrite = True, makes the mask identical to input mask
        overwrite = False, appends the current mask to match
        if you want to make a temporal mask for  a condition
        do it your self with logical
        ie.
        DY.build_mask()
        DY.mask[DY.data > limit] = 0
        will temporarily mask out all data  
        over the limit
        """
        if type(self.mask) == bool:
            self.build_mask()
        if (type(mask[0,0])!=bool) and (type(mask[0,0])!=np.bool_):
            print('mask needs to be binary, not',type(mask[0,0]))
            return
        if type(points) == bool:
            points = np.arange(0,self.n_t,dtype=int)
        if overwrite:
            temp_mask = np.ones_like(self.data,dtype=bool)
            for tt in points:
                temp_mask[tt][mask==False] = False
            self.mask = temp_mask
        else:
            for tt in points:
                self.mask[tt][mask==False] = False
        
    def append(self,date,data):
        # check if there's data here already
#         if ~self.files:
#             print("nothing here, so can't append")
#             return False
        # check the new data is the correct size
        m_check,n_check = np.shape(data)
        if m_check==self.m&n_check==self.n:
        # append the data
            self.data = np.append(self.data,data,axis = 0)
#             self.data = np.append(self.data,np.expand_dims(data,axis=0),axis = 0)
            self.dates.append(date)
        # find the final entry in the yrpd
            loc = np.where(self.yrpd == self.n_t)
        # add the next entry(ies)
            if loc[1][0] + nt > self.periods - 1:
                # how many new rows
        # add new rows if needed - keep the yrmth consistent
                self.yrpd = np.ma.append(self.yrpd,
                            np.ma.masked_values(np.ones([1,self.periods]),1),axis=0)
                self.yrpd[-1,0] = self.n_t + 1
            else:
                self.yrpd[-1,loc[1][0]] = self.n_t + 1
            self.n_t += 1
            return True
        else: return True
        # adds another time slice to the data, and sorts out the yrpd array


    def cut_date(self,date1,date2):
        """
        resizes the data to cover only the date range given - inclusive
        remakes yrpd and mask files too
        """
        S1 = self.get_index(date1)
        S2 = self.get_index(date2)
        print("removing "+str(S1+(self.n_t-S2-1))+" time points")
        
        new_data = self.data[S1:S2+1,:,:]
        
        new_dates = self.dates[S1:S2+1]
        
        if type(self.mask) == bool:
            # just return the data because it's not masked
            new_mask = self.mask
        else:
            new_mask = self.mask[S1:S2+1,:,:]
            
        self.data = new_data 
        self.dates= new_dates 
        self.mask = new_mask 
        self.n_t = np.shape(new_data)[0]
        
        # rebuild yrpd
        self.nyrs = self.dates[-1].year - self.dates[0].year + 1 
        self.yrpd = np.ma.empty([self.nyrs,self.periods],dtype = int)
        self.yrpd[:,:] = -1
        
        if self.periods == 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
                mt = self.dates[tt].month - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        elif self.periods < 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
                mt = int(self.dates[tt].month*periods/12) - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        elif self.periods > 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
#                 day_of_year = (self.dates[tt] - dt.datetime(self.dates[tt].year,1,1)).days + 1
                day_of_year = self.dates[tt].timetuple().tm_yday
                mt = int(day_of_year*self.periods/366) - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        self.yrpd.mask = self.yrpd < 0
        print("New data_year, size "+str(self.n_t)+", for "+str(self.nyrs)+" years")

    def mean_series(self,periods=[],mask = False,year_set = [],time_set = [],year_form = 'index', method='mean',mult_array=False,year_select=None,moving_av = False,):
        """
        Cash in on making the effort to organise your yearly gridded data correcty.
        We can select the years and periods of interest and get back a full time series of tht data plus the dates they are defined on. The method for collecting this data can be mean, median, std and sum.
        
        Parameters
        -------
        periods: list or np.array, 1d, ints, the periods we are interested. Default is all the periods. If we are only interested in the januarys for a month data_year, we set periods = [0]
        mask: logical, this will only select data from cell within the data_year.mask array if set
        year_set: list or tuple, year limit numbers [year_start_index, year_end_index], inclusive
            default value of [] gives all years
            yeat_set = [n,n] will consider the year of index n only
        
        time_set: list or tuple, time limit numbers [time_start_index, time_end_index], inclusive
            default value of [] allows all time points to be considered
            time_set = [1,37] will consider the data from data array with time index between 1,37 inclusive.
        method: string, the data combinining method to turn seleected data into a map.
            default = 'mean' uses np.nanmean
             'median' uses np.nanmedian
             'std' uses np.nanstd
             'sum' uses np.nansum
             'yearsum' uses np.nansum
                 note that yearsum method divides result by number of years selected
                 this result can be dodgy if there are empty years of data selected
        mult_array: option, default = False, applies a constant scaling array to the data_year data slices before processing. For example multiply your arrays by a grid area array whilst calculating the sum. Very powerful
        moving_av: collect a moving average of the data for smoothing, 
            default = False, set to an int, to create a mving average of this many years.
        Returns
        -------
        dlist: list, full of date_time.date_time, size (tn)
        numpy.array: a time series Size (tn) where tn very much depends on the options, with max of data_year.n_t
        """
        # check if there's data here already
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
        if np.shape(periods)[0] == 0:
            periods = np.arange(0,self.periods)
        # using all periods within yrpd including empty
        out_list = []
        d_list = []
        dprev = self.dates[0]
        if year_select is not None:
            year_use = [y for y in year_select if (y>=y0 and y<=yE)]
        else:
            year_use = [y for y in range(y0,yE+1)]
        for y in year_use:
            for p in periods:
                if self.yrpd[y,p]>=t0 and self.yrpd[y,p]<=tE:
                    if self.yrpd.mask[y,p]:
                        out_list.append(np.nan)
                        d_list.append(dprev)
                    else:
                        tp = self.yrpd[y,p]
                        temp = self.data[tp].copy()
                        temp[mask[tp]==False] = np.nan
                        if type(mult_array) == np.ndarray:
                            temp = temp*mult_array
                        if method=='mean':
                            out_list.append(np.nanmean(temp))
                        if method=='median':
                            out_list.append(np.nanmedian(temp))
                        if method=='std':
                            out_list.append(np.nanstd(temp))
                        if method=='sum':
                            out_list.append(np.nansum(temp))
                        d_list.append(self.dates[tp])
                        dprev = self.dates[tp]
        if type(moving_av) == bool:
            return d_list, out_list
        else:
            N = moving_av
            temp_padded = np.pad(np.array(out_list), 
                                 (N//2, N-1-N//2), mode='edge')
            temp_array = np.convolve(temp_padded, 
                                     np.ones((N,))/N, mode='valid')
            return d_list, [t for t in temp_array]


    def centile_series(self,centiles,mask = False,year_set = [],time_set = [],year_form = 'index',method='mean',mult_array=False):
        """
        Mask needs to be 1 for true 0 for false
        """
        # check if there's data here already
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            
            # using all periods within yrpd including empty
            out_list = []
            d_list = []
            dprev = self.dates[0]
            for y in range(y0,yE+1):
                for p in range(self.periods):
                    if self.yrpd.mask[y,p]:
                        out_list.append(np.nan)
                        d_list.append(dprev)
                    else:
                        tp = self.yrpd[y,p]
                        temp = self.data[tp]
                        temp[mask[tp]==False] = np.nan
                        if type(mult_array) == np.ndarray:
                            temp = temp*mult_array
                        out_list.append(np.nanpercentile(temp,centiles))
                        d_list.append(self.dates[tp])
                        dprev = self.dates[tp]
        return d_list, out_list


    def clim_mean(self,mask = False,year_set = [],time_set = [],year_form = 'index',method='mean',
                  moving_av = False,full_time=False,mult_array=False,first_p = 0):
        """
        This is the money method, cash in on making the effort to organise your yearly gridded data correcty.
        We can select the years and periods of interest and get back a map of that data. The method for collecting this data can be mean, median, std and sum.
        
        Parameters
        -------
        mask: logical, this will only select data from cell within the data_year.mask array if set
        year_set: list or tuple, year limit numbers [year_start_index, year_end_index], inclusive
            default value of [] gives all years
            yeat_set = [n,n] will consider the year of index n only
            
        year_select: list of year indexes, only these years included.
            default value of [] gives all years
            yeat_set = [1,2,4,5] will consider the years with index [1,2,4,5] only
        
        time_set: list or tuple, time limit numbers [time_start_index, time_end_index], inclusive
            default value of [] allows all time points to be considered
            time_set = [1,37] will consider the data from data array with time index between 1,37 inclusive.
        method: string, the data combinining method to turn seleected data into a map.
            default = 'mean' uses np.nanmean
             'median' uses np.nanmedian
             'std' uses np.nanstd
             'sum' uses np.nansum
             'yearsum' uses np.nansum
                 note that yearsum method divides result by number of years selected
                 this result can be dodgy if there are empty years of data selected
        full_time: option, default = False, makes the method return a value for all periods even if that period never has data. Empty periods have a value of np.nan
        moving_av: option, default = False, applies a moving average to resultant series. Implemented for noisey daily data. setting moving_av =  5 makes a 5 window moving average. Uses a combo of numpy.pad and numpy.convolve to create a even box moving average
        mult_array: option, default = False, applies a constant scaling array to the data_year data slices before processing. For example multiply your arrays by a grid area array whilst calculating the sum. Very powerful
        first_p: option, default = 0, first period point of the year. Setting this to be an integer other than one shifts the end reusult by that many points, cyclicly. Useful if you wish to present a seasonal cycle that doesn't start on the first day of the year, October to September the following year for example.
        Returns
        -------
        numpy.array, a time series Size (tn) where tn very much depends on the options, but tends to be equal to data_year.periods
        """ 
        # check if there's data here already
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            
            # only use the periods that have data
            pmask = np.sum(self.yrpd.mask ,axis = 0)<self.nyrs
            p_use = np.sum(pmask)
            if full_time:
                # locating array to fill at the end
                l = 0
                ploc = []
                for p in pmask:
                    if p:
                        ploc.append(l)
                        l+=1
                    else:
                        ploc.append(None)
            temp_array = np.empty([p_use])
            pind=0
            p_count = np.empty([p_use])
            for mn,t in enumerate(pmask):
                if t:
#             for mn in range(self.periods):
                    idx = self.yrpd[y0:yE+1,mn].compressed()
                    if np.shape(idx)[0] == 0:
                        temp_array[pind] = np.nan
                        pind+=1
                    else:
                        t_mn = np.sum((idx>=t0)&(idx<=tE))
                        temp      = np.empty([t_mn,self.m,self.n])
                        temp_mask = np.empty([t_mn,self.m,self.n],dtype=bool)
                        if type(mult_array) == np.ndarray:
                            temp[:,:,:]      = [self.data[i]*mult_array for i in idx if i>=t0 and i<=tE]
                        else:
                            temp[:,:,:]      = [self.data[i] for i in idx if i>=t0 and i<=tE]
                        p_count[pind] = np.sum([1 for i in idx if i>=t0 and i<=tE])
                        temp_mask[:,:,:] = [mask[i] for i in idx if i>=t0 and i<=tE]
                        temp[temp_mask==False] = np.nan
                        if method=='mean':
                            temp_array[pind] = np.nanmean(temp)
                        elif method=='median':
                            temp_array[pind] = np.nanmedian(temp)
                        elif method=='std':
                            temp_array[pind] = np.nanstd(temp)
                        elif (method=='ysum' or method=='sum'):
                            temp_array[pind] = np.nansum(temp)
                        pind+=1
            if method=='ysum':
                ### then we need to normalize by no.of time points
                temp_array = temp_array/p_count
            if full_time:
                out_array = []
                for p in ploc:
                    if type(p) == int:
                        out_array.append(temp_array[p])
                    else: 
                        out_array.append(np.nan)
                temp_array = out_array
            if type(moving_av) == bool:
                if first_p > 0:
                    return np.roll(temp_array,-first_p)
                else:
                    return temp_array
            else:  # calc moving average
                N = moving_av
                temp_padded = np.pad(temp_array, 
                                     (N//2, N-1-N//2), mode='edge')
                temp_array = np.convolve(temp_padded, 
                                         np.ones((N,))/N, mode='valid')
                if first_p > 0:
                    return np.roll(temp_array,-first_p)
                else:
                    return temp_array

        
    def clim_centile(self,centiles,mask = False,year_set = [],time_set = [],year_form = 'index',method='mean',mult_array=False,first_p = 0):
        
        """
        Mask needs to be 1 for true 0 for false
        """
        # check if there's data here already
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            
            # only use the periods that have data
            pmask = np.sum(self.yrpd.mask ,axis = 0)<self.nyrs
            p_use = np.sum(pmask)
            if full_time:
                # locating array to fill at the end
                l = 0
                ploc = []
                for p in pmask:
                    if p:
                        ploc.append(l)
                        l+=1
                    else:
                        ploc.append(None)
            temp_array = np.empty([p_use,len(centiles)])
            pind=0
            for mn,t in enumerate(pmask):
                if t:
#             for mn in range(self.periods):
                    idx = self.yrpd[y0:yE+1,mn].compressed()
                    t_mn = np.sum((idx>=t0)&(idx<=tE))
                    temp      = np.empty([t_mn,self.m,self.n])
                    temp_mask = np.empty([t_mn,self.m,self.n],dtype=bool)
                    if type(mult_array) == np.ndarray:
                        temp[:,:,:]      = [self.data[i]*mult_array for i in idx if i>=t0 and i<=tE]
                    else:
                        temp[:,:,:]      = [self.data[i] for i in idx if i>=t0 and i<=tE]
                    temp_mask[:,:,:] = [mask[i] for i in idx if i>=t0 and i<=tE]
                    temp[temp_mask==False] = np.nan
                    temp_array[pind,:] = np.nanpercentile(temp,centiles)
                    pind+=1
            if full_time:
                out_array = []
                for p in ploc:
                    if type(p) == int:
                        out_array.append(temp_array[p])
                    else: 
                        out_array.append(np.nan)
                temp_array = out_array
            if first_p > 0:
                return np.roll(temp_array,-first_p)
            else:
                return temp_array



    def clim_map(self,periods=[],mask = False,year_set = [],year_select=None,time_set = [],year_form = 'index',method = 'mean',calc_mask = False):
        """
        This is the money method, cash in on making the effort to organise your yearly gridded data correcty.
        We can select the years and periods of interest and get back a map of that data. The method for collecting this data can be mean, median, std and sum.
        
        Parameters
        -------
        periods: list or numpy.array (1d) of integers
            periods is the list of period no.s to use in the map
        ie. periods = [0,1,2] with give the map of average over
        the first three months of a monthly data_year
            default value is [], which will give all periods
        
        mask: logical, this will only select data from cell within the data_year.mask array if set
            note this will consider all data which is unmasked from any of the data time points considered
            default value mask = False, unmasked
        calc_mask: logical, this will only select data from cells within the data_year.mask array if set
            note this will only select data unmasked the data time points considered. This result can well be different to the mask option if the data mask is highly variable.
            default value calc_mask = False, unmasked
        year_set: list or tuple, year limit numbers [year_start_index, year_end_index], inclusive
            default value of [] gives all years
            yeat_set = [n,n] will consider the year of index n only
            
        year_select: list of year indexes, only these years included.
            default value of [] gives all years
            yeat_set = [1,2,4,5] will consider the years with index [1,2,4,5] only
        
        time_set: list or tuple, time limit numbers [time_start_index, time_end_index], inclusive
            default value of [] allows all time points to be considered
            time_set = [1,37] will consider the data from data array with time index between 1,37 inclusive.
        method: string, the data combinining method to turn seleected data into a map.
            default = 'mean' uses np.nanmean
             'median' uses np.nanmedian
             'std' uses np.nanstd
             'sum' uses np.nansum
             'yearsum' uses np.nansum
                 note that yearsum method divides result by number of years selected
                 this result can be dodgy if there are empty years of data selected
                
        Returns
        -------
        numpy.array, the map. Size (data_year.m, data_year.n)
        """ 
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
        if year_select is not None:
            year_select = [y for y in year_select if (y>=y0 and y<=yE)]
            idx = [[self.yrpd[y,mn] for mn in periods for y in year_select
                             if not self.yrpd.mask[y,mn]]]
            if method == 'ysum':
                ysum = len(year_select)
            else:
                ysum  = 1.0
        else:
            if len(periods) == 0: periods = np.arange(0,self.periods)
            idx = [self.yrpd[y0:yE+1,mn].compressed() for mn in periods]
            if method == 'ysum':
                ysum = yE - y0 +1
            else:
                ysum  = 1.0
        temp_mask = np.sum([mask[j,:,:]
                        for i in idx for j in i if j>=t0 and j<=tE],axis = 0)
        if method=='mean' and calc_mask:
            temp = np.nanmean([self[j] for i in idx for j in i if j>=t0 and j<=tE],axis = 0)
        elif method=='mean':
            temp = np.nanmean([self.data[j] for i in idx for j in i if j>=t0 and j<=tE],axis = 0)
        elif method=='median' and calc_mask:
            temp = np.nanmedian([self[j] for i in idx for j in i if j>=t0 and j<=tE],axis = 0)
        elif method=='median':
            temp = np.nanmedian([self.data[j] for i in idx for j in i if j>=t0 and j<=tE],axis = 0)
        elif method=='std' and calc_mask:
            temp = np.nanstd([self[j] for i in idx for j in i if j>=t0 and j<=tE],axis = 0)
        elif method=='std':
            temp = np.nanstd([self.data[j] for i in idx for j in i if j>=t0 and j<=tE],axis = 0)
        elif (method=='sum' or method=='ysum') and calc_mask:
            temp = np.nansum([self[j] for i in idx for j in i if j>=t0 and j<=tE],axis = 0)/ysum
        elif (method=='sum' or method=='ysum'):
            temp = np.nansum([self.data[j] for i in idx for j in i if j>=t0 and j<=tE],axis = 0)/ysum
        if np.sum(temp_mask<1)>1:
            temp[np.where(temp_mask<1)] = np.nan
        return temp

    def centile_map(self,centile,periods=[],mask = False,year_set = [],year_select=None,time_set = [],year_form = 'index',calc_mask = False):
        """
        This is the money method, cash in on making the effort to organise your yearly gridded data correcty.
        We can select the years and periods of interest and get back a map of that data. 
        
        centile is the centiles you're interested in, ie. 50th, 90th, this takes a single number currently
        
        Parameters
        -------
        periods: list or numpy.array (1d) of integers
            periods is the list of period no.s to use in the map
        ie. periods = [0,1,2] with give the map of average over
        the first three months of a monthly data_year
            default value is [], which will give all periods
        
        mask: logical, this will only select data from cell within the data_year.mask array if set
            note this will consider all data which is unmasked from any of the data time points considered
            default value mask = False, unmasked
        calc_mask: logical, this will only select data from cells within the data_year.mask array if set
            note this will only select data unmasked the data time points considered. This result can well be different to the mask option if the data mask is highly variable.
            default value calc_mask = False, unmasked
        year_set: list or tuple, year limit numbers [year_start_index, year_end_index], inclusive
            default value of [] gives all years
            yeat_set = [n,n] will consider the year of index n only
            
        year_select: list of year indexes, only these years included.
            default value of [] gives all years
            yeat_set = [1,2,4,5] will consider the years with index [1,2,4,5] only
        
        time_set: list or tuple, time limit numbers [time_start_index, time_end_index], inclusive
            default value of [] allows all time points to be considered
            time_set = [1,37] will consider the data from data array with time index between 1,37 inclusive.
                
        Returns
        -------
        numpy.array, the map. Size (data_year.m, data_year.n)
        """ 
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
        if year_select is not None:
            year_select = [y for y in year_select if (y>=y0 and y<=yE)]
            idx = [[self.yrpd[y,mn] for mn in periods for y in year_select
                             if not self.yrpd.mask[y,mn]]]
        else:
            if len(periods) == 0: periods = np.arange(0,self.periods)
            idx = [self.yrpd[y0:yE+1,mn].compressed() for mn in periods]
        temp_mask = np.sum([mask[j,:,:]
                        for i in idx for j in i if j>=t0 and j<=tE],axis = 0)
        if calc_mask:
            temp = np.nanpercentile(
                [self[j] for i in idx for j in i if j>=t0 and j<=tE],
                centile,axis = 0)
        else:
            temp = np.nanpercentile(
                [self.data[j] for i in idx for j in i if j>=t0 and j<=tE],
                centile,axis = 0)
        temp[temp_mask==False] = np.nan
        return temp    


    def trend_map(self,p,mask = False,year_set = [],year_select=None,time_set = [],year_form = 'index',calc_mask = False,return_detrended=False,t_axis = 'year',t_vec_min = 2,verbos=False):
        """
        This is the money method, cash in on making the effort to organise your yearly gridded data correcty.
        We can select the years and periods of interest and get back linear trend of that data map
        Adapted from code by W. Gregory
        
        -------
        period: single integer
            periods is the list of period no.s to use in the map
        ie. periods = [0,1,2] with give the map of average over
        the first three months of a monthly data_year
            default value is [], which will give all periods
        
        mask: logical, this will only select data from cell within the data_year.mask array if set
            note this will consider all data which is unmasked from any of the data time points considered
            default value mask = False, unmasked
        calc_mask: logical, this will only select data from cells within the data_year.mask array if set
            note this will only select data unmasked the data time points considered. This result can well be different to the mask option if the data mask is highly variable.
            default value calc_mask = False, unmasked
        year_set: list or tuple, year limit numbers [year_start_index, year_end_index], inclusive
            default value of [] gives all years
            yeat_set = [n,n] will consider the year of index n only
        time_set: list or tuple, time limit numbers [time_start_index, time_end_index], inclusive
            default value of [] allows all time points to be considered
            time_set = [1,37] will consider the data from data array with time index between 1,37 inclusive.
        t_axis - the xaxis we're computing the trend for. Default is year (only option)
        t_vec_min - the minimum no of. time points we need to compute a trend. Default = 2
                
        Returns
        -------
        numpy.array, the trend map. Size (data_year.m, data_year.n)
        if return_detrended = True
        numpy.array, the trend map. Size (data_year.m, data_year.n)
        numpy.array, the detrend data. Size (n_points,data_year.m, data_year.n), npoints depends on data limit options, max self.nyrs
        list of datetimes, the date points of the slices of the detrended data, size n_points as before
        """ 
        import itertools
        from scipy import stats
        
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
        if year_select is not None:
            year_select = [y for y in year_select if (y>=y0 and y<=yE)]
            idx = [self.yrpd[y,p] for y in year_select
                             if not self.yrpd.mask[y,mn]]
        else:
            idx = self.yrpd[y0:yE+1,p].compressed()
        trend = np.zeros((self.m,self.n,2))*np.nan
        if return_detrended:
            detrended = np.zeros((len(idx),self.m,self.n))*np.nan
            dtvec = [self.dates[d] for d in idx],
        if verbos:
            count_apr = int(np.sum(np.nanmean(self.mask[idx],axis=0)))
            tr_no = 1
            print('Appox. '+str(count_apr)+' trend points to calculate')
            count_apr_o4 = int(count_apr/4)
            verb_count = 0
        if t_axis == 'year':
            tvec  = np.array([self.dates[d].year for d in idx])
        for i,j in itertools.product(range(self.m),range(self.n)):
            if np.array([mask[ii,i,j] for ii in idx]).any():
                if calc_mask:
                    data_lr = np.array([self[d,i,j]      for d in idx])
                else:
                    data_lr = np.array([self.data[d,i,j] for d in idx])
                dmsk = np.isfinite(data_lr)
                if dmsk.sum()>2:
                    tvec_use = tvec[dmsk]
                    data_lr = data_lr[dmsk]
                    trendT, interceptT, r_valsT, probT, stderrT =stats.linregress(
                        tvec_use,data_lr)
                    trend[i,j,0] = trendT
                    trend[i,j,1] = interceptT
                else:
                    trend[i,j,0] = np.nan
                    trend[i,j,1] = np.nan
                if return_detrended:
                    lineT = (trendT*tvec) + interceptT
                    detrended[:,i,j]=data[:,i,j]-lineT
                if verbos:
                    if np.mod(tr_no,count_apr_o4) == 0:
                        print('Update '+str(verb_count)+', point '+str(tr_no))
                        verb_count += 1
                    tr_no += 1
        if return_detrended:
            return trend,detrended,dtvec
        else:
            return trend

    def year_centile(self,centiles,periods=[],mask = False,year_set = [],time_set = [],year_form = 'index',mult_array=False):
        """
        Similar to clim_centile method, but give it a list of periods and it will
        av over the periods and show how they change over the years
        """
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            if len(periods) == 0: periods = np.arange(0,self.periods)
            idx = [self.yrpd[y0:yE+1,mn].compressed() for mn in periods]
            pu,yu = np.shape(idx)
            # only use the periods that have data
            pu,yu = np.shape(idx)
            temp_array = np.empty([yu,len(centiles)])
# idx
            for yy in range(yu):
                temp_mask = np.sum([mask[i[yy]] for i in idx if i[yy]>=t0 and i[yy]<=tE],
                        axis = 0)
                temp = np.nanmean(
                        [self.data[i[yy]] for i in idx if i[yy]>=t0 and i[yy]<=tE],
                        axis = 0)
                if type(mult_array) == np.ndarray:
                    temp = temp*mult_array
                temp[temp_mask==False] = np.nan
                temp_array[yy,:] = np.nanpercentile(temp,centiles)
            return temp_array
        

        
    def year_mean(self,periods=[],mask = False,year_set = [],time_set = [],year_form = 'index',method='mean',mult_array=False):
        """
        Mask needs to be 1 for true 0 for false
        """
        # check if there's data here already
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            if len(periods) == 0: periods = np.arange(0,self.periods)
            idx = [self.yrpd[y0:yE+1,mn].compressed() for mn in periods]
            pu,yu = np.shape(idx)
            # only use the periods that have data
            pu,yu = np.shape(idx)
            temp_array = np.empty([yu])
# idx
            for yy in range(yu):
                temp_mask = np.sum([mask[i[yy]] for i in idx if i[yy]>=t0 and i[yy]<=tE],
                        axis = 0)
                temp = np.nanmean(
                        [self.data[i[yy]] for i in idx if i[yy]>=t0 and i[yy]<=tE],
                        axis = 0)
                if type(mult_array) == np.ndarray:
                    temp = temp*mult_array
                temp[temp_mask==False] = np.nan
                if method=='mean':
                    temp_array[yy] = np.nanmean(temp)
                elif method=='median':
                    temp_array[yy] = np.nanmedian(temp)
            return temp_array
        
    

    def ravel(self,mask = False,remove_nan=False,periods=[],year_set = [],year_form = 'index',time_set = []):
        """
        Mask needs to be 1 for true 0 for false
        """
        # check if there's data here already
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            if len(periods) == 0: periods = np.arange(0,self.periods)
            temp_mask = np.zeros([self.n_t,self.m,self.n],dtype=bool)
#             pmask = np.sum(self.yrpd.mask ,axis = 0)<self.nyrs
#             p_use = np.sum(pmask)
#             temp_array = np.empty([p_use,len(centiles)])
#             pind=0
#             for mn,t in enumerate(pmask):
#                 if t:
            idx = [self.yrpd[y0:yE+1,mn].compressed() for mn in periods]
#             print(idx)
            for i in idx:
                for j in i:
                    if j>=t0 and j<=tE:
                        temp_mask[j]=mask[j]
            if remove_nan:
                temp_mask[np.isnan(self.data)] = False
            return self.data[temp_mask]



    def save(self,filename):
        """saves all the DY in a npz file"""
        if not self.saved:
            temp_time = [self.dates[t].toordinal() for t in range(self.n_t)]
            np.savez_compressed(filename,
                DY_data = self.data,
                DY_mask = self.mask,
                DY_yrpd = self.yrpd,
                DY_n_t  = self.n_t,
                DY_m = self.m,
                DY_n = self.n,
                DY_nyrs = self.nyrs ,
                DY_periods = self.periods ,
                DY_dates = temp_time)
            self.saved = True

    def save_nc(self,filename,d0,data_attr, collect = 'all',tdim='days',save_mask = False,mask=False,
                grid = None,grid_info=False,
                extra_data = [],save_yrpd=False,add_attr=[],
                description='default data_year ',data_name = 'dy_data',
                      year_set = [],time_set = [],year_form = 'index',verbos=True):
        import grid_set as gs
        """saves all the DY in an netcdf file
        Works as default with the dy format, all data saved in a single array
        use the collect option for individual files per date, or per year
        
        parameters
        ----------
        filename: string, name of the written netCDF
        
        d0: datetime.datetime, time start point for the time dimension
        
        data_attr: dict: Attribute information for the variable. For example 
                        {'Description':'data description extra blurb',
                         'Units':'meters'}
        
        Options:
        description = 'whatever'
            writes a description
        data_name = 'whatever'
            changes the name of the main data variable
        data_attr = ['att1_name','the attr']
            the data needs info before you can write it to file
        add_attr = [['att1_name','the attr'],['att2_name','the attr'],.....]
            puts in attributes to the nc file
        extra_data = [gs.grid_set,'label',attr(dict)]
            write another field to the nc file, gs must have identical dimensions
            'label' will be the name of filed in the nc_file
        grid = grid_set.grid_set, add the lon and lat info for a grid_set class
            defualt = False
        grid_info = logical, adds additional gird spacing and orientation information 
            from the grid_set.grid_set class
            defualt = False
        """
        if year_form == 'year':
            if len(time_set) == 2:
                time_set = self.get_index(time_set[0],time_set[1])
            if len(year_set) == 2:
                year_set = self.get_years_index(year_set)
        mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
        
        if collect == 'all': ## SINGLE FILE - filename
            short_file = filename.split('.nc')[0]
            files = [short_file+'.nc']
            y0use = [y0]
            yEuse = [yE]
            t0use = [t0]
            tEuse = [tE]
        elif collect == 'years': ### find all the years
            ### we use the y0,yE options to force it
            short_file = filename.split('.nc')[0]
            y0use = [y for y in range(y0,yE+1)]
            yEuse = [y for y in range(y0,yE+1)]
            t0use = [t0]*((yE-y0)+1)
            tEuse = [tE]*((yE-y0)+1)
            ####
            files = [short_file+'_'+str(self.years[y])+'.nc' 
                      for y in range(y0,yE+1)]
        elif collect == 'slice': ### find every slice
            ### we use the t0,tE options to force it
            short_file = filename.split('.nc')[0]
            y0use = [y0]*((tE-t0)+1)
            yEuse = [yE]*((tE-t0)+1)
            t0use = [t for t in range(t0,tE+1)]
            tEuse = [t for t in range(t0,tE+1)]
            files = [short_file+self.dates[t].strftime('_%Y-%m-%d')+'.nc'
                     for t in range(t0,tE+1)]
        for t0,tE,y0,yE,file in zip(t0use,tEuse,y0use,yEuse,files):
            # build new yrpd
            yrpd_cp = np.ma.empty([yE-y0+1,self.periods],dtype = int)
            yrpd_cp[:,:] = self.yrpd[y0:yE+1,:]
            yrpd_cp.mask[:] = self.yrpd.mask[y0:yE+1,:]
            yrpd_cp.mask[yrpd_cp>tE] = True
            yrpd_cp.mask[yrpd_cp<t0] = True
            t_use = yrpd_cp.compressed()
            if len(t_use) == 0: continue
            if tdim == 'days':
                temp_time = [(self.dates[t] - d0).days for t in t_use]
            if tdim == 'months':
                temp_time = [(self.dates[t] - d0).months for t in t_use]
    #             ydiff =d1.year-d0.year
    #             d1.month-d0.month + ydiff*12
            NC_f = Dataset(file, 'w', format='NETCDF4')
            NC_f.description = description

            NC_f.createDimension('time', np.shape(t_use)[0])
            NC_f.createDimension('x', self.m)
            NC_f.createDimension('y', self.n)
            if save_yrpd:
                NC_f.createDimension('periods', self.periods)
                NC_f.createDimension('nyrs', np.shape(yrpd_cp)[0])
                # to save:
                # A,swh,t0
            DY_data = NC_f.createVariable(data_name, 'f4', ('time','x','y'))
            DY_data.setncatts(data_attr)
            if save_mask:
                DY_mask = NC_f.createVariable('mask', 'i1', ('time','x','y'))
                DY_mask.setncatts({'long_name':'data array binary mask',
                                  'units':'0 - no data, 1 - data present'})
            DY_time = NC_f.createVariable('time', 'f4', ('time',))
            if tdim == 'days':
                DY_time.setncatts({'Time_dimension':'Days since '+
                                      d0.strftime('%Y-%m-%d')})
            if tdim == 'months':
                DY_time.setncatts({'Time_dimension':'Months since '+
                                      d0.strftime('%Y-%m')})
            if tdim == 'seconds':
                DY_time.setncatts({'Time_dimension':'days since '+
                                      d0.strftime('%Y-%m-%dT%H:%M:%S')})
            if save_yrpd:
                DY_yrpd = NC_f.createVariable('dy_yrpd', 'i8', ('nyrs','periods'))
            e_d = []
            for extra_d in extra_data:
                if type(extra_d[0])==vec_data_year:
                    e_d.append([NC_f.createVariable(extra_d[1]+'_x', 'f4', ('time','x','y')),
                                NC_f.createVariable(extra_d[1]+'_y', 'f4', ('time','x','y'))])
                    if verbos: print('Adding vector data')
                    e_d[-1][0].setncatts(extra_d[2])
                    e_d[-1][0].setncatts({'component':'x_component'})
                    e_d[-1][1].setncatts(extra_d[2])
                    e_d[-1][1].setncatts({'component':'y_component'})
                else:
                    e_d.append(NC_f.createVariable(extra_d[1], 'f4', ('time','x','y')))
                    if verbos: print('Adding scalar data')
                    e_d[-1].setncatts(extra_d[2])
    #         DY_yrpd_mask = NC_f.createVariable('dy_yrpd_mask', 'i1', ('nyrs','periods'))

            # save a grid too
            if grid is not None:
                if verbos: print('Adding grid points')
                lons = NC_f.createVariable('lons', 'f4', ('x','y'))
                lons.setncatts({'long_name':'Longitude',
                                'units':'degrees East'})
                lats = NC_f.createVariable('lats', 'f4', ('x','y'))
                lats.setncatts({'long_name':'Latitude',
                                'units':'degrees North'})
                if grid_info:
                    if verbos: print('Adding extended grid information')
                    xdist = NC_f.createVariable('xdist', 'f4', ('x','y'))
                    xdist.setncatts({'long_name':'grid cell spacing, x dimension',
                                'units':'meters'})
                    ydist = NC_f.createVariable('ydist', 'f4', ('x','y'))
                    ydist.setncatts({'long_name':'grid cell spacing, y dimension',
                                'units':'meters'})
                    ang_c = NC_f.createVariable('ang_c', 'f4', ('x','y'))
                    ang_c.setncatts({'long_name':'grid cell orientation, angle between y axis and north, cosine component',
                                'units':'cos(local orientation)'})
                    ang_s = NC_f.createVariable('ang_s', 'f4', ('x','y'))
                    ang_s.setncatts({'long_name':'grid cell orientation, angle between y axis and north, sine component',
                                'units':'sin(local orientation)'})

            # attributes
            for att in add_attr:
                NC_f.setncattr_string(att[0],att[1])
            # Time format attribute
            if tdim == 'days':
                NC_f.setncattr_string('Time_dimension','Days since '+
                                      d0.strftime('%Y-%m-%d'))
            if tdim == 'months':
                NC_f.setncattr_string('Time_dimension','Months since '+
                                      d0.strftime('%Y-%m'))
            if tdim == 'seconds':
                NC_f.setncattr_string('Time_dimension','days since '+
                                      d0.strftime('%Y-%m-%dT%H:%M:%S'))

            # fill variables
            # if masking make the array a np.ma array with the correct mask
            DY_data[:] = [self.data[t,:,:] for t in t_use]
            if save_mask:
                DY_mask[:] = [self.mask[t,:,:] for t in t_use]
            if save_yrpd:
                DY_yrpd[:] = yrpd_cp
            # lat,lon,time(time in timestamp format)
            DY_time[:] = temp_time
            # save a grid too
            if grid is not None:
                lons[:] = grid.lons
                lats[:] = grid.lats
                if verbos: print('Adding grid information')
                if grid_info:
                    xdist[:] = grid.xdist
                    ydist[:] = grid.ydist
                    ang_c[:] = grid.ang_c
                    ang_s[:] = grid.ang_s
            # now the extra data
            for n,extra_d in enumerate(extra_data):
                if type(e_d[n]) == list:
                    ## two vec components
                    e_d[n][0][:] = [extra_d[0].x[t,:,:] for t in t_use]
                    e_d[n][1][:] = [extra_d[0].y[t,:,:] for t in t_use]
                    if verbos: print('Filling vector data')
                else:
                    e_d[n][:] = [extra_d[0].data[t,:,:] for t in t_use]
                    if verbos: print('Filling scalar data')
            NC_f.close()
        

            
#     def save_nc_slice(self,filename_ext,mask = False,
#                       add_attr=[],description='default data_year ',data_name = 'dy_data',
#                       year_set = [],time_set = []):
#         """saves all the DY in an netcdf file
#         Works as default with the dy format, all data saved in a single array
#         use save_nc for bulk save
#         mask,year_set,time_set work as default
        
#         Options:
#         description = 'whatever'
#             writes a description
#         data_name = 'whatever'
#             changes the name of the main data variable
#         add_attr = [['att1_name','the attr'],['att2_name','the attr'],.....]
#             puts in attributes to the nc file
#         """
#         mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
        
#         temp_time = [self.dates[t].toordinal() for t in range(self.n_t)]
        
#         # loop over the required times
#         for .... 
#             # create_name of file from dates
#             # date in filename depends on periods
#             if periods <= 12:
#                 filename = (filename_ext+
#                 self.dates[tt].strftime('%Y%m%')+'.nc')
#             elif periods <= 366:
#                 filename = (filename_ext+
#                 self.dates[tt].strftime('%Y%m%d%')+'.nc')
#             NC_f = Dataset(filename, 'w', format='NETCDF4')
#             NC_f.description = description

#             NC_f.createDimension('x', self.m)
#             NC_f.createDimension('y', self.n)
#             NC_f.createDimension('periods', self.periods)
#             NC_f.createDimension('nyrs', self.nyrs)
#                 # to save:
#                 # A,swh,t0
#             DY_data = NC_f.createVariable(data_name, 'f4', ('x','y'))
#             if self.mask:
#                 DY_mask = NC_f.createVariable('mask', 'i8', ('x','y'))

#             # extra date attribute
#             NC_f.setncattr_string('data_date',self.dates[tt].strftime('%Y%m%d%'))
#             # attributes
#             for att in add_attr:
#                 NC_f.setncattr_string(att[0],att[1])

#             # fill variables
#             # if masking make the array a np.ma array with the correct mask
#             DY_data[:] = self.data[tt,:,:]
#             if self.mask:
#                 DY_mask[:] = self.mask
#             # lat,lon,time(time in timestamp format)



# and now this will do the same but with vectors
class vec_data_year:
    # this will init with an input array
    # take the datetime time dims
    # construct a yrpd masked array to allow indexing
    def __init__(self,datax,datay,dates,periods = 12):
        """
        data map all tidy like
        initialise with how many data points per year
        default is periods = 12 (monthly)
        the data is the data - just links the object to methods 
        # the date are a list of datetimes, seperated by months
        it's up to the user to supply a list of datetimes, spaced by the correct periods
        ----------
        the datetimes and periods need to comply, if you want 3 monthly seasosn
        make sure the datetime lists work for this
        likewise if you have daily data, give the correct periods,
        or data every 10 days, etc etc
        """
        # check if the data timeis the same shape as the date
        # get range of years from dates list
        # 
        
        [n_t,m,n] = datax.shape
        self.n_t = n_t
        self.m = m
        self.n = n
        self.saved = False
        self.nyrs = dates[-1].year - dates[0].year + 1 
        self.yrpd = np.ma.empty([self.nyrs,periods],dtype = int)
        self.mask = False
        self.yrpd[:,:] = -1
        self.periods = periods
        self.dates = dates
        self.x = datax
        self.y = datay
        if periods == 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
                mt = self.dates[tt].month - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        elif periods < 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
                mt = int(self.dates[tt].month*periods/12) - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        elif periods > 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
                day_of_year = self.dates[tt].timetuple().tm_yday
#                 print(day_of_year)
                mt = int(day_of_year*periods/366) - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        self.yrpd.mask = self.yrpd < 0
        self.files = True

    def get_index(self,dt1,dt2=None):
        """
        if one datetime is given, the nearest index is returned
        if two are given, an array bracketed by the two nearest is returned
        output is time index of data array, for use in time_set for example
        """
        if dt2 is None:
            time_diff = np.abs([d - dt1 for d in self.dates])
            return time_diff.argmin(0)
        else:
            time_diff = np.abs([d - dt1 for d in self.dates])
            p1 = time_diff.argmin(0)
            time_diff = np.abs([d - dt2 for d in self.dates])
            p2 = time_diff.argmin(0)
            return p1,p2
        

    def get_periods(self,dt1,dt2=None,wrap_years = True):
        """
        if one datetime is given, the nearest periods is returned
        if two are given, an array bracketed by the two nearest is returned
        out put is periods, for use in periods for clim_map for example
        default is wrap years = True, so if the dates are on seperate years the periods that are between are selected, if more than a year seperates them, all years are selected.
        """
        if dt2 is None:
            time_diff = np.abs([d - dt1 for d in self.dates])
            x,y = np.where(self.yrpd == time_diff.argmin(0))
            return y[0]
        else:
            time_diff = np.abs([d - dt1 for d in self.dates])
            x1,y1 = np.where(self.yrpd == time_diff.argmin(0))
            time_diff = np.abs([d - dt2 for d in self.dates])
            x2,y2 = np.where(self.yrpd == time_diff.argmin(0))
            if x1==x2:
                return np.arange(y1,y2+1)
            elif x2==x1+1:
                if wrap_years:
                    return np.concatenate((np.arange(0,y2+1),
                                       np.arange(y1,self.periods)))
                else:
                    return np.arange(y1,y2+1)
            else:
                return np.arange(0,self.periods)


    def cut_date(self,date1,date2):
        """
        resizes the data to cover only the date range given - inclusive
        remakes yrpd and mask files too
        """
        S1 = self.get_index(date1)
        S2 = self.get_index(date2)
        print("removing "+str(S1+(self.n_t-S2-1))+" time points")
        
        new_datax = self.x[S1:S2+1,:,:]
        new_datay = self.y[S1:S2+1,:,:]
        
        new_dates = self.dates[S1:S2+1]
        
        if type(self.mask) == bool:
            # just return the data because it's not masked
            new_mask = self.mask
        else:
            new_mask = self.mask[S1:S2+1,:,:]
            
        self.x = new_datax
        self.y = new_datay
        self.dates= new_dates 
        self.mask = new_mask 
        self.n_t = np.shape(new_datax)[0]
        
        # rebuild yrpd
        self.nyrs = self.dates[-1].year - self.dates[0].year + 1 
        self.yrpd = np.ma.empty([self.nyrs,self.periods],dtype = int)
        self.yrpd[:,:] = -1
        
        if self.periods == 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
                mt = self.dates[tt].month - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        elif self.periods < 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
                mt = int(self.dates[tt].month*periods/12) - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        elif self.periods > 12:
            for tt in range(self.n_t):
                yr = self.dates[tt].year - self.dates[0].year
#                 day_of_year = (self.dates[tt] - dt.datetime(self.dates[tt].year,1,1)).days + 1
                day_of_year = self.dates[tt].timetuple().tm_yday
                mt = int(day_of_year*self.periods/366) - 1
    #             print(yr,mt,tt)
                self.yrpd[yr,mt] = tt
        self.yrpd.mask = self.yrpd < 0
        print("New vec_data_year, size "+str(self.n_t)+", for "+str(self.nyrs)+" years")

        
    def __getitem__(self,indx):
        if (type(indx) == int) or (type(indx) == np.int64):
            t_p = indx
            m = slice(None)
            n = slice(None)
        else:
            t_p,m,n = indx
        if type(self.mask) == bool:
            # just return the data because it's not masked
            return self.x[t_p,m,n], self.y[t_p,m,n]
        else:
            temp_x = self.x[t_p,m,n]
            temp_y = self.y[t_p,m,n]
            temp_mask = self.mask[t_p,m,n]
        if type(temp_x) == np.ndarray:
            temp_x[temp_mask==False] = np.nan
            temp_y[temp_mask==False] = np.nan
            return temp_x, temp_y
        elif temp_mask:
            return temp_x, temp_y
        else:
            return np.nan, np.nan
 
    def mag(self,indx):
        """
        crude getitem for magnitudes
        indx = [3,item,thing] to get the bits you want
        we cannot pass : 
        instead of ':' write 'slice(None)'
        instead of 'a:b' write 'slice(a,b,None)'
        """
        if (type(indx) == int) or (type(indx) == np.int64):
            t_p = indx
            m = slice(None)
            n = slice(None)
        else:
            t_p,m,n = indx
        if type(self.mask) == bool:
            # just return the data because it's not masked
            return np.hypot(self.x[t_p,m,n], self.y[t_p,m,n])
        else:
            temp_x = self.x[t_p,m,n]
            temp_y = self.y[t_p,m,n]
            temp_mask = self.mask[t_p,m,n]
        if type(temp_x) == np.ndarray:
            temp_x[temp_mask==False] = np.nan
            temp_y[temp_mask==False] = np.nan
            return np.hypot(temp_x, temp_y)
        elif temp_mask:
            return np.hypot(temp_x, temp_y)
        else:
            return np.nan


    def mean_series(self,mask = False,year_set = [],time_set = [],year_select=None,method='mean',mult_array=False,magnitude= False):
        """
        Mask needs to be 1 for true 0 for false
        """
        # check if there's data here already
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            
            # using all periods within yrpd including empty
            out_listx = []
            out_listy = []
            d_list = []
            dprev = self.dates[0]
            if year_select is not None:
                year_use = [y for y in year_select if (y>=y0 and y<=yE)]
            else:
                year_use = [y for y in range(y0,yE+1)]
            if method == 'ysum':
                ysum = len(year_select)
            else:
                ysum  = 1.0
            for y in year_use:
                for p in range(self.periods):
                    if self.yrpd.mask[y,p]:
                        out_listx.append(np.nan)
                        out_listy.append(np.nan)
                        d_list.append(dprev)
                    else:
                        tp = self.yrpd[y,p]
                        tempx = self.x[tp]
                        tempx[mask[tp]==False] = np.nan
                        tempy = self.y[tp]
                        tempy[mask[tp]==False] = np.nan
                        if type(mult_array) == np.ndarray:
                            tempx = tempx*mult_array
                            tempy = tempy*mult_array
                        if magnitude:
                            if method=='mean':
                                out_listx.append(np.nanmean(
                                    np.hypot(tempx,tempy)))
                            if method=='median':
                                out_listx.append(np.nanmedian(
                                    np.hypot(tempx,tempy)))
                            if method=='std':
                                out_listx.append(np.nanstd(
                                    np.hypot(tempx,tempy)))
                        else:
                            if method=='mean':
                                out_listx.append(np.nanmean(tempx))
                                out_listy.append(np.nanmean(tempy))
                            if method=='median':
                                out_listx.append(np.nanmedian(tempx))
                                out_listy.append(np.nanmedian(tempy))
                            if method=='std':
                                out_listx.append(np.nanstd(tempx))
                                out_listy.append(np.nanstd(tempy))
                        d_list.append(self.dates[tp])
                        dprev = self.dates[tp]
        if magnitude:
            return d_list, out_listx
        else:
            return d_list, out_listx, out_listy


    def centile_series(self,centiles,mask = False,year_set = [],time_set = [],method='mean',mult_array=False):
        """
        Mask needs to be 1 for true 0 for false
        """
        # check if there's data here already
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            
            # using all periods within yrpd including empty
            out_listx = []
            out_listy = []
            d_list = []
            dprev = self.dates[0]
            for y in range(y0,yE+1):
                for p in range(self.periods):
                    if self.yrpd.mask[y,p]:
                        out_list.append(np.nan)
                        d_list.append(dprev)
                    else:
                        tp = self.yrpd[y,p]
                        tempx = self.x[tp]
                        tempx[mask[tp]==False] = np.nan
                        tempy = self.y[tp]
                        tempy[mask[tp]==False] = np.nan
                        if type(mult_array) == np.ndarray:
                            tempx = tempx*mult_array
                            tempy = tempy*mult_array
                        if magnitude:
                            out_listx.append(np.nanpercentile(
                                np.hypot(tempx,tempy),centiles))
                        else:
                            out_listx.append(np.nanpercentile(tempx,centiles))
                            out_listy.append(np.nanpercentile(tempy,centiles))
                        d_list.append(self.dates[tp])
                        dprev = self.dates[tp]
        if magnitude:
            return d_list, out_listx
        else:
            return d_list, out_listx, out_listy



    def print_date(self,t,string='auto',year_only=False):
        """
        Quickly return a date lable from a given data_year time point
        return format can be overidden by setting string to datetime string format
        otherwise it is 'auto'
        year_only = true overides everything and just gives the year
        """
        # simply get a date string for a time point
        if year_only: # year_only overides
            str_option = '%Y'
        elif string=='auto':
            # auto generate the strftime option from no. of periods
            # if periods = 4 then year + JFM etc...
            # Add this option later use yrpd to find quarter 
                # manually set JFM etc
            # if periods < 12 then months only
            if self.periods <= 12:
                str_option = '%Y-%m'
            elif self.periods <= 366:
                str_option = '%Y-%m-%d'
            # longer then days too
            else:
                str_option = '%Y-%m-%d-T%H'
        else:
            str_option = string
        return self.dates[t].strftime(str_option)
 


    def build_mask(self):
        """
        fills the mask array with ones if it isn't there yet
        """
        if type(self.mask) == bool:
            self.mask = np.ones(self.x.shape,dtype=bool)
        


    def build_static_mask(self,mask,points = False,overwrite=False):
        """
        Makes a satic 2d mask for all data
        or only the time points listed in points option
        mask is 1/True for good, nan/False bad
        overwrite = True, makes the mask identical to input mask
        overwrite = False, apends the current mask to match
        if you want to make a temporal mask for  a condition
        do it your self with logical
        ie.
        DY.build_mask()
        DY.mask[DY.data > limit] = 0
        will temporarily mask out all data  
        over the limit
        """
        if type(self.mask) == bool:
            self.build_mask()
        if (type(mask[0,0])!=bool) and (type(mask[0,0])!=np.bool_):
            print('mask needs to be binary, not',type(mask[0,0]))
            return
        if type(points) == bool:
            points = np.arange(0,self.n_t,dtype=int)
        if overwrite:
            temp_mask = np.ones_like(self.x,dtype=bool)
            for tt in points:
                temp_mask[tt][mask==False] = False
            self.mask = temp_mask
        else:
            for tt in points:
                self.mask[tt][mask==False] = False
        
    def append(self,date,datax,datay):
        # check if there's data here already
        if ~self.files:
            print("nothing here, so can't append")
            return False
        # check the new data is the correct size
        m_check,n_check = np.shape(datax)
        if m_check==self.m&n_check==self.n:
        # append the data
            self.x = np.append(self.x,np.expand_dims(datax,axis=0),axis = 0)
            self.y = np.append(self.y,np.expand_dims(datay,axis=0),axis = 0)
            self.dates.append(date)
        # find the final entry in the yrpd
            loc = np.where(self.yrpd == self.n_t)
        # add the next entry(ies)
            if loc[1][0] == self.periods - 1:
        # add new rows if needed - keep the yrmth consistent
                self.yrpd = np.ma.append(self.yrpd,
                            np.ma.masked_values(np.ones([1,self.periods]),1),axis=0)
                self.yrpd[-1,0] = self.n_t + 1
            else:
                self.yrpd[-1,loc[1][0]] = self.n_t + 1
            self.n_t += 1
            return True
        else: return True
        # adds another time slice to the data, and sorts out the yrpd array

    def clim_map(self,periods,mask = False,magnitude = False,year_set = [],time_set = []):
        """
        periods is the list of period no.s to use in the map
        ie. periods = [0,1,2] with give the map of average over
        the first three months of a monthly data_year
        setting magnitude = True, takes the average of the vector
        hypot, rather than the average of each component
        """ 
        # check if there's data here already
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            if len(periods) == 0: periods = np.arange(0,self.periods)
            idx = [self.yrpd[y0:yE+1,mn].compressed() for mn in periods]
            temp_mask = np.sum([mask[j,:,:]
                            for i in idx for j in i if j>=t0 and j<=tE],axis = 0)
            if magnitude:
                temp_x = np.nanmean(
                [np.hypot(self.x[j],self.y[j]) for i in idx for j in i if j>=t0 and j<=tE],
                            axis = 0)
                temp_x[temp_mask==False] = np.nan
                return temp_x
            else:
                temp_x = np.nanmean(
                    [self.x[j] for i in idx for j in i if j>=t0 and j<=tE],
                            axis = 0)
                temp_y = np.nanmean(
                    [self.y[j] for i in idx for j in i if j>=t0 and j<=tE],
                            axis = 0)
                temp_x[temp_mask==False] = np.nan
                temp_y[temp_mask==False] = np.nan
                return temp_x,temp_y

            # check if there's data here already

    def clim_mean(self,mask = False,magnitude = False,year_set = [],time_set = [],method='mean',first_p = 0):
        """
        Mask needs to be 1 for true 0 for false
        """
        # check if there's data here already
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            temp_x = np.empty([self.periods])
            temp_y = np.empty([self.periods])
            if magnitude:
                for mn in range(self.periods):
                    idx = self.yrpd[y0:yE+1,mn].compressed()
                    t_mn = np.sum((idx>=t0)&(idx<=tE))
                    temp_x1   = np.empty([t_mn,self.m,self.n])
                    temp_y1   = np.empty([t_mn,self.m,self.n])
                    temp_mask = np.empty([t_mn,self.m,self.n],dtype=bool)
                    temp_x1[:,:,:]   = [self.x[i]    for i in idx if i>=t0 and i<=tE]
                    temp_y1[:,:,:]   = [self.y[i]    for i in idx if i>=t0 and i<=tE]
                    temp_mask[:,:,:] = [self.mask[i] for i in idx if i>=t0 and i<=tE]
                    temp = np.hypot(temp_x1,temp_y1)
                    temp[temp_mask==False] = np.nan
                    if method=='mean':
                        temp_x[mn] = np.nanmean(temp)
                    elif method=='median':
                        temp_x[mn] = np.nanmedian(temp)
                if first_p > 0:
                    return np.roll(temp_x,-first_p)
                else:
                    return temp_x
            else:
                for mn in range(self.periods):
                    idx = self.yrpd[y0:yE+1,mn].compressed()
                    t_mn = np.sum((idx>=t0)&(idx<=tE))
                    temp      = np.empty([t_mn,self.m,self.n])
                    temp_mask = np.empty([t_mn,self.m,self.n],dtype=bool)
                    temp[:,:,:]      = [self.x[i]    for i in idx if i>=t0 and i<=tE]
                    temp_mask[:,:,:] = [self.mask[i] for i in idx if i>=t0 and i<=tE]
                    temp[temp_mask==False] = np.nan
                    if method=='mean':
                        temp_x[mn] = np.nanmean(temp)
                    elif method=='median':
                        temp_x[mn] = np.nanmedian(temp)
                    temp[:,:,:]      = [self.y[i]    for i in idx if i>=t0 and i<=tE]
                    temp[temp_mask==False] = np.nan
                    if method=='mean':
                        temp_y[mn] = np.nanmean(temp)
                    elif method=='median':
                        temp_y[mn] = np.nanmedian(temp)
                return temp_x,temp_y
                if first_p > 0:
                    return np.roll(temp_x,-first_p),np.roll(temp_y,-first_p)
                else:
                    return temp_x,temp_y

    def ravel(self,mask = False,magnitude = False,
              periods=[],year_set = [],time_set = []):
        """
        Mask needs to be 1 for true 0 for false
        """
        # check if there's data here already
        if self.files:
            mask,y0,yE,t0,tE = get_range_mask(self,mask,year_set,time_set)
            temp_mask = np.zeros([self.n_t,self.m,self.n],dtype=bool)
            idx = [self.yrpd[y0:yE+1,mn].compressed() for mn in periods]
#             print(idx)
#             print(idx)
            for i in idx:
                for j in i:
                    if j>=t0 and j<=tE:
                        temp_mask[j]=mask[j]
            if magnitude:
                return np.hypot(self.x[temp_mask],self.y[temp_mask])
            else:
                return self.x[temp_mask],self.y[temp_mask]

    def save(self,filename):
        """saves all the DY in a npz file"""
        if not self.saved:
            temp_time = [self.dates[t].toordinal() for t in range(self.n_t)]
            np.savez_compressed(filename,
                DY_x = self.x,
                DY_y = self.y,
                DY_mask = self.mask,
                DY_yrpd = self.yrpd,
                DY_n_t  = self.n_t,
                DY_m = self.m,
                DY_n = self.n,
                DY_nyrs = self.nyrs ,
                DY_periods = self.periods ,
                DY_dates = temp_time)
            self.saved = True







def load_data_year(filename):
    """
    give a save filename with a saved data year
    creates a data year from what was saved
    """
    try:
        npzfile =  np.load(filename)
        data = npzfile["DY_data"] 
    except KeyError:
        print('file in incorrect format')
        pass
    else:
        temp_time = npzfile["DY_dates"] 
        dates = [dt.datetime.fromordinal(t) for t in temp_time]
        DY_out = data_year(data,dates)
        mask = npzfile["DY_mask"]
        if np.shape(mask) == ():
            DY_out.mask = False
        else:
            DY_out.mask = mask
        DY_out.periods = npzfile["DY_periods"] 
        DY_out.nyrs = npzfile["DY_nyrs"] 
        DY_out.yrpd = np.ma.empty([DY_out.nyrs,DY_out.periods],dtype = int)
        temp_yrpd = npzfile["DY_yrpd"] 
        npzfile.close()
        DY_out.yrpd[:,:] = temp_yrpd[:,:]
        DY_out.yrpd.mask = DY_out.yrpd < 0
#         DY.n_t = npzfile["DY_n_t"]  
#         self.m = npzfile["DY_m"] = m
#         self.n = npzfile["DY_n"] = n
        DY_out.saved = True
        return DY_out


def load_vec_data_year(filename):
    """
    give a save filename with a saved vector data year
    creates a vector data year from what was saved
    """
    try:
        npzfile =  np.load(filename)
        x = npzfile["DY_x"] 
        y = npzfile["DY_y"] 
    except KeyError:
        pass
    else:
        temp_time = npzfile["DY_dates"] 
        dates = [dt.datetime.fromordinal(t) for t in temp_time]
        DY_out = vec_data_year(x,y,dates)
        mask = npzfile["DY_mask"]
        if np.shape(mask) == ():
            DY_out.mask = False
        else:
            DY_out.mask = mask
        DY_out.periods = npzfile["DY_periods"] 
        DY_out.nyrs = npzfile["DY_nyrs"] 
        DY_out.yrpd = np.ma.empty([DY_out.nyrs,DY_out.periods],dtype = int)
        temp_yrpd = npzfile["DY_yrpd"] 
        npzfile.close()
        DY_out.yrpd[:,:] = temp_yrpd[:,:]
        DY_out.yrpd.mask = DY_out.yrpd < 0
#         DY.n_t = npzfile["DY_n_t"]  
#         self.m = npzfile["DY_m"] = m
#         self.n = npzfile["DY_n"] = n
        DY_out.saved = True
        return DY_out
    

def get_range_mask(dy,mask,year_set,time_set):
    """
    simplifies the preamble from before
    """
    if type(mask)==bool:
        if mask and type(dy.mask)== bool:
            print("Data year mask not set: ignoring it")
            mask = np.ones([dy.n_t,dy.m,dy.n],dtype=bool)
        elif mask:
            mask = dy.mask
        else:
            mask = np.ones([dy.n_t,dy.m,dy.n],dtype=bool)
    elif (np.shape(mask)[1] != dy.m) or (np.shape(mask)[2] != dy.n):# check mask dimension)
        print("Mask array incorrect shape, ignoring it")
        mask = np.ones([dy.n_t,dy.m,dy.n],dtype=bool)
    # build a year limit if it's wanted
    if np.size(year_set)==0:
        y0 = 0
        yE = dy.nyrs-1
    # check if months is in the correct form
    elif year_set[0]>dy.nyrs-1:
        print("year_set inconistent with data, ignoring it")
        y0 = 0
        yE = dy.nyrs-1
    else:
        y0 = year_set[0]
        yE = np.min([dy.nyrs,year_set[1]])
    # build a time limit if it's wanted
    if np.size(time_set)==0:
        t0 = 0
        tE = dy.n_t-1
    # check if months is in the correct form
    elif time_set[0]>dy.n_t:
        print("time_set inconistent with data, ignoring it")
        t0 = 0
        tE = dy.n_t-1
    else:
        t0 = time_set[0]
        tE = np.min([dy.n_t,time_set[1]])
    return mask,y0,yE,t0,tE

# def join(dyl):
#     """
#     take a list of dy and join them together into one dy
#     CAREFUL, if the dates are too misaligned this will be wierd
#     CAREFUL, if there are dulicate dates they will be over written in theyrpd
#     """
#     n_j =  np.shape(dyl)[0]
#     # find the n_yrs for the combined dy
#     y0 = 3000
#     ye = 0
#     for nj in range(n_j):
#         y0 = np.min(dyl[nj].dates[0].year,y0)
#         ye = np.max(dyl[nj].dates[0].year+dyl[nj].nyrs,ye)
#     periods = 
        
def clim_mean_combine2(dy1,dy2,op,
                       mask = False,year_set = [],time_set = [],method='mean',mult_array=False,diag = False):
    """
    performs a mean_sereis operation on two data_years,
    coaligns them, returns the clim mean
    functionalitly is as the mean_series and clim_mean methods
    Not fully tested, but uses dy1 to find dates to get dy2.
    developed with 2 similar data years in mind. Don't go mixing ones with very different ranges. 
    Only use with dy's with the same amount of periods.
    """
    ### first call a mean_series on dy1
    dd1,tt1 = dy1.mean_series(method=method,mask=mask,mult_array=mult_array,
                           year_set=year_set,time_set=time_set)
    ### use the output of the first one to call the mean_series ont he second
    d2t1 = dy2.get_index(dd1[0])
    d2t2 = dy2.get_index(dd1[-1])
    
    dd2,tt2 = dy1.mean_series(method=method,mask=mask,mult_array=mult_array,
                           time_set=[d2t1,d2t2])
    ### check alignment
    len1 = np.shape(tt1)
    len2 = np.shape(tt2)
    if len1!=len2: print('Warning dy\' misaligned')
    
    ### print dates and all that if diag = True
    if diag:
        print(dd1.strftime('%Y%m%d-'),dd2.strftime('%Y%m%d'))
    
    ### apply the opp and format the output like a clim_mean
    
