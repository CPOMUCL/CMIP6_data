import numpy as np
import grid_set as gs
from netCDF4 import Dataset

def gs_from_CMIP(m,filename,opt,grid_list=False):
    ### we will auto generate
    G = gs.grid_set(m)
    anames = [
        ['latitude','longitude'],
        ['lat','lon'],
        ['nav_lat','nav_lon'],
    ]
    G_nc = Dataset(filename)
    for an in anames:
        try:
            lon_in= G_nc.variables[an[1]][:]
            lat_in= G_nc.variables[an[0]][:]
            print('Got variables for a grid')
            print(an)
#             print(lon_in.shape,lat_in.shape)
            if opt == 1:
                lon_in= lon_in.T
                lat_in= lat_in.T
            elif opt == 2:
                lon_in,lat_in = np.meshgrid(lon_in,lat_in)
            break
        except KeyError:
            continue

    G_nc.close()
    lat_in[lat_in> 89.9] = 89.9
    lat_in[lat_in<-89.9] =-89.9

    lon_in[lon_in>180.0] -= 360.0

    G.set_grid_lon_lat(lon_in,lat_in,grid_list=grid_list)
    return G



class accumulator:
    def __init__(self,shape):
        self.count = np.zeros(shape,dtype = int)
        self.data = np.ma.masked_all(shape)
        self.data[:] = 0.0
        self.data.mask[:] = True
        
    def update(self,new_data,mask=None,count_all = True):
        new_mask = np.isnan(new_data)
        if mask is not None:
            new_mask[mask] = True
        self.data.mask[~new_mask] = False
        self.data[~new_mask] += new_data[~new_mask]
        if count_all: self.count += 1
        else: self.count += ~new_mask
    
    def clean(self):
        self.count[:] = 0
        self.data[:] = 0.0
        self.data.mask[:] = True
        
    def mean(self):
        out_array = self.data.data/self.count
        out_array[self.data.mask] = np.nan
        return out_array
        
    def total(self):
        out_array = self.data.data
        out_array[self.data.mask] = np.nan
        return out_array
    

class recent_array:
    """
    list of nlist arrays that can be updated and the mean returned
    give the shape of each array, and the no. of slices to remember
    """
    def __init__(self,shape,nlist):
        self.shape = (nlist,)+shape
        self.nlist = nlist
        self.data = np.zeros(self.shape)*np.nan
    def update(self,array):
        ### push each entry along the list by copying entries
        for n in range(self.nlist-1):
            self.data[n,:] = self.data[n+1,:] 
        self.data[self.nlist-1,:] = array[:]
    def mean(self):
        return np.nanmean(self.data,axis = 0)
    def clean(self):
        for n in range(self.nlist):
            self.data[:] = np.zeros(self.shape)*np.nan
        