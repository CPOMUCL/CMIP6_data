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