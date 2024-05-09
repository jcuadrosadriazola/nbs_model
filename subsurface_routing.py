
#import libraries
import xarray as xr
import rasterio
from rasterio.plot import show
from matplotlib import pyplot
import os
import rioxarray as rxr
import pysheds
from pysheds.grid import Grid
import numpy as np
import pandas as pd

ephemeral = str(os.environ['EPHEMERAL'])
wd_hpc=str(os.environ['PBS_O_WORKDIR'])
home=str(os.environ['HOME'])

## Load rock hydraulic conductivity
hyd_cond = xr.open_dataset(home+"/hydrodata/hydrogeo_k.nc").fillna(0) # this is in m/s. jules output is in kg/m2/s=mm/s
# for some reasons rasterization produces NA in the lower limit

jules = xr.open_dataset(ephemeral+"/jules_runs/coupled_jules_oggm_historical_00_18.nc")


## I am copying the hydraulic conductivity values across time dimension. This allows us to do array operations in xarray without dimension problems
hyd_cond0 = xr.zeros_like(jules.sub_surf_roff)

for t in range(jules.sizes["time"]):
    hyd_cond0[t,...]=hyd_cond.Band1.values
### Create deep reservoir

q_deep = xr.zeros_like(jules.sub_surf_roff)
# calculate deep percolation. Hydraulic conductivity is multiplied by 1000 (m/s to mm/s )
q_deep[...]=np.minimum(jules.sub_surf_roff,1000*hyd_cond0)

## correct shallow subsurface flow
jules.sub_surf_roff.values = np.add(jules.sub_surf_roff,-q_deep)
## Read unit hydrograph
uh = pd.read_csv(home+'/hydrodata/unit_hydro.csv').fillna(0) ## unit hydrograph for shallow. this is boris UH for now.
## unit hydrograph i,j wise for shallow subsurface flow
shallow_asnp = jules.sub_surf_roff.values
shallow_asnp = np.append(shallow_asnp,np.zeros((365,75,102)),axis=0)
shallow_asnp1 = np.zeros_like(shallow_asnp)
for i in range(jules.sizes["lon"]):
    for j in range(jules.sizes["lat"]):q
        for t in range(jules.sizes["time"]):
            if shallow_asnp[t,j,i] > 0:
                shallow_asnp1[t:t+uh["conv"].size,j,i]+= shallow_asnp[t,j,i]*uh["conv"].values
                

# this will correct jules subsurface flow to account only for shallow partition
jules.sub_surf_roff.values = shallow_asnp1[:jules.sizes["time"],...]
jules.to_netcdf(home+'/runs/shallowrouted.nc')

# V1. Note that here, we lump all the deep drainage and pass it through a unique unit hydrograph.
#       This might only be valid for KM105 station.
#Unit hydrograph for deep
# lumping flows to a daily timeseries
#q_deep_lump=q_deep.sum(dim=["lat","lon"]) # it would be better to import pixel area raster and multiply to have a water balance
## pending implementation of UH





## Implementing PCR-GLOBWB approach
# This is based in Kraijenhoff van de Leur (1958)
# 