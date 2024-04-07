
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

ephemeral = str(os.environ['EPHEMERAL'])
wd_hpc=str(os.environ['PBS_O_WORKDIR'])
home=str(os.environ['HOME'])

## Load rock hydraulic conductivity
hyd_cond = xr.open_dataset(home+"/hydrodata/hydrogeo_k.nc").fillna(0) # this is in m/s. jules output is in kg/m2/s=mm/s
# for some reasons rasterization produces NA in the lower limit

jules = xr.open_dataset(ephemeral+"/jules_runs/coupled_jules_oggm_historical_00_18.nc")

### Create deep reservoir

q_deep = np.minimum(jules.sub_surf_roff,1000*hyd_cond.Band1)

q_deep = np.minimum(1000*hyd_cond.Band1,jules.sub_surf_roff)
## this has erroneus dimensions


jules.sub_surf_roff = jules.sub_surf_roff - q_deep ## this will only contain "shallow" subsurface flow
### this 

## Partition the subsurface runoff
uh = pd.read_csv(home+'/hydrodata/unit_hydro.csv').fillna(0) ## unit hydrograph for shallow. this is boris UH for now.


np.array(0,uh["conv"].values)

Unit hydrograph model for shallow

## meta model
uh 

for x 
    for y
        for t

### Will consult with HPC support team. 
## This is an inefficient version and might be better to multiply by a time operator array with timedelta dimension
qochas_area[0,j,i].values


## Delay function for deep subsurface
# V1. Note that here, we lump all the deep drainage and pass it through a unique unit hydrograph.
#       This might only be valid for KM105 station.
Unit hydrograph for deep

Lump all together
Calculate delay


## Implementing PCR-GLOBWB approach
# This is based in Kraijenhoff van de Leur (1958)
# 