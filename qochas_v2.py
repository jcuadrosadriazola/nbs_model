#import libraries
import xarray as xr
import rasterio
from rasterio.plot import show
from matplotlib import pyplot
import os
import rioxarray as rxr
import geopandas as gpd
import numpy as np
import pandas as pd
import scipy.stats as stats
from shapely.geometry import mapping
import xesmf as xe

ephemeral = str(os.environ['EPHEMERAL'])
wd_hpc=str(os.environ['PBS_O_WORKDIR'])
home=str(os.environ['HOME'])

## jules output. change accordingly
jules = xr.open_dataset(home+"/hydrodata/jules_runs/coupled_for_nbs/rcp85_ACCESS1-3_coupled_jules_oggm_00_99.nc",decode_coords="all").sel(time=slice("2000-01-01", "2018-12-31"))

# import high res dem
dem = rxr.open_rasterio(home+'/hydrosheds/sa_dem_3s.tif').isel(band=0).rio.clip_box(minx=-74,miny=-16,maxx=-70,maxy=-12)

## hydraulic conductivity
satcon = xr.open_dataset("hydrodata/soil/jules_soil_props_2015_rosetta3_ESA_rahu_modified_v2.nc")['satcon']


## regrid dem 
threshold = 4000
regridder=xe.Regridder(dem,jules,"bilinear")
dr_out=regridder((dem>threshold).astype(float),keep_attrs=True)


### Gamma parameters fitted in R (fitdistr). 
### qochas contribution area (ha): shape=1.197986 , rate= 10.179759
### qochas area (m2/1000): shape=1.050518 ,rate=8.509921 
### qochas volume capacity (m3): shape=1.777842, rate=7.424725
### rescaling: qochas_acc max= ,qochas_area max= ,qochas_cap max= 
qochas_n = 8
qochas_prop = ["qochas_acc","qochas_area","qochas_cap"]
gamma_shapes = [1.197986,1.050518,1.777842]
gamma_rates = [10.179759,8.509921,7.424725]
rescaling = [275.5467,184568.9,83112.31] #This are the maximum of each property to re-scale. Qochas contribution area is further multiplied by 10,000 to convert from ha to m2
gamma_df =pd.DataFrame(list(zip(qochas_prop,gamma_shapes,gamma_rates)),columns=["Properties","Shape","Rate"])

## masking cells where there exists flows over the entire time
masking=jules.surf_roff.sum(dim="time")
masking=masking.where(masking>0)


## ------ Generating qochas rasters ------- ##
qochas_cap=xr.zeros_like(dr_out)

qochas_acc=xr.zeros_like(dr_out)
np.random.seed(0) ## set seed to get pseudo-random numbers (same results whenever this is run)
## allocate random fields to qochas properties
#qochas accumulation area
a=np.zeros(qochas_cap.shape)
for i in range(qochas_n):
    a+=np.random.gamma(shape=gamma_df["Shape"][0] ,scale=1/gamma_df["Rate"][0] ,size=qochas_cap.shape)
qochas_acc.values=a
qochas_acc.values=(qochas_acc*dr_out).values*rescaling[0]
qochas_acc.values=qochas_acc.where(masking>0).values*10000

qochas_area=xr.zeros_like(dr_out)
#qochas area
a=np.zeros(qochas_cap.shape)
for i in range(qochas_n):
    a+=np.random.gamma(shape=gamma_df["Shape"][1] ,scale=1/gamma_df["Rate"][1] ,size=qochas_cap.shape)
qochas_area.values=a
qochas_area.values=(qochas_area*dr_out).values*rescaling[1]
qochas_area.values=qochas_area.where(masking>0).values

qochas_cap=xr.zeros_like(dr_out)
#qochas volume capacity
a=np.zeros(qochas_cap.shape)
for i in range(qochas_n):
    a+=np.random.gamma(shape=gamma_df["Shape"][2] ,scale=1/gamma_df["Rate"][2] ,size=qochas_cap.shape)
qochas_cap.values=a
qochas_cap.values=(qochas_cap*dr_out).values*rescaling[2]
qochas_cap.values=qochas_cap.where(masking>0).values


### Zero arrays for various intermediate variables
R=xr.zeros_like(satcon[0,...])
St=xr.zeros_like(jules.surf_roff) #Storage at time step t
Qav = xr.zeros_like(satcon[0,...]) # storage volume + contributing runoff - losses (Et + drainage). 
R = xr.zeros_like(jules.surf_roff) # Potential recharge at time step t
Qin = xr.zeros_like(satcon[0,...])


### Qocha model

cell_area=15526711 #grid_area[0,j,i].values for grid area
for t in range(jules.sizes["time"]):
    R[t,...]=np.minimum(86.4*satcon[0,...]*qochas_area,St[t,...])
    Qav= (St[t,...] + 
                 qochas_acc * jules.surf_roff[t,...] * 86.4 - 
                 R[t,...] - 
                 np.minimum(jules.fao_et0[t,...] * qochas_area * 86.4,
                            St[t,...]))
    Qin = (xr.where(Qav>qochas_cap, qochas_cap-St[t,...], qochas_acc*jules.surf_roff[t,...]*86.4))
    jules.surf_roff[t,...]= (jules.surf_roff[t,...]-(Qin.values/cell_area/86.4))
    jules.sub_surf_roff[t,...] = (jules.sub_surf_roff[t,...]+R[t,...]/cell_area/86.4)
    if t<jules.sizes["time"]-1:
        St[t+1,...] = (xr.where(Qav>qochas_cap, qochas_cap,Qav))

### Add shallow UH as needed