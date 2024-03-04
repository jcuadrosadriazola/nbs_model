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

jules = rxr.open_rasterio(ephemeral+"/jules_runs/coupled_jules_oggm_historical_00_18.nc")
jules["fao_et0"].attrs['units']="kg m-2 s-1"
jules["surf_roff"].attrs['units']="kg m-2 s-1"
jules["sub_surf_roff"].attrs['units']="kg m-2 s-1"
jules["rain"].attrs['units']="kg m-2 s-1"
jules["melt"].attrs['units']="kg m-2 s-1"

# hydraulic conductivity
satcon = xr.open_dataset(ephemeral+"/data/jules_soil_props_2015_rosetta3_ESA_rahu.nc")['satcon']

# Load HydroSheds dem and accumulation
dem = rxr.open_rasterio(home+'/hydrosheds/sa_dem_3s.tif').isel(band=0).rio.clip_box(minx=-74,miny=-15,maxx=-70,maxy=-12)
acc =  rxr.open_rasterio(home+'/hydrosheds/sa_acc_3s.tif').isel(band=0).rio.clip_box(minx=-74,miny=-15,maxx=-70,maxy=-12)

# Load qochas arrays
qochas_n = rxr.open_rasterio(home+'/qochas/qochas_n.nc',masked=True)
qochas_acc = rxr.open_rasterio(home+'/qochas/qochas_acc.nc',masked=True)
qochas_area = rxr.open_rasterio(home+'/qochas/qochas_area.nc',masked=True)
qochas_cap = rxr.open_rasterio(home+'/qochas/qochas_cap.nc',masked=True)
grid_area = rxr.open_rasterio(home+'/qochas/grid_area.nc',masked=True)

#t1= time.time()
## Create variables
St=0 #Qochas starts empty # CG: is this the volume of water in the qocha at time t? yes
Qav= 0 #Dummy variable for maximum runoff available, CG: maximum Qin to the qocha
R=0 #Dummy for recharge
Ks=0 #Dummy for hydraulic conductivity

#print("a")
## xarray.size output coordinates size for each dimensions. HOWEVER, indexing can be tricky. We will follow xarray.DataArray structure, i.e. time, y, x for indexing
for i in range(jules.sizes["x"]):
  for j in range(jules.sizes["y"]):
    if qochas_n[0,j,i]>0:
      Ks=satcon[0,j,i].values
      vmax=qochas_cap[0,j,i].values
      qocha_area=qochas_area[0,j,i].values
      cell_area=grid_area[0,j,i].values
      Ac=qochas_acc[0,j,i].values

      # diagnostic purposes
      print('Running cell',i,",",j)
      #we now read the timeseries
      for t in range(jules.sizes["time"]):
        # variables
        Qs=jules.surf_roff[t,j,i]
        Qss=jules.sub_surf_roff[t,j,i]
        Et=jules.fao_et0[t,j,i]
        
        ### Water balance in the qocha will be in [m3]. we will then update JULES fluxes in [kg/m2/s]. We will introduce 
        # The recharge will equal the minimum between the maximum infiltration potential (assuming the qocha is full) and the remaining water in the qocha
        #     Remember that we assume that the water above the bottom of the qocha is fully saturated. i.e. Darcy law applies
        R = min(Ks*86.4*qocha_area, St) # Convert kg/m2/s to m3 infiltrating daily
        # With the qocha partially drained we will calculate the maximum amount of water available (Qav) at step t
        Qav = St + Ac * Qs * 86.4 - R - Et * qocha_area * 86.4 # CG modified 
        # We should consider here a better estimate of Et for open water, it should be higher than the one outputted in JULES for reference crop

        # If the water available is greater than the qocha capacity, it will overflow and the qocha will be totally filled, i.e. St=vmax
        if Qav>vmax:
          Qin=vmax-St
          St=vmax
        else:
          St=Qav
          Qin=Ac*Qs*86.4
        
        #Correct Qs and Qss
        jules.surf_roff[t,j,i]=Qs-Qin*qocha_area/cell_area/86.4
        jules.sub_surf_roff[t,j,i]=Qss+R*qocha_area/cell_area/86.4
      
      St=0
      R=0
      Qav=0

jules.to_netcdf(home+'/runs/mockrun.nc')