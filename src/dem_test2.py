# Testing topographic metrics on test DEM (see dem_test1 & demonstrating requesting pointclouds)
#%%
import sys
import os

#Add the path to the demAnalysis & Components scripts
sys.path.append('Modules/demAnalysisComponents/') 
sys.path.append('Modules/demCreationComponents/')

import dem as dg
import pointCloudCreation as pCC

from matplotlib import pyplot as plt
import numpy as np

from importlib import reload
import geopandas as gpd
from scipy.special import erf
from scipy import optimize

#%%
grid_name_dem = os.path.join('testData','dem_test1_slv')
dem= dg.demGrid([],rasterPath=grid_name_dem+'.tif')
dem[dem <0] = np.nan
hs = dem.calcHillshade()
f,axs = plt.subplots(1,1,figsize = (6,6),dpi = 200,sharex=True,sharey=True)

#Plot the hillshade
hs.plotGrid(axs = axs,cmap = 'gray',vmin = 0, vmax = 255)

#Plot the DEM as two different very transparent colormaps, just to play with overlay effects
dem.plotGrid(axs = axs,cmap = 'Blues_r',alpha = 0.1, vmin = 2340, vmax = 2750)
dem.plotGrid(axs = axs,cmap = 'Oranges',alpha = 0.1, vmin = 2500, vmax = 2950)
#Remove axis labels for clean axes
axs.set_xticklabels([])
axs.set_yticklabels([])

# %%
#Paramaters for DEM derivatives
HS_AZ = 315 #Hillshade azimuth
HS_ELV = 45 #HS elevation angle
TPI_IN_MULT = 2 #Multiplier on cell size for inner annulus radius of TPI
TPI_OUT_MULT = 4 #Multiplier on cell size for outer annulus radius of TPI
SLOPE_WIN = 1 #Window size for slope magnitude kernel
CRV_WIN = 2 #Window size for second derivative kernel

print('Loaded DEM, Building derivatives...')

#Slope magnitude
print('Building Slope magnitude...')
Smag = dem.calcFiniteSlopeMagOverWindow(SLOPE_WIN)
Smag.overwriteSaveFile()

#%%
#Visualize the result - Slope magnitude
f,axs = plt.subplots(1,1,dpi = 120)


Smag.plotGrid(axs = axs, cmap = 'Reds',vmin = 0, vmax = 1, alpha = 1)
#Remove axis labels for clean axes
axs.set_xticklabels([])
axs.set_yticklabels([])

# %%
#Curvature
print('Building Curvature...')
Crv = dem.calcFiniteLaplacianOverWindow(CRV_WIN)
Crv.overwriteSaveFile()
Crv_name = Crv._filePath

#Visualize the result - Slope magnitude
f,axs = plt.subplots(1,1,dpi = 120)

Crv.plotGrid(axs = axs, cmap = 'coolwarm',vmin = -1, vmax = 1, alpha = 1)
#Remove axis labels for clean axes
axs.set_xticklabels([])
axs.set_yticklabels([])
# %%
