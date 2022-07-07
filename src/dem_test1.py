# %%
# DEM snatch test --- code grabbed from 'demonstrating requesting pointclouds.iynb'
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

# %%
##_________________________________________________
#First, example creating a DEM a heigt above ground model - this may take a while

testExtent = ([-105.602157,-105.548255],[37.606802,37.649582]) #Define an extent to request
testEPSG = 4326 #Define the SRS of the extent
outEPSG = 32613 #Define the output SRS
print(outEPSG)

grid_name_dem = os.path.join('testData','dem_test1_slv') #Where do we want to save the data for the DEM?
#grid_name_hag = os.path.join('testData','dem_test1_slv_hag') #Where do we want to save the data for the HAG?
usgs_survey_name = 'CO_San-Luis-Valley_2011' #What is the name of the dataset we are requesting from? Found at: https://usgs.entwine.io/
cell_size = 2.0 #Pixel size of resultant raster

# %%
#### Build DEM
#First get the JSON request description
pipeline = pCC.request_build_dem_creation_pipeline(testExtent,testEPSG,outEPSG,grid_name_dem,
                                                   usgs_survey_name,cell_size = cell_size)
print("pipeline defined")

# %%
#this step takes a while
pCC.run_pipeline(pipeline) #Then execute this request
print("dem request executed")

# %%
# ISSUE 7/5/22: Running hag crashes-- Dev Env disconnects
### Build a height above ground raster (e.g., cannopy height)
#First get the JSON request description
#pipeline = pCC.request_build_hag_creation_pipeline(testExtent,testEPSG,outEPSG,grid_name_hag,usgs_survey_name,
#                                                   cell_size = cell_size)
#pCC.run_pipeline(pipeline) #Execute the request
#print("hag request executed")

# %%
#Now lets visualize what 
dem= dg.demGrid([],rasterPath=grid_name_dem+'.tif')
dem[dem <0] = np.nan
hs = dem.calcHillshade()

#hag = dg.demGrid([],rasterPath = grid_name_hag+'.tif')

# %%
#First lets load in a trace of the fault scarp above
import geopandas as gpd

pathToFaultTrace = os.path.join('testData','vectorData','SLV_Zapata_FaultScarp.shp')

fsGDF = gpd.read_file(pathToFaultTrace)
f,axs = plt.subplots(1,1,figsize = (6,6),dpi = 200,sharex=True,sharey=True)

#Plot the hillshade
hs.plotGrid(axs = axs,cmap = 'gray',vmin = 0, vmax = 255)

#Plot the DEM as two different very transparent colormaps, just to play with overlay effects
dem.plotGrid(axs = axs,cmap = 'Blues_r',alpha = 0.1, vmin = 2340, vmax = 2750)
dem.plotGrid(axs = axs,cmap = 'Oranges',alpha = 0.1, vmin = 2500, vmax = 2950)

#Plot the fault scarp
fsGDF.plot(ax = axs,color = 'k')

#Remove axis labels for clean axes
axs.set_xticklabels([])
axs.set_yticklabels([])

# %%
#Here we will request a series of orthogonal profiles to the line above
from scipy.special import erf
from scipy import optimize
from importlib import reload
reload(pCC)

nSwaths = 5 #Number of orthogonal profiles we will construct
swathPositions = np.linspace(0,1,nSwaths) #Relative positions along the line of the swaths

faultLine = fsGDF.geometry[0] #Get the line geometry

swathLength = 200.0
swathWidth = 3.0

f,axs = plt.subplots(1,nSwaths,figsize = (10,3),dpi = 90,sharex = True,sharey = True)

#Create some functions for fitting fault scarp morphology
def fsFun(L,params):
    #Solution for a scarp on an inclined surface smoothed by linear diffusion
    a,kt,b,c = params
    return a*erf((L)/(2.0*np.sqrt(kt))) + b*(L) + c

def fitFS(L,Z):
    #Optimizing for the best fitting paramaters in the above scarp profile solution
    objFun = lambda params : np.sum((Z - fsFun(L,params))**2)
    paramGuess = [1.0,10,0.0,np.mean(Z)]
    
    res = optimize.minimize(objFun,paramGuess).x
    return res


for i,relative_position in enumerate(swathPositions):
    
    
    #Execute the pipeline
    pipe = pCC.get_orthogonal_swathprofile_request_pointcloud(faultLine,
                                                             relative_position, swathLength,
                                                             swathWidth,outEPSG,
                                                             usgs_survey_name,
                                                             'profileTest',doSavePointCloud=False,
                                                            doReclassify = True)
    #Get the relevant data from the pipeline
    arrays = pipe.arrays[0]
    L = np.array(arrays['L']) - swathLength/2.0
    Z = np.array(arrays['Z'])
    C = np.array(arrays['Classification'])
    D = np.array(arrays['D'])
    X = np.array(arrays['X'])
    Y = np.array(arrays['Y'])
    
    #Which points are ground:
    isGrnd = C == 2
    
    #Fit the fault scarp
    bfParams = fitFS(L[isGrnd],Z[isGrnd])
    
    #Create a predicted FS line
    LtoFit = np.linspace(-swathLength/2.0,swathLength/2.0,50)
    Zfit = fsFun(LtoFit,bfParams)
    
    #Plot the result
    axs[i].plot(L[~isGrnd],Z[~isGrnd],'.',color = 'darkgreen',alpha =0.2)
    axs[i].plot(L[isGrnd],Z[isGrnd],'.',color = 'k',alpha =0.2)
    axs[i].plot(LtoFit,Zfit,'-r',label = 'A: {:.1f},\n kt: {:.1e},\n b: {:.2f}'.format(bfParams[0],bfParams[1],bfParams[2]))
    
    axs[i].legend(fontsize = 'x-small')