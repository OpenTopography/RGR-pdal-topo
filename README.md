# RGR-pdal-topo
RGR-pdal-topo is a containerized development environment for processing point cloud data into raster grids, computing topographic metrics from DEMs, and performing basic DEM classiciation for geological mapping.
This repository contains Python modules, scripts, and notebooks that: 
- Use PDAL to request data from The National Map and the Entwine Amazon S3 bucket that hosts USGS 3DEP lidar point clouds
- Generate DEM grids from point clouds using PDAL and GDAL-based python modules
- Compute topographic metrics on DEMs (e.g. surface roughness, slope, curvature, shaded-relief)
- Perform basic pixel classification using sci-kitlearn
- The provided example notebooks show analysis of Rio Grande rift normal faults and alluvial fan deposits.

** Setup is tested for Windows running Docker & VSCode-Remote Containers

## Required Software:
- Docker Desktop installed
- VScode
- VSCode Extensions: (Local) Docker, Dev Containers, WSL (DEV Container) Python, Pylance, isort, Docker, Jupyter

## Table of contents:
- Dockerfile -- builds a miniconda container for developing in VSCode
- devcontainer -- config files for VSCode Remote Container
- src -- folder containing DEM creation & Analysis Modules, tutorial jupyter_notebooks, test scripts, and testData
## src Table of Contents
- Modules -- folder containing the demAnalysisComponents and demCreationComponents modules
- DEMs -- folder containing the raster and vector data files used in the provided notebook
- geo7x_sangres_all_shp -- shapefile containing GPS points collected for use in the Projecting_GPS_points_to_line notebook
- profiles_for_project -- gps data used in profile notebooks
- testData -- datasets for the tutorial notebooks
- tutorial_notebooks -- notebooks for utilizing AWS point cloud and DEM creation/analysis modules
- ClusterScarp -- notebook demonstrating Kmeans clustering of fault scarp derivative layers
- FanRoughness_Stats -- notebook demonstrating computing and binning roughness values associated with geological map units
- PointCloud_Profiles -- notebook demonstrating plotting and analysis of elevation profiles across fault scarps
- Projecting_GPS_points_to_line -- notebook demonstrating how to project GPS points to a profile line shp

## Getting started:
- Clone this repo
## For VSCode Container setup:
- Launch Docker Desktop & VSCode
- In VSCode, install extensions (listed above)
- Navigate to repo folder
- Select ><, in Remote options, select 'reopen folder in container'
- Once container is launched, add installed VSCode extensions to the Dev container
- Change Python Interpreter to 'Python 3.8 ('rgr':conda) /opt/conda/bin/python'




