## Containerized Development Environment for the RGR NSF-INTERN Project
- Use PDAL to request pointclouds from The National Map and 3DEP Entwine Server
- Generate DEM grids using PDAL and GDAL
- Compute metrics on DEMs 

** Setup is tested for Windows running Docker & VSCode-Remote Containers


## Required Software:
- Docker Desktop installed
- VScode
- VSCode Extensions: (Local) Docker, Dev Containers, WSL (DEV Container) Python, Pylance, isort, Docker, Jupyter

## File List:
- Dockerfile -- builds a miniconda container for developing in VSCode
- devcontainer -- config files for VSCode Remote Container
- src -- folder containing DEM creation & Analysis Modules, tutorial jupyter_notebooks, test scripts, and testData
- .travis.yml, docker-compose.image.yml, docker-compose.yml -- configuration files for Docker.mini.rgr
- data -- folder destination for Dockerfile 'VOLUME'

## Getting started:
- Clone this repo
## For VSCode Container setup:
- Launch Docker Desktop & VSCode
- In VSCode, install extensions (listed above)
- Navigate to repo folder
- Select ><, in Remote options, select 'reopen folder in container'
- Once container is launched, add installed VSCode extensions to the Dev container
- Change Python Interpreter to 'Python 3.9 ('rgr':conda) /opt/conda/bin/python'




