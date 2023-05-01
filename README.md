## Containerized Development Environment for the RGR NSF-INTERN Project
- request pointclouds from 3DEP to create grids
- compute metrics on DEMs
- Setup is tested for Windows running Docker & VSCode-Remote Containers


## Required Software:
- Docker Desktop installed
- VScode
- VSCode Extensions: Docker, Remote Containers, Python, Jupyter

## File List:
- Dockerfile -- builds a miniconda container for developing in VSCode
- devcontainer -- config files for VSCode Remote Container
- src -- folder containing DEM creation & Analysis Modules, tutorial jupyter_notebooks, test scripts, and testData
- .travis.yml, docker-compose.image.yml, docker-compose.yml -- configuration files for Docker.mini.rgr
- data -- folder destination for Dockerfile 'VOLUME'

## Getting started:
- Clone this repo
## For VSCode Container setup (recommended):
- Launch Docker Desktop & VSCode
- In VSCode, install extensions (listed above)
- Navigate to repo folder
- Select ><, in Remote options, select 'reopen folder in container'
- Once container is launched, add installed VSCode extensions to the Dev container
- Change Python Interpreter to 'Python 3.9.13 ('base':conda) /opt/conda/bin/python'

## For Jupyter Notebook container setup:
- Run `docker-compose up --build` in terminal
- copy and paste last url into browser to launch local jupyter notebook
- Note: when re-launching container after exiting, you may need to enter the token for access (located in the launch url)

## Default Settings
- Port: `8888:8888`
- Volume: `./src:/src`
- --notebook-dir: `/src



