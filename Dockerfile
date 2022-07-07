## Directory structure:
#rgr --
#       --src
#              --jupyter_notebooks
#              --Modules
#                       --demCreationComponents
#                       --demAnalysisComponents
#              --testData
#              --dem_test1.py
#
#       --Dockerfile

FROM continuumio/miniconda3

WORKDIR /rgr
RUN conda update conda \
  && conda update --all \
  && conda config --add channels conda-forge \
  && conda config --add channels default \
  #install base level packages:
  && conda install jupyter matplotlib numpy numexpr gdal scipy shapely pyproj geopandas\
  #pdal needs super special install, otherwise won't work
  && conda install -c conda-forge python-pdal\
  && conda clean --all

VOLUME /data

EXPOSE 8888
COPY . .
