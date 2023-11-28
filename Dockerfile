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

WORKDIR /rgr/src
ENV PATH /opt/conda/envs/base/bin:$PATH

RUN conda update conda \
    && conda update --all \
    && conda config --add channels conda-forge \
    && conda config --add channels default
RUN conda install jupyter matplotlib numpy ipympl numexpr scipy shapely pyproj
RUN pip install geopandas 
RUN conda install -c conda-forge gdal
RUN conda install -c conda-forge python-pdal

RUN conda clean --all

CMD ["/bin/bash"]

