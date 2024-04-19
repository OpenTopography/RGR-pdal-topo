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

RUN conda create --name rgr python=3.8

RUN echo "source activate rgr" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

RUN conda install -n rgr jupyter matplotlib numpy ipympl numexpr scipy shapely pyproj scikit-image scikit-learn
RUN conda install -n rgr -c conda-forge gdal
RUN conda install -n rgr -c conda-forge python-pdal
RUN pip install geopandas 

RUN conda clean --all

CMD ["/bin/bash"]
