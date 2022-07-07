import setuptools

setuptools.setup(
    name="demCreationComponents",
    version="1.0.0",
    author="Sam Johnstone",
    author_email="sjohnstone@usgs.gov",
    description="DEM Creation Components",
    packages=['demCreationComponents'],
    install_requires=[
                      'numpy', 'json', 'pdal', 'gdal', 'pyproj', 'shapely',                     
                      ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)