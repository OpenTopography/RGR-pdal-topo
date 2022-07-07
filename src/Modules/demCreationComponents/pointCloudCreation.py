
import numpy as np
import json
import pdal
from osgeo import gdal
import pyproj # Python interface to PROJ (cartographic projections and coordinate transformations library)
from shapely.geometry import Point, Polygon, LineString

#TODO: don't require name of resources
'''
Update code to accomodate variable extents
https://github.com/hobu/usgs-lidar/blob/master/boundaries/resources.geojson
'''
ENTWINE_EPSG = 3857 #Expected EPSG code for Entwine point tiles
DEFAULTGRIDEXT = '.tif' #Desired output of grids

'''
Helpful install hint: On windows, installing with anaconda, my install often won't appropriately identify
two key pdal components: filters.python and readers.numpy.  Re-specifying the PDAL_PATH in anaconda seemed
to solve this for me:


conda env config vars set PDAL_DRIVER_PATH= "path-to-PDAL environment\Lib\site-packages\ bin" fixed the issue.

In that folder you should see two files, libpdal_plugin_filter_python.dll and libpdal_plugin_reader_numpy.dll

This sol'n was found here: https://github.com/PDAL/python/issues/64

You can confirm correct behavior by looking for filters.python in the list of drivers available when executing

pdal --drivers

from an anaconda prompt

'''

########################################################################################################
###
########################################################################################################

def projectPointsOntoLine(lineXY,xyPoints):
    '''

    :param lineXY:
    :type lineXY:
    :param xyPts:
    :type xyPts:
    :return:
    :rtype:
    '''

    #Preallocate space for data
    projectedDistances = np.zeros(len(xyPoints) )* np.nan
    alongProfileDistances = np.zeros_like(projectedDistances)

    # Preallocate some space to check for how far and if things are in bounds
    isProjected = np.zeros_like(projectedDistances) == 1

    L = 0

    # For each line segment
    for i in range(len(lineXY) - 1):


        # Get the squared distance along this line
        l2 = np.sum((lineXY[i] - lineXY[i + 1]) ** 2)

        # For each point
        for j, p in enumerate(xyPoints):

            # If point is not yet projected
            if not isProjected[j]:

                # Determine fraction along line segment that this point lies
                t = np.sum((p - lineXY[i]) * (lineXY[i + 1] - lineXY[i])) / l2

                # If point projection falls onto segment
                if (t >= 0) and (t <= 1):
                    # Mark point as projected
                    isProjected[j] = True

                    # Project the point
                    projPoint = lineXY[i] + t * (lineXY[i + 1] - lineXY[i])

                    # Get distance to projected point from original point
                    projectedDistances[j] = np.sqrt(np.sum((projPoint - p) ** 2))

                    # Store projected point and distance to start of line
                    alongProfileDistances[j] = L + np.sqrt(np.sum((projPoint - lineXY[i]) ** 2))

        # Calculate along line distance from the start of this line segment to start of line
        L += np.sqrt(np.sum((lineXY[i] - lineXY[i + 1]) ** 2))

    return projectedDistances, alongProfileDistances

def projectPointsOntoLine_pipeline(ins,outs):

    lineXY = pdalargs['line_coords']

    xs = ins['X']
    ys = ins['Y']
    xyPoints = [np.array([x,y]) for x,y in zip(xs,ys)]
    lineXY = [np.array([x, y]) for x, y in lineXY]

    projectedDistances, alongProfileDistances = projectPointsOntoLine(lineXY,xyPoints)

    outs['D'] = projectedDistances
    outs['L'] = alongProfileDistances

    return True


def reproject_extent(extent,in_epsg,out_epsg):
    '''

    :param extent: ([minx, maxx], [miny, maxy]).
    :param in_epsg:
    :param out_epsg:
    :return:
    '''

    fromCRS = pyproj.CRS("EPSG:{}".format(in_epsg))
    toCRS = pyproj.CRS("EPSG:{}".format(out_epsg))

    transformer = pyproj.Transformer.from_crs(fromCRS, toCRS, always_xy=True)

    LL = transformer.transform(extent[0][0],extent[1][0])
    UR = transformer.transform(extent[0][1],extent[1][1])

    extent_reproj = ([LL[0],UR[0]],[LL[1],UR[1]])

    return extent_reproj

def reprojectXYPoints(xyPoints,in_epsg,out_epsg):
    '''

    :param xyPoints:
    :type xyPoints:
    :param in_epsg:
    :type in_epsg:
    :param out_epsg:
    :type out_epsg:
    :return:
    :rtype:
    '''
    fromCRS = pyproj.CRS("EPSG:{}".format(in_epsg))
    toCRS = pyproj.CRS("EPSG:{}".format(out_epsg))

    transformer = pyproj.Transformer.from_crs(fromCRS, toCRS, always_xy=True)

    outPoints = [transformer.transform(xyPoint[0],xyPoint[1]) for xyPoint in xyPoints]

    return outPoints
########################################################################################################
###
########################################################################################################

def request_pointcloud_pipeline(extent, extent_epsg, out_epsg, usgs_survey_name, grid_name, doReclassify = False,
                                doSavePointCloud = False, pointResolution = None):
    '''

    :param extent:
    :param extent_epsg:
    :param out_epsg:
    :param usgs_survey_name:
    :return:
    '''
    if extent_epsg == ENTWINE_EPSG:
        ent_extent = extent
    else:
        ent_extent = reproject_extent(extent, extent_epsg, ENTWINE_EPSG)

    filename = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/{}/ept.json".format(usgs_survey_name)
    CRS = "EPSG:{}".format(out_epsg)

    build_pipeline = {"pipeline": [
        {
            "bounds": str(ent_extent),
            "filename": filename,
            "type": "readers.ept",
            "tag": "readdata"
        },
        {
            "limits": "Classification![7:7]",
            "type": "filters.range",
            "tag": "nonoise"
        }

    ]}

    if not(pointResolution is None):
        build_pipeline['pipeline'][0]['resolution'] = pointResolution

    build_pipeline = append_pointcloud_manipulation_pipeline(build_pipeline,extent,extent_epsg,out_epsg,filename,
                                                             grid_name,doReclassify,doSavePointCloud)

    return build_pipeline



def load_pointcloud_pipeline(extent, extent_epsg, out_epsg, filename, grid_name, doReclassify = False,
                                doSavePointCloud = False):
    '''

    :param extent:
    :param extent_epsg:
    :param out_epsg:
    :param usgs_survey_name:
    :return:
    '''
    if extent_epsg == ENTWINE_EPSG:
        ent_extent = extent
    else:
        ent_extent = reproject_extent(extent, extent_epsg, ENTWINE_EPSG)

    CRS = "EPSG:{}".format(out_epsg)

    build_pipeline = {"pipeline": [
        {
            "filename": filename,
            "type": "readers.las",
            "tag": "readdata"
        },
        {
            "limits": "Classification![7:7]",
            "type": "filters.range",
            "tag": "nonoise"
        }

    ]}

    build_pipeline = append_pointcloud_manipulation_pipeline(build_pipeline,extent,extent_epsg,out_epsg,filename,
                                                             grid_name,doReclassify,doSavePointCloud)
    return build_pipeline

def append_pointcloud_manipulation_pipeline(inPipeline, extent, extent_epsg, out_epsg, filename, grid_name, doReclassify = False,
                                doSavePointCloud = False):

    CRS = "EPSG:{}".format(out_epsg)

    if doReclassify:
        inPipeline['pipeline'].append(
            {
                "assignment": "Classification[:]=0",
                "tag": "wipeclasses",
                "type": "filters.assign"
            })

    if not(extent_epsg == out_epsg):
        inPipeline['pipeline'].append(
            {
                "out_srs": CRS,
                "tag": "reprojectUTM",
                "type": "filters.reprojection"
            })

    if doReclassify:
        inPipeline['pipeline'].append(
            {
                "tag": "groundify",
                "type": "filters.smrf"
            })

    if doSavePointCloud:
        inPipeline['pipeline'].append(
            {
                "filename": grid_name+".laz",
                "inputs": [inPipeline['pipeline'][-1]['tag']],
                "tag": "writerslas",
                "type": "writers.las"
            },
        )

    return inPipeline


def request_build_dem_creation_pipeline(extent,extent_epsg,out_epsg,grid_name, usgs_survey_name,cell_size = 1.0,
                           window_size = 6,nodatavalue = -9999, doReclassify = False,doSavePointCloud = False,
                                        pointResolution = None):
    '''

    :param extent: ([minx, maxx], [miny, maxy]).
    :param extent_epsg:
    :param grid_name:
    :param usgs_survey_name:
    :param cell_size:
    :param window_size:
    :param nodatavalue:
    :param doReclassify:
    :return:
    '''


    build_pipeline = request_pointcloud_pipeline(extent, extent_epsg, out_epsg, usgs_survey_name, grid_name,
                                                 doReclassify, doSavePointCloud,pointResolution)

    build_pipeline = append_dem_gridding(build_pipeline,grid_name,nodatavalue,cell_size, window_size)

    return build_pipeline

def append_dem_gridding(inPipeline,grid_name,nodatavalue, cell_size, window_size):

    #TODO: Perhaps I should split this further - adding a ground classification filter and then gridding
    inPipeline['pipeline'].append(
        {
            "limits": "Classification[2:2]",
            "type": "filters.range",
            "tag": "classifyground"
        })

    inPipeline['pipeline'].append(
        {
            "filename": grid_name+DEFAULTGRIDEXT,
            "gdalopts": "tiled=yes,     compress=deflate",
            "inputs": ["classifyground"],
            "nodata": nodatavalue,
            "output_type": "idw",
            "resolution": cell_size,
            "type": "writers.gdal",
            "window_size": window_size
        }
    )

    return inPipeline

def request_build_groundIntensity_creation_pipeline(extent,extent_epsg,out_epsg,grid_name, usgs_survey_name,cell_size = 1.0,
                           window_size = 6,nodatavalue = -9999, doReclassify = False,doSavePointCloud = False,
                                                    pointResolution = None):
    '''

    :param extent: ([minx, maxx], [miny, maxy]).
    :param extent_epsg:
    :param grid_name:
    :param usgs_survey_name:
    :param cell_size:
    :param window_size:
    :param nodatavalue:
    :param doReclassify:
    :return:
    '''

    build_pipeline = request_pointcloud_pipeline(extent, extent_epsg, out_epsg, usgs_survey_name, grid_name,
                                                 doReclassify, doSavePointCloud,pointResolution)

    build_pipeline = append_Intensity_gridding(build_pipeline,grid_name,nodatavalue,cell_size,window_size)

    return build_pipeline

def append_Intensity_gridding(inPipeline, grid_name, nodatavalue, cell_size, window_size):

    inPipeline['pipeline'].append(
        {
            "limits": "Classification[2:2]",
            "type": "filters.range",
            "tag": "classifyground"
        })

    inPipeline['pipeline'].append(
        {
            "filename": grid_name+DEFAULTGRIDEXT,
            "gdalopts": "tiled=yes,     compress=deflate",
            "inputs": ["classifyground"],
            "dimension": "Intensity",
            "nodata": nodatavalue,
            "output_type": "idw",
            "resolution": cell_size,
            "type": "writers.gdal",
            "window_size": window_size
        }
    )

    return inPipeline


########################################################################################################
###
########################################################################################################

def request_build_hag_creation_pipeline(extent,extent_epsg,out_epsg,grid_name, usgs_survey_name,cell_size = 1.0,
                           window_size = 6,nodatavalue = -9999, doReclassify = False, doSavePointCloud = False,
                                        pointResolution = None,rasterPath = None,doUseDelaunay = True):
    '''

    :param extent: ([minx, maxx], [miny, maxy]).
    :param extent_epsg:
    :param grid_name:
    :param usgs_survey_name:
    :param cell_size:
    :param window_size:
    :param nodatavalue:
    :param doReclassify:
    :param rasterPath:
    :param doUseDelaunay:
    :return:
    '''

    build_pipeline = request_pointcloud_pipeline(extent, extent_epsg, out_epsg, usgs_survey_name, grid_name,
                                                 doReclassify, doSavePointCloud,pointResolution)

    build_pipeline = append_hag_gridding(build_pipeline,grid_name,nodatavalue,cell_size,
                                         window_size,rasterPath,doUseDelaunay)

    return build_pipeline

def append_hag_gridding(inPipeline,grid_name,nodatavalue,cell_size,window_size,rasterPath=None,doUseDelaunay = True):
    '''

    :param inPipeline:
    :param grid_name:
    :param nodatavalue:
    :param cell_size:
    :param window_size:
    :param rasterPath:
    :param doUseDelaunay:
    :return:
    '''

    if not(rasterPath is None):
        inPipeline['pipeline'].append(
            {
                "type": "filters.hag_dem",
                "tag": "hag",
                "raster":rasterPath
            })

    elif doUseDelaunay:
        inPipeline['pipeline'].append(
            {
                "type": "filters.hag_delaunay",
                "tag": "hag"
            })
    else:
        inPipeline['pipeline'].append(
            {
                "type": "filters.hag_nn",
                "tag": "hag"
            })

    inPipeline['pipeline'].append(
        {
            "filename": grid_name+DEFAULTGRIDEXT,
            "gdalopts": "tiled=yes,     compress=deflate",
            "inputs": ["hag"],
            "dimension": "HeightAboveGround",
            "nodata": nodatavalue,
            "output_type": "idw",
            "resolution": cell_size,
            "type": "writers.gdal",
            "window_size": window_size
        }
    )

    return inPipeline

########################################################################################################
###
########################################################################################################
def get_tiled_request_pipelines(extent,extent_epsg,tileWidth,tileOverlap, baseName,build_pipeline_fun):
    '''

    :param extent: [minx, maxx], [miny, maxy]).
    :param extent_epsg:
    :param tileWidth:
    :param tileOverlap:
    :param baseName:
    :param build_pipeline_fun: function that takes input extent, extent_epsg, name, and returns a pipeline
    :return:
    '''

    if extent_epsg == ENTWINE_EPSG:
        ent_extent = extent
    else:
        ent_extent = reproject_extent(extent, extent_epsg, ENTWINE_EPSG)

    #Get the bounding coordinates for the tiles
    minX,maxX = ent_extent[0]
    minY,maxY = ent_extent[1]

    #Get the x boundaries and y boundaries of tiles closest to the tilewidth, while evenly dividing the full extent
    xBounds = np.linspace(minX,maxX,np.int((maxX-minX)/tileWidth)+1)
    yBounds = np.linspace(minY,maxY,np.int((maxY - minY)/tileWidth)+1)

    #Store a list of all the pipelines to run
    tiledPipelines = []

    #Loop through the tile grid
    for r in range(len(yBounds)-1):
        for c in range(len(xBounds)-1):
            thisExtent = ([xBounds[c], xBounds[c+1]+tileOverlap], [yBounds[r], yBounds[r+1]+tileOverlap])
            thisName = baseName+'_r{}_c{}'.format(r,c)
            thisPipeline = build_pipeline_fun(thisExtent,ENTWINE_EPSG,thisName)

            tiledPipelines.append(thisPipeline)

    return tiledPipelines


def get_orthogonal_swathprofile_request_pointcloud(lineXY,swath_relative_position,swath_length, swath_width,line_epsg,usgs_survey_name,save_name,doReclassify = False,
                                             doSavePointCloud = False, pointResolution = None):
    """[summary]

    Args:
        lineXY ([type]): [description]
        swathRelativePosition ([type]): [description]
        swathLength ([type]): [description]
        swathWidth ([type]): [description]
        line_epsg ([type]): [description]
        usgs_survey_name ([type]): [description]
        save_name ([type]): [description]
        doReclassify (bool, optional): [description]. Defaults to False.
        doSavePointCloud (bool, optional): [description]. Defaults to False.
        pointResolution ([type], optional): [description]. Defaults to None.
    """
    x,y = lineXY.xy
    x = np.array(x)
    y = np.array(y)
    
    #To determine position to place orthogonal line, calculate distances along line
    alongLineDists = np.zeros_like(x)
    alongLineDists[1:] = np.sqrt((x[1:] - x[:-1])**2 + (y[1:] - y[:-1])**2)
    absSwathDist = alongLineDists[-1]*swath_relative_position
    
    #Determine line segment that contains swath
    pastIdx = np.argwhere(alongLineDists >= absSwathDist)[0][0]
    
    if pastIdx == 0:
        print('Warning: based on the specified swath_relative_position, '+
              'the orthogonal profile will be constructed at the start of the line.')
        pastIdx = 1
    x0,x1 = x[pastIdx-1],x[pastIdx]
    y0,y1 = y[pastIdx-1],y[pastIdx] 
    
    
    #Determine the orientation of the line and something orthogonal to it
    dy,dx = y1 - y0, x1 - x0
    theta = np.arctan2(dy,dx)
    phi = theta+ np.pi/2.0
    
    #Determine the starting point of the orthogonal line
    lineLength = np.sqrt(dy**2 + dx**2)
    x_m = x0 + np.cos(theta)*lineLength*swath_relative_position
    y_m = y0 + np.sin(theta)*lineLength*swath_relative_position
    
    #Determine the start and end coordinates of the orthogonal line
    x_e = x_m + np.cos(phi)*swath_length/2.0
    y_e = y_m + np.sin(phi)*swath_length/2.0
    x_s = x_m - np.cos(phi)*swath_length/2.0
    y_s = y_m - np.sin(phi)*swath_length/2.0
    
    #Construct the new orthogonal line
    lineXY_orth = LineString([[x_s,y_s],[x_e,y_e]])
    
    #Use existing function to get the swath profile for the orthogonal line
    pipeline = get_swathprofile_request_pointcloud(lineXY_orth,swath_width,line_epsg,usgs_survey_name,save_name,doReclassify = doReclassify,
                                             doSavePointCloud = doSavePointCloud, pointResolution = pointResolution)
    
    return pipeline
    

def get_swathprofile_request_pointcloud(lineXY,swathWidth,line_epsg,usgs_survey_name,save_name,doReclassify = False,
                                             doSavePointCloud = False, pointResolution = None):

    #First buffer the line to create a polygon
    XYPolyBuffer = lineXY.buffer(swathWidth)
    x, y = lineXY.xy
    lineXY = [t for t in zip(x, y)]

    #Get the bounding coordinates of the buffered polygon
    bound_minX, bound_minY, bound_maxX,bound_maxY = XYPolyBuffer.bounds

    #Get the bounding extent of this data
    extent = ([bound_minX, bound_maxX], [bound_minY, bound_maxY])

    #Get the request pipeline
    loadPipeline = request_pointcloud_pipeline(extent,line_epsg,line_epsg,usgs_survey_name,save_name,
                                               doReclassify,doSavePointCloud=False,pointResolution=pointResolution)

    #Add a step to reproject to the coordinates of the line
    loadPipeline['pipeline'].append(
        {"type": "filters.reprojection",
         "out_srs":"EPSG:{}".format(line_epsg)
        }
    )

    #Add the projection onto the line
    loadPipeline['pipeline'].append(
        {"tag":'getswath',
         "type":'filters.python',
         "module":'pointCloudCreation',
         "script":'pointCloudCreation.py',
         "function":"projectPointsOntoLine_pipeline",
         "add_dimension":['L','D'],
         "pdalargs":{'line_coords':lineXY}

        }
    )

    # Add the projection onto the line
    loadPipeline['pipeline'].append(
        {"tag": 'trimswath',
         "type": 'filters.range',
         "limits": "D[0:{}]".format(swathWidth),
         }
    )
    #Should do saving down here
    if doSavePointCloud:
        loadPipeline['pipeline'].append(
            {
                "filename": save_name+".laz",
                "inputs": [loadPipeline['pipeline'][-1]['tag']],
                "tag": "writerslas",
                "type": "writers.las"
            },
        )

    pipeline = run_pipeline(loadPipeline,doReturnPipeline=True)

    return pipeline

def run_pipeline(pipeline_json, doReturnPipeline = False):
    '''

    :param pipeline_json:
    :return:
    '''

    pipeline = pdal.Pipeline(json.dumps(pipeline_json))

    #if pipeline.validate():
    pipeline.execute()

    #if not(doReturnPipeline):
        #pipeline = True

    return pipeline

########################################################################################################
###
########################################################################################################

def merge_warp_dems(inFileNames, outFileName,outExtent,outEPSG,pixelSize, doReturnGdalSourceResult = False,
                    resampleAlg = 'cubic', noDataValue = -9999):
    '''

    :param inFileNames: A list of all the filenames to merge
    :param outFileName: the output path to safe the file as
    :param outExtent:  ([minx, maxx], [miny, maxy]).
    :param outEPSG: EPSG code for the coordinate system of the specified output extent (also sets the output
    coordinate system)
    :param pixelSize: Dimension of the the output pixel (x and y direction) in the native units of the
    output coordinate system)
    :param doReturnGdalSourceResult: Boolean, if true returns the gdal source object for the newly created
    dataset. If False (the default) returns none and closes the connection to the newly created dataset.
    :param resampleAlg: The resampling algorithm to use in reprojecting and merging the raster. Can be
    any option allowed by GDAL. Prefered options will likely be: 'near', 'bilinear', 'cubic', 'cubicspline',
    'average'
    :param noDataValue: No data value to set for the output data.
    :return:
    '''

    wrpOptions = gdal.WarpOptions(
        outputBounds=[outExtent[0][0], outExtent[1][0], outExtent[0][1], outExtent[1][1]],
        outputBoundsSRS='EPSG:{}'.format(outEPSG),
        format='GTiff',
        xRes=pixelSize, yRes=pixelSize,
        resampleAlg=resampleAlg,
        dstSRS='EPSG:{}'.format(outEPSG),
        dstNodata=noDataValue,
        srcNodata=noDataValue
    )

    gridSource = gdal.Warp(outFileName,inFileNames, options=wrpOptions)

    if not(doReturnGdalSourceResult):
        gridSource = None

    return gridSource


#%%
import numpy as np
t = np.linspace(0,50,5)
np.argwhere(t > 30)[0][0]

#
# # testExtent = ([-107.685,-107.655],[37.42,37.46])
# # testExtent = ([-107.70,-107.62],[37.39,37.47])
# testExtent = ([-108.59, -108.57],[36.601, 36.644])
# testEPSG = 4326
# outEPSG = 32613
# grid_name = 'testData\\tileTest\\JMCLandslide'
# usgs_survey_name = 'NM_NorthWest_Navajo_A_TL_2018'
# cell_size = 2.0
# tileWidth = 1e3
# tileOverlap = 10*cell_size
#
# pointDensityRequest = 5 #Point per sq m
# pointResolution = 1.0/np.sqrt(pointDensityRequest)
#
# #Request tiled dems
# pipeline_builder = lambda ext,ext_epsg,name: request_build_dem_creation_pipeline(ext,ext_epsg,outEPSG,name,
#                                                                                  usgs_survey_name,cell_size)
#
# tiled_pipelines = get_tiled_request_pipelines(testExtent,testEPSG,tileWidth,tileOverlap,grid_name,pipeline_builder)
#
# #loop through and run the pipelines
# for i,pipeline in enumerate(tiled_pipelines):
#     thisGrid = pipeline['pipeline'][-1]['filename']
#     print('Working on {} of {}, name: {}'.format(i+1,len(tiled_pipelines),thisGrid))
#     run_pipeline(pipeline)
#
# #Merge those DEMS
# outEPSG_final = 32613
# outExtent = reproject_extent(testExtent,testEPSG,outEPSG_final)
# demNamesToMerge = [p['pipeline'][-1]['filename'] for p in tiled_pipelines]
#
# merge_warp_dems(demNamesToMerge,grid_name+'.tif',outExtent,outEPSG_final,cell_size)
#
#
#
# #Request points
# pipeline = request_pointcloud_pipeline(testExtent,testEPSG,outEPSG,usgs_survey_name,grid_name,
#                                                doReclassify=False,doSavePointCloud=True, pointResolution=pointResolution)
# run_pipeline(pipeline)
#
# newExtent = reproject_extent(testExtent,testEPSG,outEPSG)
# pipeline = load_pointcloud_pipeline(newExtent,outEPSG,outEPSG,grid_name+'.laz',grid_name)
# pipeline = append_dem_gridding(pipeline,grid_name+'.tif',-9999,cell_size,window_size=6)
#
# run_pipeline(pipeline)
#
# dem= dg.demGrid([],rasterPath=grid_name+'.tif')
# dem[dem <0] = np.nan
# hs = dem.calcHillshade()
# Smag = dem.calcSlopeMagnitude()
# hs.overwriteSaveFile()
# Smag.overwriteSaveFile()
# #
# f,axs = plt.subplots(1,1,figsize = (16,6),dpi = 200,sharex=True,sharey=True)
#
# hs.plotGrid(axs = axs,cmap = 'gray',vmin = 0, vmax = 255)
# # Smag.plotGrid(axs = axs,cmap = 'Reds',alpha = 0.5, vmin = 0, vmax = 0.75)
# %%
