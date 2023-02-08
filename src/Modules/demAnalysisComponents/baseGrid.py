import os
from time import sleep

import numpy as np
from matplotlib import pyplot as plt

import scipy.ndimage as ndi

from osgeo import gdal

from shapely.geometry import Point, Polygon, LineString

class baseGrid:

    # Kernel for eight neighbors
    __ROW_KERNEL = np.array([1, 1, 1, 0, 0, -1, -1, -1])
    __COL_KERNEL = np.array([-1, 0, 1, -1, 1, -1, 0, 1])

    def __init__(self, numpyGrid: np.ndarray, dx: float= 1.0, xllcenter : float  = 0, yllcenter : float = 0,
                 geotransform : list = None, projection = None, rasterPath = None):
        '''

        :param numpyGrid: a 2d numpy array of values
        :param dx: (default = 1) The grid spacing of that grid
        :param xllcenter: (default = 0) The x coordinate of the lower left pixel
        :param yllcenter: (default = 0) The y coordinate of the lower left pixel
        :param geotransform: (default = None). The list of 6 items used by gdal that specify the georeferencing
        and rotation of the grid
        :param projection: (default = None). The gdal projection information for this grid
        :param rasterPath: (default = None). A path to a raster that can be loaded with gdal
        '''

        #If there wasn't a georefed raster specified, load attributes based on what was provided
        if rasterPath is None:

            self._nrows, self._ncols = numpyGrid.shape
            self.shape = (self._nrows, self._ncols)
            self.grid = np.copy(numpyGrid)
            self._filePath = None

            if not (geotransform is None):
                try:
                    # Check that the geotransform has the write number of inputs
                    if len(geotransform) == 6:

                        self._dx = geotransform[1]
                        self._dy = np.abs(geotransform[5])

                        self._xllcenter = geotransform[0] + self._dx / 2.0
                        self._yllcenter = geotransform[3] - (self._dy * (self._nrows - 0.5))

                        self._geotransform = geotransform
                    else:
                        raise Exception('Incorrect format for input, check your geotransform')
                except:
                    print('Whoops, geotransform needs to have 6 items. This is what was provided: {}'.format(
                        geotransform))
                    print('WARNING: Ignoring geotransform input and proceeding with defaults.')
                    geotransform = None

            # If the geotransform specified was faulty, or if one wasn't specified
            if geotransform is None:
                self._dx = dx
                self._dy = dx

                self._xllcenter = xllcenter
                self._yllcenter = yllcenter

                xUL = self._xllcenter - self._dx / 2.0
                yUL = self._yllcenter + (self._dy * (self._nrows - 0.5))

                geotransform = (xUL, self._dx, 0, yUL, 0, -self._dy)
                self._geotransform = geotransform

            self._projection = projection

        #If there was a raster file path provided, load based on that
        else:
            if os.path.isfile(rasterPath):
                self._loadFromRaster(rasterPath)
            else:
                raise Exception('Error, it looks like you did not specifiy a valid raster path')

    def _loadFromRaster(self, filePath: str):
        '''
        Loads in a raster and populates class variables with georeferencing, grid size information)
        :param filePath:  The path to the raster data that is to be loaded
        :return:

         TODO: reset pixel size if this is a geographic grid - think about other ways to handle that...
         '''

        # Read the grid in using GDAL
        src = gdal.Open(filePath)
        self.grid = src.ReadAsArray().astype(np.float)

        # Get the size information from the dataset
        self._ncols = src.RasterXSize
        self._nrows = src.RasterYSize
        self.shape = (self._nrows, self._ncols)

        self._geotransform = src.GetGeoTransform()  # Steal the coordinate system information from the old dataset
        self._projection = src.GetProjection()  # Steal the projection from the old dataset

        self._dx = self._geotransform[1]
        self._dy = np.abs(self._geotransform[-1]) #TODO: Hmmm... can grids be flipped in different orientations? This value of _geotransform is not always negative...
        self._xllcenter = self._geotransform[0] + self._dx / 2.0

        #TODO: I've encountered a raster that isn't oriented normally.... geotransform[-1] is positive and the raster appears flipped
        self._yllcenter = self._geotransform[3] - (self._dy * (self._nrows - 0.5))

        # Get the information from the existing file needed to resave as the same type
        self._filePath = filePath
        self._gdalDriver = src.GetDriver().ShortName
        self._gdalDataType = src.GetRasterBand(1).DataType
        
        # Close the original source
        src = None
  
    
    def overwriteSaveFile(self):
        '''
        Overwrites the save file from which this grid was loaded (or previously saved) with the current data in the grid.
        :return:
        '''
        if self._filePath is None:
            print('Whoops, this grid does not have an associated file, so there is nothing to overwrite!')
            print('Try calling .saveGrid instead. ')
        else:
            self.saveGrid(self._filePath, self._gdalDriver, self._gdalDataType)

    def _copyAttributesFromGrid(self,gridToCopy):
        '''
        Copy the attributes (those that aren't the actual data grid itself) from another grid.

        TODO: Is there an alternative way to implement this, so that I do not have to list individual assignment operations?
        :param gridToCopy: A baseGrid or derived class that has the spatial referencing attributes we want to assign
        to this grid
        :return:
        '''

        self._nrows, self._ncols = gridToCopy._nrows, gridToCopy._ncols
        self.shape = (self._nrows, self._ncols)
        self._dx = gridToCopy._dx
        self._dy = gridToCopy._dy

        # Create spatial coordinates for this grid
        self._xllcenter = gridToCopy._xllcenter
        self._yllcenter = gridToCopy._yllcenter

        # Get the assigned geotransform from the grid
        self._geotransform = gridToCopy._geotransform
        self._projection = gridToCopy._projection

    def __setitem__(self, key, value):
        '''Overwrite indexing of numpy array'''
        if isinstance(key,baseGrid):
            if isinstance(value,baseGrid):
                self.grid[key.grid] = value.grid
            else:
                self.grid[key.grid] = value
        else:
            self.grid[key] = value

    def __getitem__(self, key):
        '''Overwrite indexing of numpy array'''
        if isinstance(key,baseGrid):
            return self.grid[key.grid]
        else:
            return self.grid[key]

    def __lt__(self, other):
        '''
        Overwrites '<' operator. (ThisGrid < other)

        :param other:
        :return:
        '''

        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid < other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid < other)
        else:
            raise Exception('Can only compare to other demGrids or numpy arrays')

    def __le__(self, other):
        '''
        Overwrites '<=' operator. (ThisGrid < b)

        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid <= other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid <= other)
        else:
            raise Exception('Can only compare to other demGrids or numpy arrays')

    def __eq__(self, other):
        '''
        Overwrites '==' operator. (ThisGrid == b.grid)

        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid == other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid == other)
        else:
            raise Exception('Can only compare to other demGrids or numpy arrays')

    def __ne__(self, other):
        '''
        Overwrites '!=' operator. (ThisGrid != b.grid)

        :param b:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid != other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid !=other)
        else:
            raise Exception('Can only compare to other demGrids or numpy arrays')

    def __gt__(self, other):
        '''
        Overwrites '==' operator. (ThisGrid == other.grid)

        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid > other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid > other)
        else:
            raise Exception('Can only compare to other demGrids or numpy arrays')

    def __ge__(self, other):
        '''
        Overwrites '>=' operator. (ThisGrid >= other.grid)

        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid >= other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid >= other)
        else:
            raise Exception('Can only compare to other demGrids or numpy arrays')

    def __and__(self, other):

        '''
        Overwrites '&' operator. (ThisGrid & other.grid)
        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid & other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid & other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __ior__(self,other):
        '''
        TODO: Is this and __and__ overwriting the operators I think they are?
        Overwrites '|' operator. (ThisGrid | other.grid)
        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid | other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid | other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __iand__(self, other):

        '''
        Overwrites '&' operator. (ThisGrid & other.grid)
        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid & other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid & other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __or__(self,other):
        '''
        TODO: Is this and __and__ overwriting the operators I think they are?
        Overwrites '|' operator. (ThisGrid | other.grid)
        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid | other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid | other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __and__(self, other):

        '''
        Overwrites '&' operator. (ThisGrid & other.grid)
        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid & other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid & other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __add__(self, other):

        '''
        Overwrites '+' operator. (ThisGrid + b.grid)
        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid + other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid + other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __sub__(self, other):

        '''
        Overwrites '-' operator. (ThisGrid - other.grid)
        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid - other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid - other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __mul__(self, other):
        '''
         Overwrites '*' operator. (ThisGrid * other.grid)
         :param other:
         :return:
         '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid * other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid * other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __truediv__(self, other):
        '''
         Overwrites '/' operator. (ThisGrid / other.grid)
         :param other:
         :return:
         '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid / other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid / other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __mod__(self, other):
        '''
         Overwrites '%' operator. (ThisGrid % other.grid)
         :param other:
         :return:
         '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid % other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid % other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __pow__(self, other):
        '''
         Overwrites '**' operator. (ThisGrid ** other.grid)
         :param other:
         :return:
         '''
        if isinstance(other, baseGrid):
            return self.duplicateGridWithNewArray(self.grid ** other.grid)
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            return self.duplicateGridWithNewArray(self.grid ** other)
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

    def __iadd__(self, other):

        '''
        Overwrites '+' operator. (ThisGrid + b.grid)
        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            self.grid += other.grid

        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            self.grid += other

        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

        return self

    def __isub__(self, other):

        '''
        Overwrites '-' operator. (ThisGrid - other.grid)
        :param other:
        :return:
        '''
        if isinstance(other, baseGrid):
            self.grid -= other.grid

        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            self.grid -= other

        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

        return self

    def __imul__(self, other):
        '''
         Overwrites '*' operator. (ThisGrid * other.grid)
         :param other:
         :return:
         '''
        if isinstance(other, baseGrid):
            self.grid *= other.grid
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            self.grid *= other
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

        return self

    def __idiv__(self, other):
        '''
         Overwrites '/' operator. (ThisGrid / other.grid)
         :param other:
         :return:
         '''
        if isinstance(other, baseGrid):
            self.grid /= other.grid
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            self.grid /= other
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

        return self

    def __imod__(self, other):
        '''
         Overwrites '%' operator. (ThisGrid % other.grid)
         :param other:
         :return:
         '''
        if isinstance(other, baseGrid):
            self.grid %= other.grid
        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
           self.grid %= other
        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

        return self

    def __ipow__(self, other):
        '''
         Overwrites '**' operator. (ThisGrid ** other.grid)
         :param other:
         :return:
         '''
        if isinstance(other, baseGrid):
            self.grid **= other.grid

        elif isinstance(other, np.ndarray) or isinstance(other, int) or isinstance(other, float):
            self.grid **= other

        else:
            raise Exception('Can only compute with demGrids or numpy arrays')

        return self

    def __abs__(self):
        '''
        Overwrite absolute value function
        :return:
        '''
        return self.duplicateGridWithNewArray(np.abs(self.grid))

    def _isRowsColsInBounds(self,rows : np.ndarray,cols : np.ndarray):
        '''
        Return a boolean numpy array for the row, column pairs that are in the bounds of the grid
        :param rows: numpy array of rows in the grid to index
        :param cols: numpy array of columns in the grid to index
        :return: inBounds, boolean array of which row, column pairs are in the bounds of this grid
        '''
        # Trim off anything outside the bounds
        return (cols >= 0) & (cols < self._ncols) & (rows >= 0) & (rows < self._nrows)

    def _getWrappedRowsCols(self,rows : np.ndarray, cols : np.ndarray, equivalentEdges = False):
        '''

        In some cases it may be desirable to be able to index kernels to wrap from one end of the grid to the other.
        That is, if we want something that is two columns to the right of the right hand side of the grid, we will get
        the second item from the left hand side. This function provides that indexing.

        :param rows: row indices of the grid to wrap
        :param cols: column indices of the grid to wrap
        :param equivalentEdges: Boolean. If true, considers that the (for example) left and right sides of the grid
        are equivalent, so that an index that is 1 column to the right hand side will be the second column from the
        left (i.e., column index 1). If false, an index 1 column to the right hand side would be column index 0.
        :return:
        '''

        # Hmm... to do another minus 1 (as if edges are equivalent)? Maybe make an option
        if equivalentEdges:
            incrementor = 1
        else:
            incrementor = 0

        #Wrap negatives back to positive values on opposite side of grid
        rows[rows < 0] = self._nrows + rows[rows < 0] - incrementor
        cols[cols < 0] = self._ncols + cols[cols < 0] - incrementor

        #Wrap overshoots back to begenning of array
        rows[rows >= self._nrows] = rows[rows >= self._nrows] - self._nrows + incrementor
        cols[cols >= self._ncols] = cols[cols >= self._ncols] - self._ncols + incrementor

        return rows, cols

    def median(self):
        '''
        Get the median value of the entire grid, accounting for any Nan values used to mask out areas
        :return: np.nanmedian(self.grid) - the mean value of the entire grid
        '''

        return np.nanmedian(self.grid)

    def mean(self):
        '''
        Get the mean value of the entire grid, accounting for any Nan values used to mask out areas
        :return: np.nanmean(self.grid) - the mean value of the entire grid
        '''

        return np.nanmean(self.grid)

    def std(self):
        '''
        Get the standard deviation value of the entire grid, accounting for and NaN values used to mask out areas
        :return: np.nanstd(self.grid) - the mean value of the entire grid
        '''

        return np.nanstd(self.grid)

    def percentile(self, percents: list):
        '''
        Get the requested percentile values of the grid, accounting for and NaN values used to mask out areas
        :param percents: iterable (e.g., list) of percent values you want from the grid
        :return: np.percentile(self.grid,percents) - list of percentile values from the distribution of grid values
        '''

        return np.nanpercentile(self.grid,percents)

    def min(self):
        '''
        Get the minimum value of the entire grid, accounting for and NaN values used to mask out areas
        :return: np.nanmin(self.grid) - the minimum value of the entire grid
        '''

        return np.nanmin(self.grid)

    def max(self):
        '''
        Get the maximum value of the entire grid, accounting for and NaN values used to mask out areas
        :return: np.nanmax(self.grid) - the maximum value of the entire grid
        '''

        return np.nanmax(self.grid)
    
    def ndv(self, filepath:str): # added by MS 12/7/22
        src = gdal.Open(filepath)
        src.ReadAsArray().astype(float)
        ndv = src.GetRasterBand(1).GetNoDataValue()
        #print("NoData Value: ", ndv)
        #Close file
        src = None
        return ndv
    
    def relief(self):
        '''
        Get the difference between the maximum and minimum elevations within the grid.

        :return:  np.nanmax(self.grid) - np.nanmin(self.grid), the total relief within this grid
        '''

        return np.nanmax(self.grid) - np.nanmin(self.grid)

    def sumSquaredResiduals(self, comparisonGrid):
        '''
        Computes the sum of the squared errors between this grid and a comparison grid
        :param comparisonGrid: Another baseGrid or derived grid instance or numpy array to compare to this grid. It is
        critical that these grids are the same shape and cover the same extent. No effort will be made to co-locate
        mis-aligned grids.
        :return: L2Norm (float), the sum of the squared difference between this grid and comparisonGrid.
        '''

        L2Norm = None

        if isinstance(comparisonGrid, baseGrid):
            L2Norm = np.sum((self.grid - comparisonGrid.grid)**2)

        elif isinstance(comparisonGrid, np.ndarray):
            L2Norm = np.sum((self.grid - comparisonGrid)**2)

        return L2Norm

    def getGridCoords(self):
        '''
        Get the x and y coordinates of this grid as arrays corresponding to the columns and rows of the grid cells.
        :return: xArray, yArray: 1d numpy arrays of the x and y coordinates of this grid
        '''
        # Get grid size
        nx, ny = self._ncols, self._nrows

        # Create arrays of the x and y coordinates of each pixel (the axes)
        xcoordinates = np.array([x * self._dx + self._xllcenter for x in range(nx)])
        # Flip the ys so that the first row corresponds to the first entry of this array
        ycoordinates = np.array([y * self._dy + self._yllcenter for y in range(ny)][::-1])
        return xcoordinates, ycoordinates

    def getGridCoordinatesMeshgrid(self):
        '''
        Get grids of the X and Y coordinates for this grid
        :return: X,Y - np.meshgrid(xCoordinates, yCoordinates); grids of the size of this one with all corresponding
        X and Y coordinates
        '''

        x, y = self.getGridCoords()
        X, Y = np.meshgrid(x, y)

        return X, Y

    def getRowColFromXY(self, xy):
        '''
        TODO: Update this for iterables

        Get the nearest row,column index to the specified xy point.

        :param xy: a single x,y pair that lies within this grid.
        :return: row, col . Integers of the row, column index of the nearest grid cell to these points
        '''

        row = self._nrows - np.int(np.round((xy[1] - self._yllcenter) / self._dy)) - 1
        col = np.int(np.round((xy[0] - self._xllcenter) / self._dx))  # Because GDAL counts up from the bottom, note this requires zero indexing

        return row, col

    def getXYFromRowCol(self,rowCol):
        '''
        TODO: Update this for iterables

        :param rowCol:
        :return:
        '''
        y = (self.shape[0] - (rowCol[0] + 1)) * self._dy + self._yllcenter
        x = self._xllcenter + (self._dx * rowCol[1])

        return x,y

    def getGridValueAtPoint(self,xy):
        '''
        TODO: Update this for iterables

        Get the value of the grid at the nearest row,column index to the specified xy point.
        :param xy:
        :return: value, grid(x,y). The value at the nearest row,column index to this point.
        '''

        rowCol = self.getRowColFromXY(xy)

        return self.grid[rowCol]

    def getGridExtent(self):
        '''Return the bounding extent of the grid

        :return (minX, maxX, minY, maxY) - the bounding coordinates of this raster
        '''

        #TODO: I encountered a raster that appeared flipped... how do I deal with this

        return (self._xllcenter, self._xllcenter+(self._dx*(self._ncols-1)), self._yllcenter, self._yllcenter+(self._dy*(self._nrows-1)))

    def getNeighborIndices(self, row, col):
        '''
        Get the indices of the neighboring cells according to this flow grid search kernel.
        :param row: the row index of the current cell
        :param col: the column index of the current cell
        :return: neighborRows,neighborCols,neighbordDists : numpy ndarrays of the row, column, and distance
        '''
        # Find all the surrounding indices
        outRows = self.__ROW_KERNEL + row
        outCols = self.__COL_KERNEL + col
        dxs = np.sqrt(((outRows - row)*self._dx)**2 + ((outCols - col)*self._dy)**2)

        # Determine which indices are out of bounds
        inBounds = self._isRowsColsInBounds(outRows,outCols)

        return outRows[inBounds], outCols[inBounds],dxs[inBounds]

    def createMaskFromGeoDataFrame(self,gdf,name:str = ''):
        '''

        :param gdf: geopandas geodataframe of polygons
        :param name: optional name to add to the filepath of the mask
        :return: maskedGrid: an instance of demGrid with just this value masked
        '''


        #TODO: Check for whether this is a polygon, could alternative not use geodataframes but just the list of polygon
        g = [geo for geo in gdf.geometry]

        #Preallocate space for if this is a masked grid
        isInBoundsGrid = np.zeros_like(self.grid)

        # For each polygon in this dataframe
        for j in range(len(g)):
            try:
                thisXYbounds = g[j].exterior.coords.xy
            except:
                # Exception catches multi-geometry
                print('Encountered ' + g[j].geom_type + ' at ' + str(j) )
                thisXYbounds = g[j][0].exterior.coords.xy

            #Covert the polygon to indices w/i that polygon
            rIdx,cIdx = self.getRowColIndicesBoundByPolygon(thisXYbounds)

            isInBoundsGrid[rIdx, cIdx] = 1

        return self.duplicateGridWithNewArray(isInBoundsGrid==1,'_'+name+'Mask')

    def getRowColIndicesBoundByPolygon(self,XYPolygon):

        # Revisit original roughness image, look at roughness w/i a band
        from skimage.draw import polygon

        # For each coordinate in the boundary of the polygon
        rs, cs = np.zeros(len(XYPolygon[0])), np.zeros(
            len(XYPolygon[0]))  # Preallocate space for row,colum of boundary
        # Calculate this coordinates closest pixel position
        for k in range(len(rs)):
            rs[k], cs[k] = self.getRowColFromXY((XYPolygon[0][k], XYPolygon[1][k]))

        # convert to polygon
        rIdx, cIdx = polygon(rs, cs)

        #Trim out any bad indices
        inBounds = (rIdx > 0) & (rIdx < self._nrows) & (cIdx > 0) & (cIdx < self._ncols)

        if np.any(~inBounds):
            rIdx = rIdx[inBounds]
            cIdx = cIdx[inBounds]

        return rIdx,cIdx

    def getProfile(self,lineXY,swathWidth:float = None, doTrimBeyondSwath = True):


        if swathWidth is None:
            swathWidth = np.sqrt(self._dx**2 + self._dy**2)

        #First create mask around profile to limit which data points to check
        #Buffer line by swath width
        XYPolyBuffer = LineString(lineXY).buffer(swathWidth).exterior.coords.xy

        #Get mask of buffer
        rIdx,cIdx = self.getRowColIndicesBoundByPolygon(XYPolyBuffer)

        xCoords,yCoords = self.getGridCoords()

        #Extract points w/i buffer
        xyPoints = [np.array([xCoords[cIdx[i]],yCoords[rIdx[i]]]) for i in range(len(rIdx))]
        vals = self.grid[rIdx,cIdx]

        projectedDistances = np.zeros_like(vals) * np.nan
        alongProfileDistances = np.zeros_like(vals)

        #Preallocate some space to check for how far and if things are in bounds
        isProjected = np.zeros_like(vals) == 1
        isInBounds = np.zeros_like(vals) == 1

        L = 0
        #For each line segment
        for i in range(len(lineXY)-1):

            #Get the squared distance along this line
            l2 = np.sum((lineXY[i] - lineXY[i+1]) ** 2)

            #For each point
            for j,p in enumerate(xyPoints):

                #If point is not yet projected
                if not isProjected[j]:

                    #Project points onto line

                    # Determine fraction along line segment that this point lies
                    t = np.sum((p - lineXY[i]) * (lineXY[i+1] - lineXY[i])) / l2

                    #If point projection falls onto segment
                    if (t >= 0) and (t <= 1):
                        # Mark point as projected
                        isProjected[j] = True

                        #Project the point
                        projPoint = lineXY[i] + t * (lineXY[i+1] - lineXY[i])

                        #Get distance to projected point from original point
                        projectedDistances[j] = np.sqrt(np.sum((projPoint - p)**2))

                        # Store projected point and distance to start of line
                        alongProfileDistances[j] = L + np.sqrt(np.sum((projPoint - lineXY[i]) ** 2))

                        #If distance is within swath width
                        isInBounds[j] = projectedDistances[j] <= swathWidth

            # Calculate along line distance from the start of this line segment to start of line
            L+= np.sqrt(np.sum((lineXY[i] - lineXY[i+1])**2))

        # If requested, trim the profile to what is in bounds
        if doTrimBeyondSwath:
            alongProfileDistances = alongProfileDistances[isInBounds]
            vals = vals[isInBounds]
            projectedDistances = projectedDistances[isInBounds]

        return alongProfileDistances,vals,projectedDistances

    def returnMaskedGrid(self,goodValueMask,returnAsArray: bool = False, newSuffix: str = '_masked'):
        '''

        :param goodValueMask:
        :param returnAsArray:
        :param newSuffix:
        :return:
        '''

        #Copy the grid
        gridCopy = np.copy(self.grid)

        #Set non-masked values to nan
        if isinstance(goodValueMask,baseGrid):
            gridCopy[~goodValueMask.grid] = np.nan
        else:
            gridCopy[goodValueMask] = np.nan

        #Default to returning this as a new grid instance unless requested otherwise
        if not(returnAsArray):
            gridCopy = self.duplicateGridWithNewArray(gridCopy,newSuffix = newSuffix)

        return gridCopy



    def calcMedianValuesBinnedByAnotherGrid(self,grid, bins = 20, percRange = [2.5, 97.5]):
        '''

        :param grid:
        :param bins:
        :param percRange:
        :return:
        '''

        #If a number of bins was provided (rather than bin edges) we need to generate the bin edges
        if isinstance(bins,int):
            #Get bin edges in a way that tries to space data reasonably well, even if it spans multiple orders of mag.
            bins = self._getWellSpacedBinEdgesForGrid(grid,bins)

        xMidpoints = (bins[:-1] + bins[1:])/2
        medianVals = np.zeros_like(xMidpoints)*np.nan
        percentileVals = np.zeros((2, len(xMidpoints)))

        for i in range(len(xMidpoints)):
            theseIndices = (grid >= bins[i]) & (grid < bins[i + 1])
            data_i = self[theseIndices]
            if len(data_i) > 3:
                medianVals[i] = np.nanmedian(data_i)
                percentileVals[:, i] = np.nanpercentile(data_i, percRange)

        return xMidpoints, medianVals,percentileVals

    def _getWellSpacedBinEdgesForGrid(self, grid, nBins):
        '''

        :param grid:
        :return:
        '''
        # If there are any negative or zero values, we can't necessarily identify a log-scaled range
        if np.any(grid <= 0):
            # If there is a three order of mag variation between extreme values, make a 'sym log' type scale
            # This is a log scale on either side of 0
            if np.abs(grid.max() / grid.min()) > 1e3:
                negVals = grid < 0
                posVals = grid > 0

                # TODO: Test all this - yikes! Could probably make a simpler suite of cases...
                # How do we split up the bins between positive and negative values? Perhaps in a way that tries to
                # Divide thins up evenly?
                if len(grid[negVals])>0: #Can't use np.any, because I don't know if this grid is a numpy array or a baseGrid
                    # Determine how many negative bins we should have
                    negOOM = np.log10(-np.min(grid[negVals]))  # How many o.o.m. below zero do negatives extend?

                    # If there are any positive values
                    if np.any(posVals):
                        posOOM = np.log10(np.max(grid[posVals]))  # How many o.o.m. above zero do positives extend?
                    # If not, there are 0 positive value bins
                    else:
                        posOOM = 0

                    negBins = int(nBins*(negOOM / (negOOM + posOOM)))

                    # Find log spacing of values less than zero
                    negLog = np.logspace(np.log10(np.min(-grid[negVals])), np.log10(np.max(-grid[negVals])),
                                         negBins)
                else:
                    negBins = 0
                    negLog = np.array([])

                if len(grid[posVals]) > 0:
                    # Find log spacing of values greater than zero
                    posLog = np.logspace(np.log10(np.min(grid[posVals])), np.log10(np.max(grid[posVals])),
                                         nBins - negBins)
                else:
                    posLog = []

                # Stick those two segments together
                bins = np.hstack((-1.0*negLog, posLog))
            else:
                bins = np.linspace(grid.min(), grid.max(), nBins)

        # If all values positive, and there is a 3 order of magnitude variability use a log scale
        elif (grid.max() / grid.min()) > 1e3:
            bins = np.logspace(np.log10(grid.min()), np.log10(grid.max()), nBins + 1)

        # Default to a linear spacing:
        else:
            bins = np.linspace(grid.min(), grid.max(), nBins)

        return bins

    def plotGrid(self,axs = None, logTransform = False, **kwargs):
        '''This is a wrapper for matplotlib plt.imshow that plots this grid

         :return axs (the axis that this was plotted on)'''

        if axs is None:
            f,axs = plt.subplots(1,1)
           
        if logTransform:
            axs.imshow(np.log10(self.grid), extent=self.getGridExtent(), **kwargs)
        else:
            axs.imshow(self.grid, extent=self.getGridExtent(), **kwargs)
        
        
        
        return axs

    def plotThisGridAgainstAnotherGrid(self, xAxisGrid, plotFractionOfObservations: float = 1, axs = None,
                                       xAxisLabel = 'Label me', yAxisLabel ='Label me', **kwargs):
        '''

        Plots the values of this grid on the y axis against the values of another grid on the x axis.
        Optionally sub-samples grid values from binned intervals based on plotFractionOfObservations. When this
        is between 0 and 1, we will only grab a random fraction of datapoints to plot - simplifying visualization. This
        random data points will be spread out within 20 bins on the x-axis, so that low data-density proportions of the
        plot still get some data shown.

        :param grid: the grid whose values should go on the x-axis. Can be a numpy ndarray or an instance of baseGrid
        :param fractionToSample: What proportion of data points should be plotted?
        :param axs: the matplotlib.pyplot axis to plot on. If None is specified, one will be created for this plot
        :return: axs - the matplotlib.pyplot axis of the plot
        '''

        # We'll divide everything into a number of bins, and do our sub-sampling from that so show data more uniformly
        samplingBins = 20

        if axs is None:
            f, axs = plt.subplots(1, 1)

        if (plotFractionOfObservations > 0) & (plotFractionOfObservations < 1):
            # Figure out about how many observations would occur in each bin
            nObservationsPerBin = (self._nrows * self._ncols * plotFractionOfObservations) / samplingBins

            x_bins = self._getWellSpacedBinEdgesForGrid(xAxisGrid,samplingBins)

            # Store data in lists
            allXs = []
            allYs = []

            # Loop through each of the bins
            for i in range(len(x_bins) - 1):
                theseIndices = (xAxisGrid >= x_bins[i]) & (xAxisGrid < x_bins[i + 1])

                # Pull out these indices
                x_i = xAxisGrid[theseIndices].flatten()
                y_i = self[theseIndices].flatten()

                # If there are more than the desired operations per bin, subsample
                if len(x_i) > nObservationsPerBin:
                    # Get a random subset of these indices
                    plotIndices = np.arange(len(x_i))  # Create an index for all observations
                    np.random.shuffle(plotIndices)  # Shuffle the index in place
                    plotIndices = plotIndices[:int(nObservationsPerBin)]
                    x_i = x_i[plotIndices]
                    y_i = y_i[plotIndices]

                allXs.append(x_i)
                allYs.append(y_i)

            XtoPlot = np.hstack(allXs)
            YtoPlot = np.hstack(allYs)
            axs.plot(XtoPlot, YtoPlot, '.', **kwargs)

        elif plotFractionOfObservations == 1:
            #Plot these grids, checking if we are dealing with a numpy grid or a baseGrid
            if isinstance(xAxisGrid,baseGrid):
                axs.plot(xAxisGrid.grid.flatten(), self.grid.flatten(), '.', **kwargs)
            else:
                axs.plot(xAxisGrid.flatten(), self.grid.flatten(),'.',**kwargs)
        else:
            print('Whoops, plotFractionOfObservations must be between 0 and 1, not plotting anything...')


        #If this axis hasn't been labelled yet, or something other than the default was specified
        if (axs.xaxis.label.get_text() == '') or not(xAxisLabel == 'Label me'):
            axs.set_xlabel(xAxisLabel)
        if (axs.yaxis.label.get_text() == '') or not(yAxisLabel == 'Label me'):
            axs.set_ylabel(yAxisLabel)

        return axs

    def plotMedianValuesBinnedByOtherGrid(self, grid, bins = 20, percRange = [2.5, 97.5],axs=None,
                                   yAxisLabel: str = 'Label me',xAxisLabel: str = 'Label me', **kwargs):
        '''
        Plots the median values of this grid (y-axis) calculated within bins of a specified grid (x-axis), along with
        'errorbars' to highlight the spread of the distribution of values in this grid.

        :param grid: The grid to use as the x-axis values that will be binned.
        :param bins: Either a number of bins (integer) or the values of bin edges. If a number of bins is specified,
        this will assume that these bins should be evenly spaced in log space.
        :param percRange: The upper and lower percentile bounds of y-axis data to show as error bars
        :param axs: The matplotlib.pyplot axis to do any plotting on. If None, a new axis will be created.
        :param kwargs: Any additional optional arguments to pass to matplotlib.pyplot.errorbar
        :return:
        '''

        #Get the binned mean values
        xMids, medVals, percVals = self.calcMedianValuesBinnedByAnotherGrid(grid,bins = bins, percRange = percRange)

        #Transform the percentile range to be relative 'errors' to the median value
        percErrs = np.zeros_like(percVals)
        for i in range(len(medVals)):
            percErrs[:,i] = np.abs(medVals[i] - percVals[:,i])

        if axs is None:
            f,axs = plt.subplots(1,1)

        axs.errorbar(xMids,medVals,yerr = percErrs,fmt = 'o',**kwargs)

        #If this axis hasn't been labelled yet, or something other than the default was specified
        if (axs.xaxis.label.get_text() == '') or not(xAxisLabel == 'Label me'):
            axs.set_xlabel(xAxisLabel)
        if (axs.yaxis.label.get_text() == '') or not(yAxisLabel == 'Label me'):
            axs.set_ylabel(yAxisLabel)

        return axs

    def plotHistogram(self,axs = None, **kwargs):
        '''
        Plots a histogram of this grids values, passes all non np.nan grid values to matplotlib.pyplot.hist

        :param axs: the matplotlib.pyplot axis to plot on
        :param kwargs: Any optional keyword arguments to pass to pyplot.histogram
        :return: axs the axis that was plotted on
        '''

        if axs is None:
            f,axs = plt.subplots(1,1)

        goodData = self.grid[~np.isnan(self.grid)].flatten()

        axs.hist(goodData,**kwargs)

        return axs

    def updateFileNameWithSuffix(self, suffix: str):
        '''

        Adds a suffix string ('filePath' + suffix + 'fileExtension') to the existing file path for this grid.

        :param suffix: the suffix to be added to this file name
        :return:
        '''


        #If there is an existing filepath
        if not(self._filePath is None):
            #Split the file path at its extension
            preExt,Ext = os.path.splitext(self._filePath)

            #Insert the extension in between the file name and extension
            self._filePath = preExt+suffix+Ext
        else:
            print('Whoops, no file path exists for this grid. Doing nothing.')

    def saveGrid(self, outfilePath: str, GDALDRIVERNAME : str = 'GTiff', GDALDATATYPE : int = gdal.GDT_Float32):
        '''
        Save as a georeferenced grid using gdal.  Function writes the data in the grid into a georeferenced dataset
        of type specified by GDALDRIVERNAME, a string, options here: http://www.gdal.org/formats_list.html

        :param outfilePath: the path to save the file as
        :param GDALDRIVERNAME: The name of the gdal drive to use to save the data (default: save as geotiff)
        :param GDALDATATYPE: The datatype to use (e.g., precision) in saving the raster (default: gdal.GDT_Float32)
        '''


        # Initialize new data
        drvr = gdal.GetDriverByName(GDALDRIVERNAME)  # Get the desired driver
        outRaster = drvr.Create(outfilePath, self._ncols, self._nrows, 1, GDALDATATYPE)  # Open the file

        # Write geographic information
        outRaster.SetGeoTransform(self._geotransform)  # Steal the coordinate system from the old dataset

        if not(self._projection is None):
            outRaster.SetProjection(self._projection)  # Steal the Projections from the old dataset

        # Write the array
        outRaster.GetRasterBand(1).WriteArray(self.grid)  # Writes my array to the raster

        #We may want to quickly resave over this file after some operations, store the necessary information to do so
        self._filePath = outfilePath
        self._gdalDriver = outRaster.GetDriver().ShortName
        self._gdalDataType = outRaster.GetRasterBand(1).DataType

        #Close out the memory for this raster
        outRaster = None

    def loadDerivedGrid(self,derivativeSuffix: str):
        '''
        For a grid derived from this one, that has a suffix to the file name that is known, we can load it based on
        this grids file name.
        :param derivativeSuffix:
        :return:
        '''

        #If there is an existing filepath
        if not(self._filePath is None):
            #Split the file path at its extension
            preExt,Ext = os.path.splitext(self._filePath)

            #Insert the extension in between the file name and extension
            derivFilePath = preExt+derivativeSuffix+Ext

            #Load this other grid
            loadedGrid = baseGrid([],rasterPath=derivFilePath)

        else:
            print('Whoops, no file path exists for this grid. Doing nothing.')
            loadedGrid = None

        return loadedGrid

    def duplicateGridWithNewArray(self, newGrid : np.ndarray, newSuffix: str = None):
        '''
        Creates a copy of these grid using the data from the specified grid
        :param newGrid: a numpy ndarray of the same size as this grid
        :param newSuffix: a string to append to the existing file name (if present)
        :return: returns a new grid object with the data from this array
        '''

        if newGrid.shape == self.shape:
            newGridObj = baseGrid(newGrid,geotransform=self._geotransform,projection=self._projection)

            #If this grid was loaded from a file, update file name information for the return grid with a new suffix
            if not(self._filePath is None) and isinstance(newSuffix,str):
                newGridObj._filePath = self._filePath
                newGridObj._gdalDriver = self._gdalDriver
                newGridObj._gdalDataType = self._gdalDataType
                newGridObj.updateFileNameWithSuffix(newSuffix)

        else:
            raise Exception('Whoops! Can not duplicate this grid with the specified array because it is the wrong shape')

        return newGridObj

    def _findMaskBoundaries(self, goodDataMask):
        '''

        TODO: Do I really want this to be a function of this class? Perhaps just a general function...

        This function seeks to find the cells at the edge of a mask. This could be the boundaries of a non-rectangular
        DEM stored within a rectangular grid (e.g., there being a lot of no-data values around the boundary of the grid)

        :param goodDataMask: A 2d numpy array (expected to be the same size as this grid)
        :return:
        '''

        # Function to find the cells at the edge of a dem. Dem is a ny x nx array, but may be largely padded
        # by nans. Determines where the edge of the real data is. Does this by finding the maximum value within a 3x3 kernel,
        # if that is not zero, but the original data at the corresponding location is, then that is an edge cell

        kern = np.ones((3, 3))

        if isinstance(goodDataMask,baseGrid) and (self.shape == goodDataMask.shape):
            mx = ndi.maximum_filter(goodDataMask.grid, size=None, footprint=kern)
            mn = ndi.minimum_filter(goodDataMask.grid, size=None, footprint=kern)
        elif isinstance(goodDataMask,np.ndarray) and (self.shape == goodDataMask.shape):
            mx = ndi.maximum_filter(goodDataMask, size=None, footprint=kern)
            mn = ndi.minimum_filter(goodDataMask, size=None, footprint=kern)
        else:
            raise Exception('Error, goodDataMask must be either a baseGrid instance or a numpy array, ' +
                            'and must be the same size as this grid.')

        edges = np.ones_like(goodDataMask)
        edges[1:-1, 1:-1] = 0

        # Return edge rows and columns as a tuple
        return np.where((((mx == 1) & (mn == 0) & goodDataMask)) | (
                goodDataMask & edges))