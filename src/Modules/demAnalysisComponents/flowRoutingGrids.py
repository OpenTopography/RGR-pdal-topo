'''

'''

import os
import sys
import numpy as np
from osgeo import gdal
from matplotlib import pyplot as plt

#Import local related components
try:
    from demAnalysisComponents.baseGrid import baseGrid
    from demAnalysisComponents.dem import demGrid
    from demAnalysisComponents.stablePriorityQueue import priorityQueue
    from demAnalysisComponents.FileSuffixDictionary import FILESUFFIXDICT
except:
    from baseGrid import baseGrid
    from dem import demGrid
    from stablePriorityQueue import priorityQueue
    from FileSuffixDictionary import FILESUFFIXDICT

#TODO: I had this as a baseGrid inherited class, presumably because I wanted to minimize the number of methods
#That could be calculated... but why not be able to calculate multiple derivatives with just one grid available?
class flowRoutingGrids(demGrid):
    '''
    TODO: Create a class for flow routing operations, that hosts sub-grids within it (e.g., flow directions, etc)

    TODO: How do I want to actually manage each grid? As each its own demgrid instance or as each their own numpy array?
    #There are differences in overhead here... but how drastic are these differences really...

    current plan - make this a 'baseGrid'. Give it an attribute 'grid' that is the original DEM elevation, ensure
    it has all other appropriate attributes.



    Should: Calculate new flow routing operations and save them, unless there are already appropriate grids present.

    '''

    #Perhaps I should place all these in some other file thats easy to access?
    #These could either by just a python file that everything loads in that sets some global variables...
    __AREA_SUFFIX = '_area'
    __FLOWDIR_SUFFIX = '_fd'
    __FILLEDGRID_SUFFIX = '_fill'
    __MAXFLOWLENGTH_SUFFIX = '_maxL'
    __MEANFLOWDIR_SUFFIX = '_meanDir'
    __BASINID_SUFFIX = '_basinIDs'
    __CHIGRID_SUFFIX = '_chi'
    __ORDER_SUFFIX = '_order'

    # Search kernel for D8 flow routing, the relative indices of each of the 8 points surrounding a pixel
    # ArcGIS convection
    # | grid indices |           |Downstream codes|
    # |i-1,j-1  i-1,j  i-1,j+1|  |32 64 128|
    # |i,j-1     i,j     i,j+1|  |16  X  1 |
    # |i+1,j-1  i+1,j  i+1,j+1|  |8   4  2 |
    __ROW_KERNEL = np.array([1, 1, 1, 0, 0, -1, -1, -1])
    __COL_KERNEL = np.array([-1, 0, 1, -1, 1, -1, 0, 1])

    # What are the flow direction codes necessary for those to be upstream of the current cell?
    __US_FLOW_CODES = np.array([128, 64, 32, 1, 16, 2, 4, 8])

    # What are the flow direction codes necessary for those to be downstream of the current cell?
    __DS_FLOW_CODES = np.array([8, 4, 2, 16, 1, 32, 64, 128])

    def __init__(self, baseGridInstance = None, filePath:str = None, forceReCalculation = False,
                 aggradationSlope: float = None, isAlreadyFilled = False, closedAreaGrid = None):
        '''

        :param baseGridInstance:
        :param filePath:
        :param forceReCalculation:
        :param aggradationSlope:
        '''

        #Store what the specified aggradation slope is (if this is None it will be calculated later)
        self._aggradationSlope = aggradationSlope
        self._filePath = filePath

        #Do we need to fill pits in this grid, or should they already be filled
        self.arePitsFilled = isAlreadyFilled

        #TODO: Need to test various loading/creation options, had to make an addition here...
        if not(baseGridInstance is None) and (filePath is None):
            self._createFromExistingDemGrid(baseGridInstance, closedAreaGrid)

        elif not(filePath is None):
            if self._areGridsAlreadyCalculated(filePath) and not(forceReCalculation):
                self._loadExistingGrids(filePath)
            else:
                baseGridInstance = baseGrid([], rasterPath=filePath)
                self._createFromExistingDemGrid(baseGridInstance, closedAreaGrid)
        else:
            #TODO: - update this to note what actually went wrong
            print('Whoops, a valid baseGrid was not specified or the existing files could not be found')

    def _createFromExistingDemGrid(self, demGridInstance: baseGrid, closedAreaMask: np.ndarray = None):
        '''
        TODO: Create this function, which needs to create all appropriate sub-grids from a demGrid and save them
        with the appropriate extensions
        :param demGridInstance:
        :param maskGrid:
        :return:
        '''

        #Get the basic georeferencing information
        self._copyAttributesFromGrid(demGridInstance)

        #Get the actual grid as this grid
        self.grid = np.copy(demGridInstance.grid)

        #If the existing grid has file path information, store that
        if not(demGridInstance._filePath is None):
            self._filePath = demGridInstance._filePath
            self._gdalDriver = demGridInstance._gdalDriver
            self._gdalDataType = demGridInstance._gdalDataType

        #Preform the actual flow routing
        self.fillAndRouteFlow(aggradationSlope=self._aggradationSlope,closedAreaMask= closedAreaMask)

    def _areGridsAlreadyCalculated(self,filePath):
        '''
        TODO: Write a function that returns a boolean if grids are already calculated
        :param filePath:
        :return:
        '''

        #Seperate the path and the extension
        path,extension = os.path.splitext(filePath)

        #Filled grid? (Currently I'm not storing this...)
        isFillCalcd = os.path.isfile(path + FILESUFFIXDICT['filled_elevation'] + extension)

        #Is the flow direction calculated?
        isFDCalcd = os.path.isfile(path + FILESUFFIXDICT['flow_dir'] + extension)

        #Is the drainage area calculated
        isAreaCalcd = os.path.isfile(path + FILESUFFIXDICT['drainage_area'] + extension)

        return isFDCalcd & isAreaCalcd & isFillCalcd


    def _loadExistingGrids(self, filePath):
        '''

        :param filePath:
        :return:
        '''

        #Load in the elevation grid
        baseGridInstance = baseGrid([], rasterPath=filePath)

        #Get the basic georeferencing information
        self._copyAttributesFromGrid(baseGridInstance)

        # Get the information from the existing file needed to resave as the same type
        self._filePath = baseGridInstance._filePath
        self._gdalDriver = baseGridInstance._gdalDriver
        self._gdalDataType = baseGridInstance._gdalDataType

        #Get the actual grid as this grid
        self.grid = np.copy(baseGridInstance.grid)

        #Load in the drainage area and flow direction grids
        fdGrid = self.loadDerivedGrid(FILESUFFIXDICT['flow_dir'])
        areaGrid = self.loadDerivedGrid(FILESUFFIXDICT['drainage_area'])
        fillGrid = self.loadDerivedGrid(FILESUFFIXDICT['filled_elevation'])

        #Just store the actual numpy arrays for these
        self.flowDirGrid = fdGrid.grid
        self.areaGrid = areaGrid.grid
        self.filledGrid = fillGrid.grid #Hmm, need to store the fill grid to as it allows for top-down search sorting...

    def saveFlowGrids(self, baseDEMFilePath: str = None, GDALDRIVERNAME : str = 'GTiff', GDALDATATYPE : int = gdal.GDT_Float32):
        '''
        Save as a georeferenced grid using gdal.  Function writes the data in the grid into a georeferenced dataset
        of type specified by GDALDRIVERNAME, a string, options here: http://www.gdal.org/formats_list.html

        :param outfilePath: the path to save the file as
        :param GDALDRIVERNAME: The name of the gdal driver to use to save the data (default: save as geotiff)
        :param GDALDATATYPE: The datatype to use (e.g., precision) in saving the raster (default: gdal.GDT_Float32)
        '''

        #If there is no file path, use the existing file paths
        if baseDEMFilePath is None:
            if self._filePath is None:
                self._filePath = 'UnnamedGrid'
                print('Whoops, no existing file path found, proceeding with filename ''{}'' '.format(self._filePath))
            baseDEMFilePath = self._filePath

        #Save the dem by creating a base grid instance of it
        self.saveGrid(baseDEMFilePath, GDALDRIVERNAME, GDALDATATYPE)

        #TODO: Check on whether this is an appropriate way to do this... this creates baseGrid instances b/c duplicate grid isn't set for this class... which we don't necessarily want to do
        #Save the flow direction grid
        fdGridInst = self.duplicateGridWithNewArray(self.flowDirGrid,newSuffix=FILESUFFIXDICT['flow_dir'])
        fdGridInst.overwriteSaveFile() #This will remember the save parameters from the resave through save grid
        fdGridInst = None #Clear out the memory assigned to this duplicate instance

        #Save the drainage area grid
        areaGridInst = self.duplicateGridWithNewArray(self.areaGrid,newSuffix=FILESUFFIXDICT['drainage_area'])
        areaGridInst.overwriteSaveFile()
        areaGridInst = None

        #Save the filled elevation area grid
        #TODO: Do I actually need to store the filled grid? I was thinking I needed to do that for top-down sorting algorithms, but can't I also do that with drainage area?
        fillGridInst = self.duplicateGridWithNewArray(self.filledGrid,newSuffix=FILESUFFIXDICT['filled_elevation'])
        fillGridInst.overwriteSaveFile()
        fillGridInst = None


    def fillAndRouteFlow(self,aggradationSlope:float = None, closedAreaMask:np.ndarray = None,
                         areaWeighting:np.ndarray = None):
        '''

        :param aggradationSlope:
        :param closedAreaMask:
        :return:
        '''

        #If no aggradation slope was specified, generate one that seems like it would be a good fit
        if aggradationSlope is None:
            aggradationSlope = self._guessAggradationSlope()

        self._aggradationSlope = aggradationSlope

        # print('Performing grid filling.')
        if self.arePitsFilled:
            self.filledGrid = self.grid
        else:
            self.filledGrid = self._priorityFlood(aggSlope = self._aggradationSlope, closed = closedAreaMask)

        # print('Grid filled, calculating drainage area')
        self._calcD8Area(self.filledGrid,closedAreaMask = closedAreaMask,areaWeighting=areaWeighting )

    def plotDrainageArea(self, axs = None, minArea = None, logTransform = True, **kwargs):
        '''
        Plot the drainage area grid, masking out area below 'minArea' if specified.
        :param axs:
        :param minArea:
        :param logTransform:
        :param kwargs:
        :return: axs, the axis that was plotted on
        '''

        #If no axis was specified for plotting, create one
        if axs is None:
            f,axs = plt.subplots(1,1)

        #Here I am going to mask out area's below the minimum with transparency.
        plotArea = np.copy(self.areaGrid)

        #If a minimum area was specified, need to mask out what was provided
        if not(minArea is None):
            plotArea[plotArea <= minArea] = np.nan

        #If this is log transformed, plot accordingly
        if logTransform:
            axs.imshow(np.log10(plotArea), extent=self.getGridExtent(),**kwargs)
        else:
            axs.imshow(plotArea, extent=self.getGridExtent(), **kwargs)


        return axs


    def _priorityFlood(self,aggSlope=0.0, closed=None):
        # dem is a numpy array of elevations to be flooded, aggInc is the minimum amount to increment elevations by moving upstream
        # use priority flood algorithm described in  Barnes et al., 2013
        # Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models
        # NOTE: They have another algorithm to make this more efficient, but to use that and a slope makes things more
        # complicated

        # Copy the dem so that we preserve the unfilled elevations
        filledDEM = self.grid.copy()
        openQueue = priorityQueue()  # priority queue to sort filling operation

        if closed is None:
            # Create a grid to keep track of which cells have been filled
            closed = np.zeros(self.shape)
        else:
            closed = 1.0*closed

        #Specify that things that were not closed get bad data values
        filledDEM[closed==1] = np.nan

        # Add all the edge cells to the priority queue, mark those cells as draining (not closed)
        edgeRows, edgeCols = self._findMaskBoundaries(np.logical_not(closed))

        #For each of the edge rows and cols
        for row,col in zip(edgeRows,edgeCols):
            #Mark this cell as 'closed' (noting that it drains to boundary).
            closed[row, col] = True

            # store the indices as a vector of row column, in the priority queue prioritized by the dem value
            openQueue.put(filledDEM[row, col], [row,col])

        # While there is anything left in the priority queue, continue to fill holes
        while not openQueue.isEmpty():

            #Get the highest priority (lowest elevation) cell thats currently been visited
            elevation, rowCol = openQueue.get()
            row, col = rowCol

            #Find the neighbors of this cell
            neighborRows, neighborCols, dxs = self._getNeighborIndices(row, col)

            # Look through the upstream neighbors
            for i in range(len(neighborCols)):
                #If this cell hasn't already been identified as closed
                if not closed[neighborRows[i], neighborCols[i]]:

                    # If this was a hole (lower than the cell downstream), fill it
                    if filledDEM[neighborRows[i], neighborCols[i]] <= elevation:
                        filledDEM[neighborRows[i], neighborCols[i]] = elevation + aggSlope * dxs[i]

                    #Mark the cell as closed
                    closed[neighborRows[i], neighborCols[i]] = True

                    #Add this cell to the priority queue
                    openQueue.put(filledDEM[neighborRows[i], neighborCols[i]], [neighborRows[i], neighborCols[i]])

        return filledDEM

    def _calcD8Area(self,fillGrid, closedAreaMask=None, areaWeighting=None):

        pxlArea = np.abs(self._dx * self._dy)  # area of a pixel

        # Get the sorted indices of the array in reverse order (e.g. largest first)
        #TODO: I often want to put 'nan' into grids to mask out values - but this can screw up the sorting
        idcs = fillGrid.argsort(axis=None)[::-1]

        # All pixels have at least their own area
        self.areaGrid = pxlArea * np.ones(self.shape)

        # Mask out the indexes outside of the mask
        if not closedAreaMask is None:
            allRows, allCols = np.unravel_index(idcs, self.shape)
            gdIdcs = ~closedAreaMask[allRows, allCols]
            idcs = idcs[gdIdcs]

            #Set the out-of-bounds areas to have a no-data value for drainage area
            self.areaGrid[closedAreaMask] = np.nan

        # Weight the initial drainage area
        if not areaWeighting is None:
            self.areaGrid *= areaWeighting

        #Preallocate space for the flow direction grid
        self.flowDirGrid = np.zeros(self.shape, dtype=int)

        for idx in idcs:  # Loop through all the data in sorted order
            [i, j] = np.unravel_index(idx, self.shape)  # Get the row/column indices

            if not np.isnan(self.grid[i, j]):
                iNeighbs, jNeighbs, dxs = self._getNeighborIndices(i, j)  # Find the actual indices of the neighbors

                # Assign the flow direction of the current point
                thisFD, iRel, jRel = self._assignFlowDir(iNeighbs - i, jNeighbs - j,
                                                         (fillGrid[i, j] - fillGrid[iNeighbs, jNeighbs]) / dxs)

                #So long as this cell actually drained somewhere...
                if not np.isnan(thisFD):
                    # accumulate current area, downstream area
                    self.flowDirGrid[i, j] = thisFD
                    self.areaGrid[i + iRel, j + jRel] += self.areaGrid[i, j]

    def findBoundaryOutlets(self,gridMask = None):
        '''
        Find all the row,col positions of the outlets at the margin of this grid.

        Assumes that flow routing and pit filling worked correctly, such that outlets are those points in
        the grid with flowdirections of '0'

        :param gridMask:
            boolean np.ndarray of those values that should be included in the search for outlets.
        :return:
            outletRowCols: array of row,col pairs of the indices in the grid that have flow directions of 0.
        '''

        outDraining = self.flowDirGrid == 0
        if not(gridMask is None):
            outDraining = outDraining & gridMask

        return np.argwhere(outDraining)

    def calcD8SlopeGrid(self, returnAsArray = False):
        '''
        
        
        '''

        # Preallocate space for the final answers
        slopeGrid = np.zeros(self.shape, dtype=np.float)

        # Loop through all the indices in sorted order
        for i in range(self._nrows):
            for j in range(self._ncols):
                # Get the index of the point downstream of this and the distance to that point
                [downI, downJ, inBounds] = self._getFlowToCell(i, j, self.flowDirGrid[i, j])

                # So long as we are not draining to the outskirts, proceed
                if inBounds:
                    dL = np.sqrt((self._dy*(downI-i))**2 + (self._dx*(downJ - j))**2)
                    slopeGrid[i,j] = (self.grid[i,j] - self.grid[downI,downJ])/dL

        # If it was requested, return these outputs as arrays
        if not (returnAsArray):
            slopeGrid = self.duplicateGridWithNewArray(slopeGrid,newSuffix=FILESUFFIXDICT['slope_d8'])

        return slopeGrid

    def calcChiGrid(self,A_0, theta, outlets = None, Amin = None, returnAsArray = False,
                    doUseRecursion = True):
        '''
        TODO: benchmark
        :param A_0:
        :param theta:
        :param outlets:
        :param returnAsArray:
        :return:
        '''

        if outlets is None:
            outlets = self.findBoundaryOutlets()

        if Amin is None:
            Amin = self._dx*self._dy

        #Preallocate space for this grid (can I do recursive searches that take as inputs a numpy grid?)
        chiGrid = np.zeros(self.shape)

        if doUseRecursion:
            #For each outlet, need to recursively search upstream and incrementally integrate chi
            for i,outlet in enumerate(outlets):
                self._recursiveUpstreamChiSearch(outlet,A_0,theta,chiGrid, Amin)

        else:
            for i,outlet in enumerate(outlets):
                self._nonRecursiveUpstreamChiSearch(outlet,A_0,theta,chiGrid, Amin)

        if not(returnAsArray):
            chiGrid = self.duplicateGridWithNewArray(chiGrid,newSuffix=FILESUFFIXDICT['chi_grid'])

        return chiGrid

    def calcOrderGrid(self,channelNetworkMask = None, returnAsArray = False):
        '''

        :param channelNetworkMask:
        :param returnAsArray:
        :return:
        '''
        # Get the sorted indices of drainage area from smallest to largest
        idcs = self.areaGrid.argsort(axis=None)

        # Mask out the indexes outside of the mask
        if not (channelNetworkMask is None):
            allRows, allCols = np.unravel_index(idcs, self.shape)
            gdIdcs = channelNetworkMask[allRows, allCols]
            idcs = idcs[gdIdcs]

        # Preallocate space for the final answers
        orderGrid = np.zeros(self.shape, dtype=np.float)

        # Loop through all the indices in sorted order
        for idx in idcs:
            # Get the currents row column index
            [i, j] = np.unravel_index(idx, self.shape)  # Get the row/column indices

            # Get the index of the point downstream of this and the distance to that point
            [downI, downJ, inBounds] = self._getFlowToCell(i, j, self.flowDirGrid[i, j])

            # So long as we are not draining to the outskirts, proceed
            if inBounds:

                #If there are two cells of equal order coming together here, increment them
                #TODO: Is there a wolrd where the order that cells are visited matters here? Perhaps when multiple
                #relatively high-order cells are coming into a single cell?
                if orderGrid[i,j] == orderGrid[downI,downJ]:
                    orderGrid[downI,downJ]+=1
                elif orderGrid[i,j] > orderGrid[downI,downJ]:
                    orderGrid[downI,downJ] = orderGrid[i,j]

        # If it was requested, return these outputs as arrays
        if not (returnAsArray):
            orderGrid = self.duplicateGridWithNewArray(orderGrid,newSuffix=FILESUFFIXDICT['stream_order'])

        return orderGrid

    def _nonRecursiveUpstreamChiSearch(self, outlet,A_0,theta,chiGrid, Amin):
        '''


        '''

        outletsInQueue = [outlet]



        while outletsInQueue:

            #Get the first outlet that was added
            outlet_i = outletsInQueue.pop(0)

            # Find all the cells immediately upstream of this cell
            usNeighbors, dxs = self._findUpstreamNeighbors(outlet_i[0], outlet_i[1])

            # For each of the upstream cells
            for i, usNeighbor in enumerate(usNeighbors):
                # If this isn't too small a drainage area
                if self.areaGrid[usNeighbor[0], usNeighbor[1]] >= Amin:
                    # Integrate the chi value to this upstream neighbor
                    chiGrid[usNeighbor[0], usNeighbor[1]] = chiGrid[outlet[0], outlet[1]] + \
                                                            (A_0 / self.areaGrid[usNeighbor[0],usNeighbor[1]]) ** theta * dxs[
                                                                i]
                    #Add this neighbor back into the search list
                    outletsInQueue.append(usNeighbor)

    def _recursiveUpstreamChiSearch(self,outlet,A_0,theta,chiGrid, Amin):
        '''

        :param outlet:
        :param A_0:
        :param theta:
        :param chiGrid:
        :return:
        '''

        #Find all the cells immediately upstream of this cell
        usNeighbors,dxs = self._findUpstreamNeighbors(outlet[0], outlet[1])

        #If there are any upstream cells proceed (else: recursion is over)
        if len(usNeighbors)>0:
            #For each of the upstream cells
            for i,usNeighbor in enumerate(usNeighbors):
                #If this isn't too small a drainage area
                if self.areaGrid[usNeighbor[0],usNeighbor[1]] >= Amin:

                    #Integrate the chi value to this upstream neighbor
                    chiGrid[usNeighbor[0],usNeighbor[1]] = chiGrid[outlet[0],outlet[1]] +\
                                                           (A_0/self.areaGrid[usNeighbor[0],usNeighbor[1]])**theta*dxs[i]
                    #Proceed upstream with another level of recursion
                    self._recursiveUpstreamChiSearch(usNeighbor,A_0,theta,chiGrid,Amin)



    def _assignFlowDir(self, iRel, jRel, slopes):
        '''

        :param iRel:
        :param jRel:
        :param slopes:
        :return:
        '''
        ## iRel and jRel are the relative indices from the current point to the surrounding points, the slopes of which are
        ## stored in 'slopes'

        # Search kernel for D8 flow routing, the relative indices of each of the 8 points surrounding a pixel, this is
        # ArcGIS convection
        # |i-1,j-1  i-1,j  i-1,j+1|  |32 64 128|
        # |i,j-1     i,j     i,j+1|  |16  X  1 |
        # |i+1,j-1  i+1,j  i+1,j+1|  |8   4  2 |

        idx = np.argmax(slopes)  # Find steepest surrounding slope
        iOut = iRel[idx]
        jOut = jRel[idx]  # Find the index of the steepest surrounding slope

        fd = np.nan

        #Think this will work for different flow routings that are specified as routing to single cell with array
        fdArray = self.__DS_FLOW_CODES[(iOut == self.__ROW_KERNEL) & (jOut == self.__COL_KERNEL)]

        #If this array doesn't have one item, fd will stay nan. Maximum slope also needs to be positive... (flowing away)
        if len(fdArray) == 1 & (slopes[idx] > 0):
            fd = fdArray[0]

        #Old static search way
        # if iOut == 0 and jOut == 1:
        #     fd = 1
        # elif iOut == 1 and jOut == 1:
        #     fd = 2
        # elif iOut == 1 and jOut == 0:
        #     fd = 4
        # elif iOut == 1 and jOut == -1:
        #     fd = 8
        # elif iOut == 0 and jOut == -1:
        #     fd = 16
        # elif iOut == -1 and jOut == -1:
        #     fd = 32
        # elif iOut == -1 and jOut == 0:
        #     fd = 64
        # elif iOut == -1 and jOut == 1:
        #     fd = 128

        return fd, iOut, jOut

    def _guessAggradationSlope(self):
        '''
        TODO: Think about how to actually do this... perhaps looking at grid spacing? Or average slopes in the grid?

        :return: aggradationSlope: The slope that will be used to fill closed cells in the DEM
        '''

        return 1e-7

    def _getNeighborIndices(self, row, col):
        '''
        Get the indices of the neighboring cells according to this flow grid search kernel.
        :param row: the row index of the current cell
        :param col: the column index of the current cell
        :return: neighborRows,neighborCols, neighborDistance : numpy ndarrays of the row, column, and distance
        '''
        # Find all the surrounding indices
        outRows = self.__ROW_KERNEL + row
        outCols = self.__COL_KERNEL + col
        dxs = np.sqrt(((outRows - row) * self._dx) ** 2 + ((outCols - col) * self._dy) ** 2)

        # Determine which indices are out of bounds
        inBounds = self._isRowsColsInBounds(outRows,outCols)

        return outRows[inBounds], outCols[inBounds], dxs[inBounds]

    def calculateMaxLMeanDir(self,mask = None, returnAsArray = False):
        '''
        :param flowDirectionGrid:
        :param elevationGrid:
        :param dx:
        :param noDataValue:
        :param mask:
        :return:
        '''
        #returns Two grids of the same size as those input, maxTransportLength, meanDirs, the mean direction to that point.
        #Input grids are  1) flowDirectionGrid, a grid of flow directions in ArcGIS convention
        # 2) elevationGrid, a grid of the elevations - this is used to choose the order in which we move downstream.
        # The elevation grid needs to have the pits filled and be hydrologically conditions to assure we traverse the grid
        # in the appropriate order

        ## Notes: I have an if statement in here to prevent calculations through negative elevations

        # Get the sorted indices of drainage area from smallest to largest
        idcs = self.areaGrid.argsort(axis= None)

        #Mask out the indexes outside of the mask
        if not mask is None:
            allRows,allCols = np.unravel_index(idcs,self.shape)
            gdIdcs = mask[allRows,allCols]
            idcs = idcs[gdIdcs]

        #Preallocate space for the final answers
        maxTransportLength = np.zeros(self.shape,dtype=np.float)
        delXPath = np.zeros(self.shape, dtype=np.float)
        delYPath = np.zeros(self.shape,dtype=np.float)

        #Loop through all the indices in sorted order
        for idx in idcs:
            #Get the currents row column index
            [i, j] = np.unravel_index(idx, self.shape)  # Get the row/column indices

            #Get the index of the point downstream of this and the distance to that point
            [downI,downJ,inBounds] = self._getFlowToCell(i, j,self.flowDirGrid[i,j])

            #So long as we are not draining to the outskirts, proceed
            if inBounds and not(np.isnan(self.grid[downI,downJ])):

                #How far do we move to the next cell?
                newDist = np.sqrt((self._dy*(downI-i))**2 + (self._dx*(downJ - j))**2)
                newL = maxTransportLength[i,j] + newDist

                #Keep track of the direction, length, and number of cells of the longest flow path
                if maxTransportLength[downI,downJ] < newL:
                    maxTransportLength[downI,downJ] = newL
                    delXPath[downI,downJ] = delXPath[i,j] + self._dx*(downJ - j)
                    delYPath[downI,downJ] = delYPath[i,j] + self._dy*(i - downI)

        meanDirs = np.arctan2(delYPath,delXPath)

        #If this cell doesn't have flow into it, don't count it as having a direction
        meanDirs[maxTransportLength == 0] = np.nan

        #If it was requested, return these outputs as arrays
        if not(returnAsArray):
            maxTransportLength = self.duplicateGridWithNewArray(maxTransportLength,
                                                                newSuffix=FILESUFFIXDICT['max_flow_length'])
            meanDirs = self.duplicateGridWithNewArray(meanDirs,newSuffix=FILESUFFIXDICT['mean_flow_direction'])

        return maxTransportLength, meanDirs

    def _getFlowToCell(self, row, col, flowDirCode):
        '''
        TODO: NEED TO REFACTOR AND TEST.
        TODO: THINK ABOUT BETTER WAY TO ACCOMODATE DIFFERENT FLOW ROUTINGS
        :param row:
        :param col:
        :return:
        '''
        # Function to get the indices of the cell that is drained to based on the flow direction specified in fd

        iOut = None
        jOut = None
        isGood = False

        if flowDirCode == 1 and col + 1 < self._ncols:
            iOut = row
            jOut = col + 1
        elif flowDirCode == 2 and row + 1 < self._nrows and col + 1 < self._ncols:
            iOut = row + 1
            jOut = col + 1
        elif flowDirCode == 4 and row + 1 < self._nrows:
            iOut = row + 1
            jOut = col
        elif flowDirCode == 8 and row + 1 < self._nrows and col - 1 >= 0:
            iOut = row + 1
            jOut = col - 1
        elif flowDirCode == 16 and col - 1 >= 0:
            iOut = row
            jOut = col - 1
        elif flowDirCode == 32 and row - 1 >= 0 and col - 1 >= 0:
            iOut = row - 1
            jOut = col - 1
        elif flowDirCode == 64 and row - 1 >= 0:
            iOut = row - 1
            jOut = col
        elif flowDirCode == 128 and row - 1 >= 0 and col + 1 < self._ncols:
            iOut = row - 1
            jOut = col + 1

        if not (iOut is None):
            isGood = True

        return iOut, jOut, isGood

    def returnMaskedGrid(self,goodValueMask,returnAsArray: bool = False, newSuffix: str = '_masked'):
        '''

        :param goodValueMask:
        :param returnAsArray:
        :param newSuffix:
        :return:
        '''
        #TODO: This method should override basegrid version of this function, otherwise when masking weird behavior
        #may arise.
        print('Need to implement this.')
        return None

    def findNonOverlappingBasins_areaRange(self,minArea,maxArea, returnAsArray = False,doUseRecursion = False,
                                           gridMask = None):
        '''

        :param minArea:
        :param maxArea:


        :return basinIDgrid:
            A grid of the same size as this grid, with values corresponding to the index of the outlet *in the original
            array*.
        :return unqOutlets:
            A list of row,col pairs of outlet coordinates. This is those outlets in the input list that did not
            represent sub-basins of a larger basin.
        '''

        if gridMask is None:
            gridMask = np.ones_like(self.areaGrid)==1

        #Get all the cells in this area range
        inAreaRange = (self.areaGrid >= minArea) & (self.areaGrid <= maxArea) & gridMask

        #Find the index of all of these cells
        outletIdcs = np.argwhere(inAreaRange)

        #get the actual indices and grid
        basinGrid,outletRowCols = self.findNonOverlappingBasins(outletIdcs,returnAsArray=True,
                                                                doUseRecursion = doUseRecursion)

        #Re-number the basin ID grid to correspond to these indices
        basinIDcopy = np.ones_like(basinGrid)*np.nan

        #Loop through the outlet rows and columns
        for i, outlet in enumerate(outletRowCols):

            #Get the row,col of this outlet
            r,c = outlet

            #Find all the positions in the original grid that match this outlets value
            thisBasinID = basinGrid[r,c]
            thisBasinMask = basinGrid == thisBasinID

            basinIDcopy[thisBasinMask] = i

        #Convert to a grid, unless requested otherwise
        if not(returnAsArray):
            basinIDcopy = self.duplicateGridWithNewArray(basinIDcopy,newSuffix=self.__BASINID_SUFFIX)

        return basinIDcopy,outletRowCols

    def findNonOverlappingBasins(self,outletRowCols, returnAsArray = False, doUseRecursion = False):
        '''
        Returns a grid, the same shape as this one, with cells upstream of each outlet coded by the index of the
        outlet in the specified input list of outlets.

        TODO: Should I make this flexible to accept x,y locations as well?
        :param outlets:
            list of row,col pairs of outlet coordinates

        :return basinIDgrid:
            A grid of the same size as this grid, with values corresponding to the index of the outlet *in the original
            array*.
        :return unqOutlets:
            A list of row,col pairs of outlet coordinates. This is those outlets in the input list that did not
            represent sub-basins of a larger basin.
        '''

        unqOutlets = [] #Store a list of those row,col indices that actually generate non-overlapping basins

        #Get all the areas for the basins
        outletAreas = np.array([self.areaGrid[o[0],o[1]] for o in outletRowCols])

        #Get the sorted order of these areas
        srtOutlets = np.argsort(outletAreas)[::-1] #Sort from largest to smallest

        #Preallocate space to store the array of numbers
        basinIDgrid = np.ones(self.shape)*np.nan #Set entire grid to nan to start

        #Loop through all the outlets in order from largest area to smallest
        for idx in srtOutlets:

            #Get the row,col of this point
            r, c = outletRowCols[idx]

            #First, make sure this outlet hasn't already been assigned as part of another basin
            if np.isnan(basinIDgrid[r,c]):

                #Find the upstream cells
                basinIndices = self.findBasinIndices((r,c),doUseRecursion=doUseRecursion)

                #Assign all these upstream cells this basin index
                basinIDgrid[basinIndices[:,0],basinIndices[:,1]] = idx

                #Add this outlet to the list of unqiue outlets
                unqOutlets.append([r,c])

        #Unless return requested as an array, convert to a grid
        if not(returnAsArray):
            basinIDgrid = self.duplicateGridWithNewArray(basinIDgrid,newSuffix=self.__BASINID_SUFFIX)

        return basinIDgrid, unqOutlets

    def findBasinIndices(self,outletRowCol, doUseRecursion = False):
        '''

        :param outletRowCol:
        :return:
        '''
        upstreamCells = [[outletRowCol[0], outletRowCol[1]]]  # Create a 2-dimensional python list

        if doUseRecursion:
            self.__findUpstreamCells(outletRowCol, upstreamCells)
        else:
            self.__findUpstreamCells_nonRecursive(outletRowCol, upstreamCells)

        return np.array(upstreamCells)

    def __findUpstreamCells(self, outletRowCol, upstreamCells):
        '''
        TODO: Refactor this function and test
        :param upstreamCells:
        :param outletRowCol:
        :return:
        '''
        # Recursive upstream search for area within basins
        upstreamNeighbors = self._findUpstreamNeighbors(outletRowCol[0], outletRowCol[1])[0]

        if len(upstreamNeighbors) > 0:
            for neighbor in upstreamNeighbors:
                # TODO: This is outputting a weird format... need to think/ worrk on this Yeah no shit - though perhaps I fixed?... i was doing upstreamCells += [[neighbor[0], neighbor[1]]
                upstreamCells += [[neighbor[0], neighbor[1]]]
                self.__findUpstreamCells(neighbor, upstreamCells)

    def __findUpstreamCells_nonRecursive(self, outletRowCol, upstreamCells):
        '''
        TODO: Refactor this function and test
        :param upstreamCells:
        :param outletRowCol:
        :return:
        '''

        outletsInQueue = [outletRowCol]

        while outletsInQueue:

            # Get the first outlet that was added
            outlet_i = outletsInQueue.pop(0)

            # Find all the cells immediately upstream of this cell
            usNeighbors, dxs = self._findUpstreamNeighbors(outlet_i[0], outlet_i[1])

            # For each of the upstream cells
            for i, neighbor in enumerate(usNeighbors):
                    # Add this neighbor back into the search list
                    upstreamCells += [[neighbor[0], neighbor[1]]]
                    outletsInQueue.append(neighbor)

    def findGreatestAreaPathUpstream(self, outletRowCol, Amin=None):
        '''
        TODO: Refactor this function

        :param outletRowCol:
        :param Amin:
        :return:
        '''

        #If no area was specified, go until we run out of upstream cells
        if Amin is None:
            Amin = 0

        # Initialize output array for rows and columns
        pathRowCol = [[outletRowCol[0], outletRowCol[1]]]  # This is a python list, which is easy to grow

        # Upstream search for path up largest area from an outlet
        upstreamCells = np.array(self._findUpstreamNeighbors(outletRowCol[0], outletRowCol[1])[0])

        # While there is still area upstream, keep moving upstream
        while len(upstreamCells) > 0:
            # What are the drainage areas of neighboring nodes
            # What are the areas of the upstream cells
            neighbAreas = self.areaGrid[upstreamCells[:, 0], upstreamCells[:, 1]]

            # Which of the upstream cells has the largest area
            largestAreaIdx = np.argmax(neighbAreas)

            # What are the indices of the upstream neighbor with the largest area
            nextRowCol = upstreamCells[largestAreaIdx]
            pathRowCol.append([nextRowCol[0], nextRowCol[1]])  # Append this item to the python list

            # Find the next upstream cells
            upstreamCells = np.array(self._findUpstreamNeighbors(nextRowCol[0], nextRowCol[1])[0])

            #If we have gotten below the minimum area, empty the list to end the loop
            if neighbAreas[largestAreaIdx] <= Amin:
                upstreamCells = []


        return pathRowCol

    def _findUpstreamNeighbors(self, row, col):
        '''
        TODO: Refactor this function
        :param row:
        :param col:
        :return:
        '''

        rowKern = self.__ROW_KERNEL+row
        colKern = self.__COL_KERNEL+col
        dxs = np.sqrt(((rowKern - row)*self._dx)**2 + ((colKern - col)*self._dy)**2)

        # What is within the grid bounds?
        inBounds = self._isRowsColsInBounds(rowKern,colKern)

        # Trim codes
        rowKern = rowKern[inBounds]
        colKern = colKern[inBounds]
        dxs = dxs[inBounds]
        usFlowCodes = self.__US_FLOW_CODES[inBounds]

        # Find areas where flow codes match
        isUpstream = self.flowDirGrid[rowKern, colKern] == usFlowCodes

        # Index out those areas
        ## TODO: THERE IS CERTAINLY A MORE EFFICIENT WAY TO PASS/PRODUCE THIS OUTPUT
        rowKern = rowKern[isUpstream]
        colKern = colKern[isUpstream]
        dxs = dxs[isUpstream]
        outRowCols = np.vstack((rowKern, colKern)).T.tolist()

        return outRowCols, dxs

    def calcAreaBinnedmedianValues(self,grid: baseGrid, bins = 20, percRange = [2.5, 97.5]):
        '''
        Calculates the median values of the specified grid (which must be a baseGrid type or derivative, e.g., demGrid).

        Values are calculated within the bins specified by bins. If bins is an integer number of bins (as opposed to
        a list of bin edges), a set of bin edges that span the available drainage area range will be calculated such
        that they are evenly spaced in log space.

        :param grid: The grid values that we want to bin as a function of drainage area
        :param bins: the number of bins (integer), or actual bin values (list, array), that we want to calculate binned
        averages is. If a number of bins is specified, a range of values within even log-spacing is defined.
        :param percRange: The lower and upper percentile values of data within each bin
        :return:
            areaMidPoints (numpy ndarray) - the midpoints of the drainage area bins
            medianVals (numpy ndarray) - the median values of the specified grid within the computed bins
            percentileVals (numpy ndarray) - the upper and lower percentile bounds to extract from each bin of data
        '''

        if isinstance(grid,baseGrid):
            areaMidpoints, medianVals,percentileVals = grid.calcMedianValuesBinnedByAnotherGrid(self.areaGrid,
                                                                                                bins=bins,
                                                                                                percRange=percRange)
        else:
            areaMidpoints, medianVals, percentileVals = None,None,None
            raise Exception('Whoops, must supply an instance of baseGrid. Returning nothing')

        return areaMidpoints, medianVals,percentileVals

    def plotAreaBinnedMedianValues(self, grid: baseGrid, bins = 20, percRange = [2.5, 97.5],axs=None,
                                   yAxisLabel: str = 'Label me', **kwargs):
        '''

        Plot the values of another grid as medians within bins of drainage area, with error bars corresponding to
        the percentile ranges within those bins.

        This is a wrapper for calculations that first compute binned averages, and then plot using
        matplotlib.pyplot.errorbar

        :param grid: the grid whose values we want to look at as a function of drainage area
        :param bins: the number of bins (integer), or actual bin values (list, array), that we want to calculate binned
        averages is. If a number of bins is specified, a range of values within even log-spacing is defined.
        :param percRange: The lower and upper percentile values of data within each bin that should be plotted as
        error bars
        :param axs: The axis to plot on. Default is None, in which case a new axis will be created.
        :param yAxisLabel: The text to display on the y axis label
        :param kwargs: other plotting arguments to pass to matplotlib.pyplot.errorbar
        :return: axs - the axis of the plot
        '''


        if isinstance(grid,baseGrid):
            axs = grid.plotMedianValuesBinnedByOtherGrid(self.areaGrid, bins=bins, percRange=percRange, axs=axs,
                                                         yAxisLabel=yAxisLabel, xAxisLabel=r'Area (m^2)', **kwargs)
        else:
            axs = None
            raise Exception('Whoops, must supply an instance of baseGrid. Returning nothing.')

        return axs

    def plotGridAgainstArea(self, grid: baseGrid, plotFractionOfObservations: float = 1, axs = None,
                            yAxisLabel = 'Label me', **kwargs):
        '''

        Plots the values of the specified grid on the y axis against the drainage area of this grid on the x axis.
        Optionally sub-samples grid values from binned intervals based on plotFractionOfObservations. When this
        is between 0 and 1, we will only grab a random fraction of datapoints to plot - simplifying visualization. This
        random data points will be spread out within 20 bins on the x-axis, so that low data-density proportions of the
        plot still get some data shown.

        :param grid: the grid whose values should go on the y-axis. Can be a numpy ndarray or an instance of baseGrid
        :param fractionToSample: What proportion of data points should be plotted?
        :param axs: the matplotlib.pyplot axis to plot on. If None is specified, one will be created for this plot
        :return: axs - the matplotlib.pyplot axis of the plot
        '''

        if isinstance(grid,baseGrid):
            axs = grid.plotThisGridAgainstAnotherGrid(self.areaGrid,plotFractionOfObservations=plotFractionOfObservations,
                                                  axs = axs,xAxisLabel=r'Area (m^2)', yAxisLabel=yAxisLabel,**kwargs)
        else:
            print('Whoops, must supply a baseGrid instance. Not plotting anything.')

        return axs