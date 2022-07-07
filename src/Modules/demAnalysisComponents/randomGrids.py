'''
Functions for generating random grids
'''

try:
    from demAnalysisComponents.baseGrid import baseGrid
    from demAnalysisComponents.flowRoutingGrids import flowRoutingGrids
    from demAnalysisComponents.dem import demGrid
    from demAnalysisComponents.fftGrid import fftGrid
    from demAnalysisComponents.stablePriorityQueue import priorityQueue
except:
    from baseGrid import baseGrid
    from dem import demGrid
    from flowRoutingGrids import flowRoutingGrids
    import fftGrid
    from stablePriorityQueue import priorityQueue

from scipy import interpolate
import numpy as np
from matplotlib import pyplot as plt

class randomGrid(baseGrid):
    _filePath = None

    def __init__(self):
        pass

    def _calculateRandomInstanceArray(self):
        pass

    def getPermutation(self):
        '''

        :param gridToMatch:
        :return:
        '''

        randArray = self._calculateRandomInstanceArray()

        return demGrid(randArray, dx=self._dx,xllcenter=self._xllcenter,
                       yllcenter=self._yllcenter,geotransform=self._geotransform,projection = self._projection)

    def getGridEnsemble(self,nPermutations: int):
        '''

        :param nPermutations:
        :return:
        '''

        self.nGridPermutations = nPermutations
        self.ensembleGrids = np.zeros((self._nrows,self._nCols,nPermutations))

        for i in range(self.nGridPermutations):
            self.ensembleGrids[:,:,i] = self._calculateRandomInstanceArray()

        return self.ensembleGrids

class proceduralGrid_diamondSquare(randomGrid):

    def __init__(self,gridToMatch: baseGrid, roughness: float, startingNoiseScale : float= None,  randomSeed = None,
                 matchMeanElevationAndRelief = False):
        '''

        :param gridToMatch:
        :param roughness:
        :param startingNoiseScale:
        :param randomSeed:
        :param matchMeanElevationAndRelief:
        '''

        self.H = roughness
        self._randomState = np.random.RandomState(randomSeed)

        # Get the attributes from this instance
        self._copyAttributesFromGrid(gridToMatch)

        # Get values for the corners of this grid
        # self.origGrid = np.copy(gridToMatch.grid)

        #Get the starting scale of the random grid values - this will be used to populate the initial values of the alg.
        if startingNoiseScale is None:
            self._startingScale = 1.0*np.std(gridToMatch.grid)
        else:
            self._startingScale = startingNoiseScale

        #If requested, keep track of the original grids mean elevation and relief so that we can make sure that the
        # procedurally generated grids have the same attributes
        self._doMatchMeanElevAndRelief = matchMeanElevationAndRelief
        if matchMeanElevationAndRelief:
                self._matchMeanElev = gridToMatch.mean()
                self._matchRelief = gridToMatch.relief()


    def _getRandomValue(self,scale: float):
        '''

        :return:
        '''

        # return 2.0*scale*(self._randomState.random_sample() - 0.5)
        return scale * self._randomState.randn()

    def _calculateRandomInstanceArray(self):
        '''

        :param gridToMatch:
        :return:
        '''

        # Initialize grid with random values of same roughness scale, unlike traditional diamond-square, we
        # wont require a grid of size (2^n +1), will instead rely on wrapping indexing
        self.grid = self._randomState.randn(self._nrows,self._ncols)*self._startingScale

        # Starting width needs to be of size (2^n +1),
        width = int((2**np.ceil(np.log2(np.max(self.shape) - 1)) + 1))

        scale = 1.0*self._startingScale

        # self.iteration = 0

        #Preform the iterations
        while width > 1:
            # print('Iteration {}'.format(self.iteration))
            self._diamondSquareIteration(width,scale)
            width = int(width/2) #Make sure width stays an int
            scale/=2.0**self.H
            # self.iteration+=1

        #If we reduced the resolution of the grid, we need to re-interpolate to the real grid resolution
        # if self._doReduceResolution:
        #
        #     #TODO: Perhaps this should instead be a method of all dems? Resampling to a different grid?
        #     interpolateSubsampledGrid = interpolate.interp2d(self._X_s.flatten(),self._Y_s.flatten(), self.grid, kind = 'linear')
        #     self.grid = interpolateSubsampledGrid(self._X_0.flatten(),self._Y_0.flatten())

        #If we need to match the basic scale of the input grid, rescale the mean and relief to the input grid
        if self._doMatchMeanElevAndRelief:
            self.grid = (self.grid - self.mean())*self._matchRelief/self.relief() + self._matchMeanElev

        return self.grid

    def _diamondSquareIteration(self,width: int,scale: float):
        '''

        :param row:
        :param col:
        :param width:
        :param scale:
        :return:
        '''

        halfStep = int(width/2)

        ##Square iterations
        for r_i in range(halfStep,self._nrows,width):
            for c_i in range(halfStep,self._ncols,width):
                if (r_i < self._nrows) & (c_i < self._ncols):
                    self._squareUpdate(r_i, c_i, halfStep, scale)
                    # print('row, col - value: {},{} - {}'.format(r_i, c_i, self.grid[r_i, c_i]))

        ##diamond iterations (more complicated)

        #Need to keep track of if we are in an even or an odd row to offset diamond placement
        row = 0

        for r_i in range(0,self._nrows,halfStep):

            #On odd rows, start at the first column, on even rows, offset start by a half width
            if row % 2 == 0:
                for c_i in range(halfStep,self._ncols,width):
                    self._diamondUpdate(r_i, c_i, halfStep, scale)
                    # print('row, col - value: {},{} - {}'.format(r_i, c_i, self.grid[r_i, c_i]))
            else:
                for c_i in range(0,self._ncols,width):
                    self._diamondUpdate(r_i, c_i, halfStep, scale)
                    # print('row, col - value: {},{} - {}'.format(r_i, c_i, self.grid[r_i, c_i]))

            row += 1

    def _diamondUpdate(self, row: int, col: int, halfWidth: int, scale: float):
        '''

        :return:
        '''

        #Get the diamond kernel surrounding this row and column
        baseKernRows = row + np.array([-1, 0, 0, 1],dtype=np.int)*halfWidth
        baseKernCols = col + np.array([0, -1, 1, 0], dtype=np.int)*halfWidth

        #update the grid based on these values
        self._updateBasedOnRowsCols(row,col,baseKernRows,baseKernCols,scale)

    def _squareUpdate(self, row: int, col: int, halfWidth: int, scale: float):
        '''

        :return:
        '''

        #Get the square kernel surrounding this row and column
        baseKernRows = row + np.array([-1, -1, 1, 1], dtype=np.int)*halfWidth
        baseKernCols = col + np.array([-1, 1, -1, 1], dtype= np.int)*halfWidth

        # update the grid based on these values
        self._updateBasedOnRowsCols(row,col,baseKernRows,baseKernCols,scale)

    def _updateBasedOnRowsCols(self,row: int,col: int,rows: np.ndarray,cols: np.ndarray, scale: float):
        '''
        Update the row,col position of this grid based on the average of grid values at the locations in rows,cols
        and a random value

        :param row: row location of grid to update
        :param col: col location of grid to update
        :param rows: array of rows to average
        :param cols: array of cols to average
        :return:
        '''
        # Wrap any indices back to the other side of the grid
        inBounds = self._isRowsColsInBounds(rows,cols)
        rows = rows[inBounds]
        cols = cols[inBounds]

        # rows, cols = self._getWrappedRowsCols(rows, cols, equivalentEdges= True)
        self.grid[row,col] = np.mean(self.grid[rows,cols]) + self._getRandomValue(scale)


#####################################################################################################################
## Random flow routing grids
## These are procedural routines for generating random drainage networks in gridded data
#####################################################################################################################

class proceduralGrid_randomPriorityFilling(randomGrid):

    def __init__(self,gridToMatch: flowRoutingGrids, fillSlope: float, doScalePriorityByElevation :bool = False, randomSeed = None,
                 matchMeanElevationAndRelief = False, outlets = None, gridMask = None):
        '''
        This is a procedural terrain generation routine that modifies the priority flood algorithm described in:

        Barnes, R., Lehman, C., and Mulla, D., 2014, Priority-flood: An optimal depression-filling and
         watershed-labeling algorithm for digital elevation models: Computers and Geosciences, v. 62, p. 117–127,
          doi:10.1016/j.cageo.2013.04.024.

        Specifically, the algorithm fills terrain upstream from outlets, with channels allowed to 'grow' upstream one
        cell at a time based on randomly assigned priorities.  This approach was applied and described in:

            Johnstone, S.A., Finnegan, N.J., and Hilley, G.E., 2017, Weak bedrock allows north-south elongation of
             channels in semi-arid landscapes: Earth and Planetary Science Letters, v. 478, p. 150–158,
              doi:10.1016/j.epsl.2017.08.037.


        with an addition made in that work to iteratively correct large elevation discrepencies on either side of a
        divide produced by flow paths with different integrated drainage area paths.


        :param gridToMatch:
            A flowRoutingGrids instance to serve as the template for the random grid we will create. Provides the
            shape information, georeferencing, and optionally the mean elevation and relief

        :param fillSlope:
            The (constant) slope of all channel segments in the grid. The main use for this is to produce a
            grid that flows out to the specified outlets, so that other derivatives of the drainage network
            can be derived.

        :param doScalePriorityByElevation:
            A boolean that specifies whether the priority with which a network will propogate upstream scales
            with the current elevation of that point in the network. If False (the default) the priority with
            which filling occurs is purely random, which can lead to large elevation discrepencies across internal
            drainage divides. If True, filling is still random but the random priority is scaled based on how far
            up the network a cell is, this acts to minimize the elevation gaps in internal divides by limiting the
            likelihood that large channels curve back around on themselves (as opposed to having long channels tend
            to orient themselves parallel to the catchment axis).

        :param randomSeed (optional):
            A seed for numpy.random.RandomState . Default is None.

        :param matchMeanElevationAndRelief (optional):
            Whether to rescale the topography generated in order to match the mean elevation and relief of the input
            grid.

        :param outlets:
            A numpy array of shape n, 2. Specifies the n row,column coordinates of outlets to search upstream from.
            Generated flow will be out to these coordinates.

        :param gridMask:
            A boolean mask that limits the extent of the generated grid to within True values. For example,
            if a list of one basin coordinate was provided this could specify the basin extent to generate a random
            grid within.

        '''

        #Store parameters that will be necessary for later calculations
        self._aggSlope = fillSlope
        self._randomState = np.random.RandomState(randomSeed)

        # Get the attributes from this instance
        self._copyAttributesFromGrid(gridToMatch)

        # If there is not a mask specified, make the whole grid in bounds
        if isinstance(gridMask,baseGrid):
            gridMask = gridMask.grid == 1
        elif gridMask is None:
            gridMask = np.ones(gridToMatch.shape) == 1

        #Store the mask for this grid
        self._gridMask = gridMask

        # If there are not outlets specified, get all the outlets that drain to boundaries
        if outlets is None:
            outlets = gridToMatch.findBoundaryOutlets(gridMask)
        #If only a single outlet was specified, make it iterable
        elif isinstance(outlets,np.ndarray): #Is this a numpy array
            if outlets.shape == (2,):
                outlets = np.array([outlets])
        elif ~isinstance(outlets,np.ndarray) & len(outlets) == 2:
            outlets = np.array([np.array(outlets)])

        #Assign the outlets
        self._outlets = outlets

        #Set the priority function
        if doScalePriorityByElevation:
            #Set the priority getting function to account for current elevations
            self._getPriority = self.__randomPriority_elevationScaled
        else:
            #Set the priority getting function to just get random priorities
            self._getPriority = self.__randomPriority

        # If requested, keep track of the original grids mean elevation and relief so that we can make sure that the
        # procedurally generated grids have the same attributes
        self._doMatchMeanElevAndRelief = matchMeanElevationAndRelief
        if matchMeanElevationAndRelief:
            self._matchMeanElev = gridToMatch.mean()
            self._matchRelief = gridToMatch.relief()

    def getFlowGridPermutation(self):
        '''
        Create a random elevation model based on this random-filling procedural terrain generation defined here,
        then turn that into a model of flow routing.

        :return: flowGridPerm: An instance of flowRoutingGrids generated based on a random elevation model
        generated by the procedural terrain generation routine specified in this grid
        '''

        # Get a random elevation grid based on the procedure defined by this class
        gridPermutation = self.getPermutation()

        # Turn this random grid into a flow grid
        flowGridPerm = flowRoutingGrids(gridPermutation, forceReCalculation=True, aggradationSlope=self._aggSlope/100,
                                        isAlreadyFilled = True, closedAreaGrid= ~self._gridMask)

        return flowGridPerm


    def _calculateRandomInstanceArray(self):
        '''
        Preforms the random-priority upstream filling routine to produce a random grid
        :return:
        '''

        #Create a numpy array as the DEM
        self.grid = np.zeros(self.shape)

        if not self._gridMask is None:
            self.grid[~self._gridMask] = np.nan

        openQueue = priorityQueue()  # priority queue to sort filling operation

        #Mark all cells outside of this mask as closed
        closed = np.logical_not(self._gridMask == 1)

        # Add all the outlets to the priority queue, mark those cells as draining (not closed)
        for i in range(len(self._outlets)):
            #Get the outlet row col position
            row,col = self._outlets[i,0],self._outlets[i,1]

            # Mark this cell as 'closed' (noting that it drains to boundary).
            closed[row, col] = True

            #Get a random priority for this cell
            priority = self._randomState.rand(1)[0]

            #Store the indices as a vector of row column, in the priority queue prioritized by the dem value
            openQueue.put(priority, [row, col])

        # While there is anything left in the priority queue, continue to fill holes
        while not openQueue.isEmpty():

            # Get the highest priority (lowest elevation) cell thats currently been visited
            priority, rowCol = openQueue.get()
            row, col = rowCol

            # Find the neighbors of this cell
            neighborRows, neighborCols, dxs = self.getNeighborIndices(row, col)

            # Go through each of the neighbors, fill it as necessary, and add approriate items back to the queue
            self._fillGridElevations(neighborRows,neighborCols,row,col,dxs,closed, openQueue)

        # If we need to match the basic scale of the input grid, rescale the mean and relief to the input grid
        if self._doMatchMeanElevAndRelief:
            self.grid = (self.grid - self.mean()) * self._matchRelief / self.relief() + self._matchMeanElev

        return self.grid

    def _getSlope(self,row,col):
        '''

        :param row:
        :param col:
        :return:
        '''

        return self._aggSlope

    def _fillGridElevations(self, neighborRows, neighborCols, dsRow, dsCol, dsDistances, isClosed, openQueue):

        # Shuffle the order of the rows and columns
        shflIdx = np.arange(len(neighborRows))
        self._randomState.shuffle(shflIdx)
        neighborRows = neighborRows[shflIdx]
        neighborCols = neighborCols[shflIdx]
        dxs = dsDistances[shflIdx]

        # Look through the neighbors
        for i in range(len(neighborCols)):

            # If this cell hasn't already been identified as closed
            if not isClosed[neighborRows[i], neighborCols[i]]:

                #Assigns its elevation based on the elevation of the downstream cell and slope and distance to that cell
                self.grid[neighborRows[i], neighborCols[i]] = self.grid[dsRow, dsCol] + self._aggSlope*dsDistances[i]

                ### Set priority to be random, rather than a function of elevation
                priority = self._getPriority(neighborRows[i], neighborCols[i])

                ''' How do we ensure that boundary cells remain boundaries (e.g. low drainage area)? 
                If this cell is on the boundary (i.e., its neighbors contain out of bounds values), then
                give it a very high priority, ensure that its visited, but making it unlikely to be an
                intermediatte step in the flow.'''
                neighneighRows, neighneighCols, nnDX = self.getNeighborIndices(neighborRows[i], neighborCols[i])
                if np.any(~self._gridMask[neighneighRows, neighneighCols]):
                    priority = np.inf

                # So long as the priority wasn't nan, which specifies the termination of a flow path
                if ~np.isnan(priority):
                    # Add the item back to the queue with a random priority
                    openQueue.put(priority, [neighborRows[i], neighborCols[i]])

                # Mark the cell as closed
                isClosed[neighborRows[i], neighborCols[i]] = True


    def __randomPriority(self, row, col):
        '''
        This is the function that determines the priority with which a cell at the specified location will be
        visited in the queue that controls the filling. Redefining this function can modify how 'randomness' is
        assigned in the random flooding routine, for example by giving priority based on some aspect of the position
        of the cell (e.g., cells closer to the basin center have higher priority)
        :param row:
        :param col:
        :return:
        '''

        return self._randomState.rand(1)[0]

    def __randomPriority_elevationScaled(self, row, col):
        '''
        This is the function that determines the priority with which a cell at the specified location will be
        visited in the queue that controls the filling. Redefining this function can modify how 'randomness' is
        assigned in the random flooding routine, for example by giving priority based on some aspect of the position
        of the cell (e.g., cells closer to the basin center have higher priority).

        This version of the function scales the priority by the current elevation of the cell, so that low elevation
        cells will *tend* to be visited later.

        :param row:
        :param col:
        :return:
        '''

        return self._randomState.rand(1)[0]*self.grid[row,col]

    def getMaxLengthMeanDirectionPermutationGrids(self,nPermutations):
        '''

        :param nPermutations:
        :return:
        '''

        #Preallocate space for these arrays
        lengthGridPermutations = np.zeros((self.shape[0],self.shape[1],nPermutations))
        thetaGridPermutations = np.zeros_like(lengthGridPermutations)

        #For each permutation requested
        for i in range(nPermutations):
            #Get a random permutation of a flow grid
            flowGridPerm_i = self.getFlowGridPermutation()

            #Use that flow grid permutation to get maxL, meandir grids
            lengthGridPermutations[:,:,i], thetaGridPermutations[:,:,i] = flowGridPerm_i.calculateMaxLMeanDir(
                                                                        mask = self._gridMask, returnAsArray= True)


        return  lengthGridPermutations, thetaGridPermutations

    def lengthThetaRadialDensityPlot(self, referenceFlowGrid: flowRoutingGrids, lengthGrid = None, thetaGrid = None,
                                 lengthGridPermutations = None, thetaGridPermutations = None, nPermutations = 100,
                                 minL = None, maxL = None, Lbins = 20, dTheta = 3.1415/20,
                                 thetaWin = None, axs = None, **kwargs):


        #If the length, theta data aren't provided, calculate them
        if (lengthGrid is None) or (thetaGrid is None):
            lengthGrid, thetaGrid = referenceFlowGrid.calculateMaxLMeanDir(mask = self._gridMask, returnAsArray=True)
        else:
            # The third dimension of the permutations should be the number of permutations
            nPermutations = lengthGridPermutations.shape[2]

        #If permutations of length or theta grids aren't provided, calculate them
        if (lengthGridPermutations is None) or (thetaGridPermutations is None):
            lengthGridPermutations,thetaGridPermutations = self.getMaxLengthMeanDirectionPermutationGrids(nPermutations)


        #If no window was specified for the radial count
        if thetaWin is None:
            thetaWin = dTheta / 2.0

        # Construct the bins that the length counts will be in
        if isinstance(Lbins, int):
            if maxL is None:
                maxL = lengthGrid.max()
            if minL is None:
                minL = np.min(lengthGrid[lengthGrid > 0])

            # Defaulting to log spacing here...
            lengthBins = np.logspace(np.log10(minL), np.log10(maxL), Lbins + 1)

        # Find the midpoints of each grid coordinate
        lengthMids = (lengthBins[:-1] + lengthBins[1:]) / 2

        # Construct the bins that we will do counts in
        thetaMids = np.arange(-np.pi, np.pi + dTheta, dTheta)

        #Get the radial histogram for the actual grid
        binCount = self._radialLengthThetaHistogram(lengthGrid,thetaGrid,lengthBins, lengthMids,
                                                    thetaMids, thetaWin)

        #Get the radial histogram for the grid permutations
        binCountPerms = np.zeros((binCount.shape[0],binCount.shape[1],nPermutations))

        for i in range(nPermutations):
            binCountPerms[:,:,i] = self._radialLengthThetaHistogram(lengthGridPermutations[:,:,i],
                                                                    thetaGridPermutations[:,:,i],
                                                                    lengthBins, lengthMids, thetaMids, thetaWin)

        #TODO: Isn't there a way to one-line this? Something about these arrays isn't liking what *should* work:
        #mp.sum(binCount > binCountPerms, axis = 2)/np.float(nPermutations)
        binCountPercentiles = np.zeros_like(binCount)

        for i in range(nPermutations):
            binCountPercentiles+= 1.0*(binCount > binCountPerms[:,:,i])

        binCountPercentiles/= nPermutations

        # If no plotting axis was specified
        if axs is None:
            f = plt.figure()
            axs = plt.subplot(111, projection='polar')

        axs.pcolormesh(thetaMids, lengthMids, binCountPercentiles, **kwargs)

        return axs, binCount, binCountPerms, binCountPercentiles, lengthGridPermutations, thetaGridPermutations

    def _radialLengthThetaHistogram(self,lengthGrid,thetaGrid,lengthBins, lengthMids, thetaMids, thetaWin):

        # Preallocate space for the results grids
        binCount = np.zeros(
            (len(lengthMids), len(thetaMids)))  # A count of the number of pixels in each length-theta bin
        lengthCount = np.zeros_like(
            binCount)  # A count of the number of pixels at each length bin, used for normalization

        # Loop throught the length bins - alternatively, could perhaps do this with np.histogram2d, but would need
        #to get clever b/c of 'wrapped' nature of radial coordinates.
        for i in range(len(lengthMids)):
            theseLengths = (lengthGrid >= lengthBins[i]) & (lengthGrid < lengthBins[i + 1])

            for j in range(len(thetaMids)):
                thisThetaMid = thetaMids[j]

                # Find windows around this midpoint
                minTheta = thisThetaMid - thetaWin
                maxTheta = thisThetaMid + thetaWin

                # Need to determine whether these are appropriate
                if minTheta < -np.pi:
                    minTheta = (2.0 * np.pi) + minTheta
                    theseDirs = (thetaGrid >= minTheta) | (thetaGrid < maxTheta)
                elif maxTheta > np.pi:
                    maxTheta = -2.0 * np.pi + maxTheta
                    theseDirs = (thetaGrid >= minTheta) | (thetaGrid < maxTheta)

                else:
                    theseDirs = (thetaGrid >= minTheta) & (thetaGrid < maxTheta)

                binCount[i, j] = np.sum(theseDirs & theseLengths)
            lengthCount[i, :] = np.sum(binCount[i, :])

        return binCount


class proceduralGrid_longWavelengthTopoRandomFilling(proceduralGrid_randomPriorityFilling):


    def __init__(self,gridToMatch: flowRoutingGrids, fillSlope: float, minWavelength: float,
                 sigma_wavelength:float, randomSeed = None,
                 matchMeanElevationAndRelief = False, outlets = None, gridMask = None):
        '''
        This is a procedural terrain generation routine that modifies the priority flood algorithm described in:

        Barnes, R., Lehman, C., and Mulla, D., 2014, Priority-flood: An optimal depression-filling and
         watershed-labeling algorithm for digital elevation models: Computers and Geosciences, v. 62, p. 117–127,
          doi:10.1016/j.cageo.2013.04.024.

        Specifically, in this version flooding progresses upstream based on priorities that are assigned to be random.
        Those random priorities are selected from a uniform distribution, the upper bound of which is based on a
        filtered version of an original DEM. The allows there to be some variability in the resulting structure
        of drainages, but enforces the general long-wavelength topographic structure imposed by internal and external
        drainage divides.  Specifically, the maximum value of a random priority is given by the low-pass filtered
        version of the original topography, rescaled so that the filtered topography has a minimum elevation of 0
        and a relief of 1.

        :param flowGridToMatch:
            A flowRoutingGrids instance to serve as the template for the random grid we will create. Provides the
            shape information, georeferencing, and optionally the mean elevation and relief

        :param fillSlope:
            The (constant) slope of all channel segments in the grid. The main use for this is to produce a
            grid that flows out to the specified outlets, so that other derivatives of the drainage network
            can be derived.

        :param minWavelength:
           The wavelength used to preform low-pass filtering. This is the minimum wavelength allowed to pass
           through the filter without reduction.

        :param sigma_wavelength:
            The standard deviation of a gaussian window used to taper the low pass filter. See fftGrid.addFilter_lowPass
            for further details.

        :param randomSeed (optional):
            A seed for numpy.random.RandomState . Default is None.

        :param matchMeanElevationAndRelief (optional):
            Whether to rescale the topography generated in order to match the mean elevation and relief of the input
            grid.

        :param outlets:
            A numpy array of shape n, 2. Specifies the n row,column coordinates of outlets to search upstream from.
            Generated flow will be out to these coordinates.

        :param gridMask:
            A boolean mask that limits the extent of the generated grid to within True values. For example,
            if a list of one basin coordinate was provided this could specify the basin extent to generate a random
            grid within.

        '''

        #Store parameters that will be necessary for later calculations
        self._aggSlope = fillSlope
        self._randomState = np.random.RandomState(randomSeed)

        # Get the attributes from this instance
        self._copyAttributesFromGrid(gridToMatch)

        # If there is not a mask specified, make the whole grid in bounds
        if isinstance(gridMask,baseGrid):
            gridMask = gridMask.grid == 1
        elif gridMask is None:
            gridMask = np.ones(gridToMatch.shape) == 1

        #Store the mask for this grid
        self._gridMask = gridMask

        # If there are not outlets specified, get all the outlets that drain to boundaries
        if outlets is None:
            outlets = gridToMatch.findBoundaryOutlets(gridMask)
        #If only a single outlet was specified, make it iterable
        elif isinstance(outlets,np.ndarray): #Is this a numpy array
            if outlets.shape == (2,):
                outlets = np.array([outlets])
        elif ~isinstance(outlets,np.ndarray) & len(outlets) == 2:
            outlets = np.array([np.array(outlets)])

        #Assign the outlets
        self._outlets = outlets

        '''In this version of the grid we are going to scale priorities based on the long-wavelength topography
        of the grid to match.  Get that long wavelength topography and store its value'''
        self._priortyGrid = self._getPriorityGrid(gridToMatch,minWavelength,sigma_wavelength)

        # If requested, keep track of the original grids mean elevation and relief so that we can make sure that the
        # procedurally generated grids have the same attributes
        self._doMatchMeanElevAndRelief = matchMeanElevationAndRelief
        if matchMeanElevationAndRelief:
            self._matchMeanElev = gridToMatch.mean()
            self._matchRelief = gridToMatch.relief()


    def _getPriorityGrid(self,gridToMatch: flowRoutingGrids, minWavelength: float, sigma_wavelength: float):
        '''

        Filter the topography based on the provided minWavelength - this is a low-pass filter, returning the
        aspects of the topography defined by wavelengths greater than the provided min-wavelength.

        :param self:
        :param gridToMatch:
        :param minWavelength:
        :return:
        '''

        #Preform a fourier transform and filter the topography
        specGrid = fftGrid.fftGrid(gridToMatch,detrendGrid=True)
        specGrid.addFilter_lowPass(minimumWavelength=minWavelength,sigma_wavelength=sigma_wavelength)
        specGrid.applyFilters()

        #Get the inverse transform of the topograph
        pgrid = specGrid.inverseTransform(returnAsArray=True)

        #Rescale the data so that everything varies from 0 to 1
        # otherwise, areas with large mean elevation (but little relief at this wavelength) will have minor changes
        # in priority
        pgrid = (pgrid - np.min(pgrid))/(np.max(pgrid) - np.min(pgrid))

        return pgrid

    def _getPriority(self,row,col):
        '''

        :param row:
        :param col:
        :return:
        '''
        return self._priortyGrid[row,col]*self._randomState.rand(1)[0]

class proceduralGrid_chiBalancedRandomFilling(proceduralGrid_randomPriorityFilling):

    def __init__(self, gridToMatch: flowRoutingGrids, fillKsn: float,fillTheta:float,randomSeed=None,
                 matchMeanElevationAndRelief=False, outlets=None, gridMask=None):
        '''
        This is a procedural terrain generation routine that modifies the priority flood algorithm described in:

        Barnes, R., Lehman, C., and Mulla, D., 2014, Priority-flood: An optimal depression-filling and
         watershed-labeling algorithm for digital elevation models: Computers and Geosciences, v. 62, p. 117–127,
          doi:10.1016/j.cageo.2013.04.024.

        Specifically, the algorithm fills terrain upstream from outlets, with channels allowed to 'grow' upstream one


        :param gridToMatch:
            A baseGrid instance to serve as the template for the random grid we will create. Provides the
            shape information, georeferencing, and optionally the mean elevation and relief

        :param fillKsn:

        :param fillTheta:
            The concavity of the channel, this is forced to be negative (i.e., the negative of the absolute value
            is used)

        :param doScalePriorityByElevation:
            A boolean that specifies whether the priority with which a network will propogate upstream scales
            with the current elevation of that point in the network. If False (the default) the priority with
            which filling occurs is purely random, which can lead to large elevation discrepencies across internal
            drainage divides. If True, filling is still random but the random priority is scaled based on how far
            up the network a cell is, this acts to minimize the elevation gaps in internal divides by limiting the
            likelihood that large channels curve back around on themselves (as opposed to having long channels tend
            to orient themselves parallel to the catchment axis).

        :param randomSeed (optional):
            A seed for numpy.random.RandomState . Default is None.

        :param matchMeanElevationAndRelief (optional):
            Whether to rescale the topography generated in order to match the mean elevation and relief of the input
            grid.

        :param outlets:
            A numpy array of shape n, 2. Specifies the n row,column coordinates of outlets to search upstream from.
            Generated flow will be out to these coordinates.

        :param gridMask:
            A boolean mask that limits the extent of the generated grid to within True values. For example,
            if a list of one basin coordinate was provided this could specify the basin extent to generate a random
            grid within.

        '''

        # Store parameters that will be necessary for later calculations
        self.ksn = fillKsn
        self.theta = -np.abs(fillTheta)

        #Set the random state
        self._randomState = np.random.RandomState(randomSeed)

        # Get the attributes from this instance
        self._copyAttributesFromGrid(gridToMatch)

        # If there is not a mask specified, make the whole grid in bounds
        if isinstance(gridMask, baseGrid):
            gridMask = gridMask.grid == 1
        elif gridMask is None:
            gridMask = np.ones(gridToMatch.shape) == 1

        # Store the mask for this grid
        self._gridMask = gridMask

        # If there are not outlets specified, get all the outlets that drain to boundaries
        if outlets is None:
            outlets = gridToMatch.findBoundaryOutlets(gridMask)
        # If only a single outlet was specified, make it iterable
        elif isinstance(outlets, np.ndarray):  # Is this a numpy array
            if outlets.shape == (2,):
                outlets = np.array([outlets])
        elif ~isinstance(outlets, np.ndarray) & len(outlets) == 2:
            outlets = np.array([np.array(outlets)])

        # Assign the outlets
        self._outlets = outlets

        # Keep a record of what the starting area of the outlets is
        self._outletAreas = np.array([gridToMatch.areaGrid[outlets[i,0],outlets[i,1]] for i in range(len(outlets))])

        # Keep track of what the smallest slope in the grid should be
        self._aggSlope = self.ksn*np.max(self._outletAreas)**self.theta

        #Create a dictionary that maps a downstream area to the observed values of area upstream of that
        self._areaMappingDict = self.__getAreaUSAreaMapping(gridToMatch)

        # If requested, keep track of the original grids mean elevation and relief so that we can make sure that the
        # procedurally generated grids have the same attributes
        self._doMatchMeanElevAndRelief = matchMeanElevationAndRelief
        if matchMeanElevationAndRelief:
            self._matchMeanElev = gridToMatch.mean()
            self._matchRelief = gridToMatch.relief()

    def __getAreaUSAreaMapping(self,referenceFlowGrid: flowRoutingGrids):
        '''
        This populates a dictionary that makes it easy to map a drainage area to its possible upstream areas.

        :param referenceFlowGrid:
        :return:
        '''

        #Initiate a dict that will relate a drainage area to a list of all the drainage areas that drain to that area
        areaMappingDict = {}

        #Loop through each of the outlets
        for i in range(len(self._outlets)):

            #Preform a recursive upstream search from each outlet, adding to the area mapping dict along the way
            self.__searchUpstreamPopulateAreaDictionary(self._outlets[i,0],self._outlets[i,1],
                                                        areaMappingDict,referenceFlowGrid)

        return areaMappingDict

    def __searchUpstreamPopulateAreaDictionary(self, row: int,col: int, areaMappingDict: dict,flowGrid: flowRoutingGrids):
        '''

        :param row:
        :param col:
        :param upstreamCells:
        :param flowGrid:
        :return:
        '''

        # Find cells immediately upstream (e.g., neighbors draining to this cell)
        upstreamNeighbors = flowGrid._findUpstreamNeighbors(row,col)[0]

        thisArea = flowGrid.areaGrid[row,col]

        # If this key isn't in the dictionary yet, add it as a key with an empty list
        if not (thisArea in areaMappingDict):
            areaMappingDict[thisArea] = []

        if len(upstreamNeighbors) > 0:

            # #Create a list to store these neighbors
            # thisList = [flowGrid.areaGrid[neighbor[0],neighbor[1]] for neighbor in upstreamNeighbors]
            # areaMappingDict[thisArea].append(thisList)

            #For each upstream neighbor
            for neighbor in upstreamNeighbors:

                neighborArea = flowGrid.areaGrid[neighbor[0],neighbor[1]]

                ## This dict mapped an area to all upstream areas in that cell
                #Add this neighbor's area to the list of possible areas upstream of this cell - but don't add divides
                #As those should naturally emerge
                if neighborArea > self._dx*self._dy:
                    areaMappingDict[thisArea].append(neighborArea)

                #Continue the recusrive search
                self.__searchUpstreamPopulateAreaDictionary(neighbor[0], neighbor[1], areaMappingDict,flowGrid)


    def getPermutation(self):
        '''

        :param gridToMatch:
        :return:
        '''

        #Reset which areas have already been visited
        self._unusedAreas = np.ones(self.shape)*self._gridMask
        self.areaGrid = np.ones(self.shape)*self._dx*self._dx

        #Assign the outlet coordinates to be the drainage area of the actual catchment
        self.areaGrid[self._outlets[:,0],self._outlets[:,1]] = self._outletAreas

        self.areaGrid[~self._gridMask] = np.nan

        randArray = self._calculateRandomInstanceArray()

        return demGrid(randArray, dx=self._dx,xllcenter=self._xllcenter,
                       yllcenter=self._yllcenter,geotransform=self._geotransform,projection = self._projection)

    def _getPriority(self,row,col):
        '''
        Set the priority to be based on the elevation
        :param row:
        :param col:
        :return:
        '''

        # if self.areaGrid[row,col] == (self._dx*self._dy):
        #     priority = np.inf
        # else:
        priority = self._randomState.random(1)[0]*(1.0/self.ksn)*self.areaGrid[row,col]**-self.theta

        return priority

    def _getSlope(self,row,col):
        '''

        :param row:
        :param col:
        :return:
        '''

        return self.ksn*self.areaGrid[row,col]**self.theta

    def _fillGridElevations(self, neighborRows, neighborCols, dsRow, dsCol, dsDistances, isClosed, openQueue):

        # Shuffle the order of the rows and columns
        shflIdx = np.arange(len(neighborRows))
        self._randomState.shuffle(shflIdx)
        neighborRows = neighborRows[shflIdx]
        neighborCols = neighborCols[shflIdx]
        dxs = dsDistances[shflIdx]

        # First, set the area grid for this cell to be a random choice of the areas upstream of cells with the same
        # area as the ds cell here
        thisdsArea = self.areaGrid[dsRow, dsCol]

        # How much drainage area is remaining?
        remainingArea = 1.0 * thisdsArea

        # How many neighbors are open
        openNeighbors = np.sum(~isClosed)

        # How many neighbors have we visited
        visitedNeighbors = 0

        # What are the available areas?
        theseAreas = np.array(self._areaMappingDict[thisdsArea])

        # Look through the neighbors
        for i in range(len(neighborCols)):

            # If this cell hasn't already been identified as closed
            if not isClosed[neighborRows[i], neighborCols[i]]:

                # Increment the visited neighbors
                visitedNeighbors += 1
                # If this is the last neighbor to visit
                if visitedNeighbors == openNeighbors:
                    # Assign it the remaining area (restricted to values actually observed before)
                    self.areaGrid[neighborRows[i], neighborCols[i]] = remainingArea
                elif not (thisdsArea in self._areaMappingDict) or (len(theseAreas) == 0) \
                        or (self._dx * self._dy * (openNeighbors - visitedNeighbors) >= remainingArea):
                    self.areaGrid[neighborRows[i], neighborCols[i]] = self._dx * self._dy
                else:
                    self.areaGrid[neighborRows[i], neighborCols[i]] = self._randomState.choice(theseAreas)

                remainingArea -= self.areaGrid[neighborRows[i], neighborCols[i]]

                # Now update the elevation of this cell based on the slope expected for a grid cell with this area
                self.grid[neighborRows[i], neighborCols[i]] = self.grid[dsRow, dsCol] + self._getSlope(
                                                                neighborRows[i],neighborCols[i]) * dsDistances[i]

                ### Set priority to be random, rather than a function of elevation
                priority = self._getPriority(neighborRows[i], neighborCols[i])

                ''' How do we ensure that boundary cells remain boundaries (e.g. low drainage area)? 
                If this cell is on the boundary (i.e., its neighbors contain out of bounds values), then
                give it a very high priority, ensure that its visited, but making it unlikely to be an
                intermediatte step in the flow.'''
                neighneighRows, neighneighCols, nnDX = self.getNeighborIndices(neighborRows[i], neighborCols[i])
                if np.any(~self._gridMask[neighneighRows, neighneighCols]):
                    priority = np.inf

                # So long as the priority wasn't nan, which specifies the termination of a flow path
                if ~np.isnan(priority):
                    # Add the item back to the queue with a random priority
                    openQueue.put(priority, [neighborRows[i], neighborCols[i]])

                # Mark the cell as closed
                isClosed[neighborRows[i], neighborCols[i]] = True