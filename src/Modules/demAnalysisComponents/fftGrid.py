
import numpy as np
from matplotlib import pyplot as plt

from scipy import optimize

#Import local related components
try:
    import demAnalysisComponents.randomGrids as randGrid
    from demAnalysisComponents.baseGrid import baseGrid
    from demAnalysisComponents.dem import demGrid
except:
    import randomGrids as randGrid
    from baseGrid import baseGrid
    from dem import demGrid


class fftGrid(baseGrid):

    def __init__(self, demGridInstance : baseGrid, applyHannWindow = False, detrendGrid = True):
        '''
        TODO: Need to check on best way to handle repeat applicaton of filters
        TODO: Hmmm... do I actually want this to be a subclass of demGrid? What function is it actually inheriting that is of use, what is it inheriting that is not?
        :param demGridInstance:
        '''

        #Get the grid dimensions from the demGridInstance
        self._nrows, self._ncols = demGridInstance._nrows, demGridInstance._ncols
        self.shape = (self._nrows, self._ncols)

        #Get the grid reference coordinates
        self._xllcenter = demGridInstance._xllcenter
        self._yllcenter = demGridInstance._yllcenter

        #Get the coordinate size
        self._dx = demGridInstance._dx
        self._dy = demGridInstance._dy

        #Get the projection and transformation infromation from the provided grid
        self._geotransform = demGridInstance._geotransform
        self._projection = demGridInstance._projection

        #Get the sampling frequencies
        self.fsx = 1.0/demGridInstance._dx
        self.fsy = 1.0/demGridInstance._dy

        #Get the radial fourier and radial coordinates
        self._getFFT_coords()
        self._calcL2d()

        #Detrend the grid, if requested
        if detrendGrid:
            self._detrendGrid = demGridInstance.calcFittedPlane(returnAsArray=True)
        else:
            self._detrendGrid = np.zeros(self.shape)

        self._appliedHannWindow = applyHannWindow
        if applyHannWindow:
            self._windowMultiplier = self._getHannWindowNormalizer()
        else:
            #Apply window to grid to taper values at edges
            self._windowMultiplier = np.ones_like(self._detrendGrid)

        #Preform the transform
        self.grid = np.fft.fft2((demGridInstance.grid - self._detrendGrid)*self._windowMultiplier)

        #Calculate the spectral power
        self.spec_power = np.abs(self.grid)**2/(self._ncols*self._nrows*np.sum(self._windowMultiplier.flatten()**2))

        #Prallocate a list to save the filters in. TODO: Should this be a dictionary so that we save the names of the filters?
        self.filters = []

    def _getHannWindowNormalizer(self):
        '''
        Implementation of Hann window, from Perron et al., 2008, Press et al., 2007
        :return: W, the window to apply to taper grid boundaries to 0
        '''

        #Get midpoints of grid (max value of filter)
        a = (self._nrows -1)/2
        b = (self._ncols - 1)/2

        #Get grid coordinates
        n,m = np.meshgrid(np.arange(self._ncols),np.arange(self._nrows))

        #Radial angular coordinate from grid center
        theta = np.arctan((n-b)/(m-a))

        #Radial distance coordinate from grid center
        r = np.sqrt((m-a)**2 + (n-b)**2)

        #Normalizing radial distance
        rp = np.sqrt((a ** 2 * b ** 2) / (b ** 2 * np.cos(theta) ** 2 + a ** 2 * np.sin(theta) ** 2))

        #Preallocate space for window
        w = np.zeros_like(r)

        #Compute non-zero components of window
        nonZeroIdx = r<=rp
        w[nonZeroIdx] = 0.5 * (1. + np.cos(np.pi * r[nonZeroIdx] / rp[nonZeroIdx]))

        return w

    def applyFilters(self):
        '''

        :return:
        '''


        totalFilt = np.ones_like(self.grid)
        for f in self.filters:
            totalFilt*=f

        #Calculate filter normalization
        powerNorm = (self._ncols*self._nrows)*np.sum(totalFilt.flatten()**2)

        #Apply the filter
        self.grid*=totalFilt

        #Recalculate the power
        # self.spec_power = np.abs(self.grid)**2/powerNorm


    def addFilter_LinearDiffusion(self, kt : float):
        '''

        :param kt:  Morphologic age (m^2) to diffuse the topography for
        :return:
        '''

        #Calculate the diffusion filter
        diffFilter = np.exp((-4.0 * (np.pi ** 2) * kt * (self.L ** -2)))

        #Append this filter to the list of filters to aply
        self.filters.append(diffFilter)

    def inverseTransform(self,returnAsArray = False):
        '''

        :return: demGridInstance - the back transformed grid, after any fourier transform coefficients have been applied
        '''

        #Preform the inverse transform and retrend the grid
        demGridInstance = (np.real(np.fft.ifft2(self.grid))/self._windowMultiplier) + self._detrendGrid
        if not(returnAsArray):
            #Create an instance of demGrid with this data
            demGridInstance = demGrid(demGridInstance,geotransform=self._geotransform, projection= self._projection)

        return demGridInstance

    def _getFFT_coords(self):
        '''

        :return:
        '''

        self._dx = self.fsx / self._ncols  # size of frequency bins in x directions
        self._dy = self.fsy / self._nrows  # size of frequency bins in y direction

        # # Determine frequency bins
        # Fx = np.arange(self._ncols) * self._dx
        # Fy = np.arange(self._nrows)* self._dy
        #
        # # Reorder frequency bins so that values larger than fs/2 (Nyquist frequency)
        # # are negative
        # Fx[Fx >= self.fsx / 2.0] = Fx[Fx >= self.fsx / 2.0] - self.fsx
        # Fy[Fy >= self.fsy / 2.0] = Fy[Fy >= self.fsy / 2.0] - self.fsy

        # TODO: Double check that the flipping of the y is correct (think so due to image ordering?)
        self._Fxcoords = np.fft.fftfreq(self._ncols,1.0/self.fsx)
        self._Fycoords = np.fft.ifftshift(np.fft.fftshift(np.fft.fftfreq(self._nrows,1.0/self.fsy))[::-1])

    def addFilter_lowPass(self, minimumWavelength, sigma_wavelength = 0):
        '''
        preform a low pass filtering of the grid, removing all data from wavelengths less
        than the cutOffWavelength
        :param maximumWavelength: the maximum wavelength to allow to pass through to the filtered product
        :param sigma_wavelength: Optional, if specified smooths the edge of the filter with a gaussian (e.g., Perron et al., 2008)
        :return: Nothing, preforms filtering in place
        '''

        if sigma_wavelength == 0:
            lpFilter = 1.0*(self.L > minimumWavelength)
        else:
            lpFilter = np.ones_like(self.grid)
            lpFilter[self.L > minimumWavelength] = 1
            lpFilter[self.L <= minimumWavelength] = np.exp(-(self.L[self.L <= minimumWavelength] - minimumWavelength)**2/(2.0*sigma_wavelength**2))

        self.filters.append(lpFilter)

    def addFilter_highPass(self, maximumWavelength, sigma_wavelength = 0):
        '''
         preform a high pass filtering of the grid, removing all data from wavelengths greater
        than the minimumWavelength
        :param minimumWavelength: the minimum wavelength to allow to pass through to the filtered product
        :param sigma_wavelength: Optional, if specified smooths the edge of the filter with a gaussian (e.g., Perron et al., 2008)
        :return: Nothing, preforms filtering in place
        '''

        if sigma_wavelength == 0:
            hpFilter = 1.0 * (self.L < maximumWavelength)
        else:
            hpFilter = np.ones_like(self.grid)
            hpFilter[self.L < maximumWavelength] = 1
            hpFilter[self.L >= maximumWavelength] = np.exp(
                -(self.L[self.L >= maximumWavelength] - maximumWavelength) ** 2 / (2.0 * sigma_wavelength ** 2))

        self.filters.append(hpFilter)

    def addFilter_bandPass(self,minimumWavelength : float, maximumWavelength : float, sigma_wavelength : float = 0):
        '''
        preform a band pass filtering on the data, removing all data from wavelengths greater
        than the minimumWavelength
        :param minimumWavelength: the minimum wavelength to allow to pass through to the filtered product
        :param maximumWavelength: the maximum wavelength to allow to pass through to the filtered product
        :param sigma_wavelength: Optional, if specified smooths the edge of the filter with a gaussian (e.g., Perron et al., 2008)
        :return: Nothing, preforms filtering in place
        '''

        if sigma_wavelength == 0:
            bpFilter = 1.0 * ((self.L < maximumWavelength) & (self.L >minimumWavelength))
        else:
            meanWavelength = (minimumWavelength + maximumWavelength)/2.0
            bpFilter =  np.exp(-(self.L - meanWavelength) ** 2 / (2.0 * sigma_wavelength ** 2))

        self.filters.append(bpFilter)

    def _calcL2d(self):
        FX, FY = np.meshgrid(self._Fxcoords, self._Fycoords)
        self.L = 1.0 / np.sqrt(FX*FX + FY*FY)

    def getGridExtent(self):
        '''Return the bounding extent of the grid

        :return (minX, maxX, minY, maxY)'''

        #TODO: Double check that the new retrun based on array positions is correct
        # return (np.min(self._Fxcoords), np.max(self._Fxcoords), np.min(self._Fycoords),np.max(self._Fycoords))
        return (np.fft.fftshift(self._Fxcoords)[0], np.fft.fftshift(self._Fxcoords)[-1],
                np.fft.fftshift(self._Fycoords)[-1], np.fft.fftshift(self._Fycoords)[0])

    def calcWavelengthBinnedMedianPower(self, nBinsForMedian = 20, percRange = [2.5, 97.5]):
        '''

        :param nBinsForMedian: How many log bined median powers do we want to calculate? Divides these evenly (in log space)
        amongst available wavelengths
        :param percRange: will also return the upper and lower percentile bounds within each bin, this should be
        a list of two values to define those upper and lower limits. These are returned expressed as deviations from
        the median (e.g., abs(p_50 - p_2.5)), this allows them to be used with plt.errorbar most easily.
        :return: L_hat: np.ndarray (1 x n) of midpoint wavelengths for each bin
                p_hat: np.ndarray (1 x n) of median spectral power within each bin
                p_perc_err: np.ndarray(2 x n) of deviation from median of requested percentile bounds
        '''

        L_bins = np.logspace(np.log10(np.nanmin(self.L[~np.isinf(self.L)])),
                             np.log10(np.nanmax(self.L[~np.isinf(self.L)])), nBinsForMedian + 1)
        L_hat = (L_bins[1:] + L_bins[:-1]) / 2.0
        p_hat = np.zeros_like(L_hat) * np.nan
        p_perc = np.zeros((2, len(L_hat))) * np.nan

        for i in range(len(L_hat)):
            theseIndices = (self.L >= L_bins[i]) & (self.L < L_bins[i + 1])
            data_i = self.spec_power[theseIndices].flatten()
            if len(data_i) > 3:
                p_hat[i] = np.nanmedian(data_i)
                p_perc[:, i] = np.percentile(data_i, percRange)

        return L_hat, p_hat, p_perc

    def plotGrid(self,axs = None, logTransform = False, **kwargs):
        '''This is a wrapper for matplotlib plt.imshow that plots this grid

         :return axs (the axis that this was plotted on)'''

        if axs is None:
            f,axs = plt.subplots(1,1)

        if logTransform:
            axs.imshow(np.log10(np.fft.fftshift(self.spec_power)), extent=self.getGridExtent(), **kwargs)
        else:
            axs.imshow(np.fft.fftshift(self.spec_power), extent=self.getGridExtent(), **kwargs)

        return axs

    def calcPiecewiseRegressionCoefficients(self):
        '''

        :return:
        '''

        #Take an initial guess at the parameters of the regression
        #by assuming a middle scaling break and a single fit to all data
        goodData = ~np.isnan(self.L) & ~np.isnan(self.spec_power) & ~np.isinf(self.L) & ~np.isinf(self.spec_power)
        LtoFit = self.L[goodData].flatten()
        powerToFit = self.spec_power[goodData].flatten()

        p_guess = np.polyfit(np.log10(LtoFit), np.log10(powerToFit),1)
        coef_guesses = [10**p_guess[1], p_guess[0], np.mean(LtoFit),p_guess[0]]

        #Create a function to minimize
        objFun = lambda coeffs: np.sum((np.log10(powerToFit) -
                                        np.log10(self.calcPiecewiseRegressionPredictions(LtoFit,coeffs)))**2)

        #Optimize to find the minimum L2norm
        return optimize.minimize(objFun,coef_guesses,method = 'Nelder-Mead').x

    def calcPiecewiseRegressionPredictions(self,Ls, coeffs):
        '''

        :param Ls: The wavelengths to calculate the piecewise, powerlaw regression at
        :param coeffs: the coefficients of the regression [A_l, b_l, L_scalingBreak, b_r]
        :return:
        '''
        A_l,b_l,L_scalingBreak,b_r = coeffs
        predPower = np.zeros_like(Ls)
        leftOfScalingBreak = Ls <= L_scalingBreak
        powerAtScalingBreak = A_l*L_scalingBreak**b_l
        predPower[leftOfScalingBreak] = A_l*Ls[leftOfScalingBreak]**b_l
        predPower[~leftOfScalingBreak] = powerAtScalingBreak * (Ls[~leftOfScalingBreak]/L_scalingBreak)**b_r

        return predPower

    def calcWavelengthScalingBreak(self, piecewiseRegressionCoefficients = None):
        '''
        Get the wavelength of the scaling break observed in spectral power. Determined as the point where
        the two power-law fits to spectral power as a function of wavelength join.
        :return: L_scaling break: the wavelgnth where spectral power vs wavelength plots change slope
        '''

        if piecewiseRegressionCoefficients is None:
            piecewiseRegressionCoefficients = self.calcPiecewiseRegressionCoefficients()

        return piecewiseRegressionCoefficients[2]

    def plotMedianPowerWavelength(self, axs = None, nBinsForMedian = 20, percRange = [2.5, 97.5], **kwargs):
        '''

        :param axs:
        :param nBinsForMedian:
        :param kwargs:
        :return:
        '''


        if axs is None:
            f,axs = plt.subplots(1,1)

        L_hat, p_hat, p_perc = self.calcWavelengthBinnedMedianPower(nBinsForMedian = nBinsForMedian, percRange = percRange)

        #Convert the percentiles of the data in the bins to deviations from the central value
        p_perc_err = np.zeros_like(p_perc)
        for i in range(len(p_hat)):
            p_perc_err[:,i] = np.abs(p_hat[i] - p_perc[:,i])

        axs.errorbar(L_hat,p_hat,yerr = p_perc_err,fmt = 'o',**kwargs)
        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.set_xlabel('Wavelength (m)')
        axs.set_ylabel(r'DFT mean-squared amplitude (m^2)')

        return axs

    def plotPowerWavelength(self,axs = None,plotFractionOfObservations = 1.0,**kwargs):
        '''
        TODO: Update this so that it subsampling keeps fractions of data in length bins?
        :param axs:
        :param plotFractionOfObservations:
        :param kwargs:
        :return:
        '''

        #Chop ranges into bin to plot a subset of observations
        samplingBins = 20 #We'll divide everything into a number of bins, and do our sub-sampling from that so show data more uniformly

        if axs is None:
            f,axs = plt.subplots(1,1)

        if (plotFractionOfObservations > 0) & (plotFractionOfObservations < 1):
            #Figure out about how many observations would occur in each bin
            nObservationsPerBin = (self._nrows*self._ncols*plotFractionOfObservations) / samplingBins

            L_bins = np.logspace(np.log10(np.nanmin(self.L[~np.isinf(self.L)])),
                                 np.log10(np.nanmax(self.L[~np.isinf(self.L)])), samplingBins + 1)

            #Store data in lists
            allLs = []
            allPwrs = []

            #Loop through each of the bins
            for i in range(len(L_bins)-1):
                theseIndices = (self.L >= L_bins[i]) & (self.L < L_bins[i + 1])

                #Pull out these indices
                L_i = self.L[theseIndices].flatten()
                P_i = self.spec_power[theseIndices].flatten()

                #If there are more than the desired operations per bin, subsample
                if len(L_i) > nObservationsPerBin:
                    #Get a random subset of these indices
                    plotIndices = np.arange(len(L_i)) #Create an index for all observations
                    np.random.shuffle(plotIndices) #Shuffle the index in place
                    plotIndices = plotIndices[:int(nObservationsPerBin)]
                    L_i = L_i[plotIndices]
                    P_i = P_i[plotIndices]

                allLs.append(L_i)
                allPwrs.append(P_i)

            LtoPlot = np.hstack(allLs)
            PtoPlot = np.hstack(allPwrs)
            axs.plot(LtoPlot,PtoPlot,'.',**kwargs)

        elif plotFractionOfObservations == 1:
            axs.plot(self.L.flatten(),self.spec_power.flatten(),'.',**kwargs)
        else:
            print('Whoops, plotFractionOfObservations must be between 0 and 1, not plotting anything...')

        axs.set_xscale('log')
        axs.set_yscale('log')
        axs.set_xlabel('Wavelength (m)')
        axs.set_ylabel(r'DFT mean-squared amplitude (m^2)')

        return axs

    def saveGrid(self):
        '''

        :return:
        '''

        print('Whoops, saving fourier transform grids is not implemented.')
        print('Try completing your back transformation and then saving.')
        return None

    def calcDiamondSquareComparisonAnomalyGrid(self,roughness = None,nPermutations :int = 100, randomSeed :int= None):
        '''
        Calculates a suite of  procedurally generated, random grids, to use as the null hypothesis for what spectral
        power distributions might look like. This method is based on Perron et al., 2008 who use it to identify
        'anomolously' high-power wavelengths (and directions) that identify key elements of the landscape structure.

        :param roughness: (0-1) the 'H' value of the diamond square algorithm, controls how rough the random topo is
        :param nPermutations: how many permutations of random grids do we want to create to generate? More the better!
        :param randomSeed: default None, provide a random seed for the generation of procedural grids
        :return: self.permutedGrids, array with rows and columns of size this grid, and pages inclue to nPermutations
        '''

        #Get the back-transformed version of the grid in order to know the scale to create random grids at
        backTransformedGrid = self.inverseTransform()

        # If the Hann window was applied, inverting the transform will result in a lot of nans (as many hann window
        # Values will be 0)
        if self._appliedHannWindow:
            # When the hann window is removed, infs result, complicating the computation of the starting noise level
            startingScale = np.std(backTransformedGrid.grid[~np.isinf(backTransformedGrid.grid)
                                                            & ~np.isnan(backTransformedGrid.grid)])
        else:
            startingScale = np.std(backTransformedGrid.grid)

        #If the roughness wasn't specified, fit it
        if roughness is None:
            roughness = self.calcBestFittingDiamondSquareRoughness(startingScale=startingScale)


        #Create an instance of a diamond-square grid generator with the appropriate properties
        dsq = randGrid.proceduralGrid_diamondSquare(backTransformedGrid, roughness=roughness,
                                                    startingNoiseScale=startingScale,
                                                    matchMeanElevationAndRelief=False,
                                                    randomSeed=randomSeed)

        #Preallocate and save the grids of random power
        self._randomPowerGrids = np.zeros((self._nrows,self._ncols,nPermutations))

        #Iterate through each of the permutations, and create a grid for that
        for i in range(nPermutations):
            grid_i = dsq.getPermutation()
            fftRI = fftGrid(grid_i,applyHannWindow=True,detrendGrid=True)
            self._randomPowerGrids[:,:,i] = fftRI.spec_power

        return self._randomPowerGrids

    def calcBestFittingDiamondSquareRoughness(self, nBinsForMedian = 20, scalingBreakWavelength = None,
                                              nRoughnessValuesToTest = 20, startingScale = None):
        '''

        :param nBinsForMedian:
        :param scalingBreakWavelength:
        :param nRoughnessValuesToTest:
        :param startingScale:
        :return:
        '''

        if scalingBreakWavelength is None:
            scalingBreakWavelength = self.calcWavelengthScalingBreak()

        #Get the binned median wavelengths and powers
        L_obs, P_obs = self.calcWavelengthBinnedMedianPower(nBinsForMedian)[:-1]

        #Find those values above the wavelength scaling break
        aboveScalingBreak = L_obs >= scalingBreakWavelength

        #Instead of optimizing, maybe just do a brute force search over a specified range?
        HvalsToTest = np.linspace(0,1,nRoughnessValuesToTest)
        L2Norms = np.zeros_like(HvalsToTest)

        backTransformedGrid = self.inverseTransform()

        if startingScale is None:
            #If the Hann window was applied, inverting the transform will result in a lot of nans (as many hann window
            #Values will be 0)
            if self._appliedHannWindow:

                #When the hann window is removed, infs result, complicating the computation of the starting noise level
                startingScale = np.std(backTransformedGrid.grid[~np.isinf(backTransformedGrid.grid)
                                       & ~np.isnan(backTransformedGrid.grid)])
            else:
                startingScale = np.std(backTransformedGrid.grid)

        for i, H_i in enumerate(HvalsToTest):
            # Create an object for creating diamond square values with this roughness value, matching this grids
            dsq = randGrid.proceduralGrid_diamondSquare(backTransformedGrid, roughness=H_i,
                                                        startingNoiseScale=startingScale,
                                                        matchMeanElevationAndRelief=False)

            # Get a single diamond square permutation from that grid
            randInst = dsq.getPermutation()

            # Preform the fft transform on that grid - no need to detrend, as should already be zero mean (aprox.)
            fftRI = fftGrid(randInst, applyHannWindow=True, detrendGrid=False)

            # Get the binned observations for this random grid
            L_pred, P_pred = fftRI.calcWavelengthBinnedMedianPower(nBinsForMedian)[:-1]

            # Make sure that the data we are looking at here are ok to fit.
            goodData = ~np.isnan(P_pred) & ~np.isnan(P_obs) & aboveScalingBreak

            L2Norms[i] = np.sum((np.log10(P_pred[goodData]) - np.log10(P_obs[goodData])) ** 2)

        return HvalsToTest[np.argmin(L2Norms)]


    def calcPermutationTestAnomalyGrid(self, nPermutations = 100):

        '''
        TODO: Implement this.
        :return:
        '''

        #Is this going to be way to memory intensive of an operation for this?
        self._randomPowerGrids = np.zeros((self._nrows,self._ncols,nPermutations))

        nBins = 300
        L_bins = np.logspace(np.log10(np.nanmin(self.L[~np.isinf(self.L)])),
                             np.log10(np.nanmax(self.L[~np.isinf(self.L)])), nBins + 1)

        #Create a dictionary to map indices of different Lengths in the grid to their length
        LmappingGrid = np.zeros_like(self.L)
        for i,L in enumerate(L_bins[:-1]):
            theseIndices = (self.L >= L_bins[i]) & (self.L < L_bins[i + 1])
            LmappingGrid[theseIndices] = i #Doing like this so that could do 'bins' of minimal variation

        #Construct each permutation as a length-dependent resampling of available values

        for i in range(nPermutations):
            dummyGrid = np.zeros_like(self.grid)  #Create a dummy grid for each iteration to preallocate memory creation
            #For each radial wavelength
            for j,L_i in enumerate(L_bins[:-1]):
                #Find all of the values of the grid at that value
                theseValues = self.spec_power[LmappingGrid==j]

                #Add these powers back in to the grid, resampling with replacement
                dummyGrid[LmappingGrid==j] = np.random.choice(theseValues, theseValues.shape,replace = True)

            self._randomPowerGrids[:,:,i] = dummyGrid #Hmmmm... is this going to make all values the same because dummyGrid continues to chage?

    def plotSpectralAnomalies(self,axs = None,logTransform = True,**kwargs):

        if axs is None:
            f, axs = plt.subplots(1, 1)

        if self._randomPowerGrids is None:
            self.calcPermutationTestAnomalyGrid()

        countGrid = np.zeros_like(self.spec_power)
        nPermutations = self._randomPowerGrids.shape[2]

        normGrid = np.median(self._randomPowerGrids,axis = 2)

        for i in range(nPermutations):
            countGrid+= 1.0*(self.spec_power > self._randomPowerGrids[:,:,0])

        countGrid=(countGrid/nPermutations)*100.0

        if logTransform:
            # axs.imshow(np.log10(np.fft.fftshift(countGrid)), extent=self.getGridExtent(), **kwargs)
            axs.imshow(np.log10(np.fft.fftshift(self.spec_power/normGrid)), extent=self.getGridExtent(), **kwargs)
        else:
            # axs.imshow(np.fft.fftshift(countGrid), extent=self.getGridExtent(), **kwargs)
            axs.imshow(np.fft.fftshift(self.spec_power/normGrid), extent=self.getGridExtent(), **kwargs)

        return axs, np.fft.fftshift(countGrid), normGrid

    def plotRadialSpectralAnomalies(self, axs = None, Ls = None, nThetaBins :int = 40,
                                    frequencySTDofSmoothingWindow: float = None, **kwargs):
        '''

        :param Ls:
        :param nThetaBins:
        :param frequencySTDofSmoothingWindow:
        :param kwargs:
        :return:
        '''

        #Below needs to be updated for this function
        if axs is None:
            f, axs = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))

        if self._randomPowerGrids is None:
            self.calcPermutationTestAnomalyGrid()

        normGrid = np.median(self._randomPowerGrids,axis = 2)
        anomalyGrid = self.spec_power/normGrid

        #Create L,theta grid as necessary
        thetas = np.linspace(0,2.0*np.pi,nThetaBins)

        if Ls is None:
            #Figure out the max length based on the length in the shortest direction
            lmax = np.min([self.shape[i]*d for i,d in enumerate([1./self.fsx, 1./self.fsy])])

            #Create the spacing to be linear in frequency space
            Ls = 1.0/np.linspace(1./lmax,np.max([self.fsx,self.fsy]),  int(np.min(self.shape)/2))

        #Get L theta as a grid
        self.Lgrid,self.Thetagrid = np.meshgrid(Ls,thetas)

        #Interpolate the anomaly grid onto this radial coordinate system
        self._length_theta_anomalyGrid = self._transformXYgridToLengthTheta(anomalyGrid,self.Lgrid,self.Thetagrid,
                                                                            frequencySTDofSmoothingWindow)

        #Plot the plot
        axs.pcolor(self.Thetagrid, self.Lgrid,self._length_theta_anomalyGrid, **kwargs)

        return axs

    def _transformXYgridToLengthTheta(self,gridToInterpolate: np.ndarray, L: np.ndarray,Theta: np.ndarray, frequencySTDofSmoothingWindow: float = None):
        '''
        Interpolates the provided grid ('gridToInterpolate'), which is assumed to correspond to this transform
        instances frequency coordinates (e.g., should be of the same shape as this grid), onto the radial coordinate
        system provided by L, theta (2d numpy arrays). Does this with a gaussian weighted mean with a standard deviation
        specified by 'frequencySTDofSmoothingWindow', which operates on the data in frequency space.

        :param gridToInterpolate: 2d numpy array of the same size as this grid, to interpolate to L,Theta space
        :param L: 2d numpy array of the desired L coordinates (e.g., derived from np.meshgrid)
        :param Theta: 2d numpy array of the desired radial coordinates (e.g., derived from np.meshgrid)
        :param frequencySTDofSmoothingWindow: standard deviation of gaussian weighted mean, in frequency coordinates,
        if None is specified, will default to 2* minimum of the sampling frequency
        :return:
        '''

        if frequencySTDofSmoothingWindow is None:
            sig = 1.0*np.min([self._dx, self._dx])
        else:
            sig = frequencySTDofSmoothingWindow

        #Convert the requested length, theta bins into x,y frequency bins
        freq = 1.0/L
        X_p,Y_p = freq*np.cos(Theta), freq*np.sin(Theta)

        #Get the meshgridded coordinates for this grid
        X_o, Y_o = np.meshgrid(self._Fxcoords, self._Fycoords) #Todo: need to confirm whether this should really be negative

        #Iterate through the length, theta grid - calculate values based on gaussian kernel
        length_theta_grid = np.zeros_like(X_p)

        #Loop through rows and columns of length theta spectral power grid
        for i in range(X_p.shape[0]):
            for j in range(Y_p.shape[1]):
                x_i, y_i = X_p[i,j], Y_p[i,j]
                weights = np.exp(-((X_o - x_i)**2/(2.0*sig**2) + (Y_o - y_i)**2/(2.0*sig**2)))
                length_theta_grid[i,j] = np.sum(weights*gridToInterpolate)/np.sum(weights)

        return length_theta_grid