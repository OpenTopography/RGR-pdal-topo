from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np
from scipy.optimize import fmin


try:
    from demAnalysisComponents.baseGrid import baseGrid
    from demAnalysisComponents.flowRoutingGrids import flowRoutingGrids
    from demAnalysisComponents.stablePriorityQueue import priorityQueue
except:
    from baseGrid import baseGrid
    from flowRoutingGrids import flowRoutingGrids
    from stablePriorityQueue import priorityQueue

PARAMDICT = {'drainage_area': 'A',
             'elevation': 'Z',
             'upstream_distance': 'L',
             'x_coordinate': 'x',
             'y_coordinate': 'y',
             'chi': 'Chi',
             'row_coordinate':'row',
             'col_coordinate':'col',
             'stream_order': 'order',
             'topological_order': 'topol_order',
             'has_operation_visited':'_hasVisited',
             'segment_id': 'id',
             'windowed_slope': 'S_win',
             'restored_elevation':'Zrestored',
             'downstream_direction': 'ds_direction',
             'downstream_line_intercept': 'ds_lineInt',
             'segment_length': 'seg_len',
             'downstream_slope': 'ds_slope',
             'channel_steepness': 'Ksn',
             'colinearity_score':'colinearityScore',
             'function_value':'funVal',
             'windowed_DzDchi':'DzDchi'}

'''Hmmmm... I have pigeon holed the creation of network graphs in to just be possible 
from a flow routing grid. '''

class networkNode:

    def __init__(self, row:int, col:int, L:float,parentNode =  None,
                 flowRoutingGrid: flowRoutingGrids = None,
                 Z:float = None, A:float = None,
                 x:float = None, y:float = None):
        '''
        This is a class that just stores the basic units of a linked-list drainage network.
        :param flwGrid:
        :param row:
        :param col:
        :param L:
        '''


        self.children = []  # List of all the children of this node
        self.parent = parentNode # The parent of this coordinate. This is needed if we are going to move down the network

        #List the basic position coordinates
        self.row = row
        self.col = col
        self.L = L

        #ToDO: this has some vvvveerrrryy strange behavior.... that I don't fully understand:
        '''
        I load a flowRoutingGridInstance from a file, I can check its instance and it is a flowRoutingGrids. 
        
        I build a networkGraph from it, and the first if below returns false, so all attributes have 'None' assigned.
        
        I reload the networkGraph class and things run correctly...
        
        This seems to be related to my import statements...
        
        I ran
        import dem
        import flowRoutingGrids as flw
        import networkGraph as net
        from importlib import reload
        reload(flw)
        
        commenting out the reload allowed things to work correctly....
        perhaps there is something with how reload treats these classes behind the scenes?
        '''
        if isinstance(flowRoutingGrid,flowRoutingGrids):
            self.A = flowRoutingGrid.areaGrid[row, col]
            self.Z = flowRoutingGrid.grid[row, col]
            self.x,self.y = flowRoutingGrid.getXYFromRowCol([row,col])
        else:
            self.A = A
            self.Z = Z
            self.x = x
            self.y = y

    def copyNode(self, parentNode = None, children = None, doCopyAttributes = True):
        '''

        :return:
        '''

        nodeCopy = networkNode(row = self.row,col = self.col,L = self.L, parentNode=parentNode)

        if not(children is None):
            if isinstance(children,list):
                if isinstance(children[0],networkNode):
                    nodeCopy.children = children

        #If we are copying available attributes
        if doCopyAttributes:

            #Check each of the possible node attributes
            for key in PARAMDICT.keys():
                #Get the variable code for this attribute
                attr = PARAMDICT[key]

                #If this instance of a node has this attribute, add it to the copies node
                if hasattr(self,attr):
                    setattr(nodeCopy,attr,getattr(self,attr))

        return nodeCopy


class networkGraph:
    ''' Implementation of a channel network as a linked list
    '''


    def __init__(self, flwGrid:flowRoutingGrids, outletRow:int,outletCol:int,
                 channelMask: baseGrid = None, Amin:float = 1e5, doUseRecursion = False):

        L_0 = 0.0  # Starting from outlet

        if channelMask is None:
            channelMask = flwGrid.areaGrid > Amin

        #Keep a record of the top of the channel network
        self._networkTails = []

        #Populate the channel network
        self._networkHead = self.__linkToUpstream__(flwGrid, outletRow, outletCol, L_0, channelMask,
                                                    doUseRecursion=doUseRecursion)

        #Set some default parameters
        self._isChiCalcd = False

        #Default to recursive searches or not
        self._doUseRecursion = doUseRecursion


    def __linkToUpstream__(self,flwGrid:flowRoutingGrids, row: int, col:int, L_0: float, channelMask: np.ndarray,
                           parentNode = None, doUseRecursion = None):

        if doUseRecursion is None:
            doUseRecursion = self._doUseRecursion

        if doUseRecursion:
            node = self.__recursiveLinkToUpstream__(flwGrid, row, col, L_0,channelMask)

        else:
            node = self.__nonRecursiveLinkToUpstream__(flwGrid,row,col,L_0,channelMask)#Implement this

        return node


    def __nonRecursiveLinkToUpstream__(self, flwGrid:flowRoutingGrids, row: int, col: int, L_0: float, channelMask: np.ndarray,
                           parentNode = None):
        ''' Queue-based upstream search of network
        '''

        #TODO: test

        #Create the coordinates for this node
        networkHead = networkNode(row,col,L_0,parentNode=parentNode,
                               flowRoutingGrid = flwGrid)

        #Create a queue with the active network nodes
        nodeQueue = [networkHead]


        while nodeQueue:

            #Get the current node
            thisNode = nodeQueue.pop(0)

            #Get the current distance upstream
            L_0 = thisNode.L

            # Find the cells that are upstream of this cell
            upstreamNeighbors, dxs = flwGrid._findUpstreamNeighbors(thisNode.row, thisNode.col)

            #Keep a count of how many neighbors are upstream of this
            neighborCount = 0

            #For each upstream neighbor
            for i, neighbor in enumerate(upstreamNeighbors):

                #If this neighbor was determined to be part of the network
                if channelMask[neighbor[0], neighbor[1]]:
                    #Get the upstream node
                    child_i = networkNode(neighbor[0], neighbor[1],L_0 + dxs[i],parentNode=thisNode,
                               flowRoutingGrid = flwGrid)

                    #Add it to the current nodes children
                    thisNode.children.append(child_i)

                    #Tally the number of children
                    neighborCount+=1

                    #Add it to the queue
                    nodeQueue.append(child_i)

            #If there were no upstream children, this is a termination of the network graph
            if neighborCount == 0:
                self._networkTails.append(thisNode)

        return networkHead

    def __recursiveLinkToUpstream__(self, flwGrid:flowRoutingGrids, row: int, col: int, L_0: float, channelMask: np.ndarray,
                           parentNode = None):
        ''' Recursive upstream search of network
        '''

        #Create the coordinates for this node
        thisNode = networkNode(row,col,L_0,parentNode=parentNode,
                               flowRoutingGrid = flwGrid)

        #Find the cells that are upstream of this cell
        upstreamNeighbors, dxs = flwGrid._findUpstreamNeighbors(row,col)

        #Keep a count of how many neighbors are upstream of this
        neighborCount = 0

        #For each upstream neighbor
        for i, neighbor in enumerate(upstreamNeighbors):

            #If this neighbor was determined to be part of the network
            if channelMask[neighbor[0], neighbor[1]]:
                thisNode.children.append(
                    self.__linkToUpstream__(flwGrid, neighbor[0], neighbor[1], L_0 + dxs[i],channelMask))

                neighborCount+=1

        #If there were no upstream neighbors, this is a termination of the network graph
        if neighborCount == 0:
            self._networkTails.append(thisNode)

        return thisNode

    def __upstreamSetter__(self,node,function, doUseRecursion = None):
        #Re-implement to have recursive and non-resursive solutions

        if doUseRecursion is None:
            doUseRecursion = self._doUseRecursion

        if doUseRecursion:
            self.__recursiveUpstreamSetter__(node, function)
        else:
            self.__nonRecursiveUpstreamSetter__(node, function)

    def __upstreamGetter__(self,node,attribute,output,doUseRecursion = None):
        #Re-implement to have recusrive and non-recusive solutions

        if doUseRecursion is None:
            doUseRecursion = self._doUseRecursion

        if doUseRecursion:
            self.__recursiveUpstreamGetter__(node,attribute,output)
        else:
            self.__nonRecursiveUpstreamGetter__(node,attribute,output)


    def __nonRecursiveUpstreamGetter__(self,node,attribute,output):

        nodesInQueue = [node]

        while nodesInQueue:

            # Get the first outlet that was added
            node_i = nodesInQueue.pop(0)

            output.append(getattr(node_i, attribute))

            # For each of the upstream nodes
            for i, child in enumerate(node_i.children):
                    # Add this node back into the search list
                    nodesInQueue.append(child)

    def __recursiveUpstreamGetter__(self, node, attribute, output):

        output.append(getattr(node, attribute))

        for child in node.children:
            self.__recursiveUpstreamGetter__(child, attribute, output)


    def __nonRecursiveUpstreamSetter__(self,node,function):

        nodesInQueue = [node]

        while nodesInQueue:

            # Get the first outlet that was added
            node_i = nodesInQueue.pop(0)

            # For each of the upstream nodes
            for i, child in enumerate(node_i.children):

                #Apply the function
                function(node_i,child)

                # Add this node back into the search list
                nodesInQueue.append(child)

    def __recursiveUpstreamSetter__(self, node, function):
        '''Function takes a parent node and daughter node and assigns some value,
        for example:

        '''
        for child in node.children:
            function(node, child)
            self.__recursiveUpstreamSetter__(child, function)


    def __downstreamSetter__(self,node:networkNode, function,doUseRecursion = None):

        if doUseRecursion is None:
            doUseRecursion = self._doUseRecursion

        if doUseRecursion:
            self.__recursiveDownstreamSetter__(node, function)
        else:
            self.__nonRecursiveDownstreamSetter__(node, function)


    def __nonRecursiveDownstreamSetter__(self, node: networkNode, function):

        nodesInQueue = [node]

        while nodesInQueue:

            # Get the first outlet that was added
            node_i = nodesInQueue.pop(0)

            if not (node_i.parent is None):

                #Apply the functino
                function(node_i.parent, node)

                # Add this node back into the search list
                nodesInQueue.append(node_i.parent)

    def __recursiveDownstreamSetter__(self, node: networkNode, function):
        '''
        Function that takes a daughter node, and assigns some value to is parent
        :param node:
        :param function:
        :return:
        '''

        if not(node.parent is None):
            function(node.parent,node)
            self.__recursiveDownstreamSetter__(node.parent,function)

    def calcChiValues(self, theta=0.5, A_0=1.0e6):
        ''' Perform the upstream integration of drainage area described by Perron and Royden

        '''
        self._isChiCalcd = True
        self._A_0 = A_0
        self._theta = theta

        # Assign the inititial integrated drainage area to the outlet
        self._networkHead.Chi = 0

        #Create a lambda function with this parameter set
        chiFunction = lambda parent, child: self.__singleNodeChi__(parent, child, A_0, theta)

        #Preform the upstream recursive search for setting Chi
        self.__upstreamSetter__(self._networkHead, chiFunction)

    def calcColinearityScore(self, locations: np.ndarray, orientations: np.ndarray,
                                                   orientation_sigma: float, distance_sigma: float):

        self._networkHead.colinearityScore = np.nan

        colinFun = lambda parent,child: self.__singleNodeIdentifyColinearityScore__(parent, child,
                                                   locations, orientations,
                                                   orientation_sigma, distance_sigma)

        self.__upstreamSetter__(self._networkHead,colinFun)

    def calcChildAttributeValue(self,attributeName, function):

        origAttributeName = attributeName
        if attributeName in PARAMDICT:
            attributeName = PARAMDICT[attributeName]


        if hasattr(self._networkHead,attributeName):
            recFun = lambda p,chld: setattr(chld, PARAMDICT['function_value'], function(getattr(chld, attributeName)))
            recFun([],self._networkHead)
            self.__upstreamSetter__(self._networkHead,recFun)
        else:
            raise ValueError(
                'No value named {} could be found in this network. You may not have calculated this yet.'.format(
                    origAttributeName))

    def calcSegmentDirections(self):
        '''
        TODO: Test this

        :return:
        '''

        self._networkHead.ds_direction = np.nan

        self.__upstreamSetter__(self._networkHead, self.__singleNodeSegmentDirection__)

    def calculateAttributeDifference(self,attributeName1, attributeName2,newName = 'AttrDifference',
                                     doNormalizeDifference = True):
        '''

        :param attributeName1:
        :param attributeName2:
        :param newName:
        :param doNormalizeDifference:
        :return:
        '''

        if doNormalizeDifference:
            recFun = lambda p,c: setattr(c,newName,(getattr(c,attributeName1) - getattr(c,attributeName2)
                                  )/(getattr(c,attributeName1) + getattr(c,attributeName2)))
        else:
            recFun = lambda p,c: setattr(c,newName,getattr(c,attributeName1) - getattr(c,attributeName2))

        #Apply this to network head
        recFun([],self._networkHead)

        #Apply recursively to children
        self.__upstreamSetter__(self._networkHead,recFun)

    def calcSegmentLengths(self):

        '''

        :return:
        '''

        self._networkHead.seg_len = np.nan

        self.__upstreamSetter__(self._networkHead,self.__singleNodeSegmentLength__)

    def calcSegmentIntercepts(self, interceptCoordinate = 0):
        '''
        TODO: Test this

        :return:
        '''

        self._networkHead.ds_lineInt = np.nan

        intFun = lambda parent,child : self.__singleNodeSegmentIntercept__(parent,child,interceptCoordinate)

        self.__upstreamSetter__(self._networkHead, intFun)

    def calcChannelSlopes(self):
        self._networkHead.ds_slope = np.nan

        self.__upstreamSetter__(self._networkHead,self.__singleNodeChannelSlope__)

    def calcChannelSteepness(self, theta = 0.5):
        '''

        :param theta:
        :return:
        '''
        self._networkHead.Ksn = np.nan

        setFunction = lambda parent,child: self.__singleNodeChannelSteepness__(parent,child,theta)

        self.__upstreamSetter__(self._networkHead, setFunction)

    def restoreChannelProfileUpstream(self, Ksn: float, theta: float):
        ''' Integrate a steady state slope (specified by channel steepness)
        upstream
        '''

        self._networkHead.Zrestored = self._networkHead.Z

        restoreFunction = lambda parent, child: self.__singleNodeIntegrateSlopeUpstream__(parent, child, Ksn,theta)

        self.__upstreamSetter__(self._networkHead, restoreFunction)

    def restoreChannelProfileDownstream(self, Ksn: float, theta: float):
        ''' Integrate a steady state slope (specified by channel steepness)
        downstream.

        TODO: Test this function
        '''
        #Problem: This will overwrite subsequent values in strange ways depending on the order that the recursion proceeds
        print('Hmmm... should reimplement this as a queue...')
        restoreFunction = lambda parent,child: self.__singleNodeIntegrateSlopeDownstream__(parent,child, Ksn, theta)

        for tail in self._networkTails:
            tail.Zrestored = tail.Z
            self.__downstreamSetter__(tail, restoreFunction)

    def calcUpstreamTopologicOrder(self,initialOrder = 0):



        t_order = initialOrder
        setattr(self._networkHead,PARAMDICT['topological_order'],t_order)
        openQueue = [self._networkHead]

        while openQueue:
            node_i = openQueue.pop(0)
            children = node_i.children
            if not(children is None):
                order_i = getattr(node_i, PARAMDICT['topological_order'])
                if len(children) > 1:
                    order_i +=1

                for c in children:
                    setattr(c,PARAMDICT['topological_order'],order_i)
                    openQueue.append(c)


    def calcDownstreamChannelOrder(self,initialOrder = 0):
        '''
        TODO: BLARGH.... Still not quite getting this....
        :param initialOrder:
        :type initialOrder:
        :return:
        :rtype:
        '''

        #First march upstream and add things to a queue based on topological (e.g., farthest up topology = first visit)
        openQueue = priorityQueue()
        openQueue.put(0,self._networkHead)
        nodesToVisit = [(0,self._networkHead)]


        while nodesToVisit:
            priority,node_i = nodesToVisit.pop(0)
            children = node_i.children
            if not (children is None):
                order_i = priority - 1
                for c in children:
                    openQueue.put(order_i,c)
                    nodesToVisit.append((order_i,c))

        while not openQueue.isEmpty():
            priority, child = openQueue.get()

            childOrder = getattr(child, PARAMDICT['stream_order'], None)

            if childOrder is None:
                childOrder = initialOrder
                setattr(child, PARAMDICT['stream_order'], childOrder)

            if not (child.parent is None):

                parent = child.parent

                # Could define a more efficient algorithm probably...
                parentOrder = getattr(parent, PARAMDICT['stream_order'], None)


                if parentOrder is None:
                    thisOrder = childOrder
                elif parentOrder == childOrder:
                    thisOrder = childOrder + 1
                elif parentOrder < childOrder:
                    thisOrder = childOrder
                elif parentOrder > childOrder:
                    thisOrder = parentOrder

                setattr(parent, PARAMDICT['stream_order'], thisOrder)


    def getLineSegments(self, additionalParameters = None):
        '''
        TODO: Test this function with and without additional parameters, refactor non-recursive versions of
        functions to get line segments
        :param additionalParameters:
        :return:
        '''
        lineSegments = []

        if not(additionalParameters is None):
            parameterDictionary = {}
            nParams = 0
            if not(isinstance(additionalParameters,list)):
                additionalParameters = [additionalParameters]

            #Iterate through the requested parameters, make sure they exist, and initialize a list to store values
            for p in additionalParameters:
                #If this is a plain text name representation, add the actual parameter value as a dicitonary key
                if p in PARAMDICT:
                    parameterDictionary[PARAMDICT[p]] = []
                    nParams+=1
                elif hasattr(self._networkHead,p):
                    parameterDictionary[p] = []
                    nParams+=1
                else:
                    raise ValueError('Warning: could not find any values with name {}'.format(str(p)))

            #If there werent any valid parameters, get rid of the parameter dictionary
            if nParams ==0:
                additionalParameters = None

        #If extra parameters were requested, search accordingly
        if additionalParameters is None:
            self.__recursiveUpstreamGetLineSegments__(self._networkHead, lineSegments)
            toReturn = lineSegments
        else:
            self.__recursiveUpstreamGetLineSegmentsAndExtraParams__(self._networkHead,lineSegments,parameterDictionary)
            toReturn = (lineSegments, parameterDictionary)

        return toReturn

    def lineSegmentsToGDF(self, EPSGCode:int, additionalParameters:list, filePath: str = None,
                           geoPandasSaveDriver:str = 'ESRI Shapefile'):
        '''

        :param filePath:
        :param EPSGCode:
        :param additionalParameters:
        :param fieldNames:
        :return:
        '''
        from shapely.geometry import LineString
        import geopandas as gpd

        lineSegments, parameterDictionary = self.getLineSegments(additionalParameters)

        lineGeoms = [LineString(ls) for ls in lineSegments]

        gdf = gpd.GeoDataFrame(parameterDictionary, geometry=lineGeoms,crs="EPSG:{}".format(EPSGCode))

        if not(filePath is None):
            gdf.to_file(filePath,driver=geoPandasSaveDriver)

        return gdf

    def __recursiveUpstreamGetLineSegments__(self,node,lineSegments):
        '''

        :param node:
        :param lineSegments:
        :return:
        '''

        xStart = node.x
        yStart = node.y

        for child in node.children:
            xEnd = child.x
            yEnd = child.y

            lineSegments.append(np.array([[xStart,yStart],[xEnd,yEnd]]))

            self.__recursiveUpstreamGetLineSegments__(child,lineSegments)

    def __recursiveUpstreamGetLineSegmentsAndExtraParams__(self,node,lineSegments, parameterDictionary):
        '''

        :param node:
        :param lineSegments:
        :return:
        '''

        xStart = node.x
        yStart = node.y

        for child in node.children:
            xEnd = child.x
            yEnd = child.y

            lineSegments.append(np.array([[xStart,yStart],[xEnd,yEnd]]))
            for key in parameterDictionary:
                #TODO: could have a more elegant version where I aggregated up and downstream in some way?
                parameterDictionary[key].append(getattr(child,key,np.nan))

            self.__recursiveUpstreamGetLineSegmentsAndExtraParams__(child,lineSegments, parameterDictionary)

    def assignNodesClosestValues(self,x: np.ndarray,y: np.ndarray,values: np.ndarray ,fieldname: str,
                                 maxDist:float = None, maxDistVal = np.nan):
        '''

        :param x: 1-D Array of x coordinates
        :param y: 1-D Array of y coordinates
        :param values: 1-D array of values
        :param fieldname: name of values to assign, this string is used as the attribute name for the assignment
        of a value to each node
        :param maxDist: maximum allowed distance. Points farther away than this distance will be assigned the value
        specified by maxDistVal. If None, no maximum distance is specified
        :param: maxDistVal: value assigned to network nodes whose nearest value is greater than maxDist away
        :return: None, Adds a two new parameters to each node in the graph- 'fieldname' - the string supplied to specify
        this parameter, and 'fieldname'+'_dist' - an attribute specifying the distance to that field
        '''

        #Flatten the function to one of two inputs
        recFun = lambda p,c : self.__asignClosestValueSingleNode__(p,c,x,y,values,fieldname,maxDist,maxDistVal)

        #Assign values to the network head
        recFun([],self._networkHead)

        #Recurse up through the children and apply this functino
        self.__upstreamSetter__(self._networkHead,recFun)

    def __asignClosestValueSingleNode__(self, parentNode: networkNode, child: networkNode,
                                  x: np.ndarray,y: np.ndarray,values: np.ndarray ,fieldname: str,
                                        maxDist:float = None, maxDistVal = np.nan ):
        '''

        :param parentNode:
        :param child:
        :param x: 1-D Array of x coordinates
        :param y: 1-D Array of y coordinates
        :param values: 1-D array of values
        :param fieldname: name of values to assign, this string is used as the attribute name for the assignment
        of a value to each node
        :param maxDist: maximum allowed distance. Points farther away than this distance will be assigned the value
        specified by maxDistVal. If None, no maximum distance is specified
        :param: maxDistVal: value assigned to network nodes whose nearest value is greater than maxDist away
        :return: None, Adds a two new parameters to this child node - 'fieldname' - the string supplied to specify
        this parameter, and 'fieldname'+'_dist' - an attribute specifying the distance to that field
        '''

        dists = np.sqrt((x-child.x)**2 + (y-child.y)**2)
        minDistIdx = np.argmin(dists)
        value = values[minDistIdx]
        dist = dists[minDistIdx]

        setattr(child,fieldname+'_dist',dist)
        setattr(child, fieldname, value)
        if not(maxDist is None):
            if dist > maxDist:
                setattr(child,fieldname,maxDistVal)


    def flattenNetwork(self, attribute):
        '''

        :param attribute:
        :return:
        '''

        #If the plain-text version of the value was specified, get the variable name version
        if attribute.lower() in PARAMDICT:
            attribute = PARAMDICT[attribute.lower()]

        output = []

        self.__upstreamGetter__(self._networkHead, attribute, output)

        return np.array(output)

    def __singleNodeIdentifyTail__(self,parentNode:networkNode, child: networkNode):
        '''

        :param parentNode:
        :param child:
        :return:
        '''

        if len(child.children)==0:
            self._networkTails.append(child)


    def __singleNodeIdentifyColinearityScore__(self, parent: networkNode, child: networkNode,
                                                locations: np.ndarray, orientations: np.ndarray,
                                                orientation_sigma: float, distance_sigma:float):
        ''' Hmmm....

        #TODO:
        As of 06/12 I had the idea to do something like a weighted dot product, segments along strike
        got high weights, each point was has its dot product calculated. WRT a unit vector version of
        the current node.

        I theoretically like this approach - as it would seem to provide information both about what is nearly
        along strike to a network and what is in line with that network, but it doesn't seem like it is working
        correctly. I suspect that part of the problem is my distance function, which I know to be a bit incorrect
        (should do projected distance). But there seems to be more of a problem than that...

        Another thing I would like to try is just color coding vectors based on orientations, or a probability
        distribution derived from fracture orientations (could do color and transparency with probability of
        being within a fracture orientation and class of fracture orientation)?

        Update 6/16: Spent some more time experimenting with this, and am not making progress. Should consider
        revisiting to explore a different approach.

        The two things of interest would still seem to be: where are the channel segments along strike of this one,
        and which of these is oriented in a similar direction. Perhaps I am just over interpreting the test image
        that I have, but the approach(es) I have below don't seem to do a good job of defining this. One problem
        is that the distance weighting may not allow enough uncertainty in the orientations of network segments.
        Another proble is perhaps that the dot product (e.g., cos(relAngle) may handle orientations well enough.


        '''


        if parent.x == child.x:
            child.colinearityScore = np.nan
        else:
            #Fitting a line and looking at overlap
            # m = (parent.y - child.y) / (parent.x - child.x)
            # int = child.y - m * child.x
            # Spacing between these nodes
            # yPreds =  m*locations[:,0] +int
            # dists = np.sqrt((yPreds - locations[:,1])**2)
            #
            # angle_dist = np.sqrt((np.arctan(m) - orientations)**2)

            # child.colinearityScore = np.sum((dists < distance_sigma) & (angle_dist < orientation_sigma))
            # child.colinearityScore = 1.0/np.sum(dists)


            #Orientation of this line segments
            theta = np.arctan2(parent.y - child.y,parent.x - child.x)

            # Orientation of vector from this point to others
            phi = np.arctan2(locations[:,1] - child.y, locations[:,0] - child.x)

            #Length of vector to other points
            H = np.sqrt((locations[:,0] - child.x)**2 +(locations[:,1] - child.y)**2)

            # yPreds = H * np.sin(theta) + child.y
            # xPreds = H * np.cos(theta) + child.x

            yPreds = H * np.sin(np.arctan(np.tan(theta))) + child.y
            xPreds = H * np.cos(np.arctan(np.tan(theta))) + child.x

            yPreds1 = H * np.sin(np.arctan(np.tan(theta))) + child.y
            xPreds1 = H * np.cos(np.arctan(np.tan(theta))) + child.x

            xObs = H*np.sin(np.arctan(np.tan(phi))) + child.y
            yObs = H*np.cos(np.arctan(np.tan(phi))) + child.x

            #TODO: Redo this to actual projected distance
            if parent.x == child.x:
                dists = locations[:,0]
            elif parent.y == child.y:
                dists = locations[:,1]
            else:
                m = (parent.y - child.y) / (parent.x - child.x)
                int = child.y - m * child.x
                a = -(parent.y - child.y)
                b = (parent.x - child.x)
                c =  -b*int

                dists = np.abs(a*locations[:,0] + b*locations[:,1] + c)/np.sqrt(a**2 + b**2)

            dists = np.abs(H*np.cos((np.pi/2) - np.arctan(np.tan(phi)) + np.arctan(np.tan(theta))))
            # dists = np.sqrt((xPreds - locations[:, 0]) ** 2 + (yPreds - locations[:, 1]) ** 2)
            # dists = np.sqrt((child.x - locations[:, 0]) ** 2 + (child.y - locations[:, 1]) ** 2)
            # dists = np.sqrt((xPreds - xObs) ** 2 + (yPreds - yObs) ** 2)
            # distance_sigma = np.abs(H*np.sin(orientation_sigma))

            # dists, distance_sigma = np.sqrt((np.arctan(np.tan(phi)) - np.arctan(np.tan(theta)))**2), orientation_sigma

            # weights = (1.0/np.sqrt(2.0*np.pi*distance_sigma**2))*np.exp(-0.5*(dists/distance_sigma)**2)
            weights = np.exp(-0.5*(dists/distance_sigma)**2)
            weights[np.isnan(weights)] = 0

            # print('Min: {:.2e}, Max: {:.2e}'.format(np.min(dists),np.max(dists)))
            # Get the angle between the lines, removing that atan2 directions are 4 quadrant
            relativeAngle = np.abs(np.arctan(np.tan(orientations)) - np.arctan(np.tan(theta)))
            relativeAngle[relativeAngle > np.pi/2] = np.pi - relativeAngle[relativeAngle > np.pi/2]
            child.colinearityScore = np.sum(weights*np.exp(-relativeAngle/(2.0*orientation_sigma**2)))
            # child.colinearityScore = np.sum(weights * np.cos(relativeAngle))
            # child.colinearityScore = np.sum(weights*(relativeAngle < orientation_sigma))



    def __singleNodeIntegrateSlopeUpstream__(self, parent, child, Ksn, theta):
        ''' Integrate the slopes upstream
        '''

        # Spacing between these nodes
        dx = child.L - parent.L

        # Slope between these nodes, using upstream area
        S = Ksn * child.A ** -theta

        child.Zrestored = parent.Zrestored + S * dx

    def __singleNodeSegmentLength__(self,parent:networkNode, child:networkNode):
        '''

        :param parent:
        :param child:
        :return:
        '''

        child.seg_len = child.L - parent.L

    def __singleNodeSegmentDirection__(self, parent: networkNode, child: networkNode):
        '''

        :param parent:
        :param child:
        :return:
        '''

        dx = parent.x - child.x
        dy = parent.y - child.y

        child.ds_direction = np.arctan2(dy,dx)

    def __singleNodeSegmentIntercept__(self, parent: networkNode, child: networkNode, interceptCoordinate: float):
        '''

        :param parent:
        :param child:
        :return:
        '''

        #TODO: This try/except catches infinite slopes (e.g., parent.x = child.x), is there a better way to handle this?
        try:
            m = (parent.y - child.y)/(parent.x - child.x)
            child.ds_lineInt = child.y - m*(child.x - interceptCoordinate)
        except:
            child.ds_lineInt = np.inf

    def __singleNodeChannelSlope__(self,parent: networkNode, child: networkNode):
        '''

        :param parent:
        :param child:
        :return:
        '''

        child.ds_slope = (parent.Z - child.Z)/(parent.L - child.L)

    def __singleNodeChannelSteepness__(self, parent: networkNode, child: networkNode, theta: float):
        '''

        :param parent:
        :param child:
        :param theta:
        :return:
        '''
        S = (parent.Z - child.Z) / (parent.L - child.L)

        child.Ksn = S*child.A**theta

    def __singleNodeIntegrateSlopeDownstream__(self, parent,child, Ksn, theta):
        ''' Integrate the slopes upstream
        '''

        # Spacing between these nodes
        dx = child.L - parent.L

        # Slope between these nodes, using upstream area
        S = Ksn * child.A ** -theta

        # Downstream integration of slope
        parent.Zrestored = child.Zrestored - S * dx

    def __singleNodeChi__(self, parent, child, A_0, theta):
        '''Sum the integrated drainage area from the parent to the child
        '''
        dx = child.L - parent.L  # Difference in along profile length
        child.Chi = parent.Chi + dx * (A_0 / child.A) ** theta

    def plotMapRepresentation(self,axs = None,**kwargs):
        '''

        :param axs:
        :param kwargs:
        :return:
        '''

        if axs is None:
            f,axs = plt.subplots(1,1)

        lineSegments = self.getLineSegments()

        for seg in lineSegments:
            axs.plot(seg[:,0],seg[:,1],'-',**kwargs)

        return axs

    def plotColorizedMapRepresentation(self, axs=None,colorizeParameter = 'drainage_area', colormapname = 'Blues',
                                       vmin = None, vmax = None, doLogTransformValues = False,
                                       doTransformToPercentile = False,**kwargs):
        '''

        :param axs:
        :param colorizeParameter:
        :param colormapname:
        :param kwargs:
        :return:
        '''

        if axs is None:
            f, axs = plt.subplots(1, 1)

        if colorizeParameter in PARAMDICT:
            colorizeParameter = PARAMDICT[colorizeParameter]

        if hasattr(self._networkHead,colorizeParameter):
            lineSegments, additionalParamDict = self.getLineSegments([colorizeParameter])
            paramValues = np.array(additionalParamDict[colorizeParameter])

            # if doTransformToPercentile:
            #     percValues = np.zeros_like(paramValues)
            #     for i in range(len(percValues)):
            #         percValues[i] = np.sum(paramValues[i] > paramValues)/float(len(paramValues))
            #     paramValues = percValues
            #     vmin = 0
            #     vmax = 1

            if vmin is None:
                vmin = np.min(paramValues)

            if vmax is None:
                vmax = np.max(paramValues)


            if doLogTransformValues:
                paramValues = np.log10(paramValues)
                vmin = np.log10(vmin)
                vmax = np.log10(vmax)


            cmap = cm.get_cmap(colormapname)
            normedParamValues = (paramValues - vmin)/(vmax-vmin)
            colors = [cmap(v) for v in normedParamValues]

            for i,seg in enumerate(lineSegments):
                axs.plot(seg[:, 0], seg[:, 1], '-',color = colors[i], **kwargs)


        else:
            raise ValueError('No attribute named {} found in this instanct of networkGraph.'.format(colorizeParameter))

        return axs

    def plotValues(self,xValueName, yValueName, axs = None,  **kwargs):

        #If not plotting axis was specified
        if axs is None:
            f,axs = plt.subplots(1,1)

        axs.plot(self.flattenNetwork(xValueName), self.flattenNetwork(yValueName), '.', **kwargs)

        return axs

    def plotScatter(self,colorValue = 0,sizeValue = 14,axs = None,**kwargs):
        '''

        :param xValueName:
        :param yValueName:
        :param colorValueName:
        :param sizeValue:
        :param axs:
        :param kwargs:
        :return:
        '''

        if axs is None:
            f,axs = plt.subplots(1,1)

        if (colorValue.lower() in PARAMDICT):
            colorValue = self.flattenNetwork(colorValue)
        elif hasattr(self,colorValue):
            colorValue = self.flattenNetwork(colorValue)

        if (sizeValue.lower() in PARAMDICT):
            sizeValue = self.flattenNetwork(sizeValue)
        elif hasattr(self,sizeValue):
            sizeValue = self.flattenNetwork(sizeValue)

        axs.scatter(self.flattenNetwork('x_coordinate'),self.flattenNetwork('y_coordinate'),
                    s= sizeValue,c = colorValue,**kwargs)

        return axs


    def plotChiProfiles(self, axs = None, doAdjustOutletElevationToZero = True, A_0: float = None,
                        theta:float = None,**kwargs):

        #If values of A_0 and theta are provided, calculate Chi
        if not(A_0 is None) and not(theta is None):
            self.calcChiValues(theta = theta, A_0 = A_0)

        #If Chi has been calculated proceed
        if self._isChiCalcd:
            Chi = self.flattenNetwork('chi')
            Z = self.flattenNetwork('Z')

            #Adjust the plot so that the outlet elevation is zero (for comparing catchments)
            if doAdjustOutletElevationToZero:
                Z -= Z[0]

            #If no plotting axis was specified, create one
            if axs is None:
                f,axs = plt.subplots(1,1)

            axs.plot(Chi, Z, '.', **kwargs)
            axs.set_xlabel(r'$\chi$', fontsize=14)
            axs.set_ylabel(r'Elevation', fontsize=12)

        #If chi has not been calculated notify the user
        else:
            raise Exception('Whoopsy: must provide A_0 and theta parameters to calculate Chi during the call to '+
                            '.plotChiProfiles, or calculate Chi before hand by calling .calcChiValues.')

        return axs

    ################################################################################################################
    ######## Functions for creating derivatives from network graphs
    ################################################################################################################
    def getWindowedSlopeAreaArrays(self, winHeight):
        ''' return x,y,S,A the coordinates and slopes, areas measured in the specified vertical interval
        '''

        xs = []
        ys = []
        S = []
        A = []

        # Keep track of the x,y coordinates for the current segment
        theseX = [self._networkHead.x]
        theseY = [self._networkHead.y]
        theseA = [self._networkHead.A]

        self._searchUpstreamFromPointForSlopeArea(self._networkHead, self._networkHead.Z, self._networkHead.L, winHeight, xs, ys, S, A,
                                                  theseX, theseY, theseA)

        return np.array(xs), np.array(ys), np.array(S), np.array(A)

    def _searchUpstreamFromPointForSlopeArea(self, node, Z0, L0, winHeight, xs, ys, S, A, travX, travY, travA):
        ''' helper function for recursive upstreamsearch
        '''

        if (node.Z - Z0) >= winHeight:
            travX.append(node.x)
            travY.append(node.y)
            travA.append(node.A)

            # Add this point to list
            S.append((node.Z - Z0) / (node.L - L0))
            A.append(np.mean(np.array(travA)))

            midPoint = int(len(travX) / 2)
            xs.append(travX[midPoint])
            ys.append(travY[midPoint])

            # update the Z0 and coordinate list
            travX = [node.x]
            travY = [node.y]
            travA = [node.A]

            Z0 = node.Z
            L0 = node.L

        for child in node.children:
            # Copy the traverse for this trib
            newX = list(travX)
            newY = list(travY)
            newA = list(travA)

            # Extend the traverse
            newX.append(child.x)
            newY.append(child.y)
            newA.append(child.A)

            self._searchUpstreamFromPointForSlopeArea(child, Z0, L0, winHeight, xs, ys, S, A, newX, newY, newA)

    ##############################################################################################################
    ### Functions for dissolving networks
    ##############################################################################################################

    def dissolveNetwork_PreserveTributaryJunctions(self, doPreserveChannelHeads = True):

        dissFun= lambda dsNode,currNode,newParent : self.__disolveFunction_tributaryJunctions__(dsNode,currNode,
                                                                                                newParent,
                                                                                                doPreserveChannelHeads)

        newNetwork = dissolvedNetworkGraph(self,dissFun)

        return newNetwork

    def __disolveFunction_tributaryJunctions__(self,originalDSNode: networkNode, originalCurrentNode:networkNode,
                                             newDSParentNode: networkNode, doPreserveChannelHeads: bool):

        nChildren = len(originalCurrentNode.children)

        if nChildren >1 or ((nChildren == 0) and doPreserveChannelHeads):
            doAddToNetwork = True
            newNode = originalCurrentNode.copyNode(parentNode=newDSParentNode)
        else:
            doAddToNetwork = False
            newNode = None

        return doAddToNetwork, newNode

    def dissolveNetwork_PreserveTributaryJunctionsLimitLength(self,maxSegLength :float, doPreserveChannelHeads=True):

        dissFun = lambda dsNd, crNd, nwPrnt: self.__disolveFunction_segmentLengthAndJunctions__(dsNd, crNd,nwPrnt,
                                                                                                doPreserveChannelHeads,
                                                                                                maxSegLength)
        newNetwork = dissolvedNetworkGraph(self, dissFun)

        return newNetwork

    def __disolveFunction_segmentLengthAndJunctions__(self, originalDSNode: networkNode,
                                                      originalCurrentNode: networkNode,
                                                      newDSParentNode: networkNode,
                                                      doPreserveChannelHeads: bool, segmentLength: float):

        nChildren = len(originalCurrentNode.children)
        L_i = originalCurrentNode.L - originalDSNode.L

        if (L_i > segmentLength) or (nChildren >1) or ((nChildren == 0) and doPreserveChannelHeads):
            doAddToNetwork = True
            newNode = originalCurrentNode.copyNode(parentNode=newDSParentNode)
        else:
            doAddToNetwork = False
            newNode = None

        return doAddToNetwork, newNode

    def dissolveNetwork_calculateWindowedSlope(self,windowHeight: float):
        '''

        :param windowHeight:
        :return:
        '''

        dissFun = lambda dsNd, crNd, nwPrnt: self.__disolveFunction_windowedDownstreamSlope__(dsNd, crNd,nwPrnt,
                                                                                              windowHeight)
        newNetwork = dissolvedNetworkGraph(self, dissFun)

        # The head doesn't get set because the window applies calculations to the upstream node, set it now
        setattr(newNetwork._networkHead, PARAMDICT['windowed_slope'], np.nan)

        return newNetwork

    def __disolveFunction_windowedDownstreamSlope__(self, originalDSNode: networkNode,
                                                      originalCurrentNode: networkNode,
                                                      newDSParentNode: networkNode, windowHeight:float):

        nChildren = len(originalCurrentNode.children)
        Z_i = originalCurrentNode.Z - originalDSNode.Z

        if (Z_i > windowHeight) or (nChildren >1):
            doAddToNetwork = True
            newNode = originalCurrentNode.copyNode(parentNode=newDSParentNode)
            if (Z_i > windowHeight):
                setattr(newNode,PARAMDICT['windowed_slope'],Z_i/(originalCurrentNode.L - originalDSNode.L))
            else:
                setattr(newNode,PARAMDICT['windowed_slope'],np.nan)
        else:
            doAddToNetwork = False
            newNode = None

        return doAddToNetwork, newNode


    def dissolveNetwork_calculateWindowed_DzDChi(self,windowHeight: float):
        '''

        :param windowHeight:
        :return:
        '''

        if self._isChiCalcd:
            dissFun = lambda dsNd, crNd, nwPrnt: self.__disolveFunction_windowedDownstream_DzDChi__(dsNd, crNd,nwPrnt,
                                                                                                  windowHeight)
            newNetwork = dissolvedNetworkGraph(self, dissFun)

            #The head doesn't get set because the window applies calculations to the upstream node, set it now
            setattr(newNetwork._networkHead, PARAMDICT['windowed_DzDchi'], np.nan)

        else:
            print('Whoopsy, must calculate Chi first. Returning nothing.')
            newNetwork = None

        return newNetwork

    def __disolveFunction_windowedDownstream_DzDChi__(self, originalDSNode: networkNode,
                                                      originalCurrentNode: networkNode,
                                                      newDSParentNode: networkNode, windowHeight:float):

        nChildren = len(originalCurrentNode.children)
        Z_i = originalCurrentNode.Z - originalDSNode.Z

        if (Z_i > windowHeight) or (nChildren >1):
            doAddToNetwork = True
            newNode = originalCurrentNode.copyNode(parentNode=newDSParentNode)
            if (Z_i > windowHeight):
                setattr(newNode,PARAMDICT['windowed_DzDchi'],Z_i/(originalCurrentNode.Chi - originalDSNode.Chi))
            else:
                setattr(newNode,PARAMDICT['windowed_DzDchi'],np.nan)
        else:
            doAddToNetwork = False
            newNode = None

        return doAddToNetwork, newNode


class dissolvedNetworkGraph(networkGraph):
    #TODO: Perhaps a better way to do this would be to have a standard network graph have multiple 'init' methods?
    #Unless I want to keep some record of the initial grid here, it seems strange to have this added functionality.

    def __init__(self, origNetwork: networkGraph, networkDissolveFunction, doUseRecursion = None):
        '''

        :param origNetwork:
        :param networkDissolveFunction:
        '''

        if (doUseRecursion is None) or not isinstance(doUseRecursion,bool):
            doUseRecursion = origNetwork._doUseRecursion

        #Keep track of whether to use recursion
        self._doUseRecursion = doUseRecursion

        # Cope the head of the network from the master network
        self._networkHead = origNetwork._networkHead.copyNode()

        self.__searchUpstreamDissolveNetwork__(origNetwork, networkDissolveFunction, doUseRecursion)


        # Keep a record of the top of the channel network
        self._networkTails = []
        #Search upstream in the network and add the terminal nodes as the tails
        self.__upstreamSetter__(self._networkHead,self.__singleNodeIdentifyTail__)

        # Set some default parameters
        self._isChiCalcd = False


    def __searchUpstreamDissolveNetwork__(self,origNetwork: networkGraph, networkDissolveFunction, doUseRecursion):

        if doUseRecursion:
            self.__recursiveSearchUpstreamDissolveNetwork__(origNetwork._networkHead, origNetwork._networkHead,
                                                   self._networkHead, networkDissolveFunction)
        else:
            self.__nonRecursiveSearchUpstreamDissolveNetwork__(origNetwork._networkHead, origNetwork._networkHead,
                                                   self._networkHead, networkDissolveFunction)

    def __nonRecursiveSearchUpstreamDissolveNetwork__(self,origDSNode: networkNode,origCurrentNode: networkNode,
                                          dsNewNode: networkNode, networkDissolveFunction):



        ####HMMMM.... how exactly to do this.... do I need to maintain two queues? One w/ the active nodes being searched, the second with the new nodes

        newNodesInQueue = [(dsNewNode, origDSNode)]
        while newNodesInQueue:

            #Get the downstream node for this dissolved graph
            newDSNode_i, origDSNode_i = newNodesInQueue.pop(0)

            activeSearchQueue  = [origDSNode_i]

            #While there are nodes to search under this downstream node
            while activeSearchQueue:

                activeSearchNode = activeSearchQueue.pop(0)

                #For each child in the
                for child in activeSearchNode.children:

                    #Check if this is to be added to the network
                    doAddToNetwork,newNode = networkDissolveFunction(origDSNode_i,child,newDSNode_i)

                    #If we found a node in the dissolved network
                    if doAddToNetwork:
                        newDSNode_i.children.append(newNode)

                        #Add these nodes to the master queue
                        newNodesInQueue.append((newNode,child))

                    #If we didn't find a node in the dissolved network, continue looking down that chain of the graph
                    else:
                        activeSearchQueue.append(child)




    def __recursiveSearchUpstreamDissolveNetwork__(self,origDSNode: networkNode,origCurrentNode: networkNode,
                                          dsNewNode: networkNode, networkDissolveFunction):

        #Does the node pass the defined criteria?
        doAddToNetwork, newNode = networkDissolveFunction(origDSNode, origCurrentNode, dsNewNode)

        #If yes, update nodes accordingly
        if doAddToNetwork:
            dsNewNode.children.append(newNode)
            origDSNode =origCurrentNode
            dsNewNode = newNode

        # For each upstream neighbor
        for child in origCurrentNode.children:

            self.__recursiveSearchUpstreamDissolveNetwork__(
                origDSNode,child,dsNewNode,networkDissolveFunction)
