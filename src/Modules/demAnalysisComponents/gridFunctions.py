'''
Functions that perhaps don't below within a grid, but could just be generally accesible... currently I'm not sure
what I would put in here. Functions for loading / merging? Those should probably be their own thing.
'''

import numpy as np
from scipy.optimize import minimize


def predictTwoPartSplitPowerlaw(X,params):
    '''
    Predicts the y values of a two-segment powerlaw function. Intended to be used to 
    predict the slope-area relationship (or similar) of data that includes hillslope and channel segments.
    
    X: The independent variable (e.g., drainage area)
    params: 4-item array of paramaters in the split powerlaw [X_split, m_L, b_L, m_r]. Where X_split
    is the location where the two power-law segments are joined and transition. m_L is the slope (exponent)
    of the power-law fit to the fraction of data where X <= X_split, and b_L is the pre-exponential term to that data.
    m_R is the exponent to the data where X > X_split.
    
    
    return: Y_pred. The prediction of the 
    
    '''
    X_split, m_L,b_L,m_R = params
    Y_pred = np.zeros_like(X)*np.nan
    isLeft = X <= X_split #Left part of the function
    Y_pred[isLeft] = b_L*X[isLeft]**m_L
    Y_pred[~isLeft] = (b_L*X_split**m_L)*(X[~isLeft]/X_split)**m_R
    
    return Y_pred

def getBestFitParamsSplitPowerlaw(X,Y,X_split_guess = None, doLogTransformResiduals = True):
    '''
    Find the best fitting parameters that define a two-segment power-law.
    
    '''
    
    #Initial guesses for split, other terms
    if X_split_guess is None:
        X_split_guess = 10**np.mean(np.log10(X))
    
    #Fit powerlaws to left and right segments for initial guesses
    p_hs = np.polyfit(np.log10(X[(X < X_split_guess)]),np.log10(Y[(X < X_split_guess)]),1)
    p_chn = np.polyfit(np.log10(X[X > X_split_guess]),np.log10(Y[X > X_split_guess]),1)

    paramGuess = [X_split_guess, p_hs[0], 10**p_hs[1],p_chn[0]]
    print(paramGuess)

    if doLogTransformResiduals:
        objFun = lambda params: np.sum((np.log10(predictTwoPartSplitPowerlaw(X,params)) - np.log10(Y))**2)
    else:
        objFun = lambda params: np.sum((predictTwoPartSplitPowerlaw(X,params) - Y)**2)
    
    bf_params = minimize(objFun,paramGuess,method = 'nelder-mead').x
    
    return bf_params