#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 11:56:26 2019

@author: Brice Loose

Python Implementation of the Backfit algorithm described by 
Wood, Simon. (2017). Generalized Additive Models: an introduction with R 
(CRC Press).

# USAGE:  ycal, pp, mod, alpha = GAMcal(X,y,lam)
#
#
# INPUT:  (works with numpy arrays and pandas dataframes, but must be consistent
#           across data types)
#   X = N x d, Environmental inputs or corelates; this variable should be standardized.
#   y = N x 1, instrument response
#   lam = penalty
# OUTPUT:
#   ycal = reconstructed series of y from the model. 
#   pp = polynomial fit to mass 18, not presently implemented.
#   mod = list() data type of Cubic spline fits to each column in X.
#   alpha = mean of y.  This should be zero if already standardized. 
#
# DISCLAIMER:
#    This is provided "as is" without warranty of any kind.  
#=========================================================================
"""
def GAMcal(X,y,lam):
    from sklearn.linear_model import Ridge
    import numpy as np
    from patsy import dmatrix
    
    
    ok = True; alpha = y.mean(); f = X*0; m = 2; rss0 = 0;  hdri = f.columns.values
    g = f.copy(); mod = [None]*len(hdri); Ri = [None]*len(hdri)
    
    mod = [None]*len(hdri); pred = [None]*len(hdri); Ri = [None]*len(hdri);
    
    while ok:
        for k in range(0,len(X.columns)):
            idk = (list(range(0,k,1))+list(range(k+1,len(X.columns),1)))
            ep = y - f[hdri[idk]].sum(axis=1) - alpha

            #if hdri[k] == u'MASS( 18 )':
            #    pp = np.polyfit(X[hdri[k]],ep,1)
            #    newy = np.polyval(pp,X[hdri[k]])
            #else:
                # Generate natural cubic spline basis with df+4 knots.
            bas = dmatrix("cr(train,df=10)-1", {"train": X[hdri[k]]})
            Ri[k] = Ridge(alpha=lam,normalize=True)
            mod[k] = Ri[k].fit(bas,ep)
            newy = mod[k].predict(dmatrix("cr(valid, df=10)-1", {"valid": X[hdri[k]]}, return_type='dataframe'))

            g[hdri[k]] = newy
            
            
            
            f = g
        rss = ((y-f.sum(axis=1))**2).sum()
       # print(rss)
        if abs(rss-rss0) < 1e-6*rss:
            ok = False
        rss0 = rss
    
    #Reconstruct is:
    ycal = alpha + f.sum(axis=1)
    
    if 'pp' not in locals():
        # myVar no exists.
        pp = np.array([1.,0.])
    
    return ycal, pp, mod, alpha

"""
Created on Wed May 29 11:56:26 2019

@author: Brice Loose

Python Implementation of the Backfit algorithm described by 
Wood, Simon. (2017). Generalized Additive Models: an introduction with R 
(CRC Press).

# USAGE:  yr, alpha = GAMapply(X,mod,pp,hdr,Y).
#
#
# INPUT:  (works with pandas dataframes only)
#   X = N x d dataframe, Environmental inputs or corelates; this variable 
#   should be standardized.
#   mod = list() data type of Cubic spline fits to each column in X. 
#   Y = N x 1, instrument response
#   pp = polynomial fit to mass 18, not presently implemented.
#   hdr = column names to subindex. 
#
# OUTPUT:
#   yr = reconstructed series of y from the model. 
#   alpha = mean of y.  This should be zero if already standardized. 
#
# DISCLAIMER:
#    This is provided "as is" without warranty of any kind.  
#=========================================================================
"""    

def GAMapply(X,b,pp,hdr,Y):
    # Version _rm tries to use a running mean from the data itself, instead of
    # alpha.
    # Need to add df option
    from patsy import dmatrix
    import numpy as np
    #%matplotlib qt

    Yp = X*0;

    for k in range(0,len(X.columns)):        
        #if hdr[k] == u'MASS( 18 )':
        #    newy = np.polyval(pp,X[hdr[k]])
        #else:
        newy = b[k].predict(dmatrix("cr(valid, df=10)-1", {"valid": X[hdr[k]]}, return_type='dataframe'))
        
        Yp[hdr[k]] = newy

    alpha = Y.rolling(window=500).mean()
    alpha.rename(columns={alpha.name: 'ybar'},inplace=True)
    yr = alpha + Yp.sum(axis=1)

    return yr, alpha

        
        
