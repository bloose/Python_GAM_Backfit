# Python_GAM_Backfit
A Python implementation of the Backfit algorithm for fitting to Generalized Additive Models

Python Implementation of the Backfit algorithm described by 
Wood, Simon. (2017). Generalized Additive Models: an introduction with R 
(CRC Press).

USAGE:  yr, alpha = GAMapply(X,mod,pp,hdr,Y).


 INPUT:  (works with pandas dataframes only)
   X = N x d dataframe, Environmental inputs or corelates; this variable 
   should be standardized.
   mod = list() data type of Cubic spline fits to each column in X. 
   Y = N x 1, instrument response
   pp = polynomial fit to mass 18, not presently implemented.
   hdr = column names to subindex. 

 OUTPUT:
   yr = reconstructed series of y from the model. 
   alpha = mean of y.  This should be zero if already standardized. 

 DISCLAIMER:
  Provided "as is" without warranty of any kind.  
