# poll4pop_python
Python code for poll4pop model


Adapted from the original R versions on https://github.com/yclough/poll4pop/tree/v1.1.1 
These functions are analagous to the R versions but are designed to operate on numpy array land cover / edge inputs rather than rasters. This allows the code to be generic to different GIS interfaces (ArcGIS, Q) but requires the user to pre-generate numpy arrays from rasters. E.g. see https://pro.arcgis.com/en/pro-app/latest/arcpy/functions/rastertonumpyarray-function.htm for relevant function. 

Users must currently also pre-generate python dictionaries and lists containing other input parameters in the format required for functions. 

Pre-requirements:
lc          a landcover array [n x m], integer type where each value corresponds to a specific land use code
lccodes		  a list of integer codes corresponding to all the land use codes (e.g [0,1,2,...,45,46])
edges       a list of edge floating point arrays, where each array [n x m] corresponds to a specific land use code and the values are the length of that type of edge feature
edgecodes   a list of integer codes corresponding to the land use codes of the edge arrays (e.g. [5,11,12])
AttDict     a dictionary [k: f1, f2, f3, N] of floral (3 consecutive seasons f1, f2, f3) and nesting attractiveness (N) parameters where the keys (k) are the land cover codes
FlorDict    a dictionary [k: f1, f2, f3] of floral cover (3 consecutive seasons f1, f2, f3 parameters where the keys (k) are the land cover codes
LfnDict     a dictionary [k: w] of edge widths (w) where the keys (k) are the land cover codes
species     list of species codes being considered as integers (e.g. [1,2,3,4])
gr_list     a list of arrays of growth parameters where rows are queens and workers and columns are a, b, max, and pw; one for each species
pw_list     a list of proportion of workers, one for each species 
av_list     a list of nest density parameters, one for each species
dist_list   a list of arrays of distance parameters where first row is foraging and second is nesting for each species
numpy
scipy


Files

CompFlorNest.py contains 
compFlorNest function: Takes lc, lccodes, edges, edgecodes, AttDict, FlorDict, LfnDict. Returns arrays of floral resource value (by season), nesting resource value, and proportion of value attributable to lc and edge features. 

KLP.py contains:
kernCalc function: Takes mean dispersal distance (from dist_list), cell size, maximum size of kernel, decay cutoff. Returns dispersal kernel (as array)
latForDisp function: Takes array of foragers, array of resources, dispersal kernel (from kernCalc). Returns distance weighted resources, visitation per cell, resources per forager as arrays. 
pollPopGrowth function: Takes resources, lognormal cdf parameter inputs and returns closing population as array

RunPoll_SolSoc_3S.py contains
set of functions which run the poll4pop model for three particular pollinator types: 
- solitary bees active in seasons 1 and 2
- solitary bees active in season 3 only 
- social bees active in seasons 1, 2 and 3 with workers generated at the end of season 1 and 2 that forage in seasons 2 and 3 respectively
The model is intended to be run for consecutive years until the population converges. For each pollinator type there is a function to initialise the first year and a function to re-apply for consecutive years until convergences is achieved. 

For more detail see R package documentation (e.g. runpoll_3seasons) and the following papers:

Häussler J, Sahlin U, Baey C, Smith HG, Clough Y (2017) Predicting pollinator capital and pollination service responses to enhancing floral and nesting resources. Ecology and Evolution, 7: 1898-1908. http://dx.doi.org/10.1002/ece3.2765

Gardner, E., Breeze, T.D., Clough, Y., Baldock, K.C.R., Campbell, A., Garratt, M., Gillespie, M.A.K., Kunin, W.E., Mckerchar, M., Memmott, J. and Potts, S.G. (2020). Reliably predicting pollinator abundance: challenges of calibrating process-based ecological models. Methods in Ecology and Evolution (accepted 8 July 2020)

Image, M., Gardner, E., Clough, Y., Smith, H. G., Baldock, K. C. R., Campbell, A., … Breeze, T. D. (2021). Does agri-environment scheme participation in England increase pollinator populations and crop pollination services? Agriculture, Ecosystems & Environment (accepted 5 Nov 2021)

https://zenodo.org/badge/latestdoi/427044957
