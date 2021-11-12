#ComputeFlorNest (array version for Python 3)

#Takes a landcover array and a set of edge arrays, and assigns parameters for floral and nesting resources
#Works on one species only and for one scenario for 3 seasons

#Takes as inputs
    #   lc          a landcover array, integer type where each value corresponds to a specific land use code
	#	lccodes		a list of integer codes corresponding to all the land use codes (e.g [0,1,2,...,45,46])
    #   edges       a list of edge arrays, floating point type where each array has an identifer corresponding
                    #to a specific land use code and the values represent the length of that type of edge feature in the cell
    #   edgecodes   a list of integer codes corresponding to the land use codes of the edge arrays (e.g. [5,11,12])
    #   AttDict     a dictionary of floral and nesting attractiveness parameters where the keys are the land cover codes and the values are flor1, flor2, flor3, nest attractiveness
    #   FlorDict    a dictionary of floral cover parameters where the keys are the land cover codes and the values are flor1, flor2, flor3 floral covers
    #   LfnDict     a dictionary of edge widths where the keys are the land cover codes and the values are the edgewidths 

#Returns as outputs:
    #   flor        a list of arrays of floral values where each is identified by the season (flor1, flor2 etc.)
    #   nest        an array of nesting values for that species
    #   flor_lc    a list of arrays showing the proportion of floral value attributable the lc feature in each cell (one for each season)
    #   nest_lc    an array showing the proportion of nesting value attributable to the lc feature in each cell
    #   flor_ed    a list of (list of arrays showing the proportion of floral value attributable to the edge feature in each cell (one one for each edge)) (one for each season)
    #   nest_ed    a list of arrays showing the proportion of nesting value attributable to the edge feature in each cell (one for each edge feature)  

#Libraries
    #   numpy as np



def compFlorNest(lc, lccodes, edges, edgecodes, AttDict, FlorDict, LfnDict, cellsize):

    tiledims = lc.shape
    a_keys = np.array(lccodes)
	idx = np.digitize(lc.flatten(), a_keys, right = True)
	
    af1_vals = np.array([AttDict[i][0] for i in lccodes])
    af2_vals = np.array([AttDict[i][1] for i in lccodes])
    af3_vals = np.array([AttDict[i][2] for i in lccodes])
    an_vals = np.array([AttDict[i][3] for i in lccodes])
	f1_vals = np.array([FlorDict[i][0] for i in lccodes])
	f2_vals = np.array([FlorDict[i][1] for i in lccodes])
	f3_vals = np.array([FlorDict[i][2] for i in lccodes])
    
    lcat_f1 = af1_vals[idx].reshape(tiledims).astype('float32')
    lcat_f2 = af2_vals[idx].reshape(tiledims).astype('float32')
    lcat_f3 = af3_vals[idx].reshape(tiledims).astype('float32')
    lcat_n = an_vals[idx].reshape(tiledims).astype('float32')
	lcfc_f1 = f1_vals[idx].reshape(tiledims).astype('float32')
	lcfc_f2 = f2_vals[idx].reshape(tiledims).astype('float32')
	lcfc_f3 = f3_vals[idx].reshape(tiledims).astype('float32')
	
    r_edge = list(range(len(edgecodes)))

    #Multiply attract and florcov
    lcatfc_f1 = lcat_f1*lcfc_f1
    lcatfc_f2 = lcat_f2*lcfc_f2
    lcatfc_f3 = lcat_f3*lcfc_f3

    #Create the layers
    f1edge_cons = [AttDict[i][0] * FlorDict[i][0] for i in edgecodes]
    f2edge_cons = [AttDict[i][1] * FlorDict[i][1] for i in edgecodes]
    f3edge_cons = [AttDict[i][2] * FlorDict[i][2] for i in edgecodes]
    nestedge_cons = [AttDict[i][3] for i in edgecodes]

    #Create proportional edge rasters (store prop edge info)
       
    cellarea = cellsize*cellsize
    propedges = [(edges[t] * LfnDict[edgecodes[t]][0] / cellarea) for t in r_edge] 

    #Convert land cover to proportional values
    tot_pe = sum(propedges)
    tot_pl = np.maximum(1 - tot_pe, 0)   

    lcatp_f1 = (lcatfc_f1 * tot_pl).astype('float32')
    lcatp_f2 = (lcatfc_f2 * tot_pl).astype('float32')
    lcatp_f3 = (lcatfc_f3 * tot_pl).astype('float32')
    lcatp_n = (lcat_n * tot_pl).astype('float32')
    
    del lcatfc_f1, lcatfc_f2, lcatfc_f3, lcat_n, tot_pe, lcat_f1, lcat_f2, lcat_f3    #reduce memory burden

    flor1_list = [(propedges[t] * f1edge_cons[t]).astype('float32') for t in r_edge]
    flor2_list = [(propedges[t] * f2edge_cons[t]).astype('float32') for t in r_edge]
    flor3_list = [(propedges[t] * f3edge_cons[t]).astype('float32') for t in r_edge]
    nest_list = [(propedges[t] * nestedge_cons[t]).astype('float32') for t in r_edge]

    flor1 = sum(flor1_list + [lcatp_f1])
    flor2 = sum(flor2_list + [lcatp_f2])
    flor3 = sum(flor3_list + [lcatp_f3])
    nest = sum(nest_list + [lcatp_n])

    flor_lc = [lcatp_f1/flor1, lcatp_f2/flor2, lcatp_f3/flor3]
    flor_ed = [[flor1_list[t]/flor1 for t in r_edge], [flor2_list[t]/flor2 for t in r_edge], [flor3_list[t]/flor3 for t in r_edge]]
    nest_lc = lcatp_n/nest
    nest_ed = [nest_list[t]/nest for t in r_edge]


    return [[flor1, flor2, flor3], nest, flor_lc, nest_lc, flor_ed, nest_ed]




