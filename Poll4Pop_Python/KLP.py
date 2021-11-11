#Combined script for KernCalc, LatForDisp and PollPopGrowth to save loading them individually

#Libraries needed
    #   numpy as np
	#   scipy 
	#   scipy -> stats
	#   scipy -> signal 

#1. Exponential kernel calculation

#Exponential kernel computation to be used as an input for the latfordisp function
#Takes as inputs:
    #   beta        the mean of the exponential distribtion (ie the mean dispersal distance) - as a float!
    #   cellsize    the size of the array cell in m - as a float!
    #   maxsize     the maximum size of the resulting matrix - as a float!
    #   decaycut    the quantile used to cut the tails of the kernel


def kernCalc(beta, cellsize, maxsize, decaycut):
    radius = min(maxsize, round(stats.expon.ppf(decaycut, scale = beta/cellsize)))   #works out where the radius is that gives you the decaycut you specify
    nrow = (radius * 2) + 1   #ensures that the array produced has odd dimensions so there is a single central focal cell
    ncol = (radius * 2) + 1
    matnrow = int(nrow) #converts to integer to satisfy range() function. 
    matncol = int(ncol)
    matnrow_r = list(range(matnrow))
    matncol_r = list(range(matncol))
    idim = np.array(matnrow_r*matnrow)  #an array of i's organised 0,1,2,0,1,2 etc.
    jdim = np.array([[j] * matncol for j in matncol_r]).reshape(matncol*matncol, )  #an array of j's organised 0,0,0,1,1,1
    reach = np.sqrt(np.square(idim - radius) + np.square(jdim - radius)) <= radius  #boolean that is used to cut off cells outside the reach of the decay
    decay0 = np.exp(-(cellsize/beta) * np.sqrt(np.square(idim - radius) + np.square(jdim - radius)))  #result of decay function before it is cut off at the radius
    decay1 = (reach*decay0)  #cuts off the decay function at the radius
    decay = (decay1/sum(decay1)).reshape(matnrow, matncol)  #pro-rata to a 0-1 scale and convert to a n*n array
    return decay   #returns the dispersal array which can then be applied to latfordisp. In R version also returns the decaycut, but that seems unnecessary seeing as we have to specify it anyway in the main function


#2. Foraging / Dispersal model 

#Applies a spatial kernel weighted by avail resources so that foragers appear in cells that are both close by and have higher resources

#Takes as inputs:
    #   N        an array of foragers
    #   alpha    an array of resources
    #   decay    the dispersal kernel produced by kernCalc

#Returns as outputs:
    #   dwr     distance weighted resources: an array of resources weighted by the kernel
    #   vpc     visitation per cell: an array of by cell visitation rates
    #   rpf     resources per forager:  an array of distances weighted resources but only in cells where there are nesting queens
    


def latForDisp(N, alpha, decay):
    alpha = np.nan_to_num(alpha)
    N = np.nan_to_num(N)
    dwr = signal.fftconvolve(alpha, decay, mode = 'same')
    dwr[dwr < 1e-15] = 0
    relr = N*0
    relr = relr.astype(float)
    relr[dwr > 0] = N[dwr > 0] / dwr[dwr > 0]
    relr[dwr < 0] = 0
    dwrelr = signal.fftconvolve(relr, decay, mode = 'same')
    dwrelr[dwrelr < 1e-15] = 0
    vpc = dwrelr * alpha
    rpf = dwr*(N>0)
    return [rpf, vpc, relr]

#3. Pollinator growth function

#Applies a lognorm cdf to starting population 

#Takes as inputs:
    #   R        resources as an array
    #   a        the inflection point (i.e median)
    #   b        the slope of the growth
    #   max      the asymptote

#Returns as outputs:
    #   newN    the closing population as an array


def pollPopGrowth(R, a, b, max):
    mu = a
    sigma = np.sqrt(np.log(0.5 + 0.5*np.sqrt(1 + 4 * (np.square(b) / np.square(a)))))
    newN = max * stats.lognorm.cdf(R, sigma, loc=0, scale = mu)
    newN = np.nan_to_num(newN)
    return newN


