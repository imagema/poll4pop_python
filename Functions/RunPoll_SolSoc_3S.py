#Script for the runpoll functions for solitary and social bees
#Creates functions for...
#...solitaries foraging in early and late spring only (seasons 1 and 2)
#...solitaries foraging in summer only (season 3)
#...social bees foraging in early and late spring and summer with workers generated at the end of early spring and the end of late spring
#...separate functions depending on whether we are in Year 1 or Year N of the run

#Converts input of nesting and floral resources for 3 seasons into an output of floral visitation, number of workers, number of reproductive females produced and surviving

#Libraries needed
    #   numpy as np
	#   scipy 
	#   scipy -> stats
   

#Takes as inputs:
    #   nest     the nesting resources as an array
    #   fi        the floral resources (season i) as an array for the seasons active
    #   av        average nest density
    #   dist_n    nesting distance for the guild 
    #   dist_f    foraging distance for the guild 
    #   growth    array of growth factors for the guild (a, b, max) where row1 = reproductive females, row2 = workers (only first row used for solitaries)
    #   cellsize  dimension of cell in metres
    #   decaycut  threshold for decay function (normally 0.99)

#Returns as outputs:
    #   vpc_fi       the floral vistation rate per cell as an array (season i) for the seasons active
    #   rpf_fi      the resources collected per cell as an array (season i) for the seasons active
    #   np_1      the array of workers produced at end of season 1 (ie foraging in season 2)   - not returned for solitaries
    #   np_2      the array of workers produced at end of season 2 (ie foraging in season 3)   - not returned for solitaries
    #   mp_e      the array of reproductive females at the end of final season
    #   mp_0      the array of reproductive females surviving the winter - which is the input to the next year. 

#Prerequisite functions that must be loaded
    #   KernCalc
    #   LatForDisp
    #   PollPopGrowth



#Solitaries	
#active in Early and Late Spring (f1 and f2)

def runPoll3S_SolP1_Y1(nest, f1, f2, av, dist_n, dist_f, growth, cellsize, decaycut):

    #Step 1 - Create dispersal kernels
    kern_f = kernCalc(dist_f, cellsize, 0.5*len(f1), decaycut)
    kern_n = kernCalc(dist_n, cellsize, 0.5*len(nest), decaycut)

    #Step 2 - Initialise population

    mp_1 = nest*av*cellsize*cellsize/10000
    mp_1 = stats.poisson.rvs(mp_1)

    #Step 3 - Run for foraging period 1 and 2 
    lfd_f = latForDisp(mp_1, f1, kern_f)
    rpf_f1 = lfd_f[0]
    vpc_f1 = lfd_f[1]
    lfd_f = latForDisp(mp_1, f2, kern_f)
    rpf_f2 = lfd_f[0]
    vpc_f2 = lfd_f[1]
    
    mp_e = mp_1*pollPopGrowth(rpf_f1+rpf_f2, growth[0][0], growth[0][1], growth[0][2])
    lfd_n = latForDisp(mp_e, nest, kern_n)
    mp_0 = np.minimum(lfd_n[1], nest*av*cellsize*cellsize/10000)

    return [vpc_f1, vpc_f2, rpf_f1, rpf_f2, mp_e, mp_0]


def runPoll3S_SolP1_YN(nest, f1, f2, mp_0, av, dist_n, dist_f, growth, cellsize, decaycut):

    #Step 1 - Create dispersal kernels
    kern_f = kernCalc(dist_f, cellsize, 0.5*len(f1), decaycut)
    kern_n = kernCalc(dist_n, cellsize, 0.5*len(nest), decaycut)

    #Step 2 - Initialise population

    mp_1 = mp_0

    #Step 3 - Run for foraging period 1 and 2
    lfd_f = latForDisp(mp_1, f1, kern_f)
    rpf_f1 = lfd_f[0]
    vpc_f1 = lfd_f[1]
    lfd_f = latForDisp(mp_1, f2, kern_f)
    rpf_f2 = lfd_f[0]
    vpc_f2 = lfd_f[1]

    mp_e = mp_1*pollPopGrowth(rpf_f1+rpf_f2, growth[0][0], growth[0][1], growth[0][2])
    lfd_n = latForDisp(mp_e, nest, kern_n)
    mp_0 = np.minimum(lfd_n[1], nest*av*cellsize*cellsize/10000)

    return [vpc_f1, vpc_f2, rpf_f1, rpf_f2, mp_e, mp_0]


#active in Summer (f3)


def runPoll3S_SolP3_Y1(nest, f3, av, dist_n, dist_f, growth, cellsize, decaycut):

    #Step 1 - Create dispersal kernels
    kern_f = kernCalc(dist_f, cellsize, 0.5*len(f3), decaycut)
    kern_n = kernCalc(dist_n, cellsize, 0.5*len(nest), decaycut)

    #Step 2 - Initialise population

    mp_1 = nest*av*cellsize*cellsize/10000
    mp_1 = stats.poisson.rvs(mp_1)

    #Step 3 - Run for foraging period 1 and 2 
    lfd_f = latForDisp(mp_1, f3, kern_f)
    rpf_f3 = lfd_f[0]
    vpc_f3 = lfd_f[1]

    
    mp_e = mp_1*pollPopGrowth(rpf_f3, growth[0][0], growth[0][1], growth[0][2])
    lfd_n = latForDisp(mp_e, nest, kern_n)
    mp_0 = np.minimum(lfd_n[1], nest*av*cellsize*cellsize/10000)

    return [vpc_f3, rpf_f3, mp_e, mp_0]


def runPoll3S_SolP3_YN(nest, f3, mp_0, av, dist_n, dist_f, growth, cellsize, decaycut):

    #Step 1 - Create dispersal kernels
    kern_f = kernCalc(dist_f, cellsize, 0.5*len(f3), decaycut)
    kern_n = kernCalc(dist_n, cellsize, 0.5*len(nest), decaycut)

    #Step 2 - Initialise population

    mp_1 = mp_0

    #Step 3 - Run for foraging period 1 and 2
    lfd_f = latForDisp(mp_1, f3, kern_f)
    rpf_f3 = lfd_f[0]
    vpc_f3 = lfd_f[1]

    mp_e = mp_1*pollPopGrowth(rpf_f3, growth[0][0], growth[0][1], growth[0][2])
    lfd_n = latForDisp(mp_e, nest, kern_n)
    mp_0 = np.minimum(lfd_n[1], nest*av*cellsize*cellsize/10000)

    return [vpc_f3, rpf_f3, mp_e, mp_0]

#Social bees active in Early Spring, Late Spring and Summer

def runPoll3S_SocP123_Y1(nest, f1, f2, f3, av, dist_n, dist_f, growth, pw, cellsize, decaycut):

    #Step 1 - Create dispersal kernels
    kern_f = kernCalc(dist_f, cellsize, 0.5*len(f1), decaycut)
    kern_n = kernCalc(dist_n, cellsize, 0.5*len(nest), decaycut)

    #Step 2 - Initialise population

    mp_1 = nest*av*cellsize*cellsize/10000
    mp_1 = stats.poisson.rvs(mp_1)

    #Step 3 - Queens foraging in season 1
    lfd_f1 = latForDisp(mp_1, f1, kern_f)
    rpf_f1 = lfd_f1[0]
    vpc_f1 = lfd_f1[1]  

    #Step 4 - Workers produced at end of season 1
    np_1 = mp_1*pollPopGrowth(rpf_f1, growth[1][0], growth[1][1], growth[1][2])

    #Step 5 - Workers foraging in season 2
    lfd_f2 = latForDisp(pw*np_1, f2, kern_f)
    rpf_f2 = lfd_f2[0]
    vpc_f2 = lfd_f2[1]
    
    #Step 6 - Workers produced at end of season 2
    np_2 = np_1 + mp_1*pollPopGrowth(rpf_f2*pw*np_1/mp_1, growth[1][0], growth[1][1], growth[1][2])

    #Step 7 - Workers foraging in season 3
    lfd_f3 = latForDisp(pw*np_2, f3, kern_f)
    rpf_f3 = lfd_f3[0]
    vpc_f3 = lfd_f3[1]

    #Step 8 - Queens produced at the end of season 2 and then surviving the winter 
    mp_e = mp_1*pollPopGrowth(rpf_f3*pw*np_2/mp_1, growth[0][0], growth[0][1], growth[0][2])
    lfd_n = latForDisp(mp_e, nest, kern_n)
    mp_0 = np.minimum(lfd_n[1], nest*av*cellsize*cellsize/10000)

    

    return [vpc_f1, vpc_f2, vpc_f3, rpf_f1, rpf_f2, rpf_f3, np_1, np_2, mp_e, mp_0]


def runPoll3S_SocP123_YN(nest, f1, f2, f3, mp_0, av, dist_n, dist_f, growth, pw, cellsize, decaycut):

    #Step 1 - Create dispersal kernels
    kern_f = kernCalc(dist_f, cellsize, 0.5*len(f1), decaycut)
    kern_n = kernCalc(dist_n, cellsize, 0.5*len(nest), decaycut)

    #Step 2 - Initialise population

    mp_1 = mp_0

    #Step 3 - Queens foraging in season 1
    lfd_f1 = latForDisp(mp_1, f1, kern_f)
    rpf_f1 = lfd_f1[0]
    vpc_f1 = lfd_f1[1]  

    #Step 4 - Workers produced at end of season 1
    np_1 = mp_1*pollPopGrowth(rpf_f1, growth[1][0], growth[1][1], growth[1][2])

    #Step 5 - Workers foraging in season 2
    lfd_f2 = latForDisp(pw*np_1, f2, kern_f)
    rpf_f2 = lfd_f2[0]
    vpc_f2 = lfd_f2[1]
    
    #Step 6 - Workers produced at end of season 2
    np_2 = np_1 + mp_1*pollPopGrowth(rpf_f2*pw*np_1/mp_1, growth[1][0], growth[1][1], growth[1][2])

    #Step 7 - Workers foraging in season 3
    lfd_f3 = latForDisp(pw*np_2, f3, kern_f)
    rpf_f3 = lfd_f3[0]
    vpc_f3 = lfd_f3[1]

    #Step 8 - Queens produced at the end of season 2 and then surviving the winter 
    mp_e = mp_1*pollPopGrowth(rpf_f3*pw*np_2/mp_1, growth[0][0], growth[0][1], growth[0][2])
    lfd_n = latForDisp(mp_e, nest, kern_n)
    mp_0 = np.minimum(lfd_n[1], nest*av*cellsize*cellsize/10000)

    return [vpc_f1, vpc_f2, vpc_f3, rpf_f1, rpf_f2, rpf_f3, np_1, np_2, mp_e, mp_0]


    
    







