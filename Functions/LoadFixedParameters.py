#Produce fixed parameters for model

#Takes a set of input tables and returns parameters as python dictionaries or lists to facilitate lookup in subsequent functions
#Works only for 2 seasons, as required by our model

#Takes as inputs
    #   species     list of species codes being considered as integers (e.g. [1,2,3,4])
    #   florcov     a csv of floral cover parameters for each land cover and edge feature by land use code
    #   lfn         a csv containing the edge widths for each edge feature by land use code
    #   growth      a csv containing the population growth parameters by species code
    #   av          a csv containing the average nest density / ha by species code
    #   distances   a csv containing the foraging and nesting dispersal distances by species code

#Returns as outputs:
    #   FlorDict    a dictionary of floral cover parameters, where key = land use code
    #   LfnDict     a dictionary of edge width parameters, where key = land use code
    #   gr_list     a list of arrays of growth parameters where rows are queens and workers and columns are a, b, max, and pw; one for each species
    #   pw_list     a list of proportion of workers, one for each species 
    #   av_list     a list of av parameters, one for each species
    #   dist_list   a list of arrays of distance parameters where first row is foraging and second is nesting for each species

#Libraries
    #   pandas, os

import pandas as pd
import os

def ecodealFixedParams(species, florcov, lfn, growth, av, distances):

    #load in csvs   
    fc_df = pd.read_csv(florcov)
    lfn_df = pd.read_csv(lfn)
    gr_df = pd.read_csv(growth)
    av_df = pd.read_csv(av)
    dist_df = pd.read_csv(distances)

    #Build dictionaries   
    FlorDict = dict([(i, [a,b,c]) for i, a,b,c in zip(fc_df.code, fc_df.Flor_Cov_P1_mean, fc_df.Flor_Cov_P2_mean, fc_df.Flor_Cov_P3_mean)])
    LfnDict = dict([(i, [a,b,c]) for i, a,b,c in zip(lfn_df.code, lfn_df.field_edge, lfn_df.on_scenario, lfn_df.off_scenario)])

    #Build other data sources
    gr_list = [None]*len(species)
    pw_list = [None]*len(species)
    av_list = [None]*len(species)
    dist_list = [None]*len(species)
    tic = 0
    for j in species:
        aq = gr_df.best_guess[(gr_df.species == j) & (gr_df.short_notation == 'aq')]
        bq = gr_df.best_guess[(gr_df.species == j) & (gr_df.short_notation == 'bq')]
        qm = gr_df.best_guess[(gr_df.species == j) & (gr_df.short_notation == 'QEmax')]
        aw = gr_df.best_guess[(gr_df.species == j) & (gr_df.short_notation == 'aw')]
        bw = gr_df.best_guess[(gr_df.species == j) & (gr_df.short_notation == 'bw')]
        wm = gr_df.best_guess[(gr_df.species == j) & (gr_df.short_notation == 'Wmax')]
        pw = gr_df.best_guess[(gr_df.species == j) & (gr_df.short_notation == 'pw')]
        avnh = av_df.best_guess[av_df.species == j]
        fd = dist_df.best_guess[(dist_df.species == j) & (dist_df.activity == 'foraging')]
        nd = dist_df.best_guess[(dist_df.species == j) & (dist_df.activity == 'nesting')]
        gr_list[tic] = [(aq.iloc[0], bq.iloc[0], qm.iloc[0]), (aw.iloc[0], bw.iloc[0], wm.iloc[0])]
        pw_list[tic] = pw.iloc[0]
        av_list[tic] = avnh.iloc[0]
        dist_list[tic] = [fd.iloc[0], nd.iloc[0]]
        tic += 1

    return [FlorDict, LfnDict, gr_list, pw_list, av_list, dist_list]




    
    
    