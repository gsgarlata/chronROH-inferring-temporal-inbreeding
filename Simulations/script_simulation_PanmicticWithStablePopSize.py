#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Gabriele Maria Sgarlata
"""

import os
import msprime
import random
import functions_msprime_sim_and_data_processing as PanmicticMsprime
import pickle 


PopSize = 500 #Ancestral population size
rbp = 1e-8 #Recombination rate per base pair per generation
mu = 1e-8 #Mutation rate per base pair per generation
region_length = (10**9) #Lenght of the simulated chromosome
SampleSize = 30 #Diploid individual sample size
n_sims = 50 #Number of simulation replicates
G_M = region_length * rbp #Length of the simulated chromosome in Morgan 
bin_size = 100 #bin size used for grouping ROH-segments based on their age.
SamplingTime_list = [0, 50,100,150,200,250,350,550,1000] #Sampling times over the course of the simulation
path_res = [] #Specify here the directory where the results will be saved.


##START: Running script##
region_length_expFormat = "{:.0e}".format(region_length)
out_dir_name = 'Panmictic_Ne' + str(PopSize) + '_SampleSize' + str(SampleSize) + '_rbp' + str(rbp) + '_GenomeLength' + region_length_expFormat

#Create the 'out_dir' if it does not exist yet
out_dir = path_res + '/' + out_dir_name
if not(os.path.exists(out_dir) and os.path.isdir(out_dir)):
    os.makedirs(out_dir)

#Loop over each simultaion replicate
for simID in range(1,(n_sims + 1),1):
    
    #Define output name of the file to save
    output_name = out_dir_name + '_sim' + str(simID)
    
    #Create a directory for each simulation replicate, if it does not exist yet
    sim_dir = out_dir + '/sim' + str(simID)
    if not(os.path.exists(sim_dir) and os.path.isdir(sim_dir)):
        os.makedirs(sim_dir)
        
    ##START: Performing simulation in msprime##
    #Define the name of the output file where to save tree sequences 
    dict_simulation_file = out_dir + '/sim' + str(simID) + '/' + output_name + '.pkl'
    
    #If the tree sequence file already exists, open the existing one otherwise carry out the simulation
    if os.path.isfile(dict_simulation_file):
        with open(dict_simulation_file, 'rb') as f:
            res_onesim = pickle.load(f)
    else:
        #Define a random seed
        random_seedAnc = random.randint(0, int(1e7))
        
        #It creates the msprime demographic model
        demog_model = msprime.Demography()
        
        #Define the population size "PopSize"
        demog_model.add_population(initial_size=PopSize)    
        
        #Run the simulation in msprime using the custom function "run_simulationIBDpanmicticMultipleSamplingTimes" 
        res_onesim = PanmicticMsprime.run_simulationIBD_ConstantPopSize(SamplingTime_list = SamplingTime_list, SampleSize = SampleSize, demography = demog_model, 
                                            recombination_rate = rbp, sequence_length = region_length, ploidy = 2,
                                            seed = random_seedAnc)
        
        #It saves the tree sequence in the "dict_simulation_file" file.
        with open(dict_simulation_file, 'wb') as f:
            pickle.dump(res_onesim, f)
    ##END: Performing simulation in msprime##


    ##START: Estimate the Froh temporal trajectory from present-day sampled individuals (as in Orkin et al., 2025) and 
    #compute the "true" Froh temporal trajectory by inferrin ROH segments from each sampled time point.
    chronoROH_filename = out_dir + '/' + out_dir_name + '_time0_chronoRoh_sim' + str(simID) + '.csv'
    
    if os.path.isfile(chronoROH_filename):
        pass
    else:
        #It infers the Froh temporal trajectory from present-day sampled individuals based on the true ROH-segments
        #distribution while estimating the ROH-segments (bcftools roh) from each sampled time point and individual.
        PanmicticMsprime.roh_segments_MultipleSamplingTimes(ts_list = res_onesim, bin_size = bin_size, recombination_rate = rbp, output_name = out_dir_name, out_dir = out_dir, mutation_rate = mu, error_rate = 30, G = G_M, simID = simID)
    ##END: Estimate the Froh temporal trajectory from present-day sampled individuals (as in Orkin et al., 2025) and the "true" Froh temporal trajectory
    
