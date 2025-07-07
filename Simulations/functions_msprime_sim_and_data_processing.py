#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Gabriele Maria Sgarlata
"""

import msprime
import numpy as np
import pandas as pd
import gzip
import os
from subprocess import Popen


##START: Functions for msprime simulations##
def run_simulationIBD_PopSizeChange(SamplingTime_list,SampleSize,demography,recombination_rate,sequence_length,ploidy,seed):
    
    #It create an empty list, where the tree sequence at each sampling time is saved.
    tree_list = []
    
    #It loops over each sampling time point
    for i, SamplingTime in enumerate(SamplingTime_list):
        
        #It defines the "sampling" information: sample size and time of sampling.
        SampleSize_Info = [msprime.SampleSet(SampleSize, time=SamplingTime)]
        
        #Perform the msprime simulation
        ts_IBD = msprime.sim_ancestry(samples = SampleSize_Info, demography=demography,model='dtwf',
                                  recombination_rate=recombination_rate,
                                  sequence_length=sequence_length, ploidy = ploidy, random_seed=seed)
        
        #Append the tree sequence to the "tree_list" object.
        tree_list.append({SamplingTime:ts_IBD})
                
    return tree_list


def run_simulationIBD_ConstantPopSize(SamplingTime_list,SampleSize,demography,model,recombination_rate,sequence_length,ploidy,seed):
    
    #It create an empty list, where the tree sequence at each sampling time is saved.
    tree_list = []
    
    #It loops over each sampling time point
    for i, SamplingTime in enumerate(SamplingTime_list):
        
        #It defines the "sampling" information: sample size and time of sampling.
        SampleSize_Info = [msprime.SampleSet(SampleSize, time=SamplingTime)]
        
        #Perform the msprime simulation
        ts_IBD = msprime.sim_ancestry(samples = SampleSize_Info, demography=demography,model='dtwf',
                                  recombination_rate=recombination_rate,
                                  sequence_length=sequence_length, ploidy = ploidy, random_seed=seed)
        
        #Append the tree sequence to the "tree_list" object.
        tree_list.append({SamplingTime:ts_IBD})
                
    return tree_list
##END: Functions for msprime simulations##


##START: Functions for ROH segments analyses##

def to_pairs_list(arr):
    """
    Converts a NumPy array into a list of pairs of arrays.

    Args:
        arr (np.ndarray): The input NumPy array.

    Returns:
        list: A list of pairs of arrays.
              Each element in the list is a tuple containing two NumPy arrays 
              representing a pair of consecutive elements from the input array.
    """
    if arr.size % 2 != 0:
        raise ValueError("Array size must be even to form pairs.")
    
    reshaped_arr = arr.reshape(-1, 2)
    pairs_list = [(reshaped_arr[i][0], reshaped_arr[i][1]) for i in range(len(reshaped_arr))]
    return pairs_list

#"ibd_segments" function computes the distribution of IBD segments between two nodes "a" and "b".
#In thsi application, we consider as the two haploid genomes of an individual as the two nodes.
def ibd_segments(ts, a, b):
    #It extract the sequence of trees from the "ts" object.
    trees_iter = ts.trees()
    #It access the first tree in the tree sequence
    tree = next(trees_iter)
    #It extract the node ID of the most recent common ancestor between the two haploid genomes of the a-b individual
    #at the first tree in the sequence.
    last_mrca = tree.mrca(a, b)
    last_left = 0
    segment_lengths = []
    segment_age = []
    #It loops over the tree sequence
    for tree in trees_iter:       
        #It extract the node ID of the most recent common ancestor between "a" and "b".
        mrca = tree.mrca(a, b)
        #If "mrca" is different from the "last_mrca" then:
        if mrca != last_mrca:
            #It extract the end position (on the right side) of the IBD segment corresponding to the previous tree (in base-pairs units).
            left = tree.interval[0]
            #It derives the length of the IBD segment of the previous tree.
            segment_lengths.append(left - last_left)
            #It extract the age of the mrca of the previous tree.
            segment_age.append(tree.time(last_mrca))
            last_mrca = mrca
            last_left = left
    #It extract the length and age of the last IBD segment.
    segment_lengths.append(ts.sequence_length - last_left)
    segment_age.append(tree.time(last_mrca))

    return [np.array(segment_lengths),np.array(segment_age)]

#"roh_segments" function computes the distribution of IBD segments within all the sampled individuals.
def roh_segments(ts, recombination_rate, output_name, out_dir, save):
    
    #It extracts the list of samples
    inds = ts.samples()
    #It return pairs of samples, corresponding to the two haploid genomes of each individual
    inds_list = to_pairs_list(inds)
    
    #It defines the length (in base pairs) corresponding to 0.01 Morgans (M) or 1 centi Morgans (cM), given the recombination rate. 
    L_001 = 0.01 / recombination_rate   
    
    res_df = pd.DataFrame()
    count = 0
    #It loops over each individual
    for indGenome in inds_list:
        count = count + 1
        #It extract each haploid genome for a given individual in object "a" and "b".
        a = indGenome[0]
        b = indGenome[1]
        #It extracts the distribution of IBD segments between the two haploid genomes of the individual.
        segment_res = ibd_segments(ts, a, b)
        #It extract the length of each IBD segment
        segment_lengths = segment_res[0]
        #It extract the age of each IBD segment
        segment_age = segment_res[1]
        #It convert IBD segment lenghts from base-pairs to cM units
        segment_cM = (np.array(segment_lengths) / L_001)
        #It create a dataframe with all this information
        d = {'ind': count, 'IBDage': segment_age, 'IBDlength_bp' : segment_lengths,
             'IBDlength_cM' : segment_cM, 'haploid_pair' : str(a) + '-' + str(b)}
        temp_df = pd.DataFrame(data=d)
        res_df = pd.concat([res_df, temp_df], ignore_index=True)
       
    file_name = '{}/{}.csv'.format(out_dir.rstrip('/'), output_name)
    if save:
        res_df.to_csv(file_name,index=False)        
    return res_df


def bin_data(data, bin_size):
    """Bins numerical data into intervals of specified size.

    Args:
        data (list or pd.Series): Numerical data to bin.
        bin_size (int or float): Size of each bin.

    Returns:
        pd.Series: Binned data, with each value replaced by its bin number.
    """

    if not isinstance(data, pd.Series):
        data = pd.Series(data)

    min_val = data.min()
    max_val = data.max()
    num_bins = int((max_val - min_val) // bin_size) + 1
    bins = [min_val + i * bin_size for i in range(num_bins + 1)]

    return pd.cut(data, bins=bins, labels=False, include_lowest=True, right=False)

#"binDataAndEstimateFroh" function groups ROH-segment in intervals based on their 'estimated' or 'true' age.
def binDataAndEstimateFroh(df, bin_size, age_type, G, recombination_rate):
    
    #It defines the length (in base pairs) corresponding to 0.01 Morgans (M) or 1 centi Morgans (cM), given the recombination rate.
    L_001 = 0.01 / recombination_rate  

    #It bins the ROH-segments based on their age and type ('estimated age' or 'true age') in interval of size "bin_size".
    binned_Age = bin_data(df[age_type], bin_size)
    #It assigns a bin interval at each ROH-segment.
    df = df.assign(bin_age = binned_Age)
    
    #It computes the average age of the ROH-segments within the same bin.
    average_Age_bin_Age = df.groupby('bin_age')[age_type].mean().reset_index()
    average_Age_bin_Age = average_Age_bin_Age.rename(columns={age_type: 'age'})
    average_Age_bin_Age = average_Age_bin_Age.drop('bin_age', axis=1)
    average_Age_bin_Age = average_Age_bin_Age.reset_index(names='bin_age')
    
    #It sums the length of the all the ROH-segments within the same bin.
    sum_lengthBP_bin_Age = df.groupby('bin_age')['IBDlength_bp'].sum().reset_index()
    sum_lengthBP_bin_Age = sum_lengthBP_bin_Age.drop('bin_age', axis=1)
    sum_lengthBP_bin_Age = sum_lengthBP_bin_Age.reset_index(names='bin_age')
    
    #It merges the two dataframes.
    merged_df = pd.merge(average_Age_bin_Age, sum_lengthBP_bin_Age, on='bin_age', how='inner')
    
    #It converts the sum of ROH-segments from base-pairs to Morgans and centi Morgans.
    merged_df = merged_df.assign(IBDlength_cM = merged_df['IBDlength_bp']/L_001)
    merged_df = merged_df.assign(IBDlength_M = merged_df['IBDlength_cM']/100)
    #It computes the Froh for each bin.
    merged_df = merged_df.assign(F_roh = merged_df['IBDlength_M']/G)
    
    merged_df = merged_df.drop('IBDlength_bp', axis=1)
    merged_df = merged_df.drop('IBDlength_cM', axis=1)
    merged_df = merged_df.drop('IBDlength_M', axis=1)
    merged_df = merged_df.drop('bin_age', axis=1)
    return merged_df



# "ComputeChronoROH" estimates the Froh temporal trajectory from the present-day ROH-segments of all sampled inviduals (one trajectory per individual).
def ComputeChronoROH(df, bin_size, recombination_rate, Gval, save, file_name):
    
    #It extracts the list of individuals contained in the "df" dataframe.
    sample_ids = df['ind'].unique()

    finalDF = pd.DataFrame()
    #It loops over individuals.
    for sample_id in sample_ids:
        #It subsets the "df" at the focal individual "sample_id".
        df_time_sampleID = df.loc[df['ind'] == sample_id]
        #It extimate the age of the ROH segment as in Orkin et al., (2025).
        df_time_sampleID = df_time_sampleID.assign(estIBDage = 100/(2 * df_time_sampleID['IBDlength_cM']))
        
        #It groups ROH-segments (based on their ESTIMATED age) in bins of size "bin_size" and compute the Froh as = sum(Lroh)/ G (where "G" is the length of the genome in M).
        df_trueAgeChronoRoh = binDataAndEstimateFroh(df = df_time_sampleID, bin_size = bin_size, age_type = 'IBDage', G = Gval, recombination_rate = recombination_rate)
        #It groups ROH-segments (based on their TRUE age) in bins of size "bin_size" and compute the Froh as = sum(Lroh)/ G (where "G" is the length of the genome in M).
        df_estAgeChronoRoh = binDataAndEstimateFroh(df = df_time_sampleID, bin_size = bin_size, age_type = 'estIBDage', G = Gval, recombination_rate = recombination_rate)
        
        df_trueAgeChronoRoh = df_trueAgeChronoRoh.assign(type = 'IBDage')
        df_estAgeChronoRoh = df_estAgeChronoRoh.assign(type = 'estIBDage')
        
        ind_finalDF = pd.concat([df_trueAgeChronoRoh, df_estAgeChronoRoh])
        ind_finalDF = ind_finalDF.assign(ind = sample_id)
        
        finalDF = pd.concat([finalDF,ind_finalDF])

    if save:
        finalDF.to_csv(file_name,index=False)        
    return finalDF

# ""roh_segments_MultipleSamplingTimes" function infers the Froh temporal trajectory from present-day sampled individuals based on the true ROH-segments
#distribution while estimating the ROH-segments (bcftools roh) from each sampled time point. The estimated ROH-segment are used to compute
#the true Froh over time.
def roh_segments_MultipleSamplingTimes(ts_list, bin_size, recombination_rate, output_name, out_dir, mutation_rate, error_rate, G, simID):
    
    #It loops over tree sequences, each corresponding to a sampling time point.
    for ts_dict in ts_list:
        
        #Get sampling times from tree sequence dictionary
        SamplingTime, tree_t = list(ts_dict.items())[0]
        
        output_tree = output_name + '_Time' + str(SamplingTime)
        
        #It extract the distribution of within-individual IBD segments (closely related to run-of-homozygosity)
        df_time = roh_segments(ts = tree_t, recombination_rate = recombination_rate, output_name = output_tree, out_dir = out_dir, save = False)
        
        #The "if-conditioning" considers that for present-day samples we need to compute the temporal Froh as in Orkin et al., (2025),
        #using either the true age of the IBD segments or using the equation "age = 1/(2 * Lroh)" as in Orkin et al., (2025).
        if SamplingTime == 0:
            out_ibd_segm_file = out_dir + '/' + output_name + '_time' + str(SamplingTime) + '_chronoRoh' + '_sim' + str(simID) + '.csv'
            ComputeChronoROH(df = df_time, bin_size = bin_size, recombination_rate = recombination_rate, Gval = G, save = True, file_name = out_ibd_segm_file)
        
        #It add mutations to the tree sequence, a step required for saving the tree sequence as a genotype matrix (vcf file)
        mts = msprime.sim_mutations(tree_t, rate=mutation_rate)
        
        msprime_vcf_dir = ('{}/sim' + str(simID)).format(out_dir)
        if not(os.path.exists(msprime_vcf_dir) and os.path.isdir(msprime_vcf_dir)):
            os.mkdir(msprime_vcf_dir)

        vcf_path = '{}/{}.vcf.gz'.format(msprime_vcf_dir.rstrip('/'), output_tree)

        with gzip.open(vcf_path, 'wt') as vcf_file:
            mts.write_vcf(vcf_file, None)    
        
        output_file = '{}/{}_roh.txt'.format(msprime_vcf_dir.rstrip('/'), output_tree)

        Popen(["bcftools","roh","--GTs-only", str(error_rate), "--rec-rate", 
                str(recombination_rate), "--output", output_file, "-Or", "--estimate-AF", "-", vcf_path])
