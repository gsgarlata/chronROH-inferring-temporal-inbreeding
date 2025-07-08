#This script is used to plot the distribution of ROH-segments across the simulated demographic models.
#For each demographic model, it plots the frequency of ROH-segments lengths (in bins given by the 'binwidth_val' parameter)
#and, within the same bin length, it groups ROH-segments by their age (in bins of about 1,000 generations).

#It sources the functions used for making the plots.
source('./function_data_analysis_ROHsegments_sims.R')

pre_pathdata = args(1) #Path of the directory containing all the simulated demographic models.  
SampleSize = 30 #Sample size of the sampled individuals
rbp = 1e-8 #Recombination rate per base pair per generation
G_bp = 1e9 #Lemgth of the simulated chromosome in base-pairs
n_sims = 50 #Number of simulation replicates
selectedPAIR = '0-1' #Individual ID: "0" and "1" refers to the node IDs of each of the haploid genomes of the same individual. 
binwidth_val = 0.0001 #bin width for the histogram
ConstantNe = 500 #Population Size in the "Constant Size" demographic model.
TimePopChange = 200 #Time of population size change in the "Population Size Change" demographic models.
AncientNe_bottleneck = 500 #Ancient Population size in the "Bottleneck" demographic model.
RecentNe_bottleneck = 100 #Recent Population size in the "Bottleneck" demographic model.
AncientNe_expansion = 500 #Ancient Population size in the "Expansion" demographic model.
RecentNe_expansion = 10000 #Recent Population size in the "Expansion" demographic model.


suffix_plotname = 'ConstantNe500_Bottleneck500_100_Expansion500_10000_G1e9'

#It plots the ROH-segment distribution for each scenario, and combine them in the same plot (one panel for each demographic model)
SimROHsegms_BinnedByAge_plots = GetSimulatedROHsegments_BinnedByAge_Across_Scenarios(pre_pathdata,n_sims,n_bin_age,binwidth_val,selectedPAIR,G = G_bp,rbp,
                                                                                     SampleSize,ConstantNe,TimePopChange,AncientNe_expansion,
                                                                                     RecentNe_expansion,AncientNe_bottleneck,RecentNe_bottleneck)

#Name of the outout plot
output_plot_binned_by_age = paste0(pre_pathdata,'/IBDsegments_distribution_by_age_',suffix_plotname)
ggsave(file = paste0(output_plot_binned_by_age,'.svg'), plot = SimIBDsegms_BinnedByAge_plots, device = 'png', width = 12, height = 5)

