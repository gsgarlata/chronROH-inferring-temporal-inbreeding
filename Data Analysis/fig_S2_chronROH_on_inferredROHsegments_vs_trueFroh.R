#This scripts creates plots of the "inferred" Froh temporal trajectory (i.e., chronROH) based on "inferred" ROH-segments and compared it with the "true" Froh temporal trajectory,
#for each simulated demographic model. The "inferred" ROH-segments are estimated from individuals sampled at the present-day (t = 0) using bcftools roh.


#It sources the functions used for making the plots.
source('./function_data_analysis_ROHsegments_sims.R')

pre_pathdata = args(1) #Path of the directory containing all the simulated demographic models.  
G_bp = 1e9 #Lemgth of the simulated chromosome in base-pairs
n_sims = 50 #Number of simulation replicates
TimePopChange = 200 #Time of population size change in the "Population Size Change" demographic models.
SamplingTime_list = c(0,50,100,150,200,250,350,550,1000) #Sampling time points
#List of the simulated demograohic scenarios.
scenario_list = c('Panmictic_Ne500_SampleSize30_rbp1e-08_GenomeLength1e+09',
                  'Panmictic_ancientNe500_PopSizeChangeTime200_recentNe100_SampleSize30_rbp1e-08_GenomeLength1e+09',
                  'Panmictic_ancientNe500_PopSizeChangeTime200_recentNe10000_SampleSize30_rbp1e-08_GenomeLength1e+09')

labels_list = c('Constant','Bottleneck','Expansion') #List of demographic model labels.
bin_size_roh = 1 #interval of the age bin to which each inferred ROH-segment is assigned.
standardize_Froh_val = TRUE #for visualization purposes, inferred and true Froh are standardized.

#It makes plots of the inferred chronROH based on inferred ROH-segments (using bcftools roh from present-day samples) and compare it 
#to the "true" temporal trajectory of Froh (based on inferred ROH-segment for each individual and sampling time point).
GetChronROH_inferredROHsegments_vs_trueFroh = GetChronROH_inferredROHsegments_vs_trueFroh(pre_path_data = pre_pathdata,G_bp,n_sims,TimePopChange,
                                                                                 SamplingTime_list,
                                                                                 scenario_list,
                                                                                 labels_list,
                                                                                 bin_size = bin_size_roh,
                                                                                 standardize_Froh = standardize_Froh_val)


if(standardize_Froh_val==TRUE){
  output_plot_chronROH = paste0(pre_pathdata,'/chronoROH_based_on_inferredROHsegments_VS_true_temporal_Froh_',suffix_plotname,'_bin_interval',bin_size_roh,'_standardised')
}else{
  output_plot_chronROH = paste0(pre_pathdata,'/chronoROH_based_on_inferredROHsegments_VS_true_temporal_Froh_',suffix_plotname,'_bin_interval',bin_size_roh)
}

ggsave(file = paste0(output_plot_chronROH,'.png'), plot = ChronROH_bcftools_vs_age_plots_roh, device = 'png', width = 15, height = 5)
