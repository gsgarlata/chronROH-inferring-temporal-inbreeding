#This scripts creates plots of the "inferred" Froh temporal trajectory (i.e., chronROH) based on "true" ROH-segments and compared it with the "true" Froh temporal trajectory,
#for each simulated demographic model. The "true" ROH-segments are extracted from individuals sampled at the present-day (t = 0).


#It sources the functions used for making the plots.
source('./function_data_analysis_ROHsegments_sims.R')

pre_pathdata = args(1) #Path of the directory containing all the simulated demographic models.  
G_bp = 1e9 #Lemgth of the simulated chromosome in base-pairs
n_sims = 50 #Number of simulation replicates
TimePopChange = 200 #Time of population size change in the "Population Size Change" demographic models.
labels_list = c('Constant','Bottleneck','Expansion') #List of demographic model labels.
#List of the simulated demograohic scenarios.
scenario_list = c('Panmictic_Ne500_SampleSize30_rbp1e-08_GenomeLength1e+09',
                  'Panmictic_ancientNe500_PopSizeChangeTime200_recentNe100_SampleSize30_rbp1e-08_GenomeLength1e+09',
                  'Panmictic_ancientNe500_PopSizeChangeTime200_recentNe10000_SampleSize30_rbp1e-08_GenomeLength1e+09')
SamplingTime_list = c(0,50,100,150,200,250,350,550,1000) #Sampling time points


suffix_plotname = 'ConstantNe500_Bottleneck500_100_Expansion500_10000_G1e9'


#It makes plots of the inferred chronROH based on true ROH-segments (extracted from the simulated tree-sequences) and compare it 
#to the "true" temporal trajectory of Froh (based on inferred ROH-segment for each individual and sampling time point).
ChronROH_trueROHsegments_vs_trueFroh_plots = GetChronROH_trueROHsegments_vs_trueFroh(pre_path_data = pre_pathdata,G_bp,n_sims,TimePopChange,SamplingTime_list,
                                                                     scenario_list,labels_list)


output_plot_chronROH = paste0(pre_pathdata,'/chronoROH_based_on_trueROHsegments_VS_true_temporal_Froh_',suffix_plotname)
ggsave(file = paste0(output_plot_chronROH,'_raw.png'), plot = ChronROH_trueROHsegments_vs_trueFroh_plots, device = 'png', width = 15, height = 5)
