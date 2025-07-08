# Analysis of the simulated ROH-segments in a panmictic population with constant size or undergoing population size change 

This directory provides the scripts to process the simulated and inferred ROH-segments across all demographic models. The scripts reproduce the three plots:

* [chronROH True ROH-segments](fig_S1_chronROH_on_trueROHsegments_vs_trueFroh.R): it plots the *chronROH* of each simulated demographic model obtained by using the *true* ROH-segments (extracted from the tree-sequences of individuals sampled at present-day (t = 0)) and estimating their age (based on their true ROH-segment length) using the equation in Orkin et al., (2025).

* [chronROH Inferred ROH-segments](fig_S2_chronROH_on_inferredROHsegments_vs_trueFroh.R): it plots the *chronROH* of each simulated demographic model obtained by using the *inferred* ROH-segments (estimated from individuals sampled at present-day (t = 0) using bcftools roh) and estimating their age (based on their inferred ROH-segment length) using the equation in Orkin et al., (2025).

* [Distribution of ROH-segments and age](fig_S3_age_distribution_roh_segments.R): it plots the the distribution of ROH-segments across the simulated demographic models. For each demographic model, it plots the frequency of ROH-segments lengths and, within the same bin length, it groups ROH-segments by their age (in bins of about 1,000 generations).

* [Auxiliary R functions](function_data_analysis_ROHsegments_sims.R): this file containes the R functions used in the three scripts above to generate the plots for the simulated data. 

