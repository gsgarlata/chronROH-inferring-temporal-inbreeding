# Simulating ROH-segments in a panmictic population with constant size or undergoing population size change 

This directory provides the scripts to simulate a panmictic population in msprime, derive ROH-segments from the simulated tree-sequences and infer ROH-segments from the corresponding VCF files using [bcftools roh](https://samtools.github.io/bcftools/howtos/roh-calling.html). In particular, it includes:

* [Constant Population Size](script_simulation_PanmicticWithStablePopSize.py): scripts for simulating a constant population size, sampling individuals at different time points, computing F<sub>ROH</sub> at each sampling time point and inferring the temporal trajectory of F<sub>ROH</sub> using present-day (t=0) samples.
  
* [Population Size Change](script_simulation_PanmicticWithPopSizeChange.py): scripts for simulating a panmictic population undergoing population size change, sampling individuals at different time points, computing F<sub>ROH</sub> at each sampling time point and inferring the temporal trajectory of F<sub>ROH</sub> using present-day (t=0) samples.

* [Auxiliary Functions](functions_msprime_sim_and_data_processing.py): collection of python functions used to simulate the demographic models, extract true ROH-segments and infer ROH-segments using bcftools roh. 
