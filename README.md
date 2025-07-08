# Assessing *chronROH* performance: a method for inferring temporal trajectory in inbreeding (based on F<sub>ROH</sub>)

Here we tested the *chronROH* approach proposed in [Orkin et al., 2025](https://www.nature.com/articles/s41559-024-02596-1) study entitled "Ecological and anthropogenic effects on the genomic diversity of lemurs in Madagascar". This repository containes two directory, one for performing simulations in [msprime](https://tskit.dev/msprime/docs/stable/intro.html#) and the other for processing the simulated data:

* [Simulations](Simulations): It contains scripts for carrying out *msprime* simulations for each of the three demographic models (Constant, Bottleneck, Expansion), extracting the true ROH-segments (present-day sampled individuals) from the tree-sequences and inferring ROH-segments using *bcftools roh* (across sampling time points).

* [Data Analysis](Data%20Analysis): It contains scripts for processing the true and inferred ROH-segments, estimating their age based on Orkin et al., (2025) approach and inferring the temporal trajectory of F<sub>ROH</sub> (i.e., chronROH).

