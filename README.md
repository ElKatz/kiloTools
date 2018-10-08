# kiloTools

Basic tools for spikesorting ephys data with kiloSort. The repo has a number of functions that take user from raw ephys data to a user friedly(-ish) class 'neuron' of size nNeurons with associated neuron-based methods. 

Functions in the kiloTools repo rely on the following repos:
- kilosrot - https://github.com/cortex-lab/KiloSort
- npy-matlab - https://github.com/kwikteam/npy-matlab
- Phy - https://github.com/kwikteam/phy
- spikes - https://github.com/cortex-lab/spikes
- sortingquality - https://github.com/cortex-lab/sortingQuality

## General organization
The toolbox takes you from raw ephys file (e.g. plx) to matlab-based sorted data. Here are the main steps & assocaited fucntions:

**conversion** - `convertRawToDat.m` converts raw files (e.g. .plx, .mpx..) to format appropriate for kilo (.dat). (may also perform some basic pre-processing to aid in spike sorting e.g. common mean subtraction & artifact removal but I find that kilo deals with idiosycracies well enough on its own so I don't use'em).  `convertRawToDat.m` also extracts a vector of time samples in real time (as opposed to sample number) and saves it. A seperate function extracts strobed info (`getStrobedFromRaw.m`).

**spike sorting** - `masterMegaFile` - includes the 3 core files from the kilosort repo (masterfile, config, chanMap). Runs kiloSort on GPU/CPU.

**port to matlab** - `getSp` ports the kiloSorted data into matlab and into the `sp` struct for easy manipulation. `sp2su` - converts the sp struct into a `su` struct of size nNeurons with neuron specific properties.


## Workflow
An example of my workflow can be found in the batch script `batch_convertAndSort.m`. In it you'll find somethling along the lines of:

### convert your raw ephys to .dat file
convert using `convertRawToDat.m`

### run kiloSort:
run `masterMegaFile.m`

### manual curation using phy
Run Phy on your kiloSorted data. Completion of this stage will update the 'cluster_groups.csv' & 'spike_clusters.npy' files.

### get data into matlab
Run `getSp`. getSp is a bastardized version of what was in the Spikes repo, and is used to extract useful info from the kiloSorted dataset. It is a single struct that has all cluster IDs, spike times, templates, and more and more. 

### convert sp into 'su' struct
run `sp2su`. (under construction). The neuron class is the end all be all. This is a struct of size nNeurons with useful methods. For example, these inlcude methods such as plot_isi, plot_acg, plot_spikeRate, plot_waveForm. type 'help neuron' for more indormation.






