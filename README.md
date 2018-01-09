# kiloTools

Basic tools for spikesorting ephys data with kiloSort. The repo has a number of functions that take user from raw ephys data to a user friedly(-ish) class 'neuron' of size nNeurons with associated neuron-based methods. 

## Organization
**conversion** - `convertPlxToRawBInary.m` - converts raw files (e.g. .plx) to binary format appropriate for kilo (.dat). Also extracts strobed event timing and values. 

**preprocessing** - files may be pre-processed to aid spike sorting. For example artifact removal, notch filtering..

**spike sorting** - `masterMegaFile` - includes the 3 core files from the kilosort repo (masterfile, config, chanMap) to facilitate fast error-free sorting. 

** port to matlab** - port the sorted data into matlab for manipulation.

## Workflow
The following functions are all in a script title 'main_workflow.m'. 
Below are the details of each.
### run conversion:
Run the conversion script appropriate to your raw data. 

- if plx - `convertPlxToRawBInary.m`

- if alphaLab - `yyyyy`

- if open ephys - `xxxxx`

These an Import Wizard that allows you to select 1 or more files for conversion. 

This saves 3 files:

(1) a [nChannels x nSamples] .dat file 

(2) a [1 x nSamples] time vector with plexon time of each sample

(3) a [2 x nEvents] vecotr of strobe times (row 1) and strobe values (row 2).

### run kiloSort:
Input your dataset directory and name into 'masterMegaFile.m'

Start the sort by running 'masterMegaFile.m'

### manual curation using phy
Run Phy on your kiloSorted data. Completion of this stage will update the 'cluster_groups.csv' & 'spike_clusters.npy' files.

### get data into matlab
Run 'getSp'. getSp is heavily based on Spikes repo to extract useful info from the kiloSorted dataset. It is a single struct that has all cluster IDs, spike times, templates, and more and more. 

### convert sp into 'neuron' class
The neuron class is the end all be all. This is a struct of size nNeurons with useful methods. For example, these inlcude methods such as plot_isi, plot_acg, plot_spikeRate, plot_waveForm. type 'help neuron' for more indormation.






