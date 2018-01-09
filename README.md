# kiloTools

Basic tools for spikesorting ephys data with kiloSort.
The repo has a number of functions that take user from raw ephys data to a user friedly(-ish) class 'neuron' of size nNeurons with associated neuron-based methods. 

## Organization
**conversion** - `convertPlxToRawBInary.m` - converts raw plx files to binary foramt appropriate for kilo (.dat).
also extracts timing information and constructs a time vector, along with the times at which strobed evets take place. 
**preprocessing** - files may be pre-processed to aid spike sorting. Options include (but not limited to) artifact removal, notch filtering, and more and more..
**spike sorting** - One file to rule them all:

masterMegaFile - includes the 3 core files from the kilosort repo (masterfile, config, chanMap) to facilitate fast error-free sorting. 

## How to run
### convert to binary:
Run the conversion script appropriate to your raw data. 
if plx - `convertPlxToRawBInary.m`
if alphaLab - `yyyyy`
if open ephys - `xxxxx`
this saves:

(1) a [nChannels x nSamples] .dat file 

(2) a [1 x nSamples] time vector with plexon time of each sample

(3) a [2 x nEvents] vecotr of strobe times (row 1) and strobe values (row 2).

### sort sort sort:
Input your dataset directory and name into 'masterMegaFile.m'

Start the sort by running 'masterMegaFile.m'





