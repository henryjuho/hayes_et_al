# Scripts from Hayes et al. 2020

This repository contains the scripts used in the analyses presented in Hayes et al. (2020) accessible here: <add_link>

Code for the simulation results for the anavar package presented in the paper's supplementary can be found in a separate repository: <https://github.com/henryjuho/hayes_et_al_simulations>.

The scripts in this repository and their uses are listed below:

* ```calculating_divergence.R```: Calculates the divergence from a multispecies alignment. Requires two input data files: 1 - multispecies alignment for Z chromosome (FASTA format), 2 - multispecies alignment for Autosomes (FASTA format). Note that to calculate divergence for different regions of the genome the input files must first be filtered to contain only the desired regions.
* ```descriptive_stats_all_snps.py```: Calculates the statistics Nucleotide Diversity, Wattersons Theta and Tajimas D. Requires two input data files: 1 - callable sites (FATSA format), 2 - SNPs (VCF format). Note that to calculate these statistics for different regions of the genome the input files must first be filtered to contain only the desired regions.
* ```site_frequency_spectrum.py```: Extracts both the folded and unfolded site frequency spectrum (SFS) from a VCF file. Note that to get the (SFS) for different regions the VCF file must first be filtered to contain only the desired regions
* ```anavar_ctrl_file_autosomes_ancestral_repeats.txt```: Control file for anavar run on autosomal data with SNPs in ancestral repeats as neutral reference.
* ```anavar_ctrl_file_autosomes.txt```: Control file for anavar run on autosomal data with 4fold degenerate SNPs as neutral reference.
* ```anavar_ctrl_file_zchrom_ancestral_repeats.txt```: Control file for anavar run on the Z chromosome data with SNPs in ancestral repeats as neutral reference.
* ```anavar_ctrl_file_zchrom.txt```: Control file for anavar run on the Z chromosome data with 4fold degenerate SNPs as neutral reference.
* ```varne_ctrl_file_model1_H2_equalmutationrate_neldermead.txt```: VarNe control file for the equal mutation rate, 2 epoch model.
* ```varne_ctrl_file_model1_H2_fixedNeratio_neldermead.txt```: VarNe control file for the fixed effective population size ratio, 2 epoch model.
* ```varne_ctrl_file_model1_H2_general_neldermead.txt```: VarNe control file for the full 2 epoch model.