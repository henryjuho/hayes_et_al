###############################################################################
########################### Calculating Divergence ############################
###############################################################################

#installing pakcages 
library(ape)

### Z chromo ###
#read fasta file
Z_alignment <- read.FASTA("/data/boa15kjh/greattit_data/alignment_Z.fasta")

#calculating divergence
Z_divergence <- dist.dna(Z_alignment, model = "K80", variance = FALSE,
         gamma = FALSE, pairwise.deletion = FALSE,
         base.freq = NULL, as.matrix = FALSE)

#computing divergence on Great Tit branch
Z_gt_fc <-  Z_divergence[1]
Z_zf_fc = Z_divergence[2]
Z_gt_zf = Z_divergence[3]

Z_gt <- (Z_gt_fc - Z_zf_fc + Z_gt_zf) / 2 

print("Z chromosome divergence =") 
print(Z_gt)



### Autosomes ###
#read fasta file
A_alignment <- read.FASTA("/data/boa15kjh/greattit_data/alignment_A.fasta")

#calculating divergence
A_divergence <- dist.dna(A_alignment, model = "K80", variance = FALSE,
         gamma = FALSE, pairwise.deletion = FALSE,
         base.freq = NULL, as.matrix = FALSE)


#computing divergence on Great Tit branch

A_gt_fc <-  A_divergence[1]
A_zf_fc = A_divergence[2]
A_gt_zf = A_divergence[3]

A_gt <- (A_gt_fc - A_zf_fc + A_gt_zf) / 2 

print("Autosome divergence =")
print(A_gt)