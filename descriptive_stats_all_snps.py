# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 13:41:55 2018

@author: Keith
"""
###############################################################################
############################ Nucleotide Diversity  ############################
###############################################################################
    
#installing packages 
import math
import vcf
from Bio import SeqIO

#defining function for pi 
def pi(n, allele_frq_list):
    no_seqs_dif = [(1.0 - raf**2 - (1.0-raf)**2) * (n/(n-1.0)) for raf in allele_frq_list]
    seq_pi = sum(no_seqs_dif)
    return seq_pi


#number of callable sites 
Z_callable_sites = 0
A_callable_sites = 0 
callable_sites = SeqIO.parse("/data/boa15kjh/greattit_data/n20callablesites.fa", "fasta")
for seq_record in callable_sites:
    if seq_record.id.startswith("chrZ"):
        Z_callable_sites = Z_callable_sites + seq_record.seq.upper().count("K") + seq_record.seq.upper().count("R")
    else:
        A_callable_sites = A_callable_sites + seq_record.seq.upper().count("K") + seq_record.seq.upper().count("R")
        

#allele frequencies 
Z_allele_freqs = []
A_allele_freqs = []
snps = vcf.Reader(open("/data/boa15kjh/greattit_data/n20snps(no_fixed_diffs).vcf", 'r'))
for record in snps:  
    if record.CHROM == "chrZ":
        AF = record.INFO["AF"] 
        for value in AF:
            Z_allele_freqs.append(value)
    else:
        AF = record.INFO["AF"] 
        for value in AF:
            A_allele_freqs.append(value)


#calculating pi
Z_Pi = pi(20, Z_allele_freqs)
A_Pi = pi(20, A_allele_freqs)


#calculating nucleotide diversity
Z_nucl_diversity = Z_Pi / Z_callable_sites
print("Z nucleotide diversity =", Z_nucl_diversity)

A_nucl_diversity = A_Pi / A_callable_sites
print("Autosome nucleotide diversity =", A_nucl_diversity)




###############################################################################    
############################# Wattersons Theta ################################
###############################################################################


#defining funtion for theta 
def theta(n, num_seg_sites):
    ak = sum(1/i for i in range(1, n))
    theta_W = num_seg_sites / ak
    return theta_W


#number of segregating sites 
Z_seg_sites = []
A_seg_sites = []
snps = vcf.Reader(open("/data/boa15kjh/greattit_data/n20snps(no_fixed_diffs).vcf", 'r'))
for record in snps:  
    if record.CHROM == "chrZ":
        Z_seg_sites.append(record.POS)
    else:
        A_seg_sites.append(record.POS)

num_Z_seg_sites = len(Z_seg_sites)
num_A_seg_sites = len(A_seg_sites)


#calculating theta 
Z_theta = theta(20, num_Z_seg_sites)     
A_theta = theta(20, num_A_seg_sites)         
        

#calculating wattersons theta 
Z_wattersons = Z_theta / Z_callable_sites
print("Z Theta W =", Z_wattersons) 

A_wattersons = A_theta / A_callable_sites
print("Autosome Theta W =", A_wattersons)  




###############################################################################    
################################# Tajimas D ###################################
###############################################################################

#defining function for Tajimas D

def tajima(n, pi, S):
    a1 = sum(1/i for i in range(1, n))
    a2 = sum(1/i**2 for i in range(1, n))
    
    b1 = (n + 1) / (3*(n - 1))
    b2 = (2*(n**2 + n + 3)) / (9*n*(n - 1))
    
    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1*n)) + (a2 / a1**2)
    
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    
    x = math.sqrt((e1*S) + (e2*S*(S - 1)))
    d = pi - (S / a1)
    D = d / x
    return D

tajimaD_Z = tajima(20, Z_Pi, num_Z_seg_sites)
print("Z Tajimas D =", tajimaD_Z)

tajimaD_A = tajima(20, A_Pi, num_A_seg_sites)
print("A Tajimas D =", tajimaD_A)

###############################################################################













