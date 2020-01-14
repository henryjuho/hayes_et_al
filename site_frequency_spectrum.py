# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 15:42:58 2018

@author: Keith
"""

###############################################################################
########################## Site Frequency Spectrum  ###########################
###############################################################################

#installing packages 
import vcf
      
      
#allele frequencies 
unfoldedallelefreqsZ = []
unfoldedallelefreqsA = []
foldedallelefreqsZ = []
foldedallelefreqsA = []

snps = vcf.Reader(open("/data/boa15kjh/greattit_data/n20codingsnps_4fold_nofixeddiff.vcf", 'r'))
for record in snps:  
    if "AA" in record.INFO:
        if record.CHROM == "chrZ":
            AF = record.INFO["AF"]
            AA = record.INFO["AA"]
            REF = record.REF
            if AA == REF:
                for value in AF:
                    unfoldedallelefreqsZ.append(value)
            else:
                for value in AF:
                    rvalue = round((1 - value), 3)
                    unfoldedallelefreqsZ.append(rvalue)
        else:
            AF = record.INFO["AF"]
            AA = record.INFO["AA"]
            REF = record.REF
            if AA == REF:
                for value in AF:
                    unfoldedallelefreqsA.append(value)
            else:
                for value in AF:
                    rvalue = round((1 - value), 3)
                    unfoldedallelefreqsA.append(rvalue)
    else:
        if record.CHROM == "chrZ":
            AF = record.INFO["AF"]   
            for value in AF:
                if value > 0.5:
                    fvalue = round((1 - value), 3)
                    foldedallelefreqsZ.append(fvalue)
                else:
                    foldedallelefreqsZ.append(value)
        else:
            AF = record.INFO["AF"]   
            for value in AF:
                if value > 0.5:
                    fvalue = round((1 - value), 3)
                    foldedallelefreqsA.append(fvalue)
                else:
                    foldedallelefreqsA.append(value)

    
    
#creating dictionaries for sfs    
foldedsfsZ = {}
for number in foldedallelefreqsZ:
    if number not in foldedsfsZ.keys():
        foldedsfsZ[number] = 0 
    foldedsfsZ[number] += 1    
print("folded sfs Z", sorted(foldedsfsZ.items()))

unfoldedsfsZ = {}
for number in unfoldedallelefreqsZ:
    if number not in unfoldedsfsZ.keys():
        unfoldedsfsZ[number] = 0 
    unfoldedsfsZ[number] += 1    
print("unfolded sfs Z", sorted(unfoldedsfsZ.items()))

foldedsfsA = {}
for number in foldedallelefreqsA:
    if number not in foldedsfsA.keys():
        foldedsfsA[number] = 0 
    foldedsfsA[number] += 1    
print("folded sfs A", sorted(foldedsfsA.items()))
         
unfoldedsfsA = {}
for number in unfoldedallelefreqsA:
    if number not in unfoldedsfsA.keys():
        unfoldedsfsA[number] = 0 
    unfoldedsfsA[number] += 1    
print("unfolded sfs A", sorted(unfoldedsfsA.items()))  
    


















