# Transgene builder algorithm 

The gene builder algorithm can be subdivided into four different subroutines, namely, codon optimization, optimization of ribosomal binding, removal of C. elegans piRNA homology and restriction sites, and addition of non-coding sequences. During the first subroutine, the user selects a codon usage table based on the codon frequencies of the top 500 highest expressed genes identified by Serizay et al. (2020). Alternatively, the user can choose to use the highest codon used across all the genes (Max expression, Codon Adaptation Index = 1). Max expression is calculated based on Redemann et al. (2011) and C. elegans Codon Adapter (https://worm.mpi-cbg.de/codons/cgi-bin/optimize.py). Based on the kind of codon usage selected, the program will randomly sample codons with probabilities equal to the frequency usage, i.e. the output sequence should have somewhat similar codon usage to the kind of gene selected. For “Max expression”, solely the codon with the highest expression will be used (no random sampling).  

 

The second subroutine promotes the use of a 5` weak mRNA structure to improve ribosomal binding (for a review see [Tuller and Zur 2015]). To do this, the app produces up to 10000 sequences using the same first 13 amino acids of the input sequence but different codons. Later, trails of four adenines are attached to their 5’ end and their folding energy is calculated with the RNAfold software of the Vienna package [Lorenz et al. 2011]. The sequence with free energy closest to zero is selected for later optimization. 

The third subroutine consists in the removal of piRNA and restriction enzyme target sites. To identify piRNA sites, we perform an alignment of all the possible 20-mers in the optimized sequence against all the piRNA type I and type II reverse complement sequences (21st nucleotides are removed) using the stringdist package of R (Van de loo 2014). We consider alignements from 16 matches as potential piRNA target sites. Once we know which stretches of DNA are potential targets, we introduce inframe synonymous mutations to get rid of them. The search and removal of piRNA sites is performed iteratively up to 100 times in order to remove inadvertent sites created by the addition of the synonymous substitutions. A similar process is used to remove restriction enzymes sites, with the exception that we use instead the “matchPattern” function of the package “Biostrings” (Pages et al. 2019) to find the exact coordinates of the sites to change. 

Finally, we add intronic sequences to the optimized sequences as they are known to improve transgene expression (Okkema et al. 1993). Introns are indicated by lower-case letters and are inserted at the splice consensus sites “AG|R”. The first intron is inserted between the 50th and 150th bp, while the others are inserted approximately 150 bp apart from each other. 

# Link to our website
https://wormbuilder.org/transgenebuilder/
 
# Software availability 

https://github.com/AmhedVargas/GeneOptimizer 
