# Transgene builder app
Shiny app that aids the construction of gene sequences with optimized expression for <i>C. elegans</i> and other nematodes.

## General description
This tool allows any user to input coding sequences and obtain transgenes optimized for <i>C. elegans expression</i>. We have incorporated different methods for codon optimization that ultimately increase the germline expression of transgenes in C. elegans. Therefore, please make sure to install all the dependencies (se further below) in case you want to run this shiny app locally.

## Structure of the repository
The shiny app core lies onto two R scripts, the "server.R" and "ui.R" codes, and a bin folder that contains scripts and programs compiled in a Unix environment. In order to execute properly this app, we recommend its installation into somewhat recent Ubuntu distributions (16.0+) or in MacOS. Additionally, there are two extra folders, one containing text files which the app reads for its proper functioning, and another that contains html elements which formats the aesthetics and visualization of the app (more information).

**Basic processing diagram**

ui.R: Loads data and a place holder for dynamic processing -> server.R: Produces a dynamic UI that changes according to the parameters set by the user and send processes to the background thanks to the future package -> ui.R: Show results

## Dependencies
### **R**

Your R environment must have the following libraries installed:

* shiny

This library is core for the app to work as it makes R into a responsive proxy-like server that can be interacted with via a web browser.

* shinythemes
 
Aesthetics to the user interface

* shinyWidgets
 
Cooler widgets than the base shiny library

* Biostrings
 
DNA string manipulation

* Cairo
 
Used for additional font support

* stringdist
 
Set of functions compiled in C that provides support to the hamming distance calculation used for piRNA search

* shinyjs
 
Makes shiny responsivle to javascript messages

* DT
 
Dataframes and Data tables manipulation

* promises and future packages

These last two packages are essential to the app as makes it responsive between multiple user loads. If you're looking of a single user version to this app, do not hesitate in contacting [me](mailto:avargas0lcg@gmail.com).

**Other programs**
### **perl**

perl is required within the commandline path running the shnny app instance. Additionally, the following libraries are required:

* BerkeleyDB
* Bio::Seq
* Math::Random

Particularly, perl is used to run the [GLO algorithm](http://104.131.81.59/) developed by Dan Dickinson. You can find its documentation [here](https://github.com/dickinson-lab/germline). Please note the file `sequence_lib_scores.db` containing scores for each 12-mer word is required to run the GLO algorithm but Github file size limitations prevent its upload. This file can be download following this [link](https://s3.eu-central-1.amazonaws.com/wormbuilder.dev/Downloads/sequence_lib_scores.db) and has to be placed on the bin folder of this repository for proper functioning.

### RNAfold

We provide a compiled version of RNA fold from the Vienna package. However, please cite and compile the original [code](https://www.tbi.univie.ac.at/RNA/) if RBS optimization is performed. 

### SequenceViewer

The code for neXtProt sequence viewer has been adapted from this [repository](https://github.com/calipho-sib/sequence-viewer). Please cite their work if the analytical mode is used.

## Deployment and implementation
You can see the app in action by going to [here](https://wormbuilder.org/transgenebuilder/).

To run the app locally, make sure you have installed R along with all its dependences. Follow by cloning this repository:

`git clone https://github.com/AmhedVargas/GeneOptimizer_2022`

Run the program via command line specifying an open port, e.g., 5100, and open a web browser to access the app.

`R -e "shiny::runApp('GeneOptimizer_2022',host="0.0.0.0",port=5100)"`

Alternatively, you can run the app in a graphical environment such as Rstudio.

## Usage
**Sequence adaptation**

## Transgene builder algorithm 

The gene builder algorithm can be subdivided into four different subroutines, namely, codon optimization, optimization of ribosomal binding, removal of C. elegans piRNA homology and restriction sites, and addition of non-coding sequences. During the first subroutine, the user selects a codon usage table based on the codon frequencies of the top 500 highest expressed genes identified by Serizay et al. (2020). Alternatively, the user can choose to use the highest codon used across all the genes (Max expression, Codon Adaptation Index = 1). Max expression is calculated based on Redemann et al. (2011) and C. elegans Codon Adapter (https://worm.mpi-cbg.de/codons/cgi-bin/optimize.py). Based on the kind of codon usage selected, the program will randomly sample codons with probabilities equal to the frequency usage, i.e. the output sequence should have somewhat similar codon usage to the kind of gene selected. For “Max expression”, solely the codon with the highest expression will be used (no random sampling).  

The second subroutine promotes the use of a 5` weak mRNA structure to improve ribosomal binding (for a review see [Tuller and Zur 2015]). To do this, the app produces up to 10000 sequences using the same first 13 amino acids of the input sequence but different codons. Later, trails of four adenines are attached to their 5’ end and their folding energy is calculated with the RNAfold software of the Vienna package [Lorenz et al. 2011]. The sequence with free energy closest to zero is selected for later optimization. 

The third subroutine consists in the removal of piRNA and restriction enzyme target sites. To identify piRNA sites, we perform an alignment of all the possible 20-mers in the optimized sequence against all the piRNA type I and type II reverse complement sequences (21st nucleotides are removed) using the stringdist package of R (Van de loo 2014). We consider alignements from 16 matches as potential piRNA target sites. Once we know which stretches of DNA are potential targets, we introduce inframe synonymous mutations to get rid of them. The search and removal of piRNA sites is performed iteratively up to 100 times in order to remove inadvertent sites created by the addition of the synonymous substitutions. A similar process is used to remove restriction enzymes sites, with the exception that we use instead the “matchPattern” function of the package “Biostrings” (Pages et al. 2019) to find the exact coordinates of the sites to change. 

Finally, we add intronic sequences to the optimized sequences as they are known to improve transgene expression (Okkema et al. 1993). Introns are indicated by lower-case letters and are inserted at the splice consensus sites “AG|R”. The first intron is inserted between the 50th and 150th bp, while the others are inserted approximately 150 bp apart from each other. 

## Troubleshoot

Please feel free to [e-mail me](mailto:amhed.velazquez@kaust.edu.sa) for any question, doubt or error in the app.

