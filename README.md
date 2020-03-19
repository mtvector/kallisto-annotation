# kallisto-annotation
**Scripts to build kallisto k-mer index databases for new species using CAT gtfs.**

Basic  kallisto indices for model organisms can be downloaded from https://github.com/pachterlab/kallisto-transcriptome-indices/releases. For those not available, I provide this repo to help build kallisto databases. Generates both cDNA and cDNA+intron indexes for transcriptome and RNA velocity analysis from kallisto[|bustools]. Elements of the pipeline are adapted from https://www.kallistobus.tools/velocity_index_tutorial.html.  

In this script introns are defined as the complement of Exons and Intergenic space, which is the simplest and most conservative definition.

CreateKallistoReference.sh is the main script to generate the kallisto reference. It requires the python scripts in this repo, which replace elements from the awk code on the kallisto tutorial, which didn't work on all gtfs. It also requires you make a new conda environment kallisto and install a number of packages before running the script:

conda create -n kallisto
conda activate kallisto
conda install -c bioconda bustools
conda install -c bioconda kallisto
conda install -c bioconda pybedtools
pip install pandas 

CreateKallistoReferenceRhemac10.sh is an example SGE runner script for CreateKallistoReference.sh.

Scripts have currently only been tested by one person on one system for several transcriptomes. Performance may vary until further feedback facilitates refinement!
