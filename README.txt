

---------------
FILES INCLUDED
---------------

- .xlsx file - Example experimental design for this experiment. Example experiment has 2 samples, each using the same 3 hashtags to multiplex three experimental conditions - WT, KO, and control cells. 

- *_library.csv - library file that gives the location of the GEX and FB libraries for each sample. Each sample needs a separate _library.csv file to run the pipeline. Two library files are included for this example experiment.
 
- hash_FB.csv - Feature reference file, if the same barcodes are used in multiplexing each sample only one feature reference file is needed. If separate barcodes are used between samples each will need its own feature reference file. 

- run.sh - example bash script to align the CITE-seq data using cell ranger. 

- hashtag.R - example R script to run Seurat demultiplexing with hashtag oligo commands and plot QC graphs. It is currently set up for the example experiment that has 2 samples that use 3 barcodes each.
 
-------------
STEPS TO RUN
-------------

1. Edit library, feature reference, bash script, and R script file(s) to reflect the experimental design, file locations, and hashtags used in your experiment. 

2. Run the run.sh bash script to perform the cellranger alignment of gene expression (GEX) and feature barcode (FB) data.

3. Run the demultiplexing R script to generate QC plots and demultiplex the hashed data.