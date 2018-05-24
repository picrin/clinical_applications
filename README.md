1. Abstract

Myotonic Dystrophy type 1 (DM1) is a widely studied genetic disorder, characterised by anticipating muscular dystrophy, myotonia, and other symptoms. We focus on various sets of genes, whose expression has already been known to be disrupted in DM1, in accordance with various proposed mechanisms. We use high-throughput techniques to evaluate these genes as DM1 biomarkers, using genetic material obtained from blood and muscle. We provide evidence, using a simple statistical model based on Partial Least Squares Regression (PLSR), that genes with disrupted Alternative Splicing can function as effective biomarkers, as long as genetic material is obtained from muscle. We expect this work to be directly applicable in the context of evaluating emerging DM1 treatments.

2. Instructions

We advise that all of this code be run on a machine with 64 GB of RAM or more, given that some parts of the pipeline can use up in excess of 32 GB of RAM. We were able to successfuly execute the entire analysis using an AWS "m5.4xlarge" EC2 instance. We found the default amount of storage, 8 GB, to be insufficient to store both the primary data and intermediate computions. We increased the amount of storage to 100 GB. We remove all networking restrictions on the instance, to allow for remote access of jupyter notebooks, which contain our pipline.

We had to apply the following shell commands to set up the machine and 

1. `sudo apt-get update`
2. `sudo apt-get upgrade`
3. `sudo apt-get install python3-pip`
4. `git clone https://github.com/picrin/clinical_applications.git`
5. `pip3 install jupyter`

All data and metadata used in this study are available in publically accessible s3 buckets: `dm1-biomarkers/CEL`, `dm1-biomarkers/annotations` respectively. These data-sets need to be copied to the root of the `clinical_applications` repository as directories `CEL` and `annotations` respectively.

Finally, `jupyter notebook --ip 0.0.0.0 --port 8888` can be issued to start the notebook server, which we can access remotely, using a DNS entry allocated for our EC2 instance and provided that communications on our chosen port is configured to be accessible through the AWS firewall.

We now run notebooks in the following order:

1. `parse_chip_data.ipynb`. This part of the pipeline is responsible for determining probeset ids, sequences of probes and probe coordinates on the chip. It produces an intermediate file with data adhering to the following schema: `probeset, x, y, sequence`.
2. `parse_csv_annotations.ipynb`. Here, we determine "genomic" metadata, i.e. chromosomomal coordinates and strandedness.
3. `unpack_CEL_files.ipynb`. Here we use our own contribution to Biopython to parse the binary CEL v4 file format, which is what all our microarray data uses.
4. `quantile_normalise.ipynb`. Here we perform quantile normalisation of our microarray data.
5. `reannotate_probeset_level`. Here we verify Affymetrix's annotation. We determine that over 1% probes are incorrectly annotated. We discard these probes. We limit our attention to probes, which belong to chromosomes chr1-chr22, X, Y and the mitochondrial DNA (M).
6. `intervaltrees.ipynb`. Here we carry out an exclusive filtering, choosing probes, which are identified by `gencode.v26lift37.annotation.gtf` as having `transcript_type` equal to `"protein_coding`,  or `gene_type` equal to  `"protein_coding"`, as well as the exon type equal to `"CDS"` (Coding sequence) or `"UTR"` (3' or 5' untranslated regions).
7. `merge_annotation_experiment.ipynb`. Here we produce a single file per human gene, as identified in GENCODE v26, with data from all patients for all probes for that gene.
8. `predictions.ipynb`. Here we run our PLSR model to predict MAL from the microarray data.
