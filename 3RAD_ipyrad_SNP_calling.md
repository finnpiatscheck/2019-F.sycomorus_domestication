# Processing of RAD sequences and SNP calling with ipyrad 0.9.90

ipyrad is a conda package for assembly, analysis and SNP calling of RAD-seq data. This package is here used to call SNPs for *Ficus petiolaris* sequences from Mexico. This work is a collaboration between, Guillermo Ibarra Manríquez, Eva Piedra Malegon, Antonio Gonzaled, Finn Piatscheck and John Nason.

--- 

All the following commands lack information about paralleling and HPC core uses. For the following, I have used 16 cores with the arguments ```-c 16 --MPI```.

---

## Data

```Ficus_3RAD_1_R1_.fq.gz  
Ficus_3RAD_1_R2_.fq.gz  
Ficus_3RAD_2_R1_.fq.gz  
Ficus_3RAD_2_R2_.fq.gz  
Ficus_3RAD_3_R1_.fq.gz  
Ficus_3RAD_3_R2_.fq.gz  
Jose_L_Ficus_P1_R1_.fq.gz  
Jose_L_Ficus_P1_R2_.fq.gz  
Jose_L_Ficus_P2_R1_.fq.gz  
Jose_L_Ficus_P2_R2_.fq.gz  
Jose_L_Ficus_P3_R1_.fq.gz  
Jose_L_Ficus_P3_R2_.fq.gz  
```
<i>NOTES: UNAM sequences are named "Ficus_3RAD_x_Rx_.fq.gz" and Nason lab sequences are named "Jose_L_Ficus_Px.x.fq.gz" and the lattest were renamed to fit ipyrad's file name expectations.</i>

## Create parameter files

``` 
ipyrad -n params-demux-Eva-IDX1.txt
```
Edit the file with ```vim``` to the following:

```
------- ipyrad params file (v.0.9.90)-------------------------------------------
demux-Eva-IDX1                     ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
./                             ## [1] [project_dir]: Project dir (made in curdir if not present)
/home/ECOFOG/finn.piatscheck/work/3RAD/Ficus_3RAD_1_R*_.fq.gz   ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
/home/ECOFOG/finn.piatscheck/work/3RAD/ipyrad_analysis/3RAD-Eva-barcodes.txt    ## [3] [barcodes_path]: Location of barcodes file
                               ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
denovo                         ## [5] [assembly_method]: Assembly method (denovo, reference)
                               ## [6] [reference_sequence]: Location of reference sequence file
pair3rad                       ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
CGG, GATCC                     ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                              ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                             ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
6                              ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                              ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                          ## [13] [maxdepth]: Max cluster depth within samples
0.88                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
1                              ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
2                              ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
30                             ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
0.05                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus
0.05                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus
45                              ## [21] [min_samples_locus]: Min # samples per locus for output
0.2                            ## [22] [max_SNPs_locus]: Max # SNPs per locus
5                              ## [23] [max_Indels_locus]: Max # of indels per locus
0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus
0, 0, 0, 0                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
G, a, g, k, m, n, p, s, l, u, t, v                        ## [27] [output_formats]: Output formats (see docs)
                               ## [28] [pop_assign_file]: Path to population assignment file
                               ## [29] [reference_as_filter]: Reads mapped to this reference are removed in step 3
```
Create one parameter file for each fastQ file to be demultiplexed. A total of 6 parameter files should be create.

Only parameters 1,2,3,4,7,8 and 15 are used here. After merging, other parameters can be changed.

## Demultiplex and merge libraries (Step 1)

### Demultiplex
```
ipyrad -p params-demux-Eva-IDX1.txt -s 1
ipyrad -p params-demux-Eva-IDX2.txt -s 1
ipyrad -p params-demux-Eva-IDX3.txt -s 1
ipyrad -p params-demux-Nason-IDX1.txt -s 1
ipyrad -p params-demux-Nason-IDX2.txt -s 1
ipyrad -p params-demux-Nason-IDX3.txt -s 1
```

### Merge
```
ipyrad -m 3RAD-combined params-demux-Eva-IDX1.txt params-demux-Eva-IDX2.txt params-demux-Eva-IDX3.txt params-demux-Nason-IDX1.txt params-demux-Nason-IDX2.txt params-demux-Nason-IDX3.txt
```

## Search for optimal parameters by branching

For the next steps, we will try different sets of parameters to get the best number of SNPs we can for downstream analyses.

The command to use is ```ipyrad -p parameter_file -b new parameter_file name```. 

It's also important to note that if we want to subset samples fom the libraries, we will similarly use the banching option and adding individual names as arguments. For this analysis, we will not subset any samples.

## Read filtering (Step 2)
```
ipyrad -p params-3RAD-combined.txt -s 2
```
For this step, the parameters to be used are: 0,1,3,7,8,9,16 and 17 and 25.

Read the summary in s2_rawedit_stats.txt with the command ```less```. We learn that a large proportion of reads is retained. Most reads are removed because of read 2 qualiy (this is expected). Reads removed by minimum length of 100bp are not too numerous but still, we could consider reducing the minimum length of reads retained to 50bp.

Next we create a branch with trimmed reads. An inspection of the FastQC outputs from the 6 raw reads FASTQ files indicate:
- reads are overall of good quality
- we could trim both R1 and R2 of 21 nucleotides at the end
- Nason lab sequences have the first few nucleotides of lower qualities (but these should be removed after demultiplexing)

*NOTES: Use FIGARO to estimate truncation length? I don't think it is possible, the program is not designed for this, it takes into account amplicon length for both reads to overlap. In RAD, amplicon length is variable.*

```
ipyrad -p params-3RAD-combined.txt -b params-3RAD-s2_trim_20
```
Then, we use ```vim``` to edit the new parameter file, and we change the parameter 25 to: "0, 130, 0, 130". Also, barcode and illumina primers removal will reduce the read length of 28 nucleotides max (length of barcodes varis), so we need to change parameter 17 (minimum length of read after trimming). We set it to a somewhat arbitrary value of 50bp.

```
ipyrad -p params-3RAD-s2_trim_20.txt -s 2 
```

If filtering was already performed, ipyrad will skip samples with reads already filtered even if trimming parameters are changed. So make ipyrad do the filtering with the new parameter file, use the argument ```-f```. This also apply to later steps.

## Cultering/Mapping reads within samples and alignment (Step 3)

We create several branches with clustering thrsholds of 0.88, 0.92 and 0.96. A threshold of 0.88 is considered low since we are working within the same species. We will compair the statistics after step 3.

```
ipyrad -p params-3RAD-combined.txt -b params-3RAD-s3_clust_th88
```
```
ipyrad -p params-3RAD-combined.txt -b params-3RAD-s3_clust_th92
```
```
ipyrad -p params-3RAD-combined.txt -b params-3RAD-s3_clust_th96
```
Do not forget to edit the parameter files with ```vim``` to change the cluster threshold values.

*NOTES: ipyrad doesn't like names with "." in them, it's better to only use "_" and "-" in branch names.*

### Run step 3 
Here only the cluster threshold is changing. So we change the parameter 14 in the three parameter files created. Other parameters necessary for this step, 0,1,5,7, remain unchanged.
```
ipyrad -p params-3RAD-s3_clust_th88.txt -s 3
ipyrad -p params-3RAD-s3_clust_th92.txt -s 3
ipyrad -p params-3RAD-s3_clust_th96.txt -s 3
```

*NOTES: On 32 cores, the process takes ~12 hours. ~14-21 hours on 16 cores depending on parameters*

### Examine the outputs

The statistics of Step 3 are in the file s3 ....

## Heterozygosity and error rate (Step 4)

This still will be important because, as specify on ipyrad's guidelines, eror rate should be on the order of 0.001 and heterozygosity should be on the order of 0.01. If values differ largely, the cluster threshold might have been misspecified.

For this step the parameters to be considered are parameters 11,12,13 and 18. We will leave parameter 13 as the default of 10,000. Parameters 11 and 12 default of 6 seam also reasonable. Parameter 18 remains at 2.

```
ipyrad -p params-3RAD-s3_clust_thxx.txt -s 4
```
We run step 4 for all branches so we can guess which parameters were best.

We get :
Cluster threshold 0.88
Mean error  [0.01057 sd=0.00988]
Mean hetero [0.02551 sd=0.01411]

Cluster threshold 0.92
Mean error  [0.00879 sd=0.00784]
Mean hetero [0.01999 sd=0.00992]

Cluster threshold 0.96
Mean error  [0.01057 sd=0.00988]
Mean hetero [0.02551 sd=0.01411]

As indicated in the ipyrad directions, error rate should be on the order of 0.001 and heterozygosity should be on the order of 0.01