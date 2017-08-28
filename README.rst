================================================
Studying the methylome of <insert cool organism>
================================================

This pipeline relies on Bowtie2 for mapping, Bismark to do methylation calls, and my scripts to process the Bismark results, so get them installed before you follow the pipeline.

Bowtie2: preferably installed via package manager ("apt install bowtie2").
Bismark: https://www.bioinformatics.babraham.ac.uk/projects/bismark/, "Download Now".

Remember to cite both tools when you write your manuscript.

python3: preferably installed via package manager ("apt install python3"). DO NOT USE PYTHON 2!!!

Pipeline, TLDR
--------------
1. Trim FASTQ files prior to mapping.
2. Use ``Bismark`` + ``Bowtie2`` for mapping.
3. Run ``deduplicate_bismark`` to remove duplicated reads (PCR artefacts).
4. Run ``bismark_methylation_extractor`` to produce .cov files used in downstream analysis.
5. Filter for methylated positions that are unlikely to arise due to chance.
6. Analysis!

Trimming
--------
Most people use Trimmomatic for this step... not me. I prefer cutadapt for finer control over trimming sequences.

Over the years, I wrote a shell script (``cutadapt_hiseq.sh``) to automate this step. Justifications for non-default cutoffs:
- -O 10: the default of -O 3 was too... easy to trim off, there's a 1/64 chance that the end of any read would look like the first three bases of the adapter. I raised it to 10 through trial-and-error: I observed that the number of reads with 10+ bp of adapter sequence is several orders above the expected number by chance.
- -q 20,20: trims both ends of the reads if they're Phred < 20.
- --trim-n: removes flanking Ns at the ends of sequences.

Anything reasonable is fine, really. This step won't dictate the end result.

Example commands

Typical usage: ``$ cutadapt-hiseq.sh ATCATG blah_R1.gz blah_R2.gz``
Ultra lazy mode when barcode is in filename: ``$ for a in *R1.fastq.gz; do b=`echo $a | grep -Po '_\w{6}_' | sed 's/_//g'`; cutadapt-hiseq.sh $b $a ${a/R1/R2}; done``

Mapping
-------
Here's the gospel that gets updated periodically:
https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf

If my advice clashes with his, trust him. RTFM.

I'm very used to writing loops to map everything in a folder, but it makes for terrible reading (exemplified by the "ultra lazy" command in the previous section). Adapt the "typical usage" situations to your own data. If you're unable to, then I guess you... can't be lazy!

1. Genome preparation

This command assumes that bowtie2 can be called directly on the command line. If it doesn't work, then look into the ``--path_to_bowtie`` command.

Genome file has sit in a folder of its own. Use symlinks within the folder to save space.

``$ ~/tools/bismark_v0.17.0/bismark_genome_preparation --verbose --bowtie2 .``

2. Running Bismark

``$ ~/tools/bismark_v0.17.0/bismark relative/path/to/prepared/genome --bowtie2 --score_min L,0,-0.6 -1 blah_R1.trim.gz -2 blah_R2.trim.gz -p 12 --bam``

Of all the commands here, the ``--score_min`` stands out as the unexplainable non-default. This exact command was suggested by the program creator in one of the forum posts when someone complained they obtained very low mapping rates for their data. This equation controls how "loose" the mapping should be--and from testing, I can vouch that this cutoff strikes a good balance between sensitivity and accuracy. ``--score_min`` is a Bowtie2 parameter, read its docs if you want to know what it means EXACTLY.

In this pipeline, ``--score_min`` is the BIGGEST determinant of the results. This is the one thing I'll optimise if my results look weird.

Deduplication
-------------
For MOST cases, you want to deduplicate your reads to get rid of PCR duplicates. The only notable exception is if you're carrying out amplicon sequencing, which will produce reads that map to the same loci, and thus be deemed duplicate.

To be pedantic, this script works by checking whether the start and the end coordinate of the mapped read coincides with any other read. If it does, it gets discarded. Other deduplication methods CAN BE different, don't assume all deduplication methods would produce the same end product.

``$ ~/tools/bismark_v0.17.0/deduplicate_bismark --bam -p blah_R1.trim_bismark_bt2_pe.bam``

Methylation calling
-------------------
Based on the deduplicated bam (or not), call # of methylated reads and unmethylated reads on a per-position basis.

``$ ~/tools/bismark_v0.17.0/bismark_methylation_extractor -p --bedGraph --counts --scaffolds --multicore 8 blah_R1.trim_bismark_bt2_pe.deduplicated.bam``

(the filenames are intentional: bismark automatically generates the name of the output files, so it gets longer and longer and longer...)

Note that the reason I use ``--bedGraph --counts`` is to get the file that ends with ``*.cov``.

This stage produces a lot of intermediate files--I keep the results (*.cov) and log files (*.png, *.M-bias.txt, *_splitting_report.txt), and delete the rest (C*.txt, *.bed). You'll get what I mean once you run it once.

Why keep *.cov and discard *.bed? Well, that's because my scripts work on the former but not the latter. Design choice, unfortunately.

From this point onwards, what you're reading is no longer Bismark--it's my stuff, with my cutoffs, with my ideas. Alter these things to taste.

Filtering for bona fide methylation
-----------------------------------
I'll talk about the theory first, then the implementation. This process is a bit convoluted, but the basic idea of a methylated position is one that:

I. Is well-covered across all treatments.
II. Is present in all treatments (similar to (1)).
III. Is methylated in all replicates of a biologically meaningful treatment.
IV. When pooled, is significantly methylated.

A more precise description of the cutoffs are described in ``filter_pos.four_criteria.py``.

What is "significant methylation"? One must first understand that there are two sources of errors that will cause the wrong methylation call.

a) Sequencing error (as reflected by the Phred score of the base). Phred 20 is 1%, Phred 30 is 0.1%.
b) Non-conversion of the unmethylated cytosine (i.e. bisulphite treatment was suppose to convert C-->T but the chemical didn't work as expected, thus the unmethylated base appears methylated). Occurs at the rate of 0.1 to 1%.

In my work, I chose an extremely conservative error rate of 1%. I prefer to deal with fewer real stuff than more stuff that might not be real.

Given a composition of x methylated and y non-methylated reads at a certain position, one can calculate the probability of the observation arising purely by chance. I wrote a script (``filter_miscalled_Cs.py``) to apply binomial theorem and calculate P(X >= x); where P(X = x) = (x+y choose x) * 0.01^x * 0.99^y, and correct the P value with B-H. This script is run AFTER filtering the original .cov files for positions that have AT LEAST one methylation event (it's pointless to deal with the 95% CpGs that are not methylated--Benjamini-Hochberg correction is carried out on this filtered subset, not the entire dataset).

For the implementation steps, the input filenames can be changed to your files of interest, but the output filenames are mandatory--``filter_pos.four_criteria.py`` uses a lot of hardcoded filenames. Using blah1.cov, blah2.cov, blah3.cov as generic inputs, run these commands in the same directory as the files:

1. Run ``tabulate_tsvs.py`` to merge the Bismark cov files into a giant table.
tabulate_tsvs.py blah1.cov blah2.cov blah3.cov -k 0 1 -c 4 5 -v > compiled_coverage.pre_filt.meth_unmeth.tsv

2. gzip-compress this giant file.
``gzip compiled_coverage.pre_filt.meth_unmeth.tsv``

3. Merge all .cov files produced by Bismark into one giant file.
``merge_bismark_cov.py blah1.cov blah2.cov blah3.cov > all.merged.cov``

4. Run ``filter_miscalled_Cs.py`` on this merged file for Criteria IV.
``filter_miscalled_Cs.py all.merged.cov > all.bona_fide_meth_pos.cov``

5. Edit lines 71--89 of ``filter_pos.four_criteria.py`` to specify which columns are replicates of a meaningful biological condition (criteria III). To disable this, just delete these lines. Column numbering starts from 0. My comments in the script and the lines of code match up well, you should be able to figure out how to modify the script even if you do not write Python.

6. Save the script, and run ``filter_pos.four_criteria.py``.

The script picks out significant positions in all of the *.cov files, producing a *.filt.cov file per *.cov file fed into the script.

Annotation of methylated positions
----------------------------------
SANITY CHECK: *.filt.cov should all have the same number of lines.
``wc -l *.filt.cov``

DO NOT PASS GO, DO NOT COLLECT $200 IF THIS IS NOT SATISFIED.

Collected your $200? Great. I suggest merging the filtered files first:
``merge_bismark_cov.py blah1.filt.cov blah2.filt.cov blah3.filt.cov > all.merged.filt.cov``

Then annotate it using (replace <generic_filenames> with real stuff):
``annotate_bismark_cov.py <genome_fasta_file.fa> <gff_file.gff3> all.merged.filt.cov > all.merged.filt.annot.cov``

It is VERY LIKELY that this won't work for you, because ``annotate_bismark_cov.py`` uses a self-written gff3 parser that works with the genomes that we assembled. You'll need to read and understand the code of ``parse_gff3.py`` and ``annotate_bismark_cov.py`` if errors appear. I'm afraid the further you swim away from the safe Bismark shores, the more rocks you'll hit!

If things work, great! As all files have the same number of lines and the same positions in the same order (trust me), you can do some magic to annotate all your individual files:
``cut -f 7- all.merged.filt.annot.cov > tmp``
``paste blah1.filt.cov tmp > blah1.filt.annot.cov``

What's next?
------------
Well, this place is a good point to let your hand go. With the *.filt.annot.cov files, you can do a lot of wonderful stuff. ``head`` or ``less`` the files to see what's inside them. If you don't understand which column stores what information, read the scripts that produced them. I (mostly) documented their functions as comments that precede the script itself.

A few analysis suggestions below:

1. PCA of all *.filt.annot.cov files to see whether related replicates have more similar methylation patterns?

2. Check genomic distribution of methylated positions using *.filt.annot.cov (are there more methylated positions in genic regions? More in exonic regions? Start of exonic regions?)

3. Start thinking about how to compare replicates to obtain differentially expressed genes/regions/etc. I'll probably discuss this in a separate document--there's so many things I could talk about this!
