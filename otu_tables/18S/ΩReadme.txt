2024-10-24 BL

Notes:
- Eventually, move the UNIFRAC distance tables to clean_data, since they are used in downstream analyses
- Eventually, clean the files in this directory that aren't needed


Files:

18S_otu_raw.txt
- The original sequence abundance file output from QIIME

18S_otu_taxonomy.txt
- The original metadata file output from QIIME

18S_UNweighted_Unifrac.tsv
- Original presence-absence based distance table, not used in my analysis and can probably be deleted later

18S_weighted_samples_UNIFRAC.tsv
- Includes all samples in composing the distance table, and the one that should be used in sample-based work from now on. 

18S_weighted_samples_unifrac_CLEAN.csv
- This file includes only those samples which were kept during the previous strategy of rarefying data and removing samples to balance their number across fields. This strategy wastes too much data and is unnecessary. This file is kept for now but will probably be deleted. 
- File was produced by filtering 18@_weighted_samples_UNIFRAC.tsv

18S_weighted_Unifrac.tsv
- Field averages produced under the old strategy of balancing samples in fields. This file will likely be deleted. Suggest renaming to 18S_prerarefy_weighted_Unifrac.tsv.
- The new field-averaged file will need a good name. These two will be hard to tell apart. 

18S_samps_4unifrac.tsv
- Field averages produced from all samples for export to LB and processing into a UNIFRAC distance matrix

18S_avg_4unifrac.tsv
- Sample data for export and UNIFRAC processing...

18S_sequences.fasta
- Phylogenetic information for UNIFRAC distance computation