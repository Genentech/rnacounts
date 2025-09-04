
# Description

A python command line tool which counts variant-containing RNA-seq read pairs in an input BAM file and ascertains transcripts with which those pairs are consistent. 

# Assumptions/assertions

* Reads are processed as pairs. If both mates in a pair overlap the variant and contain the alt allele, they will not be double-counted. However, at this time we are not requiring that both mates contain the variant if they both overlap the variant position.
* Read pairs that are in some way not proper (e.g. mate is unmapped, mate is mapped to a different chromosome) are not considered. If each read in a pair has an `XS` tag, and they disagree as to the strand of the fragment, the pair is also considered improper (even if the flag does not indicate that to be the case).
* Read pairs are defined as ref or alt allele-containing on the *sole* bases of the position(s) of the variant, and structural/sequence consistency with the variant. For instance, if a read as the reference allele at the variant position, but has an undescribed variant (or sequencing error) at a different position, its pair will still currently be described as supporting the reference allele. 
* A variant position-overlapping read's mate must currently be within 200,000 bp of the variant. If it is not found, the read is still processed (provided the pair was proper as indicated by the SAM flag), but its mate cannot be used to help infer transcript consistency. This may result in a small amount of information loss in transcripts with enormous introns.
* Ref allele-containing read pairs for insertions are those that span the reference positions flanking the insertion with no intervening bases.  Alt allele-containing read pairs are those which contain intervening sequence of the appropriate length *and* sequence.
* Ref allele containing read pairs for deletions are those that partially or completely overlap the deletion position(s) with sequence matching that of the reference. Alt allele-containing read pairs are those that span the reference positions flanking the deletion, but lack intervening sequence.
* Reads matching the ref/alt allele of SNVs can be filtered (ignored) on the base quality at the position of the variant. Currently, no quality filtration is available for indels.

# Requirements 

```
python3.10 
pysam==0.23.0
```

# Usage


## As a standalone python script

```
python3 \
  read_checker.py \
  --gtf test_data/rnacounts_test_minimal.gtf \
  --bam test_data/rnacounts_test_input_reads.bam \
  --variants test_data/rnacounts_test_input_variants.tsv \
  --outdir /path/to/output/directory \
  --prefix outfile_prefix \
  --generate_output
```

# Output 

### `outfile_prefix_rna_tally.tsv`

Variant-level counts of variant-overlapping read pairs. The fields `total_ref_counts`, `total_alt_counts`, and `total_counts` map to the RNA-seq `refDepth*`, `altDepth*`, and `totalDepth*` variables currently created by `neoag` (and stored in `*-perspective_rank-genomic_mutations.tab`) among other places. `tx_alt_counts` and refers to the subset of variant-containing read pairs whose structure is consistent with any transcript present in the input GTF with reference to the variant.  Structural consistency is established by checking read pair aligned contigs against transcript exons, and read splice junctions against transcript splice junctions.

| seqnames | start     | end       | ref | alt | total_ref_counts | total_alt_counts | tx_ref_counts | non_tx_ref_counts | tx_alt_counts | non_tx_alt_counts | other_counts | total_counts |
|-------|-----------|-----------|-----|-----|------------------|------------------|---------------|-------------------|---------------|-------------------|--------------|--------------|
| 1     | 15630382  | 15630382  | A   | G   | 59               | 100              | 59            | 0                 | 100           | 0                 | 0            | 159          |
| 1     | 26780218  | 26780220  | GCC | G   | 4                | 18               | 4             | 0                 | 18            | 0                 | 0            | 22           |
| 1     | 64178743  | 64178743  | T   | C   | 0                | 0                | 0             | 0                 | 0             | 0                 | 0            | 0            |
| 1     | 152086898 | 152086898 | T   | C   | 0                | 0                | 0             | 0                 | 0             | 0                 | 0            | 0            |
| 1     | 152215633 | 152215633 | T   | A   | 0                | 0                | 0             | 0                 | 0             | 0                 | 0            | 0            |
| 1     | 155131410 | 155131410 | A   | C   | 66               | 0                | 66            | 0                 | 0             | 0                 | 0            | 66           |
| 1     | 197429578 | 197429578 | G   | A   | 0                | 0                | 0             | 0                 | 0             | 0                 | 0            | 0            |
| 1     | 224296573 | 224296573 | T   | G   | 35               | 4                | 27            | 8                 | 3             | 1                 | 0            | 39           |
| 2     | 11566621  | 11566621  | T   | C   | 0                | 0                | 0             | 0                 | 0             | 0                 | 0            | 0            |
| 2     | 38602619  | 38602619  | G   | A   | 22               | 2                | 22            | 0                 | 2             | 0                 | 0            | 24           |
| 2     | 70230779  | 70230779  | T   | G   | 216              | 0                | 174           | 42                | 0             | 0                 | 0            | 216          |
| 2     | 79086438  | 79086438  | C   | T   | 53               | 0                | 53            | 0                 | 0             | 0                 | 0            | 53           |
| 2     | 79086439  | 79086439  | G   | T   | 54               | 0                | 54            | 0                 | 0             | 0                 | 0            | 54           |
| 2     | 79909792  | 79909792  | A   | C   | 0                | 0                | 0             | 0                 | 0             | 0                 | 0            | 0            |
| 2     | 95285038  | 95285038  | T   | G   | 5                | 0                | 0             | 5                 | 0             | 0                 | 0            | 5            |
| 2     | 96287486  | 96287486  | T   | G   | 352              | 0                | 287           | 65                | 0             | 0                 | 0            | 352          |
| 2     | 127624265 | 127624265 | T   | C   | 79               | 0                | 78            | 1                 | 0             | 0                 | 0            | 79           |
| 2     | 151643831 | 151643831 | T   | G   | 0                | 0                | 0             | 0                 | 0             | 0                 | 0            | 0            |
| 2     | 178776006 | 178776006 | T   | G   | 0                | 0                | 0             | 0                 | 0             | 0                 | 0            | 0            |
| 2     | 201647809 | 201647809 | C   | T   | 0                | 0                | 0             | 0                 | 0             | 0                 | 0            | 0            |
| 2     | 215419289 | 215419289 | C   | T   | 2218             | 0                | 2058          | 160               | 0             | 0                 | 4            | 2222         |
| 2     | 238094504 | 238094504 | C   | T   | 3                | 20               | 2             | 1                 | 17            | 3                 | 0            | 23           |

`seqnames`: chromosome name  
  
`start`: variant start position 
  
`end`: variant end position 
  
`ref`: reference allele  
  
`alt`: variant allele  
  
`total_ref_counts`: count of overlapping read pairs consistent with the reference allele    
  
`total_alt_counts`: count of overlapping read pairs consistent with the variant allele. This should equal the sum of `tx_alt_counts` and `non_tx_alt_counts`  

`tx_ref_counts`: count of overlapping read pairs consistent with the reference allele AND at least one transcript associated to the variant in the input rowranges file  

`non_tx_ref_counts`: count of overlapping read pairs consistent with the reference allele but NOT consistent with ANY transcript associated to the variant in the input rowranges file  
  
`tx_alt_counts`: count of overlapping read pairs consistent with the variant allele AND at least one transcript associated to the variant in the input rowranges file  
  
`non_tx_alt_counts`: count of overlapping read pairs consistent with the variant allele but NOT consistent with any transcript associated to the variant in the input rowranges file  
  
`other_counts`: count of overlapping read pairs consistent with neither the reference nor the given variant allele    
  
`total_counts`: count of all overlapping read pairs (this should equal the sum of `ref_counts`, `total_alt_counts`, and `other_counts`)  
  
  

### `outfile_prefix_transcript_rna_tally.tsv`

Transcript-level counts of variant-containing read pairs.  `consistent_counts` are counts of variant-containing read pairs whose aligned contigs are consistent with the structure of the transcript specified in the `transcript` field.  Structural consistency is established by checking read pair aligned contigs against transcript exons, and read splice junctions against transcript splice junctions.

| seqnames | start     | end       | ref | alt | transcript      | consistent_ref_counts | consistent_alt_counts | exclusive_ref_counts | exclusive_alt_counts |
|-------|-----------|-----------|-----|-----|-----------------|-----------------------|-----------------------|----------------------|----------------------|
| 1     | 15630382  | 15630382  | A   | G   | ENST00000480945 | 59                    | 100                   | 59                   | 100                  |
| 1     | 26780218  | 26780220  | GCC | G   | ENST00000324856 | 4                     | 18                    | 0                    | 0                    |
| 1     | 26780218  | 26780220  | GCC | G   | ENST00000374152 | 4                     | 18                    | 0                    | 0                    |
| 1     | 26780218  | 26780220  | GCC | G   | ENST00000430799 | 4                     | 18                    | 0                    | 0                    |
| 1     | 26780218  | 26780220  | GCC | G   | ENST00000457599 | 4                     | 18                    | 0                    | 0                    |
| 1     | 64178743  | 64178743  | T   | C   | ENST00000371079 | 0                     | 0                     | 0                    | 0                    |
| 1     | 64178743  | 64178743  | T   | C   | ENST00000545203 | 0                     | 0                     | 0                    | 0                    |
| 1     | 152086898 | 152086898 | T   | C   | ENST00000368806 | 0                     | 0                     | 0                    | 0                    |
| 1     | 152215633 | 152215633 | T   | A   | ENST00000368801 | 0                     | 0                     | 0                    | 0                    |
| 1     | 155131410 | 155131410 | A   | C   | ENST00000368406 | 66                    | 0                     | 0                    | 0                    |
| 1     | 155131410 | 155131410 | A   | C   | ENST00000368407 | 66                    | 0                     | 0                    | 0                    |
| 1     | 197429578 | 197429578 | G   | A   | ENST00000367397 | 0                     | 0                     | 0                    | 0                    |
| 1     | 197429578 | 197429578 | G   | A   | ENST00000367399 | 0                     | 0                     | 0                    | 0                    |
| 1     | 197429578 | 197429578 | G   | A   | ENST00000367400 | 0                     | 0                     | 0                    | 0                    |
| 1     | 197429578 | 197429578 | G   | A   | ENST00000535699 | 0                     | 0                     | 0                    | 0                    |
| 1     | 197429578 | 197429578 | G   | A   | ENST00000638467 | 0                     | 0                     | 0                    | 0                    |
| 1     | 224296573 | 224296573 | T   | G   | ENST00000281701 | 25                    | 3                     | 0                    | 0                    |
| 1     | 224296573 | 224296573 | T   | G   | ENST00000340871 | 27                    | 3                     | 0                    | 0                    |
| 1     | 224296573 | 224296573 | T   | G   | ENST00000391875 | 25                    | 3                     | 0                    | 0                    |
| 1     | 224296573 | 224296573 | T   | G   | ENST00000469075 | 25                    | 3                     | 0                    | 0                    |
| 1     | 224296573 | 224296573 | T   | G   | ENST00000469968 | 25                    | 3                     | 0                    | 0                    |
| 1     | 224296573 | 224296573 | T   | G   | ENST00000482491 | 27                    | 3                     | 0                    | 0                    |
| 2     | 11566621  | 11566621  | T   | C   | ENST00000234142 | 0                     | 0                     | 0                    | 0                    |
| 2     | 11566621  | 11566621  | T   | C   | ENST00000263834 | 0                     | 0                     | 0                    | 0                    |
| 2     | 11566621  | 11566621  | T   | C   | ENST00000381483 | 0                     | 0                     | 0                    | 0                    |
| 2     | 11566621  | 11566621  | T   | C   | ENST00000381486 | 0                     | 0                     | 0                    | 0                    |
| 2     | 11566621  | 11566621  | T   | C   | ENST00000389825 | 0                     | 0                     | 0                    | 0                    |
| 2     | 38602619  | 38602619  | G   | A   | ENST00000378915 | 22                    | 2                     | 0                    | 0                    |
| 2     | 38602619  | 38602619  | G   | A   | ENST00000409328 | 22                    | 2                     | 0                    | 0                    |
| 2     | 38602619  | 38602619  | G   | A   | ENST00000449105 | 22                    | 2                     | 0                    | 0                    |
| 2     | 38602619  | 38602619  | G   | A   | ENST00000608859 | 22                    | 2                     | 0                    | 0                    |
| 2     | 70230779  | 70230779  | T   | G   | ENST00000282574 | 159                   | 0                     | 0                    | 0                    |
| 2     | 70230779  | 70230779  | T   | G   | ENST00000415783 | 174                   | 0                     | 0                    | 0                    |
| 2     | 70230779  | 70230779  | T   | G   | ENST00000416149 | 159                   | 0                     | 0                    | 0                    |
| 2     | 70230779  | 70230779  | T   | G   | ENST00000433529 | 159                   | 0                     | 0                    | 0                    |
| 2     | 70230779  | 70230779  | T   | G   | ENST00000445587 | 174                   | 0                     | 0                    | 0                    |


`seqnames`: chromosome name  
  
`start`: variant start position  
  
`end`: variant end position  
  
`ref`: reference allele  
  
`alt`: variant allele  
  
`transcript`: name of the overlapping transcript  
  
`consistent_ref_counts`: count of overlapping reference allele-containing read pairs that are consistent with the structure of the indicated transcript  

`consistent_alt_counts`: count of overlapping variant allele-containing read pairs that are consistent with the structure of the indicated transcript  

`exclusive_ref_counts`: count of overlapping reference allele-containing read pairs that are *exclusively* consistent with the structure of the indicated transcript (with regard only to overlapping transcripts provided in the input variant file).    

`exclusive_ref_counts`: count of overlapping variant allele-containing read pairs that are *exclusively* consistent with the structure of the indicated transcript (with regard only to overlapping transcripts provided in the input variant file).  

### `outfile_prefix_joined_transcript_rna_tally.tsv`  

Provides reference and variant allele counts of read pairs with regard to the whole set of transcripts with which a group of read pairs is consistent. This results from the fact that most read pairs that are consistent with a particular transcript will not overlap a splice junction that is completely unique to it, and consequently may also be consistent with a number of other transcripts.  

| seqnames | start     | end       | ref | alt | transcripts                                                                                                                                                                                     | consistent_ref_counts | consistent_alt_counts |
|-------|-----------|-----------|-----|-----|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|-----------------------|
| 1     | 15630382  | 15630382  | A   | G   | ENST00000480945                                                                                                                                                                                 | 59                    | 100                   |
| 1     | 26780218  | 26780220  | GCC | G   | ENST00000324856,ENST00000374152,ENST00000430799,ENST00000457599                                                                                                                                 | 4                     | 18                    |
| 1     | 155131410 | 155131410 | A   | C   | ENST00000368406,ENST00000368407                                                                                                                                                                 | 66                    | 0                     |
| 1     | 224296573 | 224296573 | T   | G   | ENST00000281701,ENST00000340871,ENST00000391875,ENST00000469075,ENST00000469968,ENST00000482491                                                                                                 | 25                    | 3                     |
| 1     | 224296573 | 224296573 | T   | G   | ENST00000340871,ENST00000482491                                                                                                                                                                 | 2                     | 0                     |
| 2     | 38602619  | 38602619  | G   | A   | ENST00000378915,ENST00000409328,ENST00000449105,ENST00000608859                                                                                                                                 | 22                    | 2                     |
| 2     | 70230779  | 70230779  | T   | G   | ENST00000282574,ENST00000415783,ENST00000416149,ENST00000433529,ENST00000445587                                                                                                                 | 159                   | 0                     |
| 2     | 70230779  | 70230779  | T   | G   | ENST00000415783,ENST00000445587                                                                                                                                                                 | 15                    | 0                     |
| 2     | 79086438  | 79086438  | C   | T   | ENST00000305089                                                                                                                                                                                 | 53                    | 0                     |
| 2     | 79086439  | 79086439  | G   | T   | ENST00000305089                                                                                                                                                                                 | 54                    | 0                     |
| 2     | 96287486  | 96287486  | T   | G   | ENST00000323853                                                                                                                                                                                 | 287                   | 0                     |
| 2     | 127624265 | 127624265 | T   | C   | ENST00000409090,ENST00000409816,ENST00000428314                                                                                                                                                 | 78                    | 0                     |
| 2     | 215419289 | 215419289 | C   | T   | ENST00000323926,ENST00000336916,ENST00000354785,ENST00000356005,ENST00000357867,ENST00000359671,ENST00000421182,ENST00000426059,ENST00000432072,ENST00000443816,ENST00000446046                 | 1998                  | 0                     |
| 2     | 215419289 | 215419289 | C   | T   | ENST00000323926,ENST00000336916,ENST00000354785,ENST00000356005,ENST00000357867,ENST00000359671,ENST00000421182,ENST00000432072,ENST00000443816,ENST00000446046                                 | 60                    | 0                     |
| 2     | 238094504 | 238094504 | C   | T   | ENST00000254663                                                                                                                                                                                 | 2                     | 17                    |
| 3     | 5210182   | 5210182   | A   | T   | ENST00000256497,ENST00000445686                                                                                                                                                                 | 160                   | 0                     |
| 3     | 13354059  | 13354059  | C   | T   | ENST00000254508                                                                                                                                                                                 | 5                     | 0                     |
| 3     | 36846371  | 36846371  | C   | T   | ENST00000429976                                                                                                                                                                                 | 12                    | 0                     |
| 3     | 45946815  | 45946815  | A   | C   | ENST00000304552,ENST00000438735,ENST00000457814,ENST00000458629                                                                                                                                 | 42                    | 0                     |
| 3     | 130999605 | 130999605 | A   | G   | ENST00000328560,ENST00000359644,ENST00000422190,ENST00000428331,ENST00000504381,ENST00000504948,ENST00000505330,ENST00000507488,ENST00000508532,ENST00000510168,ENST00000513801,ENST00000533801 | 40                    | 27                    |
| 3     | 130999605 | 130999605 | A   | G   | ENST00000359644,ENST00000422190,ENST00000428331,ENST00000504381,ENST00000504948,ENST00000505330,ENST00000507488,ENST00000508532,ENST00000510168,ENST00000513801                                 | 185                   | 100                   |
| 3     | 132528223 | 132528223 | A   | G   | ENST00000260818                                                                                                                                                                                 | 41                    | 8                     |
| 4     | 8231997   | 8231997   | G   | A   | ENST00000245105                                                                                                                                                                                 | 3                     | 0                     |
| 4     | 146325995 | 146325995 | C   | T   | ENST00000335472,ENST00000507030                                                                                                                                                                 | 0                     | 4                     |
| 2     | 11566621  | 11566621  | T   | C   | ENST00000381483                                                                                                                                                                                 | 0                     | 0                     |
| 2     | 11566621  | 11566621  | T   | C   | ENST00000381486                                                                                                                                                                                 | 0                     | 0                     |
| 2     | 11566621  | 11566621  | T   | C   | ENST00000389825                                                                                                                                                                                 | 0                     | 0                     |
| 2     | 38602619  | 38602619  | G   | A   | ENST00000378915                                                                                                                                                                                 | 22                    | 2                     |
| 2     | 38602619  | 38602619  | G   | A   | ENST00000409328                                                                                                                                                                                 | 22                    | 2                     |  


`seqnames`: chromosome name  
  
`start`: variant start position  
  
`end`: variant end position  
  
`ref`: reference allele  
  
`alt`: variant allele  
  
`transcripts`: comma separated list of transcripts with which the counted group of read pairs is consistent.  
  
`consistent_ref_counts`: count of overlapping reference allele-containing read pairs that are consistent with the structure of the indicated group of transcripts. Note that implicitly the these read pairs are inconsistent with transcripts not in this group (among transcripts provided in the input file that were associated to the variant in question).  

`consistent_alt_counts`: count of overlapping variant allele-containing read pairs that are consistent with the structure of the indicated group of transcripts. Note that implicitly the these read pairs are inconsistent with transcripts not in this group (among transcripts provided in the input file that were associated to the variant in question).  
