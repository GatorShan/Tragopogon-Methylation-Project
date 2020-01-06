# Methylpy analysis -- data processing
## 1. Description
The method is adapted from [here](https://github.com/yupenghe/methylpy/blob/methylpy/tutorial/tutorial.md).

## 2. Analysis of T.dubius (3040-6-2; Pullman)
The raw data contains ~615 million paired-end 150-bp reads.

### 2.1
Script `methylpy_dataprocessing_Tdu_3040-6-2_v1.sh` was used.
  - The job ended after running ~4 days
  - Error message below:
  ```
  Exception in thread "main" htsjdk.samtools.util.RuntimeIOException: java.io.IOException: No space left on device
        at htsjdk.samtools.util.SortingCollection.spillToDisk(SortingCollection.java:247)
        at htsjdk.samtools.util.SortingCollection.add(SortingCollection.java:167)
        at picard.sam.markduplicates.MarkDuplicates.buildSortedReadEndLists(MarkDuplicates.java:524)
        at picard.sam.markduplicates.MarkDuplicates.doWork(MarkDuplicates.java:232)
        at picard.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:282)
        at picard.cmdline.PicardCommandLine.instanceMain(PicardCommandLine.java:103)
        at picard.cmdline.PicardCommandLine.main(PicardCommandLine.java:113)
  Caused by: java.io.IOException: No space left on device
  ```
Output:
  - `methylpy_dataprocessing_Tdu_3040-6-2_v1_45540539.out`
  - `methylpy_dataprocessing_Tdu_3040-6-2_v1_45540539.error`
  - **`T.dubius_3040-6-2_libA_processed_reads.bam`**, which is 19 G

Mapping results:
```
615057146 reads; of these:
  615057146 (100.00%) were paired; of these:
    207932383 (33.81%) aligned concordantly 0 times
    91810217 (14.93%) aligned concordantly exactly 1 time
    315314546 (51.27%) aligned concordantly >1 times
66.19% overall alignment rate
615057146 reads; of these:
  615057146 (100.00%) were paired; of these:
    209508325 (34.06%) aligned concordantly 0 times
    91184678 (14.83%) aligned concordantly exactly 1 time
    314364143 (51.11%) aligned concordantly >1 times
65.94% overall alignment rate
```
### 2.2
Script `methylpy_dataprocessing_Tdu_3040-6-2_v1_Picard.sh` was used for mark duplicate analysis.

Output:
  - `methylpy_dataprocessing_Tdu_3040-6-2_v1_Picard_45599907.error`
    - `MarkDuplicates	Marking 93,694,590 records as duplicates`
  - **`T.dubius_3040-6-2_libA_processed_reads_no_clonal.bam`**, which is 15 G
  - `T.dubius_3040-6-2_libA.metric`

## 3. Analysis of T. pratensis (3058-1-2; Garfield)
The raw data contains ~259 million paired-end 150-bp reads. Script `methylpy_dataprocessing_Tpr_3058-1-2_v1.sh` was used.

Output:
  - `methylpy_dataprocessing_Tpr_3058-1-2_v1_45430393.out`
  - `methylpy_dataprocessing_Tpr_3058-1-2_v1_45430393.error`
  - `T.pratensis_3058-1-2_processed_reads_no_clonal.bam`
  - `T.pratensis_3058-1-2_libA.metric`
  - `allc_T.pratensis_3058-1-2.tsv.gz.idx`
  - `allc_T.pratensis_3058-1-2.tsv.gz`

```
There are 258651069 total input read pairs
Wed Dec 25 08:40:05 2019

There are 15891937 uniquely mapping read pairs, 6.14416057179 percent remaining
Wed Dec 25 08:40:05 2019

There are 13360495 non-clonal read pairs, 5.16545129763 percent remaining
Wed Dec 25 08:40:05 2019
...
The non-conversion rate is 0.242098066597%
Wed Dec 25 10:05:32 2019
...
Done
```


