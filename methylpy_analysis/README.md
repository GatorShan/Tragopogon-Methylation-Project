# Methylpy analysis -- data processing
## 1. Description
The method is adapted from [here](https://github.com/yupenghe/methylpy/blob/methylpy/tutorial/tutorial.md).

## 2. Analysis of T.dubius (3040-6-2; Pullman)
The raw data contains 615 million paired-end 150-bp reads.

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


