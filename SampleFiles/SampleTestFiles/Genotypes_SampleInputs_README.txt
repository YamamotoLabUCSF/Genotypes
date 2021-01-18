** Test files for Genotypes.py **
01-03-2021
This directory contains files compatible with test runs of Genotypes.py


For the following user-entered variables, available test files can be used for the indicated data inputs (underscored). Other inputs must be user-defined (e.g., local output directory, path to BLASTN executable, path to BLASTDBCMD executable).

Notes regarding # of fastq files as options for script input:
--fastq_files_subset directory contains R1 and R2 fastq files for 8 samples
--an alternative fastq source directory (fastq_files) is available in the associated Zenodo repository (https://doi.org/10.5281/zenodo.3406861), and contains R1 and R2 fastq files for 384 samples;
--file output sizes and time to complete test runs are relatively small for fastq_files_subset (available for download here in GitHub), large for fastq_files (available for download from Zenodo)
===========================================================================

output_directory: ***user-defined (system-specific)***

fastq_directory: SampleTestFiles/fastq_files_subset
---------------------------------------------------

blastn_path: ***user-defined (system-specific)***

db_path: SampleTestFiles/blastn_database
----------------------------------------

db_prefix: GRCh38
-----------------

blastdbcmd_path: ***user-defined (system-specific)***

guideRNA_seq: GTGCCAGCCACATTCAGAAC
----------------------------------

extant_seq: AGAACAGGGTGTTCTG
----------------------------