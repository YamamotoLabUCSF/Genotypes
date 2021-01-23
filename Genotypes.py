#!/usr/local/bin/anaconda3/bin/python3
# Note: edit shebang line above as appropriate for your system
# NAME: Kirk Ehmsen
# FILE: Genotypes.py
# DATE: 08-30-2018/08-02-2019/01-02-2021
# DESC: This script accepts text to standard input, and returns allele definitions for samples
# from a demultiplexed NGS fastq dataset.
# USAGE: ./Genotypes.py or python3 Genotypes.py
# REPO: https://github.com/YamamotoLabUCSF/Genotypes

#############################################################################
# Background notes:
# ==============================================
# This is Genotypes.py v1.1
# ==============================================
# https://github.com/YamamotoLabUCSF/Genotypes
# v1.1/Committed 01-02-2021
# ----------------------------------------------
# For usage details, please refer to README file at GitHub and to the following manuscript:
#   Ehmsen, Knuesel, Martinez, Aridomi, Asahina, Yamamoto (2021)
# Please cite usage as:
#   Genotypes.py
#   Ehmsen, Knuesel, Martinez, Aridomi, Asahina, Yamamoto (2021)

# Operation notes:
# ==============================================
# This script accepts text to standard input, and returns allele definitions (and extrapolated genotypes) for samples
# from a demultiplexed NGS fastq dataset. Required: Python 3.7.0 or higher + 3rd-party libraries NumPy, SciPy, psutil,
# fdpf, PyPDF2; BLASTN & BLASTDBCMD (NCBI) locally installed
# What does this script do?
# 1. counts reads per well/sample (sample defined by fastq file); fastq file name provides the sample name
# 2. identifies the top 10 reads per well (in terms of read abundance), and calculates representation among reads within the well at four levels:
#   (a) raw frequency (% read type in question, relative to total reads)
#   (b) percentile (% of other read types that fall below the frequency of the read type in question)
#   (c) adjusted frequency @ 1% (% read type in question, relative to reads that occur at >1% frequency)
#   (d) adjusted frequency @ 10% (% read type in question, relative to reads that occur at >10% frequency)
# 3. aligns top 10 reads to reference genome using BLASTN (National Center for Biotechnology Information;
#    Altschul S.F. et al. (1990) "Basic local alignment search tool")
# 4. defines alleles, imputes genotypes, and returns alignments
#    -  for mutants, the alignment shows location of Cas9 cut(s) and indel(s) relative to wt,
#       if Cas9 guide sequence(s) supplied by user;
#    -  also indicates location of test sub-sequence(s) and whether sub-sequence is altered (ablated),
#       if test sub-sequence(s) supplied by user
# 5. provides overall population statistics:
#   (a) total sample # for which genotypes were inferred
#   (b) distribution of genotypes among samples (homozygous, heterozygous, etc.)
#   (c) estimated wild-type vs. mutant allele frequencies
#   (e) summary of samples and reads that either had 'no hit' in reference database provided to BLASTN,
#       or multiple hits (>1)

# Input notes:
# ==============================================
# You will be prompted for the following user-specific information (up to 8 inputs; 6 required):
#      ** required (5 directory paths):
#         * where should output files go?                                 path to output directory for output files
#         * where are input files found?                                  path to single directory containing
#                                                                           demultiplexed fastq files
#         * where is BLASTN executable found?                             path to BLASTN installation
#         * where is the reference sequence database used for alignment?  path to directory containing six files
#                                                                           that compose the reference sequence database
#                                                                           for BLASTN alignments
#                                                                           (.nhr, .nin, .nog, .nsd, .nsi, .nsg)
#         * where is BLASTDBCMD executable found?                         path to BLASTDBCMD installation
#      ** required (1 string):
#         * what prefix is shared among all six BLASTN reference sequence alignment database files?
#
#      ** optional (2 sequences: strings or lists of strings):
#         * guide RNA sequence (in DNA representation, excluding PAM sequence)
#         * sequence of interest (sub-sequence motif(s) of interest, to query whether lost or gained in allele(s); i.e., 'test sequence')

# Output notes:
# ==============================================
# This script produces 8 output files in the user-specified output directory.
# These include:
#	  1. fasta.fa
#	  2. blastn_alignments.txt (output of BLASTN operation on fasta.fa)
#     3. allele_definitions.txt (output of script operation on blastn_alignments.txt,
#        samples returned in order of processing)
#     4. allele_evidence.pdf (*optional*; output of script operation on blastn_alignments.txt,
#        plots of calculated read/allele frequencies)
#     5. genotypes.txt (output of script operation on blastn_alignments.txt,
#        samples returned in ranked order based on genotype inference)
#     6. allele_definitions.csv (tabular representation of allele data for all samples)
#	  7. population_summary.txt (output of script operation on genotypes.txt)
#     8. script_metrics.txt (summary/analysis of script operation metrics [metadata])
#
#           Directory structure under an output directory specified as 'Genotypes', for example,
#           would contain the following subdirectories and files following Genotypes.py operations:
#
#           /Genotypes
#                          `-----allele_definitions.csv
#                          `-----allele_definitions.txt
#                          `-----allele_evidence.pdf
#                          `-----blastn_alignments.txt
#                          `-----fasta.fa
#                          `-----genotypes.txt
#                          `-----population_summary.txt
#                          `-----script_metrics.txt
# 
#############################################################################

#############################################################################
# SCRIPT:

# Check for availability of Python dependencies (libraries, modules) in path
missing_dependencies_list = []

try:
    import psutil
except ImportError:
    missing_dependencies_list.append('psutil')
    
try:
    import numpy
except ImportError:
    missing_dependencies_list.append('numpy')

try:
    import scipy
except ImportError:
    missing_dependencies_list.append('scipy')
    
try:
    import matplotlib.pyplot
except ImportError:
    missing_dependencies_list.append('matplotlib.pyplot')
    
try:
    import fpdf
except ImportError:
    missing_dependencies_list.append('fpdf')
    
try:
    import PyPDF2
except ImportError:
    missing_dependencies_list.append('PyPDF2')

try:
    import pandas
except ImportError:
    missing_dependencies_list.append('pandas')
    
if len(missing_dependencies_list) > 0:
    print('ModuleNotFoundError\n')
    print('Please note, the following required Python module(s) are not found in your Python system path:')
    for i in missing_dependencies_list:
        print('   '+i)
    print('\nPlease exit the script and install these Python dependencies in your system path.')
    print("""\nGuidelines for installation of Python dependencies can be found in the README file for Genotypes.py ('System Setup')""")
    print("""    (Creation of a Python virtual environment is recommended)""")

# Import libraries, modules
# Operating system interfaces
import os

# import subprocess management
import subprocess

# Time access and conversions, Basic data and time types
import time
from datetime import datetime

# System-specific parameters and functions
import sys

# Process and system utilities
import psutil
from psutil import virtual_memory

# Gzip to read GNU zipped files
import gzip

# Low-level networking interface
import socket

# System version information
import platform

# Unix-style pathname pattern expansion
import glob

# Object-oriented filesystem paths
from pathlib import Path

# NumPy (numeric operations)
import numpy

# SciPy (for percentile) 
from scipy import stats

# Statistics
import statistics

# Container datatypes (for Counter operation)
from collections import Counter

# Decimal fixed point and floating point arithmetic
from decimal import Decimal

# Internationalization services (for use of thousands separator in numbers where appropriate)
import locale
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

# Regular expression operations
import re

# Python plotting
import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.ticker as ticker

# Python PDF tools
from fpdf import FPDF
from PyPDF2 import PdfFileMerger, PdfFileReader

# Python panel data frames
import pandas as pd

# start time
initialTime = datetime.now()

# Define 'prompts' function for coached user input
def prompts():
    """Coached prompts to collect user input"""
    # Make variables assigned in prompts() function globally available
    global output_directory
    global fastq_directory
    global blastn_path
    global db_path
    global db_prefix
    global blastdbcmd_path
    global guideRNA_seq
    global extant_seq
    global test_seq
    # 1-Specify output directory.
    print(r"""
    ---------------------------------------------
    Location of OUTPUT DIRECTORY for output files
    ---------------------------------------------
	
    The script generates 8 separate files, all in the directory you indicate here.  It is important that this directory
    either not exist prior to running the script, or if it does exist, it must be *empty* of any files with the names
    to be created below.  These files are:
    
        1. fasta.fa

        2. blastn_alignments.txt
            (output of BLASTN operation on fasta.fa)

        3. allele_definitions.txt
            (output of script operation on blastn_alignments.txt, samples returned in order of processing)
                  
        4. allele_evidence.pdf (optional)
            (output of script operation on blastn_alignments.txt, plot of calculated read/allele frequencies)

        5. genotypes.txt
            (output of script operation on blastn_alignments.txt, samples returned in order of genotype inference)

        6. population_summary.txt
            (output of script operation on genotypes.txt)
                  
        7. allele_definitions.csv
            (allele metrics (frequency representations) and definitions for each sample, in spreadsheet format)
                  
        8. script_metrics.txt 
            (summary/analysis of script operation metrics)
            
            Notes: 
            * These files do not exist before the script is run. The files are made by the script.
            * The primary data outputs for genotypes are found in:
                 allele_definitions.txt, allele_evidence.pdf, genotypes.txt & population_summary.txt
        
        
    At this prompt, indicate an absolute path to a ** directory ** that will be created by the script as the location
    for output files.  This directory should not exist yet -- it will be created as an output of this script, and will
    be populated with the file outputs of this specific instance of the script operation.

    Example: if you'd like to create a directory ('Genotypes') in an existing directory ('Illumina'), accessed
    with absolute path of '/Users/myname/Illumina/Genotypes' (Mac) or 'C:\Users\myname\Illumina\Genotypes'
    (Windows), enter '/Users/myname/Illumina/Genotypes' at the command line prompt. Replace 'myname' with the
    appropriate intervening directory identifiers. Do *not* flank your entry with quotation marks (') at the
    command line.
    
    Alternatively, simply enter a desired directory name (e.g., 'Genotypes') and run this script from
    within a directory where you'd like to create this new directory."""+'\n')
    output_directory = input(r"""    -----> Output directory name and path:  """)
    # 2-Specify the fastq files to be used for input, by indicating directory location of the file list.
    print(r"""
    ------------------------------------------------------------------------------
    Location of INPUT FILES (single directory containing demutiplexed fastq files)
    ------------------------------------------------------------------------------

    You will now be asked to enter the path to the directory containing the fastq files
    to be processed as Genotypes.py input.  

    Example: if your fastq input files are named file1.fastq, file2.fastq, etc. and are found in a directory
    named 'Sequences' with absolute path of '/Users/myname/Sequences' (Mac) or 'C:\Users\myname\Sequences' (PC),
    enter '/Users/myname/Sequences' at the command line prompt.

    When you're done entering the fastq file location, press 'enter' again to proceed in the script."""+'\n')
    fastq_directory = input(r"""    -----> Directory name and path:  """)
    # 3-Collect path to BLASTN executable.
    print("""
    -----------------------------
    Location of BLASTN EXECUTABLE
    -----------------------------

    This script uses BLASTN (NCBI) to align reads from your fastq files to a reference sequence database
    (such as a genome database or sequence database).
    Please indicate the absolute path to the BLASTN executable.
    
    Example: if your BLASTN executable is found at absolute path /Users/myname/blastn, type '/Users/myname/blastn'
    and press Enter."""+'\n')
    blastn_path = input(r"""    -----> Path to BLASTN executable:  """)
    # 4-Collect location of BLASTN database directory.
    print("""
    -----------------------------------------------
    Location of BLASTN ALIGNMENT DATABASE DIRECTORY
    -----------------------------------------------

    Because this script uses BLASTN (NCBI) to align reads from your fastq files a reference sequence database,
    an alignment reference database is needed. This reference database consists of a single directory containing
    six files (.nhr, .nin, .nog, .nsd, .nsi, .nsg), generated by the program MAKEBLASTDB (NCBI) from a custom file
    containing sequences in fasta format (or available for some genomes as downloads from NCBI).
    
    Please indicate the absolute path to the directory you are using as your reference sequence database.
    
    Example: if your reference sequence database is found at absolute path /Users/myname/database, type
    '/Users/myname/database' and press Enter."""+'\n')
    db_path = input(r"""    -----> Path to BLASTN alignment reference sequence database:  """)
    # 5-Collect prefix of BLASTN database files.
    print("""
    ------------------------------------------------
    PREFIX common to BLASTN ALIGNMENT DATABASE FILES
    ------------------------------------------------

    A BLASTN reference sequence database consists of six files in a single directory, with each of the six
    files sharing a common prefix (usually determined by the name of the fasta file provided to MAKEBLASTDB
    during database generation).
    
    Please indicate the common prefix for files of the reference sequence database."""+'\n')
    db_prefix = input(r"""    -----> Prefix for alignment reference sequence database files:  """)
    # 6-Collect path to BLASTDBCMD executable.
    print("""
    ---------------------------------
    Location of BLASTDBCMD EXECUTABLE
    ---------------------------------

    This script uses BLASTDBCMD (NCBI) to retrieve DNA sequence spans between coordinates in a reference sequence database* (such as a genome database or sequence database) (*in cases where BLASTN has split an alignment into >1 local high-scoring pair (hsp)).
    
    Please indicate the absolute path to the BLASTDBCMD executable.
    
    Example: if your BLASTDBCMD executable is found at absolute path /Users/myname/blastdbcmd,
    type '/Users/myname/blastdbcmd' and press Enter."""+'\n')
    blastdbcmd_path = input(r"""    -----> Path to BLASTDBCMD executable:  """) 
    # 7-Specify whether to include sequence(s) of interest to query in alignment outputs.
    print("""      
    -----------------------------------------------------------------
    Optional: Nucleotide sequence(s) to identify in output alignments
    -----------------------------------------------------------------
    
    Some applications of 'allele definition' and 'genotype inference' may call for identification of the presence
    or absence of a specific anticipated sub-sequence (few nucleotides), and/or for the mapping of the location of
    a sub-sequence if present in the sequence alignment ('sequence of interest').
    
    Genotypes.py allows for the optional testing of sub-sequences.
    
    If you would like to specify subsequences, type 'Yes' and press Enter.
    Otherwise, if you do not wish to specify subsequences, type 'No' and press Enter."""+'\n')
    test_seq = input("""    -----> 'Yes' or 'No' to sub-sequence specification:  """)
    if test_seq == 'Yes':
        print("""  
    --------------------------------------------------------------------------------------------------
    Nucleotide sequence(s): guide RNA annealing sites and/or test for presence/absence of sub-sequence
    --------------------------------------------------------------------------------------------------""")
        # 7a-Collect guide RNA sequence details.
        print("""
    ............................................................
    ***** guide RNA details: specify guide RNA sequence(s) *****
    
    To specify guide RNA sequence(s), enter text for each directly at the command line,
    separated by a comma ('x,y').
    
    Please specify guide RNA sequence(s) [excluding PAM]:
    
    When text entries are entered, press ‘enter’ again to proceed in the script.
    To skip text entries for these fields, simply press ‘enter’ until the next prompt appears.
	
    Examples:
      If your single guide RNA sequence is 'ATCCAGTTCTCCAGTCTCCC', enter: 'ATCCAGTTCTCCAGTCTCCC'.
      If you have two guide RNA sequences and they are 'ATCCAGTTCTCCAGTCTCCC' and 'GCGAGCTCGTGTCTGTGACG',
      enter: 'ATCCAGTTCTCCAGTCTCCC, GCGAGCTCGTGTCTGTGACG'."""+'\n')
        guideRNA_seq_temp = input(r"""    -----> guide RNA sequence(s):  """)
        guideRNA_seq = [i.strip() for i in guideRNA_seq_temp.split(',')]
        # 7b-Collect sub-sequence details (extant in allele or not?).
        print("""
    ......................................................................................
    ***** query DNA sequence(s): specify sequence(s) to test for presence or absence *****
    
    To specify query DNA sequence(s), enter text for each directly at the command line,
    separated by a comma ('x,y').
    
    Please specify short DNA sequence to test for presence vs. ablation.
    
    When text entries are entered, press ‘enter’ again to proceed in the script.
    To skip text entries for these fields, simply press ‘enter’ until the next prompt appears.
	
    Examples:
      If your single query sequence is 'TACTCAATATCGATC', enter: 'TACTCAATATCGATC'.
      If you have two query sequences and they are 'TACTCAATATCGATC' and 'CGGGAGCCCGAG', enter:
      'TACTCAATATCGATC, CGGGAGCCCGAG'."""+'\n')
        extant_seq_temp = input(r"""    -----> query DNA sequence(s):  """)
        extant_seq = [i.strip() for i in extant_seq_temp.split(',')]

# Define 'allele_output' function to report defined alleles for samples based on genotype class designation
def allele_output(genotype_class):
    """
    This function outputs allele definitions for samples belonging to a specified class of inferred genotypes
    """
    for i in genotype_class:
        file.write(i+'\n'+(18+len(imputedgenotypes_dict.get(i)[0]))*'*'+'\n'+'INFERRED GENOTYPE: '+imputedgenotypes_dict.get(i)[0]+'\n'+(18+len(imputedgenotypes_dict.get(i)[0]))*'*'+'\n\n')
        read_checklist = []
        read_abundance_checklist = []
        for n in range(1, len(imputedgenotypes_dict.get(i))):
            if imputedgenotypes_dict.get(i)[n][0].get('allele_name').split(' ')[1] == 'R1+R2':
                if 'R1+R2' not in read_checklist:
                    read_checklist.append('R1+R2')
                    file.write('\n'+3*' '+'*'+17*'~'+'*\n'+3*' '+'| READ 1 + READ 2 |\n'+3*' '+'*'+17*'~'+'*\n\n')
                else:
                    pass
                if float(imputedgenotypes_dict.get(i)[n][0].get('allele_name').split(' ')[4].split(':')[1]) < 10:
                    if 'R1dregs' not in read_abundance_checklist:
                        read_abundance_checklist.append('R1dregs')
                        file.write(3*' '+'*'+56*'~'+'*\n'+3*' '+'|  >>>>> remaining alleles occur at frequency <10% <<<<< |\n'+3*' '+'*'+56*'~'+'*\n\n')
                    else:
                        pass
            if imputedgenotypes_dict.get(i)[n][1].get('allele_type') == 'wild-type':
                file.write(3*' '+'Allele: '+imputedgenotypes_dict.get(i)[n][0].get('allele_name')+' | '+imputedgenotypes_dict.get(i)[n][1].get('allele_type')+'\n    Locus: '+imputedgenotypes_dict.get(i)[n][0].get('chr+build')+', '+imputedgenotypes_dict.get(i)[n][0].get('locusID')+' '+imputedgenotypes_dict.get(i)[n][0].get('coordinates')+'\n')
            else:
                file.write(3*' '+'Allele: '+imputedgenotypes_dict.get(i)[n][0].get('allele_name')+' | '+imputedgenotypes_dict.get(i)[n][1].get('allele_type')+', '+imputedgenotypes_dict.get(i)[n][1].get('allele_specs')+'\n    Locus: '+imputedgenotypes_dict.get(i)[n][0].get('chr+build')+', '+imputedgenotypes_dict.get(i)[n][0].get('locusID')+' '+imputedgenotypes_dict.get(i)[n][0].get('coordinates')+'\n')           
            for guide in imputedgenotypes_dict.get(i)[n][2]:
                if imputedgenotypes_dict.get(i)[n][2].get(guide) != 'None':
                    if guide in guideRNA_seq:
                        file.write('\n'+(int(imputedgenotypes_dict.get(i)[n][2].get(guide))-2)*' '+"5'-"+guide+"-3' (guide sequence)"+'\n'+(int(imputedgenotypes_dict.get(i)[n][2].get(guide))+len(guide)-3)*' '+'v')
                    elif guide in guideRNA_seq_rev:
                        file.write('\n'+(int(imputedgenotypes_dict.get(i)[n][2].get(guide))-2)*' '+"3'-"+guide+"-5' (guide sequence)"+'\n'+(int(imputedgenotypes_dict.get(i)[n][2].get(guide))+4)*' '+'v')
            file.write(imputedgenotypes_dict.get(i)[n][0].get('alignment'))
            for seq in imputedgenotypes_dict.get(i)[n][3]:
                if imputedgenotypes_dict.get(i)[n][3].get(seq) != 'None':
                    file.write('\n')
                    if seq in extant_seq:
                        file.write((1+int(imputedgenotypes_dict.get(i)[n][3].get(seq)))*' '+len(seq)*'^'+'\n'+(int(imputedgenotypes_dict.get(i)[n][3].get(seq))-2)*' '+"5'-"+seq+"-3' (sequence of interest)\n")
                    elif seq in extant_seq_rev:
                        file.write((1+int(imputedgenotypes_dict.get(i)[n][3].get(seq)))*' '+len(seq)*'^'+'\n'+(int(imputedgenotypes_dict.get(i)[n][3].get(seq))-2)*' '+"3'-"+seq+"-5' (sequence of interest)\n")
                elif imputedgenotypes_dict.get(i)[n][3].get(seq) == 'None':
                    file.write('\n')
            file.write('\n')
            file.write('\n')

# Define 'frequency_plots' function to plot sample allele frequency metrics for visualization in pdf file
def frequency_plots():
    """
    This function plots sample allele frequency metrics for visualization in a pdf file
    """
    # Make variable assigned in frequency_plots() function globally available
    global frequencyplotsDuration
    # Start the clock on plot time duration
    startTime_frequencyplots = datetime.now()
    # Assign allele_evidence.pdf file to output path
    allele_evidence_output = Path(str(output_path)+'/'+processdate+'_allele_evidence.pdf')
    # Initiate PDF file to record allele frequency plots for each sample
    pdf = FPDF(format='letter')
    pdf.add_page()
    pdf.set_font("Arial", size=20, style='B')
    pdf.ln(20)
    pdf.write(5, 'Frequency plots to support inferred genotypes')
    pdf.output(allele_evidence_output)
    # Generate allele frequency plots for each sample, based on the following principles and frequency metrics:
    # for each sample, up to 10 candidate alleles are 'ranked' based on relative read frequency in the initial fastq file
    # four frequency plots are generated, (1) raw read frequency relative to all other reads, (2) frequency relative to the
    # top 10 most abundant reads, (3) frequency relative to reads that occur at >1% raw abundance; (4) frequency relative to
    # reads that occur at >10% raw abundance
    for samplename in imputedgenotypes_dict:
        plot_name = output_directory / (samplename+'_plot.png')
        pdf_output = output_directory / (samplename+'_.pdf')
        R1R2_allele_list = []
        R1R2_allele_names = []
        R1R2_allele_frequency = []
        R1R2_allele_type = []
        R1R2_allele_specs = []
        for x in range(1, len(imputedgenotypes_dict[samplename])):
            if imputedgenotypes_dict[samplename][x][0].get('allele_name').split(' ')[1] == 'R1+R2':
                R1R2_allele_list.append(imputedgenotypes_dict[samplename][x][0].get('allele_name'))
                R1R2_allele_type.append(imputedgenotypes_dict[samplename][x][1].get('allele_type'))
                R1R2_allele_specs.append(imputedgenotypes_dict[samplename][x][1].get('allele_specs'))
        for i in range(1, len(R1R2_allele_list)+1):
            R1R2_allele_names.append(i)
        for i in R1R2_allele_list:
            R1R2_allele_frequency.append(float(i.split(' ')[4].split(':')[1]))
# 
        x1 = R1R2_allele_names
#
        N = len(R1R2_allele_names)
        width = 0.4
        spacing1 = [float(i) for i in range(1,N+1)]
#
        y1a = R1R2_allele_frequency
        while len(y1a) < N:
            y1a.append(float(0))
        label_list_a = [value for value in zip(x1,y1a,spacing1,R1R2_allele_type,R1R2_allele_specs)]
#
        R1R2_allele_frequency_top10 = []
        for i in R1R2_allele_list:
            freq1 = i.split(' ')[6].split(':')[1]
            R1R2_allele_frequency_top10.append(float(freq1) if freq1 != 'None' else 0)
        y1b = R1R2_allele_frequency_top10
        while len(y1b) < N:
            y1b.append(float(0))
        label_list_b = [value for value in zip(x1,y1b,spacing1,R1R2_allele_type,R1R2_allele_specs)]
    
#
        R1R2_allele_frequency_1 = []
        for i in R1R2_allele_list:
            freq1 = i.split(' ')[7].split(':')[1]
            R1R2_allele_frequency_1.append(float(freq1) if freq1 != 'None' else 0)
        y1c = R1R2_allele_frequency_1
        while len(y1c) < N:
            y1c.append(float(0))
        label_list_c = [value for value in zip(x1,y1c,spacing1,R1R2_allele_type,R1R2_allele_specs)]
#
        R1R2_allele_frequency_10 = []
        for i in R1R2_allele_list:
            freq1 = i.split(' ')[8].split(':')[1]
            R1R2_allele_frequency_10.append(float(freq1) if freq1 != 'None' else 0)
        y1d = R1R2_allele_frequency_10
        while len(y1d) < N:
            y1d.append(float(0))
        label_list_d = [value for value in zip(x1,y1d,spacing1,R1R2_allele_type,R1R2_allele_specs)]
#
# Plots
        fig = plt.figure(figsize=(11,7), dpi=250)
# Subplot 1
        ax1 = fig.add_subplot(141)
        ax1.grid(color='#808080', linestyle='--', linewidth=0.2, axis='x', zorder=0)
        rects1 = ax1.barh(spacing1, y1a, width, color='#87CEFA', alpha=1, edgecolor='black', linewidth=0.5, align='center', zorder=2, label='  R1\n+R2')
        ax1.set_ylabel('Allele Rank', fontsize=10, fontname='Arial', alpha=0.8)
        ax1.set_xlabel('Frequency', fontsize=10, fontname='Arial', alpha=0.8)
        ax1.set_title('Allele frequencies\n(% total reads)', fontsize=8.5, fontweight='bold', fontname='Arial')
        plt.xlim([0,110])
        plt.ylim([0.5, N+0.5])
        plt.gca().invert_yaxis()
        ax1 = plt.gca()
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(25))
        ax1.xaxis.set_minor_locator(ticker.MultipleLocator(5))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax1.spines['right'].set_visible(False)
        ax1.legend(loc = 'lower right', frameon=True, fontsize=7)
        ax1.spines['bottom'].set_color('grey')
        ax1.spines['top'].set_color('grey')
        ax1.spines['left'].set_color('grey')
        ax1.tick_params(which='both', color='grey', labelcolor='grey')
        for i in label_list_a:
            if i[1] != 0:
                if i[1] > 20:
                    if str(i[4]) != 'None':
                        ax1.text(i[1]/10, i[2], str(i[1])+'% | '+i[3], va = 'bottom', color = 'black', weight = 'bold', fontsize=7)
                        ax1.text(i[1]/10, i[2], str(i[4]), va = 'top', color = 'black', style = 'italic', fontsize=7)
                    else:
                        ax1.text(i[1]/10, i[2], str(i[1])+'% | '+i[3], va = 'center', color = 'black', weight = 'bold', fontsize=7)
                else:
                    if str(i[4]) != 'None':
                        ax1.text(i[1]+2, i[2], str(i[1])+'% | '+i[3], va = 'bottom', color = 'black', fontsize=7)
                        ax1.text(i[1]+2, i[2], str(i[4]), va = 'top', color = 'black', style = 'italic', fontsize=7)
                    else:
                        ax1.text(i[1]+2, i[2], str(i[1])+'% | '+i[3], va = 'center', color = 'black', fontsize=7)
#
# Subplot 2
        ax2 = fig.add_subplot(142)
        ax2.grid(color='#808080', linestyle='--', linewidth=0.2, axis='x', zorder=0)
        rects1 = ax2.barh(spacing1, y1b, width, color='#87CEFA', alpha=1, edgecolor='black', linewidth=0.5, align='center', zorder=2, label='  R1\n+R2')
        ax2.set_ylabel('Allele Rank', fontsize=10, fontname='Arial', alpha=0.8)
        ax2.set_xlabel('Frequency', fontsize=10, fontname='Arial', alpha=0.8)
        ax2.set_title('Allele frequencies\n(% top 10 most abundant reads)', fontsize=8.5, fontweight='bold', fontname='Arial')
        plt.xlim([0,110])
        plt.ylim([0.5, N+0.5])
        plt.gca().invert_yaxis()
        ax2 = plt.gca()
        ax2.xaxis.set_major_locator(ticker.MultipleLocator(25))
        ax2.xaxis.set_minor_locator(ticker.MultipleLocator(5))
        ax2.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax2.spines['right'].set_visible(False)
        ax2.legend(loc = 'lower right', frameon=True, fontsize=7)
        ax2.spines['bottom'].set_color('grey')
        ax2.spines['top'].set_color('grey')
        ax2.spines['left'].set_color('grey')
        ax2.tick_params(which='both', color='grey', labelcolor='grey')
        for i in label_list_b:
            if i[1] != 0:
                if i[1] > 20:
                    if str(i[4]) != 'None':
                        ax2.text(i[1]/10, i[2], str(i[1])+'% | '+i[3], va = 'bottom', color = 'black', weight = 'bold', fontsize=7)
                        ax2.text(i[1]/10, i[2], str(i[4]), va = 'top', color = 'black', style = 'italic', fontsize=7)
                    else:
                        ax2.text(i[1]/10, i[2], str(i[1])+'% | '+i[3], va = 'center', color = 'black', weight = 'bold', fontsize=7)
                else:
                    if str(i[4]) != 'None':
                        ax2.text(i[1]+2, i[2], str(i[1])+'% | '+i[3], va = 'bottom', color = 'black', fontsize=7)
                        ax2.text(i[1]+2, i[2], str(i[4]), va = 'top', color = 'black', style = 'italic', fontsize=7)
                    else:
                        ax2.text(i[1]+2, i[2], str(i[1])+'% | '+i[3], va = 'center', color = 'black', fontsize=7)
#
# Subplot 3       
        ax3 = fig.add_subplot(143)
        ax3.grid(color='#808080', linestyle='--', linewidth=0.2, axis='x', zorder=0)
        rects1 = ax3.barh(spacing1, y1c, width, color='#87CEFA', alpha=1, edgecolor='black', linewidth=0.5, align='center', zorder=2, label='  R1\n+R2')
        ax3.set_ylabel('Allele Rank', fontsize=10, fontname='Arial', alpha=0.8)
        ax3.set_xlabel('Frequency', fontsize=10, fontname='Arial', alpha=0.8)
        ax3.set_title('Allele frequencies\n(% reads adjusted for frequency >1%)', fontsize=8.5, fontweight='bold', fontname='Arial')
        plt.xlim([0,110])
        plt.ylim([0.5, N+0.5])
        plt.gca().invert_yaxis()
        ax3 = plt.gca()
        ax3.xaxis.set_major_locator(ticker.MultipleLocator(25))
        ax3.xaxis.set_minor_locator(ticker.MultipleLocator(5))
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax3.spines['right'].set_visible(False)
        ax3.legend(loc = 'lower right', frameon=True, fontsize=7)
        ax3.spines['bottom'].set_color('grey')
        ax3.spines['top'].set_color('grey')
        ax3.spines['left'].set_color('grey')
        ax3.tick_params(which='both', color='grey', labelcolor='grey')
        for i in label_list_c:
            if i[1] != 0:
                if i[1] > 20:
                    if str(i[4]) != 'None':
                        ax3.text(i[1]/10, i[2], str(i[1])+'% | '+i[3], va = 'bottom', color = 'black', weight = 'bold', fontsize=7)
                        ax3.text(i[1]/10, i[2], str(i[4]), va = 'top', color = 'black', style = 'italic', fontsize=7)
                    else:
                        ax3.text(i[1]/10, i[2], str(i[1])+'% | '+i[3], va = 'center', color = 'black', weight = 'bold', fontsize=7)
                else:
                    if str(i[4]) != 'None':
                        ax3.text(i[1]+2, i[2], str(i[1])+'% | '+i[3], va = 'bottom', color = 'black', fontsize=7)
                        ax3.text(i[1]+2, i[2], str(i[4]), va = 'top', color = 'black', style = 'italic', fontsize=7)
                    else:
                        ax3.text(i[1]+2, i[2], str(i[1])+'% | '+i[3], va = 'center', color = 'black', fontsize=7)
#
# Subplot 4
        ax4 = fig.add_subplot(144)
        ax4.grid(color='#808080', linestyle='--', linewidth=0.2, axis='x', zorder=0)
        rects1 = ax4.barh(spacing1, y1d, width, color='#87CEFA', alpha=1, edgecolor='black', linewidth=0.5, align='center', zorder=2, label='  R1\n+R2')
        ax4.set_ylabel('Allele Rank', fontsize=10, fontname='Arial', alpha=0.8)
        ax4.set_xlabel('Frequency', fontsize=10, fontname='Arial', alpha=0.8)
        ax4.set_title('Allele frequencies\n(% reads adjusted for frequency >10%)', fontsize=8.5, fontweight='bold', fontname='Arial')
        plt.xlim([0,110])
        plt.ylim([0.5, N+0.5])
        plt.gca().invert_yaxis()
        ax4 = plt.gca()
        ax4.xaxis.set_major_locator(ticker.MultipleLocator(25))
        ax4.xaxis.set_minor_locator(ticker.MultipleLocator(5))
        ax4.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax4.spines['right'].set_visible(False)
        ax4.legend(loc = 'lower right', frameon=True, fontsize=7)
        ax4.spines['bottom'].set_color('grey')
        ax4.spines['top'].set_color('grey')
        ax4.spines['left'].set_color('grey')
        ax4.tick_params(which='both', color='grey', labelcolor='grey')
        for i in label_list_d:
            if i[1] != 0:
                if i[1] > 20:
                    if str(i[4]) != 'None':
                        ax4.text(i[1]/10, i[2], str(i[1])+'% | '+i[3], va = 'bottom', color = 'black', weight = 'bold', fontsize=7)
                        ax4.text(i[1]/10, i[2], str(i[4]), va = 'top', color = 'black', style = 'italic', fontsize=7)
                    else:
                        ax4.text(i[1]/10, i[2], str(i[1])+'% | '+i[3], va = 'center', color = 'black', weight = 'bold', fontsize=7)
                else:
                    if str(i[4]) != 'None':
                        ax4.text(i[1]+2, i[2], str(i[1])+'% | '+i[3], va = 'bottom', color = 'black', fontsize=7)
                        ax4.text(i[1]+2, i[2], str(i[4]), va = 'top', color = 'black', style = 'italic', fontsize=7)
                    else:
                        ax4.text(i[1]+2, i[2], str(i[1])+'% | '+i[3], va = 'center', color = 'black', fontsize=7)
        plt.tight_layout()
#
        plt.savefig(plot_name, format='png', dpi=250)
        plt.close(fig)
#   
        pdf = FPDF('L', 'mm', (400, 250))
        pdf.add_page()
        pdf.set_font("Arial", size=8, style='BU')
        pdf.write(5, samplename)
        pdf.ln(3)
        pdf.set_font("Arial", size=7)
        pdf.write(5, 'inferred genotype: '+imputedgenotypes_dict[samplename][0].split('|')[1]+','+imputedgenotypes_dict[samplename][0].split('|')[2])
        pdf.ln(5)
        allele_count_R1R2 = 1
        read1read2_check = []
        for x in range(1, len(imputedgenotypes_dict[samplename])):
            if imputedgenotypes_dict[samplename][x][0].get('allele_name').split(' ')[1] == 'R1+R2':
                if imputedgenotypes_dict[samplename][x][0].get('allele_name').split(' ')[8].split(':')[1] == 'None':
                    pass
                elif float(imputedgenotypes_dict[samplename][x][0].get('allele_name').split(' ')[8].split(':')[1]) > 20:
                    if len(read1read2_check) == 0:
                        read1read2_check.append('R1+R2')
                        pdf.set_font("Arial", size=7, style='U')
                        pdf.write(5, 'R1+R2, sequences with >20% representation among reads:')
                        pdf.ln(3)
                    else:
                        pass
                    pdf.set_font("Arial", size=6)
                    if imputedgenotypes_dict[samplename][x][1].get('allele_specs') is not None:
                        pdf.write(5, 'Allele '+str(allele_count_R1R2)+': '+imputedgenotypes_dict[samplename][x][0].get('allele_name').split(' ')[2]+' '+imputedgenotypes_dict[samplename][x][1].get('allele_type')+', '+imputedgenotypes_dict[samplename][x][1].get('allele_specs'))
                    else:
                        pdf.write(5, 'Allele '+str(allele_count_R1R2)+': '+imputedgenotypes_dict[samplename][x][0].get('allele_name').split(' ')[2]+' '+imputedgenotypes_dict[samplename][x][1].get('allele_type'))     
                    pdf.ln(3)
                    pdf.set_font("Courier", size=6)
                    if len(imputedgenotypes_dict[samplename][x][0].get('alignment').split('\n')[1]) > 200:
                        string = imputedgenotypes_dict[samplename][x][0].get('alignment').split('\n')[1]
                        string2 = imputedgenotypes_dict[samplename][x][0].get('alignment').split('\n')[2]
                        string3 = imputedgenotypes_dict[samplename][x][0].get('alignment').split('\n')[3]
                        chunks = range(12, len(string), 100)
                        for index, chunk in enumerate(chunks):
                            if index == 0:
                                pdf.write(5, string[chunk-12:chunk+100])
                                pdf.ln(3)
                                pdf.write(5, string2[chunk-12:chunk+100])
                                pdf.ln(3)
                                pdf.write(5, string3[chunk-12:chunk+100])
                            elif 0 < index < (len(chunks)-1):
                                pdf.write(5, '\n'+4*' '+'query'+'  '+string[chunk:chunk+100])
                                pdf.ln(3)
                                pdf.write(5, 11*' '+string2[chunk:chunk+100])
                                pdf.ln(3)
                                pdf.write(5, 'reference'+'  '+string3[chunk:chunk+100])
                            else:
                                pdf.write(5, '\n'+4*' '+'query'+'  '+string[chunk:])
                                pdf.ln(3)
                                pdf.write(5, 11*' '+string2[chunk:])
                                pdf.ln(3)
                                pdf.write(5, 'reference'+'  '+string3[chunk:])
                        allele_count_R1R2 = allele_count_R1R2+1
                        pdf.ln(4)
                    else:
                        pdf.write(5, imputedgenotypes_dict[samplename][x][0].get('alignment').split('\n')[1])
                        pdf.ln(3)
                        pdf.write(5, imputedgenotypes_dict[samplename][x][0].get('alignment').split('\n')[2])
                        pdf.ln(3)
                        pdf.write(5, imputedgenotypes_dict[samplename][x][0].get('alignment').split('\n')[3])
                        allele_count_R1R2 = allele_count_R1R2+1
                        pdf.ln(4)
        pdf.ln(5)
        pdf.image(str(plot_name), x = None, y = None, w = 195, h = 0, type = '', link = '')
        pdf.output(str(pdf_output))
#    
        merger = PdfFileMerger()
#   
        filenames =[str(allele_evidence_output), str(pdf_output)]
#
        for filename in filenames:
            merger.append(PdfFileReader(filename, 'rb'))
#   
        merger.write(str(allele_evidence_output))
    # Remove sample png file and pdf file as intermediaries
        try:
            os.remove(plot_name)
        except OSError:
            pass
#    
        try:
            os.remove(pdf_output)
        except OSError:
            pass  
# Log frequency plotting time duration
    frequencyplotsDuration = str(datetime.now()- startTime_frequencyplots).split(':')[0]+' hr|'+str(datetime.now() - startTime_frequencyplots).split(':')[1]+' min|'+str(datetime.now() - startTime_frequencyplots).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_frequencyplots).split(':')[2].split('.')[1]+' microsec'

# Define 'merge' function to merge R1 & R2 reads
def merge(s1, s2):
    i = 0
    while not s2.startswith(s1[i:]):
        i += 1
    if i < len(s2):
        return s1[:i] + s2
    else:
        return 'no overlap'
      
# Define 'merge1' function to append two strings that do not overlap
def merge1(s1, s2):
    i = 0
    while not s2.startswith(s1[i:]):
        i += 1
    return s1[:i] + s2

# Define nt complement dictionary      
nt_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', '-':'-'}    
    
# Define 'convert_bytes' and 'path_size' functions to be used in data collection for script_metrics.txt
def convert_bytes(num):
    """
    This function converts bytes to convenient order of magnitude prefixes
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0

def path_size(given_path):
    """
    This function returns file or directory size
    """
    if os.path.isfile(given_path):
        file_info = os.stat(given_path)
        return convert_bytes(file_info.st_size)
    elif os.path.isdir(given_path):
        dir_info = os.stat(given_path)
        return convert_bytes(dir_info.st_size)
        
# Welcome/orient to script
print("""
    ===================================================
    Genotypes.py v1.1
    ===================================================
    https://github.com/YamamotoLabUCSF/Genotypes
    v1.1/Committed 01-02-2021
    ---------------------------------------------------
    This script accepts text to standard input, and returns allele definitions for samples from a demultiplexed
    next-generation sequencing (NGS) fastq dataset.
    
    Python3, BLASTN and BLASTDBCMD (NCBI) are required for operation.
    BLASTN and BLASTDBCMD can be downloaded and locally installed at 'https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/'.
    
    For usage details, please refer to README file at GitHub and to the following manuscript:
        Ehmsen, Knuesel, Martinez, Asahina, Aridomi, Yamamoto (2021)
    
    Please cite usage as:
        Genotypes.py
        Ehmsen, Knuesel, Martinez, Asahina, Aridomi, Yamamoto (2021)
    
    -------------------------------------------------------------------------------
    Welcome.  You will now be prompted for the following user-specific information:

              ** required:
                 * where should output files go?                        path to output directory for output files
                 * where are input files found?                         path to single directory containing
                                                                          demultiplexed fastq files
                 * where is BLASTN executable found?                    path to BLASTN installation
                 * where is the reference database used for alignment?  path to directory containing six files
                                                                          that compose the reference database
                                                                          for BLASTN alignments
                                                                          (.nhr, .nin, .nog, .nsd, .nsi, .nsg)
                 * what prefix is shared by the BLASTN alignment database files?
                 * where is BLASTDBCMD executable found?                path to BLASTDBCMD installation
                 
              ** optional:
                 * guide RNA sequence (in DNA representation, excluding PAM sequence)
                 * 'sequence of interest' (sub-sequence motif(s) of interest, to query whether lost or gained in allele(s); i.e., 'test sequence')

    """)

input("    Press Enter to continue...")

# Specify whether user input is provided at individual coached prompts or as single-list entry
user_input = input(r"""
    ---------------------------------------------------------------------
    User-specified input: choice of coached prompts vs. single-list entry
    ---------------------------------------------------------------------
    
    Values for the user-specified input indicated above can be entered at individually coached command-line prompts,
    or as a single list of variables provided in a single command-line entry without coached prompts.
    
    To proceed with input at individual command-line PROMPTS, type 'Prompt' and press Enter;
    To proceed with input provided as a single LIST in one command-line entry, type 'List' and press Enter:  """)

if user_input == 'Prompt':
    prompts()
elif user_input == 'List':
    print("""

    You specified LIST format to specify input values.
    ..............................................................................................................
    Some applications of 'allele definition' and 'genotype inference' may call for identification of the presence
    or absence of a specific anticipated sub-sequence (few nucleotides), and/or for the mapping of the location of
    a sub-sequence if present in the sequence alignment.
    
    Genotypes.py allows for the optional testing of sub-sequences.
    
    If you would like to specify subsequences, type 'Yes' and press Enter.
    Otherwise, if you do not wish to specify subsequences, type 'No' and press Enter."""+'\n')
    test_seq = input(r"""    -----> 'Yes' or 'No' to sub-sequence specification:  """)
    if test_seq == 'Yes':
        print(r"""    
        Will you specify guide RNA sequence(s) to map onto output alignments?
        """)
        user_input2 = input(r"""    -----> 'Yes' or 'No' to guide RNA specification:  """)
        print(r"""    
        Will you specify query DNA sequence(s) (sequence(s) to test for presence or absence) to map onto output
        alignments?
        """)
        user_input3 = input(r"""    -----> 'Yes' or 'No' to query sub-sequence specification:  """)
        if user_input2 == 'Yes' and user_input3 == 'Yes':
            print(r"""
    ----------------------------------
    User-specified input (list format)
    ----------------------------------
    Please paste a single list of input values directly at the command line prompt, specifying the following 8 values.
    Press 'Enter' twice to complete.

    1-Location of OUTPUT DIRECTORY for output files
    2-Location of INPUT FILES (directory containing fastq files)
    3-Location of BLASTN EXECUTABLE
    4-Location of BLASTN ALIGNMENT DATABASE DIRECTORY
    5-Prefix common to BLASTN sequence database files
    6-Location of BLASTDBCMD EXECUTABLE
    7-Optional guide RNA sequence(s) to identify in output alignments
    8-Optional sub-sequence(s) to identify in output alignments

    """)
            input_list = []
            stopword = ""
            while True:
                input_str = input()
                if input_str.strip() == stopword:
                    break
                else:
                    input_list.append(input_str)
            output_directory = input_list[0].strip()
            fastq_directory = input_list[1].strip()
            blastn_path = input_list[2].strip()
            db_path = input_list[3].strip()
            db_prefix = input_list[4].strip()
            blastdbcmd_path = input_list[5].strip()
            guideRNA_seq = [i.strip() for i in input_list[6].split(',')]
            extant_seq = [i.strip() for i in input_list[7].split(',')]
        elif user_input2 == 'Yes' and user_input3 == 'No':
            print(r"""
    ----------------------------------
    User-specified input (list format)
    ----------------------------------
    
    Please paste a single list of input values directly at the command line prompt, specifying the following 7 values.
    Press 'Enter' twice to complete.
    
    1-Location of OUTPUT DIRECTORY for output files
    2-Location of INPUT FILES (directory containing fastq files)
    3-Location of BLASTN EXECUTABLE
    4-Location of BLASTN ALIGNMENT DATABASE DIRECTORY
    5-Location of BLASTDBCMD EXECUTABLE
    6-Prefix common to BLASTN sequence database files
    7-Optional guide RNA sequence(s) to identify in output alignments
    
    """)
            input_list = []
            stopword = ""
            while True:
                input_str = input()
                if input_str.strip() == stopword:
                    break
                else:
                    input_list.append(input_str)
            output_directory = input_list[0].strip()
            fastq_directory = input_list[1].strip()
            blastn_path = input_list[2].strip()
            db_path = input_list[3].strip()
            db_prefix = input_list[4].strip()
            blastdbcmd_path = input_list[5].strip()
            guideRNA_seq = [i.strip() for i in input_list[6].split(',')]
        elif user_input2 == 'No' and user_input3 == 'Yes':
            print(r"""
    ----------------------------------
    User-specified input (list format)
    ----------------------------------
    
    Please paste a single list of input values directly at the command line prompt, specifying the following 7 values.
    Press 'Enter' twice to complete.
    
    1-Location of OUTPUT DIRECTORY for output files
    2-Location of INPUT FILES (directory containing fastq files)
    3-Location of BLASTN EXECUTABLE
    4-Location of BLASTN ALIGNMENT DATABASE DIRECTORY
    5-Prefix common to BLASTN sequence database files
    6-Location of BLASTDBCMD EXECUTABLE
    7-Optional sub-sequence(s) to identify in output alignments
    
    """)  
            input_list = []
            stopword = ""
            while True:
                input_str = input()
                if input_str.strip() == stopword:
                    break
                else:
                    input_list.append(input_str)
            output_directory = input_list[0].strip()
            fastq_directory = input_list[1].strip()
            blastn_path = input_list[2].strip()
            db_path = input_list[3].strip()
            db_prefix = input_list[4].strip()
            blastdbcmd_path = input_list[5].strip()
            extant_seq = [i.strip() for i in input_list[6].split(',')]
    elif test_seq == 'No':
        print(r"""
    ----------------------------------
    User-specified input (list format)
    ----------------------------------
    
    Please paste a single list of input values directly at the command line prompt, specifying the following 6 values.
    Press 'Enter' twice to complete.
    
    1-Location of OUTPUT DIRECTORY for output files
    2-Location of INPUT FILES (directory containing fastq files)
    3-Location of BLASTN EXECUTABLE
    4-Location of BLASTN ALIGNMENT DATABASE DIRECTORY
    5-Prefix common to BLASTN sequence database files
    6-Location of BLASTDBCMD EXECUTABLE
    
    """)  
        input_list = []
        stopword = ""
        while True:
            input_str = input()
            if input_str.strip() == stopword:
                break
            else:
                input_list.append(input_str)
        output_directory = input_list[0].strip()
        fastq_directory = input_list[1].strip()
        blastn_path = input_list[2].strip()
        db_path = input_list[3].strip()
        db_prefix = input_list[4].strip()
        blastdbcmd_path = input_list[5].strip()

# Wait to create the directories and files until after input has been reviewed and accepted.
# Convert fastq_directory input to operating system-appropriate filepath.
output_directory = Path(str(output_directory))
# Convert fastq_directory input to operating system-appropriate filepath.
fastq_directory = Path(str(fastq_directory))
# Convert blastn_path input to operating system-appropriate filepath.
blastn_path = Path(str(blastn_path))
# Convert db_path input to operating system-appropriate filepath.
db_path = Path(str(db_path))
# Convert blastdbcmd_path input to operating system-appropriate filepath.
blastdbcmd_path = Path(str(blastdbcmd_path))

# Collect fastq files from directory
myFastqFilenames = [file for file in glob.glob(str(fastq_directory)+'/*') if Path(file).suffix in [".gz",".fastq"]]

#Sort fastq file names
myFastqFilenames = sorted(myFastqFilenames)
    
# Double-check whether entries look good:
print("""
---------------------------------------------------------------
Preparation for output:
Please double-check that your inputs were recorded as expected.
---------------------------------------------------------------""")

print("""
Your OUTPUT DIRECTORY was recorded as:
""")
print(str(output_directory))

print("""
Your directory containing fastq INPUT FILES was recorded as:
""")
print(str(fastq_directory))

print("""
Please wait a moment while the following file properties are retrieved and/or
calculated across the fastq files to be processed:

* Illumina sequencing run ID(s)
* Total number of fastq files
* Total number of sequencing reads
* Size distribution of fastq files

""")

# Collect Illumina run IDs from fastq files, consolidate to unique run IDs
runIDlist = []
for sourcefile in myFastqFilenames:
    if Path(sourcefile).suffix == ".gz":
        with gzip.open(sourcefile, "rt") as f:
            runID = ":".join(f.readline().split(":",-2)[:2])
            if not runID in runIDlist:
                runIDlist.append(runID) 
    elif Path(sourcefile).suffix == ".fastq":
        with open(sourcefile, "r") as f:
            runID = ":".join(f.readline().split(":",-2)[:2])
            if not runID in runIDlist:
                runIDlist.append(runID)

# Collect total read counts for fastq files
readcount = []
for sourcefile in myFastqFilenames:
    if Path(sourcefile).suffix == ".gz":
        with gzip.open(sourcefile, "rt") as f:    
            readcount.append(int(len((f).readlines())/4))
    elif Path(sourcefile).suffix == ".fastq":
        with open(sourcefile, "r") as f:
            readcount.append(int(len((f).readlines())/4))
        
# Collect file sizes for fastq files
filesize = []
for sourcefile in myFastqFilenames:
    if Path(sourcefile).suffix == ".gz":
        with gzip.open(sourcefile, "rt") as f:
            filesize.append(round((os.path.getsize(sourcefile)/1048576),5))
    elif Path(sourcefile).suffix == ".fastq":
        filesize.append(round((os.path.getsize(sourcefile)/1048576),5))

# fastq_overview prepares summation of fastq file names, their sizes, and read counts, to be reported in script_metrics.txt    
fastq_overview = list(zip(myFastqFilenames, filesize, readcount))

print("""The following data were collected:  """)
print("    Illumina sequencing run ID(s): ")
for i in runIDlist:
    print('        '+i)

print("    # of fastq files to process: {0}".format(len(myFastqFilenames)))

print("    size distribution of fastq files to process: \n      total... "+str(round((sum(file for file in filesize))))+' MB \n      range... max: '+str(round((max(file for file in filesize)),2))+' MB; min: '+str(round((min(file for file in filesize)),5))+' MB; median: '+str(round((numpy.median([file for file in filesize])),3))+' MB; mean +/- stdev: '+str(round((numpy.mean([file for file in filesize])),3))+' +/- '+str(round((numpy.std([file for file in filesize])),3))+' MB')

print("    read distribution within fastq files to process: \n      total... "+locale.format_string("%d", sum(readcount), grouping=True)+' reads \n      range... max: '+str((max(file for file in readcount)))+' reads; min: '+str((min(file for file in readcount)))+' reads; median: '+str((numpy.median([file for file in readcount])))+' reads; mean +/- stdev: '+str(round((numpy.mean([file for file in readcount]))))+' +/- '+str(round((numpy.std([file for file in readcount]))))+' reads')

# Activate this code block if wish to print total read # for each individual fastq file
#print("fastq files to process for genotyping:")
#for sourcefile in myFastqFilenames:
#    with open(sourcefile, "r") as f:
#        reads = int(len((f).readlines())/4)
#    print('    '+sourcefile, reads)

print("""
Your BLASTN EXECUTABLE location was recorded as:
""")
print(str(blastn_path))

print("""
Your BLASTN ALIGNMENT DATABASE DIRECTORY was recorded as:
""")
print(str(db_path))

print("""
Your PREFIX common to BLASTN ALIGNMENT DATABASE FILES was recorded as:
""")
print(db_prefix)

print("""
Your BLASTDBCMD EXECUTABLE location was recorded as:
""")
print(str(blastdbcmd_path))

if test_seq == 'Yes':
    print("""
Your DNA sub-sequences were recorded as:
""")
    print(guideRNA_seq)
    print(extant_seq)
else:
    pass

check = input("""
Is this list accurately recorded? Type 'Y' or 'N': 
""")

if check == 'Y':
    pass
elif check == 'N':
    checkup = input("""
If you have corrections to make, please quit the active script and start again.
To continue in the script, type 'Continue' and press enter.
To quit the script, type 'Exit' and press enter, or press 'Ctrl+C'.  """)
    if checkup == 'Exit':
        exit(0)
    elif checkup == 'Continue':
        pass

# Present option to include ('Y') or bypass ('N') frequency plot generation (optional file output, allele_evidence.pdf)
frequency_plot_check = input("""
Genotypes.py is ready to process fastq files. Before script operations begin, please indicate whether visual
plots of allele frequencies should be rendered and delivered in an output file, allele_evidence.pdf.

Note that production of allele_evidence.pdf can require hours of processing time, although the output timing of
key text files with allele definitions and genotype inferrals (e.g., allele_definitions.txt, genotypes.txt,
allele_definitions.csv, population_summary.txt) will not be affected.

To PROCEED with script operations that INCLUDE allele_evidence.pdf, type 'Y';

To BYPASS script operations that generate allele_evidence.pdf, type 'N': 
""")

# Start the clock on script operation duration
startTime = datetime.now()
startTimestr = str(startTime).split(' ')[1].split('.')[0]
    
# Proceed to file processing
# Generate the directory and its files (to accept content later in script)
path = str(output_directory)

if not os.path.exists(path):
    os.makedirs(path)
    
output_path = Path(output_directory)

# Create output files
if frequency_plot_check == 'Y':
    filename_list = ['fasta.fa', 'blastn_alignments.txt', 'allele_definitions.txt', 'allele_evidence.pdf', 'genotypes.txt', 'population_summary.txt', 'allele_definitions.csv', 'script_metrics.txt']
elif frequency_plot_check == 'N':
    filename_list = ['fasta.fa', 'blastn_alignments.txt', 'allele_definitions.txt', 'genotypes.txt', 'population_summary.txt', 'allele_definitions.csv', 'script_metrics.txt']

# Define current date as prefix to all filenames
processdate = datetime.today().strftime("%m%d%Y")
      
for filename in filename_list:
    with open(os.path.join(path, processdate+'_'+filename), 'wb') as file:
        pass
    
# Collect RAM info
mem = virtual_memory()
ramem = mem.total/1073741824

# The file script_metrics.txt records script operation metadata (summarizes script input and performance)
# Use print redirection to write to target file, in append mode (begin script_metrics.txt)
filename = Path(str(output_path)+'/'+processdate+'_script_metrics.txt')
with open(filename, 'a') as f:
    print("""Genotypes.py: Script Metrics\nDate: """ + (datetime.today().strftime("%m/%d/%Y")) +
"""\n\nOperating system information:
    name: """ + socket.gethostname() +
'\n    platform: ' + platform.platform() +
'\n    RAM (GB): ' + str(ramem) +
'\n    physical CPU/effective CPU: ' + str(psutil.cpu_count(logical=False)) +'/'+ str(psutil.cpu_count()) +
'\n    executable: ' + psutil.Process().exe() +
"""\n\nUser-entered variables:
    output_directory: """+ str(output_directory) +
"\n    fastq_directory: "+ str(fastq_directory) +
"\n    blastn_path: "+ str(blastn_path) +
"\n    db_path: "+ str(db_path) +
"\n    db_prefix: "+ db_prefix +
"\n    blastdbcmd_path: "+ str(blastdbcmd_path), file = f)
    try:
        guideRNA_seq
    except NameError:
        print('    guideRNA_seq: not defined', file = f)
    else:
        print("    guideRNA_seq: "+ str(guideRNA_seq).strip('[]').replace("'",""), file = f)
    try:
        extant_seq
    except NameError:
        print('    extant_seq: not defined', file = f)
    else:
        print("    extant_seq: "+ str(extant_seq).strip('[]').replace("'",""), file = f)
    print("""\nfastq file information:
    Illumina sequencing run ID(s): """+ str(runIDlist).strip('[]').replace("'","") +
"\n    Number of fastq files processed: "+ str(len(myFastqFilenames)) +
"""\n    Size distribution of fastq files processed: 
        total... """ +str(round((sum(file for file in filesize))))+' MB \n        range... max: '+str(round((max(file for file in filesize)),2))+' MB; min: '+str(round((min(file for file in filesize)),5))+' MB; median: '+str(round((numpy.median([file for file in filesize])),3))+' MB; mean +/- stdev: '+str(round((numpy.mean([file for file in filesize])),3))+' +/- '+str(round((numpy.std([file for file in filesize])),3))+' MB'
"\n    Read distribution within fastq files to process: \n        total... "+locale.format_string("%d", sum(readcount), grouping=True)+' reads \n        range... max: '+str((max(file for file in readcount)))+' reads; min: '+str((min(file for file in readcount)))+' reads; median: '+str((numpy.median([file for file in readcount])))+' reads; mean +/- stdev: '+str(round((numpy.mean([file for file in readcount]))))+' +/- '+str(round((numpy.std([file for file in readcount]))))+' reads', file = f)
    print("\nfastq files processed (name, size (MB), reads): ", file = f)
    for i in (sorted(fastq_overview)):
        print("    " + str(i).strip("()").replace("'",""), file = f)
f.close()

print("""
Script is now processing the top 10 unique read sequences for each fastq file
(read identity and count [% of reads in sample]).
    
The output of this step is a fasta file (.fa) that will be created in the OUTPUT DIRECTORY you indicated.""")

# Start the clock on read count operation duration
startTime_readcount = datetime.now()

# For each fastq file (sourcefile) in fastq_directory, count top 10 most abundant read types and direct read sequence + annotation defline (sample name + frequency metrics) to fasta.fa (future alignment input)
query_input = Path(str(output_path)+'/'+processdate+'_fasta.fa')

# define Nextera adaptor sequence, in preparation to trim read 3' ends if necessary
adaptor_str = 'CTGTCTCTTATACACATCT'
adaptor_str_rev = 'AGATGTGTATAAGAGACAG'

# proceed to R1+R2 merge, unique read identification & counting
    
# Merge R1 and R2 reads, if present, as single sequence
R1_file_list = [sourcefile for sourcefile in myFastqFilenames if bool(re.split('_',os.path.basename(sourcefile))[3] == 'R1')]
R2_file_list = [sourcefile for sourcefile in myFastqFilenames if bool(re.split('_',os.path.basename(sourcefile))[3] == 'R2')]
      
# R1, R2 cluster mapping (identify fastq file pairs representing R1 & R2 data for individual clusters)
processed_files_list = []
R1_R2_map_list = []
for sourcefile in myFastqFilenames:
    if sourcefile in processed_files_list:
        pass
    else:
        testname = ''.join(re.split('_',os.path.basename(sourcefile))[0:3])
        for sourcefile1 in R1_file_list:
            if testname == ''.join(re.split('_',os.path.basename(sourcefile1))[0:3]):
                R1 = sourcefile
                if sourcefile not in processed_files_list:
                    processed_files_list.append(sourcefile)
        for sourcefile2 in R2_file_list:
            if testname == ''.join(re.split('_',os.path.basename(sourcefile2))[0:3]):
                R2 = sourcefile2
                if sourcefile2 not in processed_files_list:
                    processed_files_list.append(sourcefile2)
        R1_R2_map_list.append((R1, R2))
        
for file_pair in R1_R2_map_list:
    R1_file = file_pair[0]
    R2_file = file_pair[1]
    fastaname = re.split('_', os.path.basename(R1_file))
    cluster_sequence_R1_dict = {}
    cluster_sequence_R2_dict = {}
    cluster_sequence_R2_revcomp_dict = {}
    cluster_merged_R1_R2revcomp_dict = {}
    cluster_merged_R1_R2revcomp_dict2 = {}
    merged_read_list = []
    counter=()
    if Path(R1_file).suffix == ".gz":
        with gzip.open(R1_file, "rt") as f:
            lines_R1 = f.readlines()
    elif Path(R1_file).suffix == ".fastq":
        with open(R1_file, 'r') as f:
            lines_R1 = f.readlines()    
    for x in range(0,len(lines_R1),4):
        # trim adaptor sequence and up to 3' end of read from R1 sequence, if adaptor sequence found
        cluster_sequence_R1_dict[lines_R1[x].split(':')[5]+':'+lines_R1[x].split(':')[6].split(' ')[0]] = lines_R1[x+1].strip('\n')[:lines_R1[x+1].strip('\n').index(adaptor_str)] if adaptor_str in lines_R1[x+1].strip('\n') else lines_R1[x+1].strip('\n') 
    #cluster_IDs_list_R1 = [x.split(':')[5]+':'+x.split(':')[6].split(' ')[0] for x in lines_R1[0::4]]
    if Path(R2_file).suffix == ".gz":
        with gzip.open(R2_file, "rt") as f:
            lines_R2 = f.readlines()
    elif Path(R2_file).suffix == ".fastq":
        with open(R2_file, 'r') as f:
            lines_R2 = f.readlines()
    for x in range(0,len(lines_R2),4):
        # trim adaptor sequence and up to 3' end of read from R2 sequence, if adaptor sequence found
        cluster_sequence_R2_dict[lines_R2[x].split(':')[5]+':'+lines_R2[x].split(':')[6].split(' ')[0]] = lines_R2[x+1].strip('\n')[:lines_R2[x+1].strip('\n').index(adaptor_str)] if adaptor_str in lines_R2[x+1].strip('\n') else lines_R2[x+1].strip('\n') 
    #cluster_IDs_list_R2 = [x.split(':')[5]+':'+x.split(':')[6].split(' ')[0] for x in lines_R2[0::4]]
    for cluster in cluster_sequence_R2_dict:
        cluster_sequence_R2_revcomp_dict[cluster] = ''.join(reversed(''.join(nt_dict.get(nt) for nt in cluster_sequence_R2_dict.get(cluster))))
    for cluster in cluster_sequence_R1_dict:
        if cluster in cluster_sequence_R2_revcomp_dict:
            if merge(cluster_sequence_R1_dict.get(cluster), cluster_sequence_R2_revcomp_dict.get(cluster)) != 'no overlap':
                cluster_merged_R1_R2revcomp_dict[cluster] = merge(cluster_sequence_R1_dict.get(cluster), cluster_sequence_R2_revcomp_dict.get(cluster))
            else:
                cluster_merged_R1_R2revcomp_dict2[cluster] = merge1(cluster_sequence_R1_dict.get(cluster), cluster_sequence_R2_revcomp_dict.get(cluster))
    for cluster in cluster_merged_R1_R2revcomp_dict:
        merged_read_list.append(cluster_merged_R1_R2revcomp_dict.get(cluster))
    # create dictionary (counter) relating unique read sequence to its # of occurrences
    counter=Counter(merged_read_list)
    # assign top 10 reads by count in fastq file (sourcefile) to modified_read_list_top10
    modified_read_list_top10 = []
    for index, i in enumerate(counter.most_common(10)):
        # read frequency relative to other reads that occur at >1% raw frequency
        filtered1 = sum([x for x in counter.values() if x/(sum(counter.values())) > 0.01])
        # read frequency relative to other reads that occur at >01% raw frequency
        filtered10 = sum([x for x in counter.values() if x/(sum(counter.values())) > 0.1])
        # read raw frequency
        raw_freq = round((100*i[1]/sum(counter.values())),2)
        modified_read_list_top10.append([i[0], '['+str(i[1])+'/'+str(sum(counter.values()))+']', 'rank'+str(index+1), raw_freq, int(stats.percentileofscore([i for i in counter.values()], i[1], 'rank')), round((100*i[1]/sum([i[1] for i in counter.most_common(10)])),2), round((100*i[1]/filtered1),2) if filtered1 > 0 and raw_freq >= 1 else 'None', round((100*i[1]/filtered10),2) if filtered10 > 0 and raw_freq >= 10 else 'None'])
    # direct output in fasta format (with defline encoding sample name + frequency metrics) to fasta.fa
    with open(str(query_input), 'a+') as file:
        for i in modified_read_list_top10:
            file.write('>'+fastaname[0]+'_'+'R1+R2'+'_'+str(i[1])+'_'+i[2]+'_%totalreads:'+str(i[3])+'_percentile:'+str(i[4])+'_%top10reads:'+str(i[5])+'_%readsfilteredfor1%:'+str(i[6])+'_%readsfilteredfor10%:'+str(i[7])+'\n'+i[0]+'\n')

# Log read count time duration
readcountDuration = str(datetime.now()- startTime_readcount).split(':')[0]+' hr|'+str(datetime.now() - startTime_readcount).split(':')[1]+' min|'+str(datetime.now() - startTime_readcount).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_readcount).split(':')[2].split('.')[1]+' microsec'

print("""
Script is now using BLASTN (NCBI) to process alignments to reference sequence database.""")

# Start the clock on blastn alignments duration
startTime_alignments = datetime.now()

# Process alignments relative to reference sequence database, using blastn
# Reference database
db_input = db_path / db_prefix

# Alignment output
blast_directory = processdate+'_blastn_alignments.txt'
query_output = output_directory / blast_directory
                   
# Alignment command
cmd_align = str(blastn_path)+' -query '+str(query_input)+' -db '+str(db_input)+' -out '+str(query_output)+' -gapopen 1 -gapextend 1 -outfmt "5"'

os.system(cmd_align)

# Log alignment time duration
alignmentsDuration = str(datetime.now()- startTime_alignments).split(':')[0]+' hr|'+str(datetime.now() - startTime_alignments).split(':')[1]+' min|'+str(datetime.now() - startTime_alignments).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_alignments).split(':')[2].split('.')[1]+' microsec'

print("""
Script is now defining alleles.""")

# Start the clock on genotypes inference duration
startTime_imputation = datetime.now()

# Impute genotypes from blastn_alignments.txt file
# Import blastn alignments output as a list of strings (each string corresponds to a query alignment)
alignments_list = []
with open(str(query_output), 'r') as file:
    reader = file.read()
    for i,part in enumerate(reader.split('<Iteration_iter-num>')):
        alignments_list.append(part)
# Remove blastn header line from alignments_list
alignments_list = alignments_list[1:]

# Convert alignments_list to list of lists (i.e., each query alignment string is encapsulateed into its own sublist within alignments_list2)
alignments_list2 = [alignments_list[i:i+1] for i in range(0, len(alignments_list))]

# Identify & subset queries for which no alignments were found in reference database ('no hits found')
no_hits_list = []
for i in alignments_list2:
    if re.search('No hits found', str(i)):
        no_hits_list.append(str(i).split('<Iteration_query-def>')[1].split('</Iteration_query-def>')[0])

# Record sample names having reads with no alignment hits
no_hits_samplename_list = []
for i in no_hits_list:
    samplename = i.split('_')[0]
    if samplename not in no_hits_samplename_list:
        no_hits_samplename_list.append(samplename)

# Within each sublist of alignments_list2, split each line into an individual string, remove beginning and trailing whitespace, and recapture specified subset of alignment information in alignments_list3
alignments_list3 = []
for i in alignments_list2:
    if str(i).split('<Iteration_query-def>')[1].split('</Iteration_query-def>')[0] not in no_hits_list:
        alignments_list3.append([y.strip() for x in i for y in x.split('\n') if y.strip().startswith(('<Iteration_query-ID>', '<Iteration_query-def>', '<Hit_num>', '<Hit_id>', '<Hit_def>', '<Hit_accession>', '<Hsp_num>', '<Hsp_hit-from>', '<Hsp_hit-to>', '<Hsp_qseq>', '<Hsp_hseq>', '<Hsp_midline>'))])
    
# Identify & subset reads with >1 alignment to sequences in reference database
# Some reads with >1 alignment will be recovered to 'reconstitute' hypothesized allele (if BLASTN has split the read into multiple 'hits' or 'high-scoring pairs' (hsp's) within 1 kb)

# There are in principle at least 3 ways a read could potentially align to >1 position in reference database (1 & 2a,b below):
# (1) same sequence span aligns to >1 different locus (disparate <Hit_id>'s)
# (2) one sequence span may be split into two (ore more) different alignment matches, because of intervening gap(s) or insertion(s) that exceed ~60 bp (an apparent BLASTN gap limit)
#    (a) if the two (or more) 'split matches' align to the same <Hit_id>, but to different chromosomal coordinates of that <Hit_id>, they will be presented by BLASTN as belonging to the same <Hit_num>, but to different <Hsp_num> (Hsp=high scoring pair)
#    (b) if the two (or more) 'split matches' span different <Hit_id>'s (essentially different 'chunks' of sequence with unique names, as organized within the alignment database), they will be presented by BLASTN as belonging to different <Hit_num>
# These observations suggest that it is important to distinguish a read with alignment to >1 sequence as either one with poor genomic target resolution, vs. one that harbors sizeable deletions or insertions relative to the reference sequence
# One way to make this distinction may be to gauge how close in chromosomal space the >1 alignments are to one another, and assume their continuity if their proximity falls within a defined constraint.  Genotypes.py assigns this constraint to be 1000 bp (1kb)
# Genotypes.py attempts to reconstitute hypothesized alleles that span multiple non-overlapping hsp's (and does not attempt to reconstitute across multiple hits or for ambiguous reconstructions from overlapping hsp's)

# Organize reads with multiple hit IDs
# These reads are deprecated (not further analyzed)
multiple_alignments_hits_list = []
for i in alignments_list3:
    if len(re.findall('<Hit_num>', str(i))) > 1:
        multiple_alignments_hits_list.append(i)
        
# Organize reads with single hit, but multiple associated hsp IDs
# These reads will be processed by BLASTDBCMD in an effort to 'reconstitute' potential alleles, with high-scoring alignment pairs matched to the reference genome, but split into separate matches due to intervening non-aligning span between the alignment matches
multiple_alignments_hsp_list = []
for i in alignments_list3:
    if len(re.findall('<Hit_num>', str(i))) > 1:
        pass
    elif len(re.findall('<Hsp_num>', str(i))) > 1:
        multiple_alignments_hsp_list.append(i)

# Identify read IDs with >1 alignment to sequences in reference database (hit class)
# Note, in Genotypes.py v1.0, these reads are deprecated at this time (in step just prior)
multiple_alignments_hits_readID_list = []
for i in multiple_alignments_hits_list:
    multiple_alignments_hits_readID_list.append(i[1].split('>')[1].split('<')[0])

# Identify read IDs with >1 alignment to sequences in reference database (hsp class)
multiple_alignments_hsp_readID_list = []
for i in multiple_alignments_hsp_list:
    multiple_alignments_hsp_readID_list.append(i[1].split('>')[1].split('<')[0])

# Record sample names having reads with >1 alignment to sequences in reference database (hit class)
multiple_alignments_hits_samplename_list = []
for i in multiple_alignments_hits_readID_list:
    samplename = i.split('_')[0]
    if samplename not in multiple_alignments_hits_samplename_list:
        multiple_alignments_hits_samplename_list.append(samplename)
        
# Record sample names having reads with >1 alignment to sequences in reference database (hsp class)
multiple_alignments_hsp_samplename_list = []
for i in multiple_alignments_hsp_readID_list:
    samplename = i.split('_')[0]
    if samplename not in multiple_alignments_hsp_samplename_list:
        multiple_alignments_hsp_samplename_list.append(samplename)

# Prepare dictionary linking sample names to their reads having >1 alignment to sequences in reference database (hit class)
multiple_alignments_hits_dict = {}
for i in multiple_alignments_hits_samplename_list:
    multiple_alignments_hits_dict["{0}".format(i)] = tuple(x for x in multiple_alignments_hits_list if bool(re.search(i, x[1])))
    
# Prepare dictionary linking sample names to their reads having >1 alignment to sequences in reference database (hsp class)
multiple_alignments_hsp_dict = {}
for i in multiple_alignments_hsp_samplename_list:
    multiple_alignments_hsp_dict["{0}".format(i)] = tuple(x for x in multiple_alignments_hsp_list if bool(re.search(i, x[1])))
    
# Not all top 10 ranked alleles for every sample will be subject to a multiple-hsp search; will need to keep track of this, to regroup allele alignments for individual samples at later step (i.e., regroup and print by rank alleles that have been reconstituted by BLASTDBCMD as well as alleles that did not require reconstitution)

# For reads that have hsp's that overlap with one another, it is unclear how best to reconstitute an alignment; discard these reads from further analysis ('deprecate', but record their deprecation in population_summary.txt)

multiple_alignments_hsp_dict_valid = {}
multiple_alignments_hsp_dict_invalid = {}
for i in multiple_alignments_hsp_dict:
    valid_hsp_reads_list = []
    invalid_hsp_reads_list = []
    for x in multiple_alignments_hsp_dict.get(i):
        hsp_id_list = x[6::6]
        hsp_from_to_list = []
        hsp_sequence_identifier = x[5]
        # assign sequence identifier
        for hsp in hsp_id_list:
            hsp_index = x.index(hsp)
            hsp_from_to_list.append((x[hsp_index+1], x[hsp_index+2]))
        # check whether hsp spans overlap; if they do, discard read
        coordinate_range_subset_list = []
        for index, span in enumerate(hsp_from_to_list):
            coordinate_range_subset = sorted([int(span[0].split('>')[1].split('<')[0])]+[int(span[1].split('>')[1].split('<')[0])])
            coordinate_range_subset_list.append(coordinate_range_subset)
        coordinate_range_subset_list_expanded = [range(coordinate_range_subset[0],coordinate_range_subset[1]) for coordinate_range_subset in coordinate_range_subset_list]
        # set default overlap validity to True (no hsp overlaps)
        valid_hsp_overlap = True
        for y in coordinate_range_subset_list_expanded:
            coordinate_range_subset_list_expanded_copy = coordinate_range_subset_list_expanded.copy()
            coordinate_range_subset_list_expanded_copy.remove(y)
            for w in coordinate_range_subset_list_expanded_copy:
                if set(y).intersection(set(w)):
                    valid_hsp_overlap = False
        if valid_hsp_overlap is True:
            valid_hsp_reads_list.append(x)
        else:
            invalid_hsp_reads_list.append(x)
    multiple_alignments_hsp_dict_valid["{0}".format(i)] = valid_hsp_reads_list
    multiple_alignments_hsp_dict_invalid["{0}".format(i)] = invalid_hsp_reads_list

# Make a copy of multiple_alignments_hsp_dict_valid, removing dictionary keys with empty lists (no valid multiple hsp's: hsp(s) either span >1 kb, or overlap coordinates)
multiple_alignments_hsp_dict_valid2 = { k : v for k,v in multiple_alignments_hsp_dict_valid.items() if v}
multiple_alignments_hsp_dict_invalid2 = { k : v for k,v in multiple_alignments_hsp_dict_invalid.items() if v}
    
    
print("""
Script is now using BLASTDBCMD (NCBI) to retrieve reference sequences that span multiple high-scoring pairs identified by BLASTN.""")

# BLASTDBCMD: use BLASTN reference genome database for BLASTDBCMD sequence retrieval
# Reference database (same as for blastn in earlier alignment step)
db_input = db_path / db_prefix
    
# Now prepare alignment 'reconstitutions' with aid of BLASTDBCMD
blastdbcmd_alignments_list = []
for index, i in enumerate(multiple_alignments_hsp_dict_valid2):
    for read in multiple_alignments_hsp_dict_valid2.get(i):
        blastdbcmd_alignments_list_read_sublist = []
        intervening_bases = ''
        try:
            retrieved_sequence_updated
        except NameError:
            pass
        else:
            del retrieved_sequence_updated
        # this sublist ultimately has 10 items; 1st 5 come from multiple_alignments_hsp_dict_valid, next 5 come from blastdbcmd & Genotypes.py in upcoming code below
        blastdbcmd_alignments_list_read_sublist = [read[0]]+[read[1]]+[read[2]]+[read[3]]+[read[4]]
        coordinate_range_exceeds_1kb_threshold = False
        # number of high scoring pairs at/across single-hit accession number (sequence identifier)
        hsp_count = len(read[6::6])
        hsp_id_list = read[6::6]
        hsp_from_to_list = []
        hsp_sequence_identifier = read[5]
        # assign sequence identifier
        for hsp in hsp_id_list:
            hsp_index = read.index(hsp)
            hsp_from_to_list.append((read[hsp_index+1], read[hsp_index+2]))
        # sequence_identifier shared by hsp's
        sequence_identifier = hsp_sequence_identifier.split('>')[1].split('<')[0] 
        coordinates_list = []
        for hsp_from_to in hsp_from_to_list:
            for x in hsp_from_to:
                coordinates_list.append(int(x.split('>')[1].split('<')[0]))
        coordinate_range = str(min(set(coordinates_list)))+'-'+str(max(set(coordinates_list)))
    # Check whether the hsp's are within 1 kb of one another.  If not, deprecate as multi-mapping.  If they are within 1 kb, assume that blastn split them as separate hsp's because of intervening deletion/insertion relative to reference sequence.  Attempt to reconstitute.
        if max(set(coordinates_list))-min(set(coordinates_list)) > 1000:
            coordinate_range_exceeds_1kb_threshold = True
            # add this deprecated read to multiple_alignments_hsp_dict_invalid2
            if i in multiple_alignments_hsp_dict_invalid2:
                multiple_alignments_hsp_dict_invalid2[i].append(read)
            else:
                multiple_alignments_hsp_dict_invalid2["{0}".format(i)] = [read]
    # Extract sequence spanning multiple high-scoring pairs (hsp), (if within 1 kb), using blastdbcmd
    # blastdbcmd format is: blastdbcmd -db database -dbtype nucl -entry sequence_identifier -range coordinate_range
    # for example: blastdbcmd -db database -dbtype nucl -entry chr6 -range 20000-20500
    # Sequence retrieval (blastdbcmd) command
        if coordinate_range_exceeds_1kb_threshold is True:
            pass
        else:
            cmd_retrieve_sequence = str(blastdbcmd_path)+' -db '+str(db_input)+' -dbtype nucl -entry '+sequence_identifier+' -range '+coordinate_range+' -outfmt %s'
            retrieved_sequence = subprocess.check_output(cmd_retrieve_sequence, shell=True).decode('utf8').strip()
        # Get ready to reconstruct alignment
        # for each tuple of from-to coordinates in hsp_from_to_list, convert to integers and sort by coordinates
            coordinates_integer_list = []
            for hsp_from_to in hsp_from_to_list:
                new_tuple = tuple(sorted([int(x.split('>')[1].split('<')[0]) for x in hsp_from_to]))
                coordinates_integer_list.append(new_tuple)
            # sort coordinate spans in sequential order
            coordinates_integer_list_sorted = sorted(coordinates_integer_list)
            # take stock of whether sorting changed the relative order of the alignment blocks (hsp's), which will need to be accounted for in eventual alignment reconstruction effort
            new_index_list = []
            for index, hsp_from_to in enumerate(hsp_from_to_list):
                new_tuple = tuple(sorted([int(x.split('>')[1].split('<')[0]) for x in hsp_from_to]))
                new_index = coordinates_integer_list_sorted.index(new_tuple)
                if new_index == index:
                    pass
                else:
                    new_index_list.append(new_index)
            # address whether there needs to be a regrouping of hsp data to correspond to their order relative to an alignment reference span
            # first, if no order adjustment needed:
            if len(new_index_list) == 0:
                hsp_queries_and_midlines_list = []
                for hsp in range(len(hsp_from_to_list)):
                    intervening_bases = ''
                    cognate_read_sequence = ''
                    result = ''
                    hsp_query = read[9+(hsp*6)].split('>')[1].split('<')[0]
                    hsp_midline = read[11+(hsp*6)].split('>')[1].split('<')[0]
                    hsp_hit = read[10+(hsp*6)].split('>')[1].split('<')[0]
                    if '-' in hsp_hit:
                        hsp_hit_adjustment = 'hsp query has insertion(s) relative to hit'
                    else:
                        hsp_hit_adjustment = 'no hsp query insertion(s) relative to hit'
                    if hsp < len(hsp_from_to_list)-1:
                        gap_to_next_hsp = coordinates_integer_list_sorted[hsp+1][0] - coordinates_integer_list_sorted[hsp][1] - 1
                    else:
                        gap_to_next_hsp = 0
                    # Determine whether there is/are insertion(s) to account for in read relative to hsp's
                    hsp_crosscheck_back_to_read = []
                    hsp_sequence_identifier = read[5]
                    read_ID = '>'+read[1].split('>')[1].split('<')[0]
                    # append read_ID (as fastaname) to hsp_crosscheck_back_to_read (making it index 0)
                    hsp_crosscheck_back_to_read.append(read_ID)
                    hsp_sequences_to_join = []
                    for hsp in hsp_id_list:
                        hsp_index = read.index(hsp)
                        hsp_sequences_to_join.append((read[hsp_index+3].split('>')[1].split('<')[0]))
                    hsp_sequences_joined = ''.join(hsp_sequences_to_join)
                    # append hsp match sequences to hsp_crosscheck_back_to_read as a list (making them index 1)
                    hsp_crosscheck_back_to_read.append(hsp_sequences_to_join)
                    # refer to initial fasta.fa file (query_input) to retrieve cognate read sequence
                    with open(query_input) as file_iterator:
                        for line in file_iterator:
                            if read_ID in line:
                                cognate_read_sequence = next(file_iterator).strip()
                    # append cognate read sequence to hsp_crosscheck_back_to_read (making it index 2)
                    hsp_crosscheck_back_to_read.append(cognate_read_sequence)
                    # find any intervening read bases that were not accounted for in an hsp
                    result = re.search(hsp_crosscheck_back_to_read[1][0]+'(.*)'+hsp_crosscheck_back_to_read[1][1], hsp_crosscheck_back_to_read[2])
                    if result:
                        if hsp in range(0,len(hsp_from_to_list),2):
                            intervening_bases = result.group(1)
                            hsp_queries_and_midlines_list.append((hsp_query,intervening_bases,hsp_midline,hsp_hit,hsp_hit_adjustment,gap_to_next_hsp+len(intervening_bases)))
                            try:
                                retrieved_sequence_updated
                            except NameError:
                                retrieved_sequence_updated = ((result.span(1)[1]-result.span(1)[0])*'-').join([retrieved_sequence_updated[:result.span(1)[0]],retrieved_sequence_updated[result.span(1)[1]:]])
                                # get indices of added
                            else:
                                retrieved_sequence_updated = ((result.span(1)[1]-result.span(1)[0])*'-').join([retrieved_sequence[:result.span(1)[0]],retrieved_sequence[result.span(1)[1]:]])
                        else:
                            hsp_queries_and_midlines_list.append((hsp_query,'',hsp_midline,hsp_hit,hsp_hit_adjustment,gap_to_next_hsp)) 
                    else:
                        hsp_queries_and_midlines_list.append((hsp_query,'',hsp_midline,hsp_hit,hsp_hit_adjustment,gap_to_next_hsp))
                    # Also determine whether insertion(s) exist in hsp qseq relative to hsp hseq: if the hsp hseq sequence has any '-' in it, this means that the qseq has insertion(s) relative to the qseq
                    # Check for '-' in hseq for any of the hsp's. If any hsp hseq has '-', to accurately account for these changes in this hseq span, need to lift this qseq span entirely from the blastn file.
                    
            # alternatively, adjust hsp call series based on new index if needed:
            else:
                hsp_queries_and_midlines_list = []
                #print('reversed')
                hsp_count = 0
                for hsp in new_index_list:
                    intervening_bases = ''
                    cognate_read_sequence = ''
                    result = ''
                    hsp_query = read[9+(hsp*6)].split('>')[1].split('<')[0]
                    hsp_midline = read[11+(hsp*6)].split('>')[1].split('<')[0]
                    hsp_hit = read[10+(hsp*6)].split('>')[1].split('<')[0]
                    if '-' in hsp_hit:
                        hsp_hit_adjustment = 'hsp query has insertion(s) relative to hit'
                    else:
                        hsp_hit_adjustment = 'no hsp query insertion(s) relative to hit'
                    if hsp_count < len(hsp_from_to_list)-1:
                        gap_to_next_hsp = coordinates_integer_list_sorted[hsp_count+1][0] - coordinates_integer_list_sorted[hsp_count][1] - 1
                    else:
                        gap_to_next_hsp = 0
                    hsp_count = hsp_count+1
                    # Determine whether there is/are insertion(s) to account for in read relative to hsp's
                    hsp_crosscheck_back_to_read = []
                    hsp_sequence_identifier = read[5]
                    read_ID = '>'+read[1].split('>')[1].split('<')[0]
                    # append read_ID (as fastaname) to hsp_crosscheck_back_to_read (making it index 0)
                    hsp_crosscheck_back_to_read.append(read_ID)
                    hsp_sequences_to_join = []
                    for hsp in hsp_id_list:
                        hsp_index = read.index(hsp)
                        hsp_sequences_to_join.append((read[hsp_index+3].split('>')[1].split('<')[0]))
                    hsp_sequences_joined = ''.join(hsp_sequences_to_join)
                    # append hsp match sequences to hsp_crosscheck_back_to_read as a list (making them index 1)
                    hsp_crosscheck_back_to_read.append(hsp_sequences_to_join)
                    # refer to initial fasta.fa file (query_input) to retrieve cognate read sequence
                    with open(query_input) as file_iterator:
                        for line in file_iterator:
                            if read_ID in line:
                                cognate_read_sequence = next(file_iterator).strip()
                    # append cognate read sequence to hsp_crosscheck_back_to_read (making it index 2)
                    hsp_crosscheck_back_to_read.append(cognate_read_sequence)
                    # find any intervening read bases that were not accounted for in an hsp
                    result = re.search(hsp_crosscheck_back_to_read[1][0]+'(.*)'+hsp_crosscheck_back_to_read[1][1], hsp_crosscheck_back_to_read[2])
                    if result:
                        if hsp in range(0,len(hsp_from_to_list),2):
                            intervening_bases = result.group(1)
                            hsp_queries_and_midlines_list.append((hsp_query,intervening_bases,hsp_midline,hsp_hit,hsp_hit_adjustment,gap_to_next_hsp+len(intervening_bases)))
                            try:
                                retrieved_sequence_updated
                            except NameError:
                                retrieved_sequence_updated = ((result.span(1)[1]-result.span(1)[0])*'-').join([retrieved_sequence_updated[:result.span(1)[0]],retrieved_sequence_updated[result.span(1)[1]:]])
                                # get indices of added
                            else:
                                retrieved_sequence_updated = ((result.span(1)[1]-result.span(1)[0])*'-').join([retrieved_sequence[:result.span(1)[0]],retrieved_sequence[result.span(1)[1]:]])
                        else:
                            hsp_queries_and_midlines_list.append((hsp_query,'',hsp_midline,hsp_hit,hsp_hit_adjustment,gap_to_next_hsp)) 
                    else:
                        hsp_queries_and_midlines_list.append((hsp_query,'',hsp_midline,hsp_hit,hsp_hit_adjustment,gap_to_next_hsp))
            # retrieved_sequence
            # reconstructed alignment midline
            reconstructed_midline = ''.join([(hsp_queries_and_midlines_list[hsp][2]+(hsp_queries_and_midlines_list[hsp][5]*' ')) for hsp in range(len(hsp_from_to_list))])
            # reconstructed alignment sequence
            reconstructed_alignment_sequence = ''.join([(hsp_queries_and_midlines_list[hsp][0]+hsp_queries_and_midlines_list[hsp][1]+(hsp_queries_and_midlines_list[hsp][5]*'-')) for hsp in range(len(hsp_from_to_list))])
            # provide reconstituted 'hsp hit from-to'
            blastdbcmd_alignments_list_read_sublist.append('<Hsp_hit-from>'+coordinate_range.split('-')[0]+'</Hsp_hit-from>')
            blastdbcmd_alignments_list_read_sublist.append('<Hsp_hit-to>'+coordinate_range.split('-')[1]+'</Hsp_hit-to>')
            blastdbcmd_alignments_list_read_sublist.append('<Hsp_qseq>'+reconstructed_alignment_sequence+'</Hsp_qseq>')
            try:
                retrieved_sequence_updated
            except NameError:
                # check whether any hsp has a retrieved sequence alignment needing adjustment (in the retrieved sequence) across the hsp span
                adjustment_required = False
                for index, hsp in enumerate(hsp_queries_and_midlines_list):
                    if 'hsp query has insertion(s) relative to hit' in hsp:
                        adjustment_required = True
                if adjustment_required == True:
                    for index, hsp in enumerate(hsp_queries_and_midlines_list):
                        if 'hsp query has insertion(s) relative to hit' in hsp:
                            hit_sequence = hsp_queries_and_midlines_list[index][3].replace('-','')
                            try:
                                retrieved_sequence_updated
                            except NameError:
                                retrieved_sequence_updated = retrieved_sequence.replace(hit_sequence, hsp_queries_and_midlines_list[index][3])
                            else:
                                 retrieved_sequence_updated = retrieved_sequence_updated.replace(hit_sequence, hsp_queries_and_midlines_list[index][3])
                    blastdbcmd_alignments_list_read_sublist.append('<Hsp_hseq>'+retrieved_sequence_updated+'</Hsp_hseq>')
                elif adjustment_required == False:
                    blastdbcmd_alignments_list_read_sublist.append('<Hsp_hseq>'+retrieved_sequence+'</Hsp_hseq>')
            else:
                # check whether any hsp has a retrieved sequence alignment needing adjustment (in the retrieved sequence [already updated with any insertion(s) in the span between hsp's]) across the hsp span
                adjustment_required = False
                for index, hsp in enumerate(hsp_queries_and_midlines_list):
                    if 'hsp query has insertion(s) relative to hit' in hsp:
                        adjustment_required = True
                if adjustment_required == True:
                    for index, hsp in enumerate(hsp_queries_and_midlines_list):
                        if 'hsp query has insertion(s) relative to hit' in hsp:
                            hit_sequence = hsp_queries_and_midlines_list[index][3].replace('-','')
                            try:
                                retrieved_sequence_updated2
                            except NameError:
                                retrieved_sequence_updated2 = retrieved_sequence_updated.replace(hit_sequence, hsp_queries_and_midlines_list[index][3])
                            else:
                                 retrieved_sequence_updated2 = retrieved_sequence_updated2.replace(hit_sequence, hsp_queries_and_midlines_list[index][3])
                    blastdbcmd_alignments_list_read_sublist.append('<Hsp_hseq>'+retrieved_sequence_updated+'</Hsp_hseq>')
                elif adjustment_required == False:
                    blastdbcmd_alignments_list_read_sublist.append('<Hsp_hseq>'+retrieved_sequence_updated+'</Hsp_hseq>')
            blastdbcmd_alignments_list_read_sublist.append('<Hsp_midline>'+reconstructed_midline+'</Hsp_midline>')
            blastdbcmd_alignments_list.append(blastdbcmd_alignments_list_read_sublist)

# Take note of sample IDs with reads that were considered for multiple-hsp reconstruction, but deprecated because
# of overlapping hsp's and/or hsp span exceeding 1kb
multiple_alignments_hsp_invalid2_list = [i for i in multiple_alignments_hsp_dict_invalid2]

# Log alignment time duration
alignmentsDuration = str(datetime.now()- startTime_alignments).split(':')[0]+' hr|'+str(datetime.now() - startTime_alignments).split(':')[1]+' min|'+str(datetime.now() - startTime_alignments).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_alignments).split(':')[2].split('.')[1]+' microsec' 
###################

# Prepare alignment_list4 for reads with exclusively 1 alignment hit in reference database

# alignments_list4 = []
# for i in alignments_list3:
#     if i not in multiple_alignments_list:
#         alignments_list4.append(i)

# for reads & their alignment data in alignments_list3 that were not identified as associated with >1 hsp, trim back each sublist to 10 indices (rather than 12 [remove accession and hsp_num indices, located at indices 5 & 6 of each sublist in alignments_list3]) (make these compatible with alignments_list4)
# First exclude these reads from alignments_list4, because they will be instead delivered via blastdbcmd_alignments_list
multiple_alignments_hsp_list_to_exclude_from_alignments_list4 = []
for i in multiple_alignments_hsp_dict:
    for read in multiple_alignments_hsp_dict.get(i):
        multiple_alignments_hsp_list_to_exclude_from_alignments_list4.append(read)

alignments_list_temp = []
for i in alignments_list3:
    if i not in multiple_alignments_hits_list:
        if i not in multiple_alignments_hsp_list_to_exclude_from_alignments_list4:
            alignments_list_temp.append([i[0]]+[i[1]]+[i[2]]+[i[3]]+[i[4]]+[i[7]]+[i[8]]+[i[9]]+[i[10]]+[i[11]])
            
# Now pool alignment data delivered strictly from BLASTN (alignments_list4), with alignment data 'reconstituted' by BLASTDBCMD (blastdbcmd_alignments_list)
alignments_list4 = alignments_list_temp + blastdbcmd_alignments_list

# Among lists containing alignment data in alignments_list4, determine which queries (reads) correspond to the same sample; where querydef = i[1].split(">")[1].split("_[")[0], reads belonging to the same sample share identical querydef
# Fasta deflines encode frequency metrics for reads, based on defline format:
# sampleID_[reads/total reads]_percentile_% read abundance_% top 10 reads_% reads filtered for 1%_% reads filtered for 10%
querydef_list = []
for i in alignments_list3:
    querydef = i[1].split(">")[1].split("_")[0]
    querydef_list.append(querydef)
    
querydef_uniq_list = []
for i in querydef_list:
    if i in querydef_uniq_list:
        pass
    else:
        querydef_uniq_list.append(i)

# Prepare dictionary relating sample IDs to their associated reads ('alleles')
# Reorder ranked alleles as necessary
alignmentoutput_dict = {}
for i in querydef_uniq_list:
    temp_ranked_read_list = [x for x in alignments_list4 if bool(re.search(i, x[1]))]
    temp_ranked_read_order = []
    new_ranked_read_order = []
    for index, ranked_read in enumerate(temp_ranked_read_list):
        rank = int(ranked_read[1].split('>')[1].split('<')[0].split('_')[3].split('k')[1])
        temp_ranked_read_order.append([rank, ranked_read])
    for read in sorted(temp_ranked_read_order):
        new_ranked_read_order.append(read[1])
    alignmentoutput_dict["{0}".format(i)] = tuple(new_ranked_read_order)

#for i in querydef_uniq_list:
#    alignmentoutput_dict["{0}".format(i)] = tuple(x for x in alignments_list4 if bool(re.search(i, x[1])))
    
# Identify & subset sample ID's that do not have output alleles (empty tuple values in dictionary)
empty_sampleIDs_list = []
for i in alignmentoutput_dict:
    if bool(alignmentoutput_dict.get(i) == ()):
        empty_sampleIDs_list.append(i)

# Make a copy of alignmentoutput_dict, removing dictionary keys with empty tuple values
alignmentoutput_dict2 = { k : v for k,v in alignmentoutput_dict.items() if v}
# Alignmentoutput_dict2 is the key input dictionary for genotype inferences

#for i in alignmentoutput_dict2:
#    for x in range(0,len(alignmentoutput_dict2.get(i))):
#        print(i, alignmentoutput_dict2.get(i)[x][4].split(">")[1].split("<")[0])

print("""
Script is now extrapolating genotypes.""")

# Define nt complement dictionary
nt_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', '-':'-'}

# Prepare dictionary relating sample IDs to their associated 'alleles', allele interpretations/definitions, and inferred genotype
imputedgenotypes_dict = {}
for i in alignmentoutput_dict2:
    imputedgenotypes_dict["{0}".format(i)] = []
    imputed_genotype = []
    allele_data = ['allele_name', 'chr+build', 'locusID', 'coordinates', 'alignment']
    allele_descriptions = ['allele_type', 'allele_specs']
    guideRNA_match = ''
    extant_match = ''
    for x in range(0,len(alignmentoutput_dict2.get(i))):
        imputedgenotypes_dict[i].extend([[]])
        allele_data_x = []
        allele_data_x.append(alignmentoutput_dict2.get(i)[x][1].split(">")[1].split("<")[0].replace('_', ' '))
        allele_data_x.append(alignmentoutput_dict2.get(i)[x][4].split(">")[1].split("<")[0])
        allele_data_x.append(alignmentoutput_dict2.get(i)[x][3].split(">")[1].split("<")[0])
        allele_data_x.append(alignmentoutput_dict2.get(i)[x][5].split(">")[1].split("<")[0]+'-'+alignmentoutput_dict2.get(i)[x][6].split(">")[1].split("<")[0])
        allele_data_x.append('\n'+'    query  '+alignmentoutput_dict2.get(i)[x][7].split(">")[1].split("<")[0]+'\n'+'           '+alignmentoutput_dict2.get(i)[x][9].split(">")[1].split("<")[0]+'\n'+'reference  '+alignmentoutput_dict2.get(i)[x][8].split(">")[1].split("<")[0])
        imputedgenotypes_dict[i][x].append(((dict(zip(allele_data, allele_data_x)))))
        allele_descriptions_x = []
        if bool(re.search(' ', alignmentoutput_dict2.get(i)[x][9])):
            allele_descriptions_x.append('mutant')
            if bool(re.search('-', alignmentoutput_dict2.get(i)[x][7])) and not bool(re.search('-', alignmentoutput_dict2.get(i)[x][8])):
                allele_descriptions_x.append('likely deletion, '+str(alignmentoutput_dict2.get(i)[x][7].count('-'))+' bp')
                imputedgenotypes_dict[i][x].append(((dict(zip(allele_descriptions, allele_descriptions_x)))))
                # imputed genotype is based on alleles with frequency adjusted as relative to reads that occurred at >10% raw frequency
                if imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1] == 'None':
                    pass
                elif float(imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1].strip()) > 10:
                    imputed_genotype.append(imputedgenotypes_dict.get(i)[x][1].get('allele_specs'))
            elif bool(re.search('-', alignmentoutput_dict2.get(i)[x][8])) and not bool(re.search('-', alignmentoutput_dict2.get(i)[x][7])):
                allele_descriptions_x.append('likely insertion, '+str(alignmentoutput_dict2.get(i)[x][8].count('-'))+' bp')
                imputedgenotypes_dict[i][x].append(((dict(zip(allele_descriptions, allele_descriptions_x)))))
                if imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1] == 'None':
                    pass
                elif float(imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1]) > 10:
                    imputed_genotype.append(imputedgenotypes_dict.get(i)[x][1].get('allele_specs'))
            elif bool(re.search('-', alignmentoutput_dict2.get(i)[x][7])) and bool(re.search('-', alignmentoutput_dict2.get(i)[x][8])):
                allele_descriptions_x.append('likely complex indel')
                imputedgenotypes_dict[i][x].append(((dict(zip(allele_descriptions, allele_descriptions_x)))))
                if imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1] == 'None':
                    pass
                elif float(imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1]) > 10:
                    imputed_genotype.append(imputedgenotypes_dict.get(i)[x][1].get('allele_specs'))
            else:
                allele_descriptions_x.append('likely substitution')
                imputedgenotypes_dict[i][x].append(((dict(zip(allele_descriptions, allele_descriptions_x)))))
                if imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1] == 'None':
                    pass
                elif float(imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1]) > 10:
                    imputed_genotype.append(imputedgenotypes_dict.get(i)[x][1].get('allele_specs'))
        else:
            allele_descriptions_x.append('wild-type')
            imputedgenotypes_dict[i][x].append(((dict(zip(allele_descriptions, allele_descriptions_x)))))
            if imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1] == 'None':
                pass
            elif float(imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1]) > 10:
                imputed_genotype.append('wild-type')
        # find guide RNA sequence in reference sequence, and record location so that guide can be printed at appropriate position above reference sequence
        try:
            guideRNA_seq
        except NameError:
            imputedgenotypes_dict[i][x].extend([{}])
        else:
            guideRNA_seq_orientation_list = []
            guide_positions_list = []
            guideRNA_match = ''
            guideRNA_revcomp_match = ''
            for guide in guideRNA_seq:
                if bool(re.search('-', alignmentoutput_dict2.get(i)[x][7])):
                    if bool(re.search('-', alignmentoutput_dict2.get(i)[x][8])):
                        guideRNA_match = re.search(guide, alignmentoutput_dict2.get(i)[x][8].replace('-', '')) or re.search(guide, alignmentoutput_dict2.get(i)[x][7].replace('-', ''))
                    else:
                        guideRNA_match = re.search(guide, alignmentoutput_dict2.get(i)[x][8].replace('-', ''))
                else:
                    guideRNA_match = re.search(guide, alignmentoutput_dict2.get(i)[x][7].replace('-', '')) 
                if not guideRNA_match:
                    guide_revcomp = ''.join(reversed(''.join(nt_dict.get(nt) for nt in guide)))
                    guide_rev = ''.join(reversed(guide))
                    if bool(re.search('-', alignmentoutput_dict2.get(i)[x][7])):
                        if bool(re.search('-', alignmentoutput_dict2.get(i)[x][8])):
                            guideRNA_revcomp_match = re.search(guide_revcomp, alignmentoutput_dict2.get(i)[x][8].replace('-', '')) or re.search(guide_revcomp, alignmentoutput_dict2.get(i)[x][7].replace('-', ''))
                        else:
                            guideRNA_revcomp_match = re.search(guide_revcomp, alignmentoutput_dict2.get(i)[x][8].replace('-', '')) or re.search(guide_revcomp, alignmentoutput_dict2.get(i)[x][7].replace('-', ''))
                    if not guideRNA_revcomp_match:
                        guideRNA_seq_orientation_list.append('guide')
                        guide_positions_list.append('None')
                    else:
                        guideRNA_seq_orientation_list.append(guide_rev)
                        guide_positions_list.append(guideRNA_revcomp_match.start())
                else:
                    guideRNA_seq_orientation_list.append(guide)
                    guide_positions_list.append(guideRNA_match.start())
            imputedgenotypes_dict[i][x].append(((dict(zip(guideRNA_seq_orientation_list, guide_positions_list)))))
        try:
            extant_seq
        except NameError:
            imputedgenotypes_dict[i][x].extend([{}])
        else:
            extant_seq_orientation_list = []
            seq_positions_list = []
            extant_match = ''
            seq_revcomp_match = ''
            for seq in extant_seq:
                # if test sequence is found in forward direction in query or hit:
                extant_match = re.search(seq, alignmentoutput_dict2.get(i)[x][8]) or re.search(seq, alignmentoutput_dict2.get(i)[x][7])
                if not extant_match:
                    seq_revcomp = ''.join(reversed(''.join(nt_dict.get(nt) for nt in seq)))
                    seq_rev = ''.join(reversed(seq))
                    # if test sequence is found in reverse direction in query:
                    seq_revcomp_match = re.search(seq_revcomp, alignmentoutput_dict2.get(i)[x][8]) or re.search(seq_revcomp, alignmentoutput_dict2.get(i)[x][7])
                    if not seq_revcomp_match:
                        extant_match = re.search(seq, alignmentoutput_dict2.get(i)[x][8].replace('-', '')) or re.search(seq, alignmentoutput_dict2.get(i)[x][7].replace('-', ''))
                        if extant_match:
                            pad_match = len(re.findall('-', alignmentoutput_dict2.get(i)[x][8])) or len(re.findall('-', alignmentoutput_dict2.get(i)[x][8]))
                            extant_seq_orientation_list.append(seq)
                            seq_positions_list.append(extant_match.start()+pad_match)
                        if not extant_match:
                            seq_revcomp_match = re.search(seq_revcomp, alignmentoutput_dict2.get(i)[x][8].replace('-', '')) or re.search(seq_revcomp, alignmentoutput_dict2.get(i)[x][7].replace('-', ''))
                            if seq_revcomp_match:
                                pad_match = len(re.findall('-', alignmentoutput_dict2.get(i)[x][8])) or len(re.findall('-', alignmentoutput_dict2.get(i)[x][7]))
                            if not seq_revcomp_match:
                                extant_seq_orientation_list.append(seq)
                                seq_positions_list.append('None')
                            else:
                                extant_seq_orientation_list.append(seq_rev)     
                                seq_positions_list.append(seq_revcomp_match.start())
                    else:
                        extant_seq_orientation_list.append(seq_rev)     
                        seq_positions_list.append(seq_revcomp_match.start()) 
                else:
                    extant_seq_orientation_list.append(seq)
                    seq_positions_list.append(extant_match.start())
            imputedgenotypes_dict[i][x].append(((dict(zip(extant_seq_orientation_list, seq_positions_list))))) 
    # impute genotype based on allele(s) recorded as occurring at >10% frequency in imputed_genotypes list; deliver to index[0] position of imputedgenotypes_dict value for sampleID key
    # homozygous states
    if len(set(imputed_genotype)) == 1:
        if bool('wild-type' in imputed_genotype):
            imputedgenotypes_dict[i].insert(0, '|homozygous| wild-type (wt/wt)')
        for n in set(imputed_genotype):
            if bool(re.search('deletion', n)):
                imputedgenotypes_dict[i].insert(0, '|homozygous| deletion (delta/delta)')
            if bool(re.search('insertion', n)):
                imputedgenotypes_dict[i].insert(0, '|homozygous| insertion (++/++)')
            if bool(re.search('indel', n)):
                imputedgenotypes_dict[i].insert(0, '|homozygous| indel (indel/indel')
            if bool(re.search('substitution', n)):
                imputedgenotypes_dict[i].insert(0, '|homozygous| substitution (sub/sub)')
    #heterozygous states
    elif len(set(imputed_genotype)) == 2:
    # with wild-type allele
        if bool('wild-type' in imputed_genotype):
            for n in set(imputed_genotype):
                if re.search('deletion', n):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| deletion + wild-type (delta/wt)')
                elif re.search('insertion', n):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| insertion + wild-type (++/wt)')
                elif re.search('indel', n):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| indel + wild-type  (indel/wt)')
                elif re.search('substitution', n):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| substitution + wild-type (sub/wt)')
    #heterozygous states (no wild-type allele)
        elif not bool('wild-type' in imputed_genotype):
            genotype_impute_summary = []
            for n in set(imputed_genotype):
                if re.search('wild-type', n):
                    genotype_impute_summary.append('wild-type')
                elif re.search('deletion', n):
                    genotype_impute_summary.append('deletion')
                elif re.search('insertion', n):
                    genotype_impute_summary.append('insertion')
                elif re.search('indel', n):
                    genotype_impute_summary.append('indel')
                elif re.search('substitution', n):
                    genotype_impute_summary.append('substitution')
            if len(set(genotype_impute_summary)) == 1:
                if 'deletion' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| deletion1 + deletion2 (del1/del2)')
                elif 'insertion' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| insertion1 + insertion2 (++1/++2)')
                elif 'indel' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| indel1 + indel2 (indel1/indel2)')
                elif 'substitution' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| substitution1 + substitution2 (sub1/sub2)')
            else:
                if 'deletion' in set(genotype_impute_summary) and 'insertion' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| deletion + insertion (del/++)')
                elif 'deletion' in set(genotype_impute_summary) and 'indel' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| deletion + indel (del/indel)')
                elif 'deletion' in set(genotype_impute_summary) and 'substitution' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| deletion + substitution (del/sub)')
                elif 'insertion' in set(genotype_impute_summary) and 'indel' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| insertion + indel (++/indel)')
                elif 'insertion' in set(genotype_impute_summary) and 'substitution' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| insertion + indel (++/sub)')                    
                elif 'indel' in set(genotype_impute_summary) and 'substitution' in set(genotype_impute_summary):
                    imputedgenotypes_dict[i].insert(0, '|heterozygous| indel + substitution (indel/sub)')
    #multizygous states
    elif len(set(imputed_genotype)) > 2:
        imputed_genotype_str = '|multizygous|'
        genotype_impute_summary = []
        for n in set(imputed_genotype):
            if re.search('wild-type', n):
                genotype_impute_summary.append('wild-type')
            elif re.search('deletion', n):
                genotype_impute_summary.append('deletion')
            elif re.search('insertion', n):
                genotype_impute_summary.append('insertion')
            elif re.search('indel', n):
                genotype_impute_summary.append('indel')
            elif re.search('substitution', n):
                genotype_impute_summary.append('substitution')
        if 'wild-type' in genotype_impute_summary:
            imputed_genotype_str = imputed_genotype_str+' wild-type'
            if 'deletion' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + deletion ('+str(genotype_impute_summary.count('deletion'))+')'
            if 'insertion' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + insertion ('+str(genotype_impute_summary.count('insertion'))+')'
            if 'indel' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + indel ('+str(genotype_impute_summary.count('indel'))+')'
            if 'substitution' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + substitution ('+str(genotype_impute_summary.count('substitution'))+')'
        elif 'deletion' in genotype_impute_summary:
            imputed_genotype_str = imputed_genotype_str+' deletion ('+str(genotype_impute_summary.count('deletion'))+')'
            if 'insertion' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + insertion ('+str(genotype_impute_summary.count('insertion'))+')'
            if 'indel' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + indel ('+str(genotype_impute_summary.count('indel'))+')'
            if 'substitution' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + substitution ('+str(genotype_impute_summary.count('substitution'))+')'
        elif 'insertion' in genotype_impute_summary:
            imputed_genotype_str = imputed_genotype_str+' insertion ('+str(genotype_impute_summary.count('insertion'))+')'
            if 'indel' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + indel ('+str(genotype_impute_summary.count('indel'))+')'
            if 'substitution' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + substitution ('+str(genotype_impute_summary.count('substitution'))+')'
        elif 'indel' in genotype_impute_summary:
            imputed_genotype_str = imputed_genotype_str+' indel ('+str(genotype_impute_summary.count('indel'))+')'
            if 'substitution' in genotype_impute_summary:
                imputed_genotype_str = imputed_genotype_str+' + substitution ('+str(genotype_impute_summary.count('substitution'))+')'
        elif 'substitution' in genotype_impute_summary:
            imputed_genotype_str = imputed_genotype_str+' + substitution ('+str(genotype_impute_summary.count('substitution'))+')'   
        imputedgenotypes_dict[i].insert(0, imputed_genotype_str)
    elif not imputed_genotype:
        imputedgenotypes_dict[i].insert(0, '|unclear or multi-allelic| insufficient representation of any allele (i.e., no allele exceeds >10% of total reads when adjusted for 10% read threshold)')

# Create list containing reversed guideRNA_seq & extant_seq sequences, for use with R2 sequences
try:
    guideRNA_seq
except NameError:
    pass
else:
    guideRNA_seq_rev = [''.join(reversed(i)) for i in guideRNA_seq]
    
try:
    extant_seq
except NameError:
    pass
else:
    extant_seq_rev = [''.join(reversed(i)) for i in extant_seq]
        
# Print summaries of sample-specific allele definitions to output files; first to allele_definitions.txt, preserving sample order
allele_definitions_output = Path(str(output_path)+'/'+processdate+'_allele_definitions.txt')
with open(str(allele_definitions_output), 'a+') as file:
    file.write('Genotypes.py: Allele Definitions\nDate: ' + (datetime.today().strftime("%m/%d/%Y")) + '\n\n')
    for i in querydef_uniq_list:
        if i in imputedgenotypes_dict:
            file.write((len(i)*'=')+'\n'+i+'\n'+(len(i)*'=')+'\n')
            file.write('Inferred Genotype: '+imputedgenotypes_dict.get(i)[0]+'\n')
            # display alleles and their descriptions
            read_checklist = []
            read_abundance_checklist = []
            for n in range(1, len(imputedgenotypes_dict.get(i))):
                if imputedgenotypes_dict.get(i)[n][0].get('allele_name').split(' ')[1] == 'R1+R2':
                    if 'R1+R2' not in read_checklist:
                        read_checklist.append('R1+R2')
                        file.write('\n'+3*' '+'*'+17*'~'+'*\n'+3*' '+'| READ 1 + READ 2 |\n'+3*' '+'*'+17*'~'+'*\n\n')
                    else:
                        pass
                    if float(imputedgenotypes_dict.get(i)[n][0].get('allele_name').split(' ')[4].split(':')[1]) < 10:
                        if 'R1+R2dregs' not in read_abundance_checklist:
                            read_abundance_checklist.append('R1+R2dregs')
                            file.write(3*' '+'*'+56*'~'+'*\n'+3*' '+'|  >>>>> remaining alleles occur at frequency <10% <<<<< |\n'+3*' '+'*'+56*'~'+'*\n\n')
                        else:
                            pass
                if imputedgenotypes_dict.get(i)[n][1].get('allele_type') == 'wild-type':
                    file.write(3*' '+'Allele: '+imputedgenotypes_dict.get(i)[n][0].get('allele_name')+' | '+imputedgenotypes_dict.get(i)[n][1].get('allele_type')+'\n    Locus: '+imputedgenotypes_dict.get(i)[n][0].get('chr+build')+', '+imputedgenotypes_dict.get(i)[n][0].get('locusID')+' '+imputedgenotypes_dict.get(i)[n][0].get('coordinates')+'\n')
                else:
                    file.write(3*' '+'Allele: '+imputedgenotypes_dict.get(i)[n][0].get('allele_name')+' | '+imputedgenotypes_dict.get(i)[n][1].get('allele_type')+', '+imputedgenotypes_dict.get(i)[n][1].get('allele_specs')+'\n    Locus: '+imputedgenotypes_dict.get(i)[n][0].get('chr+build')+', '+imputedgenotypes_dict.get(i)[n][0].get('locusID')+' '+imputedgenotypes_dict.get(i)[n][0].get('coordinates')+'\n')           
                for guide in imputedgenotypes_dict.get(i)[n][2]:
                    if imputedgenotypes_dict.get(i)[n][2].get(guide) != 'None':
                        if guide in guideRNA_seq:
                            file.write('\n'+(int(imputedgenotypes_dict.get(i)[n][2].get(guide))-2)*' '+"5'-"+guide+"-3' (guide sequence)"+'\n'+(int(imputedgenotypes_dict.get(i)[n][2].get(guide))+len(guide)-3)*' '+'v')
                        elif guide in guideRNA_seq_rev:
                            file.write('\n'+(int(imputedgenotypes_dict.get(i)[n][2].get(guide))-2)*' '+"3'-"+guide+"-5' (guide sequence)"+'\n'+(int(imputedgenotypes_dict.get(i)[n][2].get(guide))+4)*' '+'v')
                file.write(imputedgenotypes_dict.get(i)[n][0].get('alignment'))
                for seq in imputedgenotypes_dict.get(i)[n][3]:
                    if imputedgenotypes_dict.get(i)[n][3].get(seq) != 'None':
                        file.write('\n')
                        if seq in extant_seq:
                            file.write((1+int(imputedgenotypes_dict.get(i)[n][3].get(seq)))*' '+len(seq)*'^'+'\n'+(int(imputedgenotypes_dict.get(i)[n][3].get(seq))-2)*' '+"5'-"+seq+"-3' (sequence of interest)\n")
                        elif seq in extant_seq_rev:
                            file.write((1+int(imputedgenotypes_dict.get(i)[n][3].get(seq)))*' '+len(seq)*'^'+'\n'+(int(imputedgenotypes_dict.get(i)[n][3].get(seq))-2)*' '+"3'-"+seq+"-5' (sequence of interest)\n")
                    elif imputedgenotypes_dict.get(i)[n][3].get(seq) == 'None':
                        file.write('\n')
                file.write('\n\n')

# Next print to genotypes.txt, using imputed genotype criteria as basis for reporting order
# First prepare lists that bin sampleIDs based on imputed genotype
# homo genotypes
imputedgenotypes_homowildtype = []
imputedgenotypes_homodeletion = []
imputedgenotypes_homoinsertion = []
imputedgenotypes_homoindel = []
imputedgenotypes_homosubstitution = []
# biallelic mutant genotypes
imputedgenotypes_biallelic_deletion = []
imputedgenotypes_biallelic_insertion = []
imputedgenotypes_biallelic_indel = []
imputedgenotypes_biallelic_substitution = []
imputedgenotypes_biallelic_other = []
# hetero genotypes (containing wt allele)
imputedgenotypes_heterodeletion = []
imputedgenotypes_heteroinsertion = []
imputedgenotypes_heteroindel = []
imputedgenotypes_heterosubstitution = []
# multizygous
imputedgenotypes_multizygous = []
# unclear
imputedgenotypes_unclear = []


for i in imputedgenotypes_dict:
    if imputedgenotypes_dict.get(i)[0] in ('|homozygous| wild-type (wt/wt)'):
        imputedgenotypes_homowildtype.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|homozygous| deletion (delta/delta)'):
        imputedgenotypes_homodeletion.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|homozygous| insertion (++/++)'):
        imputedgenotypes_homoinsertion.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|homozygous| indel (indel/indel)'):
        imputedgenotypes_homoindel.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|homozygous| substitution (sub/sub)'):
        imputedgenotypes_homosubstitution.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|heterozygous| deletion1 + deletion2 (del1/del2)'):
        imputedgenotypes_biallelic_deletion.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|heterozygous| insertion1 + insertion2 (++1/++2)'):
        imputedgenotypes_biallelic_insertion.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|heterozygous| indel1 + indel2 (indel1/indel2)'):
        imputedgenotypes_biallelic_indel.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|heterozygous| substitution1 + substitution2 (sub1/sub2)'):
        imputedgenotypes_biallelic_substitution.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|heterozygous| deletion + insertion (del/++)', '|heterozygous| deletion + indel (del/indel)', '|heterozygous| deletion + substitution (del/sub)', '|heterozygous| insertion + indel (++/indel)', '|heterozygous| insertion + indel (++/sub)', '|heterozygous| indel + substitution (indel/sub)'):
        imputedgenotypes_biallelic_other.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|heterozygous| deletion + wild-type  (delta/wt)'):
        imputedgenotypes_heterodeletion.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|heterozygous| insertion + wild-type (++/wt)'):
        imputedgenotypes_heteroinsertion.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|heterozygous| indel + wild-type (indel/wt)'):
        imputedgenotypes_heteroindel.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|heterozygous| substitution + wild-type (sub/wt)'):
        imputedgenotypes_heterosubstitution.append(i)
    elif re.search('multizygous', imputedgenotypes_dict.get(i)[0]):
        imputedgenotypes_multizygous.append(i)
    elif imputedgenotypes_dict.get(i)[0] in ('|unclear or multi-allelic| insufficient representation of any allele (i.e., no allele exceeds >10% of total reads when adjusted for 10% read threshold)'):
        imputedgenotypes_unclear.append(i)
    elif re.search('', imputedgenotypes_dict.get(i)[0]):
        imputedgenotypes_unclear.append(i)

imputedgenotypes_homowildtype.sort()
imputedgenotypes_homodeletion.sort()
imputedgenotypes_homoinsertion.sort()
imputedgenotypes_homoindel.sort()
imputedgenotypes_homosubstitution.sort()
imputedgenotypes_biallelic_deletion.sort()
imputedgenotypes_biallelic_insertion.sort()
imputedgenotypes_biallelic_indel.sort()
imputedgenotypes_biallelic_substitution.sort()
imputedgenotypes_biallelic_other.sort()
imputedgenotypes_heterodeletion.sort()
imputedgenotypes_heteroinsertion.sort()
imputedgenotypes_heteroindel.sort()
imputedgenotypes_heterosubstitution.sort()
imputedgenotypes_multizygous.sort()
imputedgenotypes_unclear.sort()
           
# Print to genotypes.txt
imputed_genotypes_output = Path(str(output_path)+'/'+processdate+'_genotypes.txt')
# samples are returned based on the following priority:
    # imputedgenotypes_homodeletion
    # imputedgenotypes_homoinsertion
    # iputedgenotypes_homoindel
    # imputedgenotypes_homosubstitution
    # imputedgenotypes_biallelic_deletion
    # imputedgenotypes_biallelic_insertion
    # imputedgenotypes_biallelic_indel
    # imputedgenotypes_biallelic_substitution
    # imputedgenotypes_biallelic_other
    # imputedgenotypes_heterodeletion
    # imputedgenotypes_heteroinsertion
    # imputedgenotypes_heteroindel
    # imputedgenotypes_heterosubstitution
    # imputedgenotypes_multizygous
    # imputedgenotypes_homowildtype
    # imputedgenotypes_unclear
# Call upon allele_output function to report alleles and inferred genotypes based on inferred genotype class
with open(str(imputed_genotypes_output), 'a+') as file:
    file.write('Genotypes.py: Genotypes\nDate: ' + (datetime.today().strftime("%m/%d/%Y")) + '\n\n')
    if len(imputedgenotypes_homodeletion) > 0:
        file.write('\n\nHOMOZYGOUS DELETION\n...................\n\n')
        allele_output(imputedgenotypes_homodeletion)
    if len(imputedgenotypes_homoinsertion) > 0:
        file.write('\n\nHOMOZYGOUS INSERTION\n....................\n\n')
        allele_output(imputedgenotypes_homoinsertion)
    if len(imputedgenotypes_homoindel) > 0:
        file.write('\n\nHOMOZYGOUS INDEL\n................\n\n')
        allele_output(imputedgenotypes_homoindel)
    if len(imputedgenotypes_homosubstitution) > 0:
        file.write('\n\nHOMOZYGOUS SUBSTITUTION\n.......................\n\n')
        allele_output(imputedgenotypes_homosubstitution)
    if len(imputedgenotypes_biallelic_deletion) > 0:
        file.write('\n\nBIALLELIC DELETION\n..................\n\n')
        allele_output(imputedgenotypes_biallelic_deletion)
    if len(imputedgenotypes_biallelic_insertion) > 0:
        file.write('\n\nBIALLELIC INSERTION\n...................\n\n')
        allele_output(imputedgenotypes_biallelic_insertion)
    if len(imputedgenotypes_biallelic_indel) > 0:
        file.write('\n\nBIALLELIC INDEL\n...............\n\n')
        allele_output(imputedgenotypes_biallelic_indel)
    if len(imputedgenotypes_biallelic_substitution) > 0:
        file.write('\n\nBIALLELIC SUBSTITUTION\n......................\n\n')
        allele_output(imputedgenotypes_biallelic_substitution)
    if len(imputedgenotypes_biallelic_other) > 0:
        file.write('\n\nBIALLELIC MUTANT (VARIOUS)\n..........................\n\n')
        allele_output(imputedgenotypes_biallelic_other)
    if len(imputedgenotypes_heterodeletion) > 0:
        file.write('\n\nHETEROZYGOUS DELETION\n.....................\n\n')
        allele_output(imputedgenotypes_heterodeletion)
    if len(imputedgenotypes_heteroinsertion) > 0:
        file.write('\n\nHETEROZYGOUS INSERTION\n......................\n\n')
        allele_output(imputedgenotypes_heteroinsertion)
    if len(imputedgenotypes_heteroindel) > 0:
        file.write('\n\nHETEROZYGOUS INDEL\n..................\n\n')
        allele_output(imputedgenotypes_heteroindel)
    if len(imputedgenotypes_heterosubstitution) > 0:
        file.write('\n\nHETEROZYGOUS SUBSTITUTION\n.........................\n\n')
        allele_output(imputedgenotypes_heterosubstitution)
    if len(imputedgenotypes_multizygous) > 0:
        file.write('\n\nMULTIZYGOUS (>2 ALLELES)\n........................\n\n')
        allele_output(imputedgenotypes_multizygous)
    if len(imputedgenotypes_homowildtype) > 0:
        file.write('\n\nHOMOZYGOUS WILD-TYPE\n....................\n\n')
        allele_output(imputedgenotypes_homowildtype)
    if len(imputedgenotypes_unclear) > 0:
        file.write('\n\nGENOTYPE UNCLEAR (e.g., UNUSUAL ALLELE FREQUENCIES)\n..................................................\n\n')
        allele_output(imputedgenotypes_unclear)

# Log allele definition & genotype inference time duration
imputationDuration = str(datetime.now()- startTime_imputation).split(':')[0]+' hr|'+str(datetime.now() - startTime_imputation).split(':')[1]+' min|'+str(datetime.now() - startTime_imputation).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_imputation).split(':')[2].split('.')[1]+' microsec'

# Start the clock on accessory file processing duration (allele_definitions.csv, population_summary.txt)
startTime_fileprocessing = datetime.now()

# Import data into pandas dataframe
imputedgenotypes_dataframe = pd.DataFrame(
    {
        "allele": [imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[0]+'_'+str(x) for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "read": [imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[1] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "sample": [imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[0] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "reads": [imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[2].split('/')[0].strip('[') for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "totalreads": [imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[2].split('/')[1].strip(']') for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "%totalreads": [imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[4].split(':')[1] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "%top10reads": [imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[6].split(':')[1] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "%readsfilteredfor>1%": [imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[7].split(':')[1] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "%readsfilteredfor>10%": [imputedgenotypes_dict.get(i)[x][0].get('allele_name').split(' ')[8].split(':')[1] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],  
        "chr": [imputedgenotypes_dict.get(i)[x][0].get('chr+build').split(',')[0] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "locusID": [imputedgenotypes_dict.get(i)[x][0].get('locusID') for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "coordinates": [imputedgenotypes_dict.get(i)[x][0].get('coordinates') for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "alignment_query": [imputedgenotypes_dict.get(i)[x][0].get('alignment').split('\n')[1] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "alignment_midline": [imputedgenotypes_dict.get(i)[x][0].get('alignment').split('\n')[2] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "alignment_hit": [imputedgenotypes_dict.get(i)[x][0].get('alignment').split('\n')[3] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "allele_type": [imputedgenotypes_dict.get(i)[x][1].get('allele_type') for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "allele_specs": [imputedgenotypes_dict.get(i)[x][1].get('allele_specs') for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))],
        "inferred_genotype": [imputedgenotypes_dict.get(i)[0] for i in imputedgenotypes_dict for x in range(1,len(imputedgenotypes_dict.get(i)))]
     }
)

# Prepare output file containing allele data in comma-separated, tabular format (allele_definitions.csv)
allele_definitions_csv_output = Path(str(output_path)+'/'+processdate+'_allele_definitions.csv')

imputedgenotypes_dataframe.to_csv(path_or_buf=allele_definitions_csv_output, sep=',')

# Prepare population summary and print to population_summary.txt
population_summary_output = Path(str(output_path)+'/'+processdate+'_population_summary.txt')

# Create list containing contents of pandas dataframe, summarizing sample-specific allele definitions and inferred genotype properties
imputedgenotypes_dataframe['sample'].unique().tolist()

# Population metrics: total sample #
total_samples = []
total_sample_count = 0
for sourcefile in myFastqFilenames:
    fastaname = re.split('_', os.path.basename(sourcefile))
    if fastaname[0] not in total_samples:
        total_samples.append(fastaname[0])
        total_sample_count = total_sample_count+1
        
# Population metrics
sample_checklist = []
genotype_checklist = []
for i in imputedgenotypes_dataframe['sample'].tolist():
    if i not in sample_checklist:
        sample_checklist.append(i)
        genotype_checklist.append(set(imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['inferred_genotype']))
    else:
        pass

# Population metrics: Genotype counts
diploid = 0
for i in genotype_checklist:
    if re.search('homozygous', str(i)) or re.search('heterozygous', str(i)):
        diploid = diploid+1
        
multiploid = 0
for i in genotype_checklist:
    if re.search('multizygous', str(i)):
        multiploid = multiploid+1

unclear = 0
for i in genotype_checklist:
    if re.search('unclear', str(i)):
        unclear = unclear+1
        
homo_wt = 0
for i in genotype_checklist:
    if re.search('homozygous', str(i)) and re.search('wild-type', str(i)):
        homo_wt = homo_wt+1
        
homo_mutant = 0
for i in genotype_checklist:
    if re.search('homozygous', str(i)) and not re.search('wild-type', str(i)):
        homo_mutant = homo_mutant+1
        
homo_deletion = 0
for i in genotype_checklist:
    if re.search('homozygous', str(i)) and re.search('deletion', str(i)):
        homo_deletion = homo_deletion+1
        
homo_insertion = 0
for i in genotype_checklist:
    if re.search('homozygous', str(i)) and re.search('insertion', str(i)):
        homo_insertion = homo_insertion+1
        
homo_substitution = 0
for i in genotype_checklist:
    if re.search('homozygous', str(i)) and re.search('substitution', str(i)):
        homo_substitution = homo_substitution+1
        
homo_indel = 0
for i in genotype_checklist:
    if re.search('homozygous', str(i)) and re.search('indel', str(i)):
        homo_indel = homo_indel+1
        
hetero_wt = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('wild-type', str(i)):
        hetero_wt  = hetero_wt+1
        
hetero_wt_deletion = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('wild-type', str(i)) and re.search('deletion', str(i)):
        hetero_wt_deletion  = hetero_wt_deletion+1
        
hetero_wt_insertion = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('wild-type', str(i)) and re.search('insertion', str(i)):
        hetero_wt_insertion = hetero_wt_insertion+1
        
hetero_wt_substitution = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('wild-type', str(i)) and re.search('substitution', str(i)):
        hetero_wt_substitution = hetero_wt_substitution+1
        
hetero_wt_indel = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('wild-type', str(i)) and re.search('indel', str(i)):
        hetero_wt_indel = hetero_wt_indel+1
        
hetero_mutant_mutant = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and not re.search('wild-type', str(i)):
        hetero_mutant_mutant = hetero_mutant_mutant+1
        
hetero_deletion_insertion = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('deletion', str(i)) and re.search('insertion', str(i)):
        hetero_deletion_insertion = hetero_deletion_insertion+1
        
hetero_deletion_substitution = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('deletion', str(i)) and re.search('substitution', str(i)):
        hetero_deletion_substitution = hetero_deletion_substitution+1

hetero_insertion_substitution = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('insertion', str(i)) and re.search('substitution', str(i)):
        hetero_insertion_substitution = hetero_insertion_substitution+1
        
hetero_deletion_indel = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('deletion', str(i)) and re.search('indel', str(i)):
        hetero_deletion_indel = hetero_deletion_indel+1
        
hetero_insertion_indel = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('insertion', str(i)) and re.search('indel', str(i)):
        hetero_insertion_indel = hetero_insertion_indel+1
        
hetero_substitution_indel = 0
for i in genotype_checklist:
    if re.search('heterozygous', str(i)) and re.search('substitution', str(i)) and re.search('indel', str(i)):
        hetero_substitution_indel = hetero_substitution_indel+1
        
# Population metrics: Allele counts
sample_checklist = []
allele_type_checklist = []
for i in imputedgenotypes_dataframe['sample'].tolist():
    sample_alleles = []
    if i not in sample_checklist:
        sample_checklist.append(i)
        R1R2_check = []
        for x in range(0, len(imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['%readsfilteredfor>10%'])):
            if imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['%readsfilteredfor>10%'].iloc[x] == 'None':
                pass
            elif float(imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['%readsfilteredfor>10%'].iloc[x]) > 10:
                if imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['read'].iloc[x] == 'R1+R2':
                    R1R2_check.append('R1+R2/'+imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['allele_type'].iloc[x])
                    sample_alleles.append(imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['allele_type'].iloc[x])
        if len(sample_alleles) != 0:
            allele_type_checklist.append(sample_alleles)
            
# Population metrics: Compile counts of total defined alleles, wild-type alleles, and mutant alleles 
wt_alleles = 0
mutant_alleles = 0

for i in allele_type_checklist:
    if len(i) == 1:
        if re.search('wild-type', str(i)):
            wt_alleles = wt_alleles+2
        if re.search('mutant', str(i)):
            mutant_alleles = mutant_alleles+2
    elif len(i) == 2:
        if re.search('wild-type', str(i)):
            wt_alleles = wt_alleles+1
        if re.findall('mutant', str(i)):
            mutant_alleles = mutant_alleles+len(re.findall('mutant', str(i)))
    elif len(i) == 3:
        if re.search('wild-type', str(i)):
            wt_alleles = wt_alleles+1
        if re.findall('mutant', str(i)):
            mutant_alleles = mutant_alleles+len(re.findall('mutant', str(i)))

total_alleles = wt_alleles+mutant_alleles

# Population metrics
sample_checklist = []
allele_specs_checklist = []
for i in imputedgenotypes_dataframe['sample'].tolist():
    sample_alleles = []
    if i not in sample_checklist:
        sample_checklist.append(i)
        R1R2_check = []
        for x in range(0, len(imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['%readsfilteredfor>10%'])):
            if imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['%readsfilteredfor>10%'].iloc[x] == 'None':
                pass
            elif imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['allele_specs'].iloc[x] is None:
                if 'wild-type' not in sample_alleles:
                    sample_alleles.append('wild-type')
            elif float(imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['%readsfilteredfor>10%'].iloc[x]) > 10:
                if imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['read'].iloc[x] == 'R1+R2':
                    if 'R1+R2/'+imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['allele_specs'].iloc[x] not in R1R2_check:
                        R1R2_check.append('R1+R2/'+imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['allele_specs'].iloc[x])
                        sample_alleles.append(imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['allele_specs'].iloc[x])
                elif imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['read'].iloc[x] == 'R1+R2':
                    if 'R1+R2'+imputedgenotypes_dataframe[imputedgenotypes_dataframe['sample'] == i]['allele_specs'].iloc[x] in R1R2_check:
                        pass
        if len(sample_alleles) != 0:
            allele_specs_checklist.append(sample_alleles)
            
# Population metrics: Compile counts of total deletion, insertion, substitution and indel alleles
deletion_alleles = 0
insertion_alleles = 0
substitution_alleles = 0
indel_alleles = 0

for i in allele_specs_checklist:
    if len(i) == 1:
        if re.search('deletion', str(i)):
            deletion_alleles = deletion_alleles+2
        elif re.search('insertion', str(i)):
            insertion_alleles = insertion_alleles+2
        elif re.search('substitution', str(i)):
            substitution_alleles = substitution_alleles+2
        elif re.search('indel', str(i)):
            indel_alleles = indel_alleles+2
    elif len(i) == 2:
        if re.findall('deletion', str(i)):
            deletion_alleles = deletion_alleles+len(re.findall('deletion', str(i)))
        if re.findall('insertion', str(i)):
            insertion_alleles = insertion_alleles+len(re.findall('insertion', str(i)))
        if re.findall('substitution', str(i)):
            substitution_alleles = substitution_alleles+len(re.findall('substitution', str(i)))
        if re.findall('indel', str(i)):
            indel_alleles = indel_alleles+len(re.findall('indel', str(i)))
    elif len(i) == 3:
        if re.findall('deletion', str(i)):
            deletion_alleles = deletion_alleles+len(re.findall('deletion', str(i)))
        if re.findall('insertion', str(i)):
            insertion_alleles = insertion_alleles+len(re.findall('insertion', str(i)))
        if re.findall('substitution', str(i)):
            substitution_alleles = substitution_alleles+len(re.findall('substitution', str(i)))
        if re.findall('indel', str(i)):
            indel_alleles = indel_alleles+len(re.findall('indel', str(i)))
            
# Population metrics: Compile list of sampleIDs for which there were no alignment hits for any of the top 10 reads
no_hits_and_hits_samplename_list = []
no_hits_for_any_top10_reads_samplename_list = []
for i in no_hits_samplename_list:
    if i in querydef_uniq_list:
        no_hits_and_hits_samplename_list.append(i)
    else:
        no_hits_for_any_top10_reads_samplename_list.append(i)
        
# Prepare population_summary.txt file
with open(str(population_summary_output), 'a+') as file:
    file.write('Genotypes.py: Population Summary\nDate: ' + (datetime.today().strftime("%m/%d/%Y")) +
"""\n\nI. Synopsis of Interpretations: Allele Definitions & Genotype Extrapolations

    (A) Sample summary
        (i) Number of samples processed: """ + str(total_sample_count) +
'\n        (ii) % samples called (genotype inferred): ' + str(len(sample_checklist)) + ' (' + str(round((100*(len(sample_checklist)/total_sample_count)),2))+'%)' +
"""\n\n    (B) Genotypes summary
        (i) % samples diploid (1-2 prominent alleles inferred): """ + str(diploid) + ' (' + str(round((100*(diploid/total_sample_count)),2))+'%)' +
'\n            (1) % homozygous wild-type (wt): ' + str(homo_wt) + ' (' + str(round((100*(homo_wt/total_sample_count)),2))+'%)' +
'\n            (2) % homozygous mutant: ' + str(homo_mutant) + ' (' + str(round((100*(homo_mutant/total_sample_count)),2))+'%)' +
'\n                -> % homozygous deletion: ' + str(homo_deletion) + ' (' + str(round((100*(homo_deletion/total_sample_count)),2))+'%)' +
'\n                -> % homozygous insertion: ' + str(homo_insertion) + ' (' + str(round((100*(homo_insertion/total_sample_count)),2))+'%)' +
'\n                -> % homozygous substitution: ' + str(homo_substitution) + ' (' + str(round((100*(homo_substitution/total_sample_count)),2))+'%)' +
'\n                -> % homozygous complex indel: ' + str(homo_indel) + ' (' + str(round((100*(homo_indel/total_sample_count)),2))+'%)' +
'\n            (3) % heterozygous (wt + mutant): ' + str(hetero_wt) + ' (' + str(round((100*(hetero_wt/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous deletion: ' + str(hetero_wt_deletion) + ' (' + str(round((100*(hetero_wt_deletion/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous insertion: ' + str(hetero_wt_insertion) + ' (' + str(round((100*(hetero_wt_insertion/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous substitution: ' + str(hetero_wt_substitution) + ' (' + str(round((100*(hetero_wt_substitution/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous complex indel: ' + str(hetero_wt_indel) + ' (' + str(round((100*(hetero_wt_indel/total_sample_count)),2))+'%)' +
'\n            (4) % heterozygous (mutant + mutant): ' + str(hetero_mutant_mutant) + ' (' + str(round((100*(hetero_mutant_mutant/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous deletion + insertion: ' + str(hetero_deletion_insertion) + ' (' + str(round((100*(hetero_deletion_insertion/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous deletion + substitution: ' + str(hetero_deletion_substitution) + ' (' + str(round((100*(hetero_deletion_substitution/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous insertion + substitution: ' + str(hetero_insertion_substitution) + ' (' + str(round((100*(hetero_insertion_substitution/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous deletion + complex indel: ' + str(hetero_deletion_indel) + ' (' + str(round((100*(hetero_deletion_indel/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous insertion + complex indel: ' + str(hetero_insertion_indel) + ' (' + str(round((100*(hetero_insertion_indel/total_sample_count)),2))+'%)' +
'\n                -> % heterozygous substitution + complex indel: ' + str(hetero_substitution_indel) + ' (' + str(round((100*(hetero_substitution_indel/total_sample_count)),2))+'%)' +
'\n        (ii) % samples multiploid (>2 prominent alleles inferred): ' + str(multiploid) + ' (' + str(round((100*(multiploid/total_sample_count)),2))+'%)' +
"""\n\n    (B) Alleles summary
        (i) % wild-type alleles: """ + str(wt_alleles) + ' (' + str(round((100*(wt_alleles/total_alleles)),2)) + '% of total alleles)' +
'\n        (ii) % mutant alleles: ' + str(mutant_alleles) + ' (' + str(round((100*(mutant_alleles/total_alleles)),2)) + '% of total alleles)' +
'\n            (1) % deletion alleles: ' + str(deletion_alleles) + ' (' + str(round((100*(deletion_alleles/total_alleles)),2)) + '% of total alleles)' +
'\n            (2) % insertion alleles: ' + str(insertion_alleles) + ' (' + str(round((100*(insertion_alleles/total_alleles)),2)) + '% of total alleles)' +
'\n            (3) % substitution alleles: ' + str(substitution_alleles) + ' (' + str(round((100*(substitution_alleles/total_alleles)),2)) + '% of total alleles)' +
'\n            (4) % complex indel alleles: ' + str(indel_alleles) + ' (' + str(round((100*(indel_alleles/total_alleles)),2)) + '% of total alleles)')
    file.write("""\n\nII. Synopsis of Reads Lost to Analysis
    Reads among the 'Top 10' reads are deprecated (not analyzed by Genotypes.py) if they fall into the following categories:
    (A) no hits, (B) multiple hits, or (C) (if >1 hsp for 1 hit) overlapping hsp's in reference database and/or hsp span (end-to-end) that exceeds 1 kb
    
    (A) Samples with reads among the 'top 10 most abundant reads', that did not map to the reference genome
        (i) For the following sample IDs ("""+str(len(no_hits_for_any_top10_reads_samplename_list))+"""), NO reads among the "top 10 most abundant reads" could be mapped to the reference genome:""")
    if len(no_hits_for_any_top10_reads_samplename_list) == 0:
        file.write('\n               None')
    else:
        for i in no_hits_for_any_top10_reads_samplename_list:
            file.write('\n             '+i)
    file.write('')
    
    file.write('\n        (ii) For the following sample IDs ('+str(len(no_hits_samplename_list))+'), the indicated reads among the "top 10 most abundant reads" did not map to the reference genome:')
    # Print this output to population summary
    if len(no_hits_samplename_list) == 0:
        file.write('\n               None')
    else:
        for i in no_hits_samplename_list:
            file.write('\n             '+i+':')
            for x in no_hits_list:
                if i == x.split('_')[0]:
                    file.write('\n               '+x)
            file.write('')
        
    file.write('\n\n    (B) Samples with reads among the "top 10 most abundant reads", that mapped to multiple loci in the reference genome' +
'\n        (i) For the following sample IDs ('+str(len(multiple_alignments_hits_samplename_list))+'), the indicated reads among the "top 10 most abundant reads" mapped to more than one locus in the reference genome:')
    if len(multiple_alignments_hits_samplename_list) > 0:
        for i in multiple_alignments_hits_samplename_list:
            file.write('\n               '+i.strip()) 
        file.write('\n\n               Details:')
        for i in multiple_alignments_hits_dict:
            file.write('\n               '+i+'\n               '+len(i)*'=')
            for x in multiple_alignments_hits_dict.get(i):
                chunked_list = []
                subchunk = 0
                for index, y in enumerate(x):   
                    if y.split('>')[0] == '<Hit_num':
                        chunked_list.append(x[subchunk:index])
                        subchunk = index
                    if index == len(x)-1:
                        chunked_list.append(x[subchunk:])
                file.write('\n               Read: '+chunked_list[0][1].split('>')[1].split('<')[0])
                for hit in chunked_list[1:]:
                    file.write('\n               Hit '+hit[0].split('>')[1].split('<')[0]+': '+hit[2].split('>')[1].split('<')[0]+': '+hit[1].split('>')[1].split('<')[0]+', '+hit[5].split('>')[1].split('<')[0]+'-'+hit[6].split('>')[1].split('<')[0])
                    file.write('\n                   '+hit[7].split('>')[1].split('<')[0])
                    file.write('\n                   '+hit[9].split('>')[1].split('<')[0])
                    file.write('\n                   '+hit[8].split('>')[1].split('<')[0])
                file.write('\n')
            file.write('\n')
    else:
        file.write('\n               None')
    file.write("""\n\n    (C) Samples with reads among the 'top 10 most abundant reads', with high-scoring pairs (hsp's) that were attempted 
        for reconstruction with BLASTDBCMD sequence retrieval across hsp span, but deprecated because the hsp's
        a) overlapped in the reference genome (potential artefact or signature of inversion, duplication and/or
        b) exceeded 1 kb span (end-to-end)
        Consult fasta.fa and/or BLASTN output if any of these reads is of further interest to examine manually:\n""")
    file.write("\n        (i) For the following sample IDs ("+str(len(multiple_alignments_hsp_samplename_list))+"), the indicated reads among the 'top 10 most abundant reads' exhibited overlapping hsp's aligned to the reference genome:")
    if len(multiple_alignments_hsp_invalid2_list) > 0:
        for i in multiple_alignments_hsp_invalid2_list:
            file.write('\n               '+i.strip())
        file.write('\n\n               Details:')   
        for i in multiple_alignments_hsp_dict_invalid2:
            file.write('\n               '+i+'\n               '+len(i)*'=')
            for x in multiple_alignments_hsp_dict_invalid2.get(i):
                chunked_list = []
                subchunk = 0
                for index, y in enumerate(x):
                    if y.split('>')[0] == '<Hsp_num':
                        chunked_list.append(x[subchunk:index])
                        subchunk = index
                    if index == len(x)-1:
                        chunked_list.append(x[subchunk:])
                file.write('\n               Read: '+chunked_list[0][1].split('>')[1].split('<')[0])
                file.write('\n               Hit: '+chunked_list[0][4].split('>')[1].split('<')[0]+': '+chunked_list[0][3].split('>')[1].split('<')[0])
                for hsp in chunked_list[1:]:
                    file.write('\n               Hsp '+hsp[0].split('>')[1].split('<')[0]+': '+hsp[1].split('>')[1].split('<')[0]+'-'+hsp[2].split('>')[1].split('<')[0])
                    file.write('\n                   '+hsp[3].split('>')[1].split('<')[0])
                    file.write('\n                   '+hsp[5].split('>')[1].split('<')[0])
                    file.write('\n                   '+hsp[4].split('>')[1].split('<')[0])
                file.write('\n')
            file.write('\n')
    else:
        file.write('\n               None')
                
# Log file processing time duration                
fileprocessingDuration = str(datetime.now()- startTime_fileprocessing).split(':')[0]+' hr|'+str(datetime.now() - startTime_fileprocessing).split(':')[1]+' min|'+str(datetime.now() - startTime_fileprocessing).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime_fileprocessing).split(':')[2].split('.')[1]+' microsec'

# Optional: print supporting evidence (frequency metrics demonstrated in plots) for allele definitions & inferred genotypes to output file, allele_evidence.pdf
if frequency_plot_check == 'Y':
    print("""
Script is now compiling evidence for extrapolated genotypes in the form of allele frequency plots.""")
    frequency_plots()
elif frequency_plot_check == 'N':
    pass

# Log script processing time duration 
processingDuration = str(datetime.now()- startTime).split(':')[0]+' hr|'+str(datetime.now() - startTime).split(':')[1]+' min|'+str(datetime.now() - startTime).split(':')[2].split('.')[0]+' sec|'+str(datetime.now() - startTime).split(':')[2].split('.')[1]+' microsec'

# Log script end time
endTime = datetime.now()
endTimestr = str(endTime).split(' ')[1].split('.')[0]     

# Assess output file set created by script
file_set = [file for file in os.listdir(output_directory) if Path(file).suffix in ('.pdf','.txt','.fa')]

# Log further script operation metrics to script_metrics.txt
filename = Path(str(output_path)+ '/'+processdate+'_script_metrics.txt')

if frequency_plot_check == 'Y':
    with open(filename, 'a') as f:
        print("""\nFile output information:
    Output directory: """ + str(output_directory) +
'\n    Total file #: ' + str(len(file_set)) +
'\n    Total file output sizes: ', file = f)
        for file in file_set:
            print('        '+file+': '+path_size(str(output_directory)+'/'+file), file = f)
        print("""\nScript operation times:
    start time: """+startTimestr+
    '\n    fasta processing time: '+readcountDuration+
    '\n    alignments processing time: '+alignmentsDuration+
    '\n    genotype inference processing time: '+imputationDuration+
    '\n    frequency plots compilation time: '+frequencyplotsDuration+
    '\n    accessory file processing time: '+fileprocessingDuration+
    '\n    total processing time: '+processingDuration+
    '\n    end time: ' + endTimestr, file = f)
    f.close()
elif frequency_plot_check == 'N':
    with open(filename, 'a') as f:
        print("""\nFile output information:
    Output directory: """ + str(output_directory) +
'\n    Total file #: ' + str(len(file_set)) +
'\n    Total file output sizes: ', file = f)
        for file in file_set:
            print('        '+file+': '+path_size(str(output_directory)+'/'+file), file = f)
        print("""\nScript operation times:
    start time: """+startTimestr+
    '\n    fasta processing time: '+readcountDuration+
    '\n    alignments processing time: '+alignmentsDuration+
    '\n    genotype inference processing time: '+imputationDuration+
    '\n    accessory file processing time: '+fileprocessingDuration+
    '\n    total processing time: '+processingDuration+
    '\n    end time: ' + endTimestr, file = f)
    f.close()
    

# End of script operations
print("""
Script has completed.  Please find output files at """+str(output_directory))

############################################################################# end