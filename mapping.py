#! /bin/python3
# Copyright (C) 2020 Mario Ernst
# This is a script that maps previously cleaned sequence reads in fastq format against a fasta reference with the assembled contigs.
# The workflow is based on the cleaning script by Mozes Blom (https://github.com/MozesBlom/NGSdata_tools/blob/master/clean_up_raw_reads.py)
# It assumes that the input files for each library are both paired and single end fastq files and performs the following steps:
#       1) indexes FASTA file containing the assembled contigs using the bwa software package
#       2) aligns the reads to the reference ref_genome using the BWA-MEM algorithm
#       3) converts the SAM files into its binary version (BAM)
#       4) removes reads that did not pair properly with the reference using SAMtools
#       5) sorts and 6) indexes the filtered BAM files using SAMtools
#       7) merges the BAM file with the single end reads and the BAM file with the paired end reads for each library using SAMtools
#       8) indexes the merged BAM files for each library using SAMtools


#########################
## Required modules
#########################
try:
        import os
        import csv
        import subprocess
        import sys
        import datetime
        import gzip
except ImportError:
        sys.exit("Cant import one of the required modules...")

#########################
# Which Analysis to run? 1 = yes, 0 = no
#########################

index_ref = 0            #[int, 0 or 1] index the reference genome
align = 0                #[int, 0 or 1] align the single end and the paired end read files to the reference genome
SAM2BAM = 0              #[int, 0 or 1] convert the SAM files to BAM files
filterBAM = 0            #[int, 0 or 1] remove the reads that have not mapped to the reference
sort_BAM = 0             #[int, 0 or 1] sort the BAM files
index_BAM = 0            #[int, 0 or 1] index the BAM files
merge_PE_and_SE = 0       #[int, 0 or 1] merge the paired and the unpaired BAM files
index_merged_BAM = 0     #[int, 0 or 1] index the merged BAM files

# ----------------------- #

#########################
## Input
#########################

# Run specific folder name, where the results will be stored [char]
results_path = '/path/to/results/folder/'
# Folder where all the sequence reads are stored [char]
reads_path = '/path/to/folder/with/reads/'
# Folder where all the fasta reference file is stored [char]
reference_path = '/path/to/reference/'
# Path to list with library file names. The list should be a column containing the names of each library in the following format: indivA_Lxxx
lib_list = '/Path/to/lib_list.txt'

#########################
## Path to dependables
#########################

pigz = '/Path/to/bin/pigz-2.3.4/pigz'               # [char] Path to compiled PIGZ executable
bwa = '/Path/to/bin/bwa-0.7.17/bwa'                 # [char] Path to compiled BWA executable
samtools = '/Path/to/bin/samtools-0.1.19/samtools'  # [char] Path to compiled samtools executable

#########################
## Parameter settings
#########################

mem = 'mem'
VIEW = 'view'
SORT = 'sort'
index = 'index'
merge = 'merge'

#########################
## Functions
#########################

# Index reference genome
def indexref(bwa, index, fasta_file):
        subprocess.call("%s %s %s" % (bwa, index, fasta_file), shell=True)

# Run BWA mem for Unpaired
def bwamemU(bwa, mem, input_reference, Uin_file, Uout_file):
        subprocess.call("%s %s %s %s > %s" % (bwa, mem, input_reference, Uin_file, Uout_file), shell=True)

# Run BWA mem for Paired
def bwamemPE(bwa, mem, reference, in_file_path_r1, in_file_path_r2, out_dir):
        subprocess.call("%s %s %s %s %s > %s" % (bwa, mem, reference, in_file_path_r1, in_file_path_r2, out_dir), shell=True)

# Convert from SAM to BAM
def SAMtoBAM(samtools, VIEW, SAM_file, BAM_file):
        subprocess.call("%s %s -bS %s > %s" % (samtools, VIEW, SAM_file, BAM_file), shell=True)

# Filter the unmapped reads
def filterunmapped(samtools, VIEW, BAM_in, BAM_out):
        subprocess.call("%s %s -F 4 -b %s > %s" % (samtools, VIEW, BAM_in, BAM_out), shell=True)

# Sort BAM
def sortBAM(samtools, SORT, BAM_in_sort, BAM_out_sort):
        subprocess.call("%s %s %s %s" % (samtools, SORT, BAM_in_sort, BAM_out_sort), shell=True)

# index the BAM file
def indexBAM(samtools, index, in_file):
        subprocess.call("%s %s %s" % (samtools, index, in_file), shell=True)

# merge the paired end and the single end file
def mergefiles(samtools, merge, out_file, in_filePE, in_fileSE):
        subprocess.call("%s %s %s %s %s" % (samtools, merge, out_file, in_filePE, in_fileSE), shell=True)

#########################
## AUTOMATED ANALYSIS
#########################

# Print the date + time, when analysis started
print('Mapping started at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
sys.stdout.flush()

# Create the output folder in case it doesn't exist yet
if os.path.isdir(results_path):
        pass
else:
        subprocess.call("mkdir '%s'" % (results_path), shell=True)

# Create a list of all the libraries to be run
# Expects libs to be in format: indivA_Lxxx. In this case the read file name would be: indivA_Lxxx.fastq.gz

libs = []
try:
        with open(lib_list, 'r') as fn:
                for line in fn:
                        lib_name = line.rsplit('\n')[0]
                        libs.append(lib_name)
except IOError:
        sys.exit("something is wrong with the path to the file with the library names")
except:
        sys.exit("Unexpected error, is the lib_list file formatted correctly?")

#index the reference genome
if (index_ref == 1):
        indexref(bwa, index, reference_path)
else:
        pass

#align the read files to the reference genome
outDir_1 = os.path.join(results_path, "1.aligned")
if (align == 1):
        if os.path.isdir(outDir_1):
                pass
        else:
                subprocess.call("mkdir '%s'" % (outDir_1), shell=True)
        indiv_done = []
        for lib in libs:
                indiv = lib.rsplit('_')[0]
                outDir_1_indiv = os.path.join(outDir_1, indiv)
                if indiv not in indiv_done:
                        subprocess.call("mkdir '%s'" % (outDir_1_indiv), shell=True)
                        indiv_done.append(indiv)
                else:
                        pass
                outDir_1_indiv_lib = os.path.join(outDir_1_indiv, lib)
                if os.path.isdir(outDir_1_indiv_lib):
                        sys.exit("path %s already exists!!" % (outDir_1_indiv_lib))
                else:
                        subprocess.call("mkdir '%s'" % (outDir_1_indiv_lib), shell=True)
                os.chdir(outDir_1_indiv_lib)
                #input
                lib_U_path_gz_in = os.path.join(reads_path, indiv, lib, (lib + '_U.fastq.gz'))
                lib_PE_R1_path_gz = os.path.join(reads_path, indiv, lib, (lib + '_R1.fastq.gz'))
                lib_PE_R2_path_gz = os.path.join(reads_path, indiv, lib, (lib + '_R2.fastq.gz'))
                #output
                lib_U_path_out = os.path.join(outDir_1_indiv_lib, (lib + '_unpaired.ali.sam'))
                lib_PE_path_out = os.path.join(outDir_1_indiv_lib, (lib+ '_paired.ali.sam'))
                #run bwa mem
                bwamemPE(bwa, mem, reference_path, lib_PE_R1_path_gz, lib_PE_R2_path_gz, lib_PE_path_out)
                bwamemU(bwa, mem, reference_path, lib_U_path_gz_in, lib_U_path_out)
else:
        pass

#convert SAM to BAM
outDir_2 = os.path.join(results_path, '2.SAMtoBAM')
if (SAM2BAM == 1):
        if os.path.isdir(outDir_2):
                pass
        else:
                subprocess.call("mkdir '%s'" % (outDir_2), shell=True)
        indiv_done = []
        for lib in libs:
                indiv = lib.rsplit('_')[0]
                outDir_2_indiv = os.path.join(outDir_2, indiv)
                if indiv not in indiv_done:
                        subprocess.call("mkdir '%s'" % (outDir_2_indiv), shell=True)
                        indiv_done.append(indiv)
                else:
                        pass
                outDir_2_indiv_lib = os.path.join(outDir_2_indiv, lib)
                if os.path.isdir(outDir_2_indiv_lib):
                        sys.exit("path %s already exists!!" % (outDir_2_indiv_lib))
                else:
                        subprocess.call("mkdir '%s'" % (outDir_2_indiv_lib), shell=True)
                os.chdir(outDir_2_indiv_lib)
                #input
                lib_U_path_SAM_in = os.path.join(outDir_1, indiv, lib, (lib + '_unpaired.ali.sam'))
                lib_PE_path_SAM_in = os.path.join(outDir_1, indiv, lib, (lib + '_paired.ali.sam'))
                #output
                lib_U_path_BAM_out = os.path.join(outDir_2_indiv_lib, (lib + '_unpaired.ali.bam'))
                lib_PE_path_BAM_out = os.path.join(outDir_2_indiv_lib, (lib+ '_paired.ali.bam'))
                #run santools
                SAMtoBAM(samtools, VIEW, lib_PE_path_SAM_in, lib_PE_path_BAM_out)
                SAMtoBAM(samtools, VIEW, lib_U_path_SAM_in, lib_U_path_BAM_out)
else:
        pass

#Get only the mapped reads in a BAM file
outDir_3 = os.path.join(results_path, "3.filterBAM")
if (filterBAM == 1):
        if os.path.isdir(outDir_3):
                pass
        else:
                subprocess.call("mkdir '%s'" % (outDir_3), shell=True)
        indiv_done = []
        for lib in libs:
                indiv = lib.rsplit('_')[0]
                outDir_3_indiv = os.path.join(outDir_3, indiv)
                if indiv not in indiv_done:
                        subprocess.call("mkdir '%s'" % (outDir_3_indiv), shell=True)
                        indiv_done.append(indiv)
                else:
                        pass
                outDir_3_indiv_lib = os.path.join(outDir_3_indiv, lib)
                if os.path.isdir(outDir_3_indiv_lib):
                        sys.exit("path %s already exists!!" % (outDir_2_indiv_lib))
                else:
                        subprocess.call("mkdir '%s'" % (outDir_3_indiv_lib), shell=True)
                os.chdir(outDir_3_indiv_lib)
                #input
                lib_U_path_BAM_in = os.path.join(outDir_2, indiv, lib, (lib + '_unpaired.ali.bam'))
                lib_PE_path_BAM_in = os.path.join(outDir_2, indiv, lib, (lib + '_paired.ali.bam'))
                #output
                lib_U_path_BAM_filtered = os.path.join(outDir_3_indiv_lib, (lib + '_unpaired.ali.mapped.bam'))
                lib_PE_path_BAM_filtered = os.path.join(outDir_3_indiv_lib, (lib + '_paired.ali.mapped.bam'))
                #remove the unmapped reads
                filterunmapped(samtools, VIEW, lib_PE_path_BAM_in, lib_PE_path_BAM_filtered)
                filterunmapped(samtools, VIEW, lib_U_path_BAM_in, lib_U_path_BAM_filtered)
else:
        pass

# sort BAM files
outDir_4 = os.path.join(results_path, "4.sorted_BAM")
if (sort_BAM == 1):
        if os.path.isdir(outDir_4):
                pass
        else:
                subprocess.call("mkdir '%s'" % (outDir_4), shell=True)
        indiv_done = []
        for lib in libs:
                indiv = lib.rsplit('_')[0]
                outDir_4_indiv = os.path.join(outDir_4, indiv)
                if indiv not in indiv_done:
                        subprocess.call("mkdir '%s'" % (outDir_4_indiv), shell=True)
                        indiv_done.append(indiv)
                else:
                        pass
                outDir_4_indiv_lib = os.path.join(outDir_4_indiv, lib)
                if os.path.isdir(outDir_4_indiv_lib):
                        sys.exit("path %s already exists!!" % (outDir_4_indiv_lib))
                else:
                        subprocess.call("mkdir '%s'" % (outDir_4_indiv_lib), shell=True)
                os.chdir(outDir_4_indiv_lib)
                #input
                U_BAM_path_in = os.path.join(outDir_3, indiv, lib, (lib + '_unpaired.ali.mapped.bam'))
                PE_BAM_path_in = os.path.join(outDir_3, indiv, lib, (lib + '_paired.ali.mapped.bam'))
                #output
                U_sorted_BAM_path_out = os.path.join(outDir_4_indiv_lib, (lib + '_unpaired_sorted'))
                PE_sorted_BAM_path_out = os.path.join(outDir_4_indiv_lib, (lib + '_paired_sorted'))
                sortBAM(samtools, SORT, U_BAM_path_in, U_sorted_BAM_path_out)
                sortBAM(samtools, SORT, PE_BAM_path_in, PE_sorted_BAM_path_out)
else:
        pass

# index BAM
if (index_BAM == 1):
        indiv_done = []
        for lib in libs:
                indiv = lib.rsplit('_')[0]
                if indiv not in indiv_done:
                        indiv_done.append(indiv)
                else:
                        pass
                #input
                U_BAM_path_in = os.path.join(outDir_4, indiv, lib, (lib + '_unpaired_sorted.bam'))
                PE_BAM_path_in = os.path.join(outDir_4, indiv, lib, (lib + '_paired_sorted.bam'))
                #index
                index = "index"
                indexBAM(samtools, index, U_BAM_path_in)
                indexBAM(samtools, index, PE_BAM_path_in)
else:
        pass

# merge paired and unpaired BAMs
outDir_5 = os.path.join(results_path, "5.merge_paired_unpaired")
if (merge_PE_and_SE == 1):
        if os.path.isdir(outDir_5):
                pass
        else:
                subprocess.call("mkdir '%s'" % (outDir_5), shell=True)
        indiv_done = []
        for lib in libs:
                indiv = lib.rsplit('_')[0]
                outDir_5_indiv = os.path.join(outDir_5, indiv)
                if indiv not in indiv_done:
                        subprocess.call("mkdir '%s'" % (outDir_5_indiv), shell=True)
                        indiv_done.append(indiv)
                else:
                        pass
                outDir_5_indiv_lib = os.path.join(outDir_5_indiv, lib)
                if os.path.isdir(outDir_5_indiv_lib):
                        sys.exit("path %s already exists!!" % (outDir_5_indiv_lib))
                else:
                        subprocess.call("mkdir '%s'" % (outDir_5_indiv_lib), shell=True)
                os.chdir(outDir_5_indiv_lib)
                PE_merge_path_in = os.path.join(outDir_4, indiv, lib, (lib + '_paired_sorted.bam'))
                U_merge_path_in = os.path.join(outDir_4, indiv, lib, (lib + '_unpaired_sorted.bam'))
                merged_path_out = os.path.join(outDir_5_indiv_lib, (lib + '_merged.bam'))
                mergefiles(samtools, merge, merged_path_out, PE_merge_path_in, U_merge_path_in)
else:
        pass

# index merged BAMs
if (index_merged_BAM == 1):
        indiv_done = []
        for lib in libs:
                indiv = lib.rsplit('_')[0]
                if indiv not in indiv_done:
                        indiv_done.append(indiv)
                else:
                        pass
                #input
                BAM_path_in = os.path.join(outDir_5, indiv, lib, (lib + '_merged.bam'))
                #index
                index = "index"
                indexBAM(samtools, index, BAM_path_in)
else:
        pass

# Print the date + time, when analysis finished
print('Mapping finished at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
sys.stdout.flush()
