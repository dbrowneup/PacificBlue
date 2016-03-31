NOTE: I have hosted this software here to easily manage the code across several machines, 
as well as make and track changes as needed. In no way do I claim ownership or authorship.

Cerulean Hybrid Genome Assembler v0.1.1

This software extends contigs assembled using short read datasets like Illumina
paired-end reads using long reads like PacBio RS long reads.

The method is fully described in:
Deshpande, V., Fung, E. D., Pham, S., & Bafna, V. (2013).
Cerulean: A hybrid assembly using high throughput short and long reads.
http://arxiv.org/abs/1307.7933

A] Requirements:
Ubuntu 12.04 (may run on other operating systems, but not tested)
Python 2.7.1 (may run on older versions, but not tested)
numpy, matplotlib libraries for Python
ABySS assembler: http://www.bcgsc.ca/platform/bioinfo/software/abyss
SMRT Analysis tookit (for BLASR): http://pacbiodevnet.com/
PBJelly: https://sourceforge.net/projects/pb-jelly/ 

B] Inputs and Pre-processing:
 i) Assembled contigs from ABySS short read assembler
 ii)Mapping of Pacbio reads to ABySS contigs using BLASR

 i) Assembly of Illumina paired-end reads:
   If the paired-end reads are stored in fastq format in the files reads1.fastq
   and reads2.fastq, then contigs may be assembled by:
   $ abyss-pe k=64 n=10 in='reads1.fastq reads2.fastq' name=<dataname>
   This will generate 2 files used for inputs to Cerulean:
   * <dataname>-contigs.fa    #This contains the contig sequences
   * <dataname>-contigs.dot   #This contains the graph structure

 ii)Mapping PacBio reads to ABySS contigs using BLASR:
   Note: sawriter and blasr are part of SMRT Analysis toolkit
   Note: You need to set the environmental variables and path:
   $ export SEYMOUR_HOME=/opt/smrtanalysis/
   $ source $SEYMOUR_HOME/etc/setup.sh
   
   Suppose PacBio reads are stored in <dataname>_pacbio.fasta
   $ sawriter <dataname>-contigs.fa
   $ blasr <dataname>_pacbio.fa <dataname>-contigs.fa -minMatch 10 \
     -minPctIdentity 70 -bestn 30 -nCandidates 30 -maxScore -500 \
     -nproc <numthreads> -noSplitSubreads \
     -out <dataname>_pacbio_contigs_mapping.fasta.m4
   
   Make sure the fasta.m4 file generated has the following format:
   qname tname qstrand tstrand score pctsimilarity tstart tend tlength \
   qstart qend qlength ncells
   The file format may be verified by adding the option -header to blasr. 

C] Execute:
 Cerulean requires that all input files are in the same directory <basedir>:
 i)   <basedir>/<dataname>-contigs.fa
 ii)  <basedir>/<dataname>-contigs.dot
 iii) <basedir>/<dataname>_pacbio_contigs_mapping.fasta.m4

 To run:
 $ python src/Cerulean.py --dataname <dataname> --basedir <basedir> \
 --nproc <numthreads>
 
 This will generate:
 i)  <basedir>_cerulean.fasta
 ii) <basedir>_cerulean.dot
 Note: The dot does not have same contigs as fasta, but intermediate graph.
 
 
D] Post-processing:
 Currently Cerulean does not include consensus sequence of PacBio reads in gaps
 The gaps may be filled using PBJelly.
 $ python $JELLYPATH/fakeQuals.py <dataname>_cerulean.fasta <dataname>_cerulean.qual
 $ python $JELLYPATH/fakeQuals.py <dataname>_pacbio.fasta <dataname>_pacbio.qual
 $ cp $JELLYPATH/lambdaExample/Protocol.xml .
 $ mkdir PBJelly
 Modify Protocol.xml as follows:
 Set <reference> to $PATH_TO_<basedir>/<dataname>_cerulean.fasta 
 Set <outputDir> to $PATH_TO_<basedir>/PBJelly
 Set <baseDir> to $PATH_TO_<basedir>
 Set <job> to <dataname>_pacbio.fasta
 Set <blasr> option -nproc <numthreads> 
 Note: PBJelly requires that the suffix be .fasta and not .fa
 Next run PBJelly:
 ($ source $JELLYPATH/exportPaths.sh)
 $ python $JELLYPATH/Jelly.py <stage> Protocol.xml
 where <stage> has to be in the order:
 setup
 mapping
 support
 extraction
 assembly
 output
 
 The assembled contigs may be view in <basedir>/PBJelly/assembly/jellyOutput.fasta

In case of any questions or errors please contact vdeshpan eng DT ucsd DT edu
