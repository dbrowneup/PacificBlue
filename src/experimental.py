# This software is Copyright 2013 The Regents of the University of
# California. All Rights Reserved.
#
# Permission to copy, modify, and distribute this software and its
# documentation for educational, research and non-profit purposes, without fee,
# and without a written agreement is hereby granted, provided that the above
# copyright notice, this paragraph and the following three paragraphs appear
# in all copies.
#
# Permission to make commercial use of this software may be obtained by
# contacting:
# Technology Transfer Office
# 9500 Gilman Drive, Mail Code 0910
# University of California
# La Jolla, CA 92093-0910
# (858) 534-5815
# invent@ucsd.edu
#
# This software program and documentation are copyrighted by The Regents of the
# University of California. The software program and documentation are supplied
# "as is", without any accompanying services from The Regents. The Regents does
# not warrant that the operation of the program will be uninterrupted or
# error-free. The end-user understands that the program was developed for
# research purposes and is advised not to rely exclusively on the program for
# any reason.
#
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO
# ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
# CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING
# OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
# EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF
# CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
# THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
# CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
# ENHANCEMENTS, OR MODIFICATIONS.

#Author: Eric Fung
#Contact: kspham@eng.ucsd.edu

# using global variable to store files names
# not a good practice and need to change later when release

illumina_graph_file = "../data/ecoli/ecoli_miseq-contigs.dot"
illumina_contig_file = "../data/ecoli/ecoli_miseq-contigs.fa"
contig_pacbiomap_file = "../data/ecoli/ecoli_pacbio_contigs_mapping.fasta.m4"
contig_refmap_file = "../data/ecoli/ecoli_contigs_ref_mapping.blastn"

pacbio_refmap_file = "../data/ecoli/ecoli_pacbio_ref_mapping.fasta.m4"



"""
illumina_graph_file = "../data/SAureusTW20/SAureusTW20_miseq-contigs.dot"
illumina_contig_file = "../data/SAureusTW20/SAureusTW20_miseq-contigs.fa"

contig_pacbiomap_file = "../data/SAureusTW20/SAureusTW20_pacbio_miseqcontig_mapping.fasta.m4"
contig_refmap_file = "../data/SAureusTW20/SAureusTW20_miseqcontigs_ref_mapping.blastn"
pacbio_refmap_file = "../data/SAureusTW20/SAureusTW20_pacbio_ref_mapping.fasta.m4"
"""
path_extend_graph_file = "../data/ecoli/path_extend_graph.dot"

path_extend_png_file = "../data/ecoli/path_extend_graph.png"
