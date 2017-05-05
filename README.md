PacificBlue Scaffolding Tool for PacBio Long Reads

Developed by Dan Browne

<<<<<<< HEAD
This tool is a rudimentary implementation of a scaffolding algorithm that utilizes BLASR alignments of PacBio long reads to a genome assembly. The algorithm analyzes each PacBio read individually, stacking the alignments in an array and calculating alignment coverage along the read. This idea was inspired by Gene Myers' work on DAZZLER that he presented at the 2016 JGI Users Meeting (https://www.youtube.com/watch?v=8vg3225G8gU). The goal of analyzing the alignment coverage along each read is to find unique alignments and exclude repetitive alignments. Thus unique alignments on a read that span unassembled gaps could provide the basis for linking two sequences together (i.e. scaffolding).
=======
This is an experimental tool still under active development. At this point in time, I do not recommend using this tool without first contacting me. If you'd like to try it out or help with development, feel free to send me an email: dbrowne at tamu dot edu

Small tweak
>>>>>>> 5e1b05751252353a0fc8cf2fdfbd8aeffa205723
