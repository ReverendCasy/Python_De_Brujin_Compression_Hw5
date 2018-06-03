# Python_De_Brujin_Compression_Hw5

## Description
De Bruijin graph assembler operates via prefix sequence to assemble nucleotide strings over a linear (or at least nearly linear) time. Here we implement custom De Bruijin assembler that can be executed from command line. In comparison to [previous version of the code](https://github.com/PreacherCasy/Python_De_Brujin_Hw4), it provides opportunity to compress the grph and get more or less informative contigs, though for now it lacks ability to resolve artifacts, such as 'bones' and 'bubbles'.

## Code structure
Basically, code comprises of two parts:
+ code body: object-oriented code comprising of three class objects for graph assembly
+ *argparse* part: integration into command line and code execution

### Code body
The body holds three classes essential for graph assembly
+ *Class **Vertex***: a class for graph vertices which stand for individual k-mers. It contains information about coincindent edges and k-mer coverage
+ *Class **Edge***: a class for graph edges representing existing reads. This class contains information about edge coverage calculated as mean between two interlinked k-mers. It also contains two functions for edge merging to operate graph compression.
+ *Class **Graph***: the largest classes representing the graph itself. It contains three functions: one for read addition (*add_read*), one for vertex coverage calculation (*calc_init_coverage*) and one for graph visualization (*visualize*) producing a pdf object either with full text (if **t=full**) or only with numeric stats (in any other case).

### Compression
Apart from two properties of **Edge** class, compression option comprises of two functions:
+ *function **merge***: allows joining of two vertices neighbouring over a redundant node. Works only if redundant node has ony one incoming and one outcoming edge;
+ *function **compress***: performs **merge** over all vertices meeting the "one in - one out" requirement;
+ *function **cut***: a function for cutting of low-coverage marginal regions.
Functions **compress** and **cut** a repeated sequentially in a cycle embodied in *argparse* ection of the code unless resulting number of vertices comes close to a chosen baseline.

### Contig obtaining
*Function **get contigs*** provides opportunity to download resulting contigs as *fasta*-flavoured entries. The function is called automaticaly; if called via **-o** flag, it writes contigs into a *fasta* file

### Argparse
Arguments for command line go as following:
+ **-i**: a path to fasta file containing reads for further assembly;
+ **-k**: a desired k (k-mer length); default is 3;
+ **-t**: if followed by 'full', draws a full graph with k-mers on vertices and read sequences on edges; if followed by any other random text, produces a stat graph with node and edge coverage
+ **-f** and **-r**: two alternative flags for assembly mode. If **-f** is stated, graph is assembled fron raw read sequences; if **-r** is chosen, graph is assemble from read reverse complements;
+ **-b**: a vertex number baseline for cutting off marginal vertices;
+ **-o**: a voluntary flag specifying the output location for compression results;
+ **function**: a none-flag variable for **t**/**r** alteration handling.

## Launch example
The following command executes graph assembly from forward reads with k equal to 55 and baseline of 15 with compression performed (see results in 3.fasta):
```{bash}
python De_Bruijin_Compression.py -i hw3_dataset.fasta -c -f -k 55 -b -o 3.fasta
```
Alternatively, this one will create a stat graph (by random text after **-t**) from reverse reads with k=10:
 ```{bash}
python De_Bruijin_Compression.py -r -t pop -k 10 -i hw3_dataset.fasta
```

## Assembly
In the testing case, we got a *fasta* file (see attachments in this repository) with short Illumina-like reads. By performing assemly with *k* equal to 55, we got larger contigs of ~140-150 bp long (see attachment in this repository). These reads were aligned against **NCBI Nucleotide** database and proved to be a part of human hepatite B virus genome.

## Acknowledgements
Eugene Bakin from Bioinformatics Institute for code backbone and testing material
Pre-existing assemblers for graph visualization function logic
