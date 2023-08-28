# Non-Coding RNA Prediction Pipeline

This a R package implementing a pipeline to predict non-coding RNA (ncRNA)
using alignments of highly similar DNA sequences.

## Background

The discovery of protein families like IscBs, TnpBs and Fanzors has
demonstrated the existence of programmable RNA-guided protein systems that
encode their guide RNA adjacent to or overlapping with their CDS. This
configuration immediately suggests a method to finding these systems. By taking
advantage of highly similar protein sequences, we can generate multiple
sequence alignments of DNA to detect sequence conservation near the CDS.
Although DNA sequence conservation does not mean some ncRNA exists, it does
offer a useful heuristic that can be used in downstream processing.

## Dependencies

* mafft (>= 7.490)
* R-scape (>= 2.0.0.q)
* ViennaRNA (>= 2.5.1)
* R (>= 4.1.0)
* R packages
	* Biostrings (>= 2.62.0)
	* changepoint (>= 2.2.3)
	* digest (>= 0.6.31)
	* ggplot2 (>= 3.3.6)
	* gggenes (>= 0.4.1)
	* jsonlite (>= 1.8.2)
	* IRanges (>= 2.32.0)
	* patchwork (>= 1.1.2)
	* svglite (>= 2.1.0)
	* xml2 (>= 1.3.0)
    * optparse (>= 1.7.3) (only needed to run the helper script)

## Pipeline

1. Parse input files
2. Make multiple sequence alignments of the 3' and 5' ends of the loci
3. Compute a numerical conservation score for each column of the alignments
4. Predict a changepoint for the conservation scores.
5. Fold the alignments with a detected changepoint
6. Generate output files

## Input

The pipeline requires two files: a FASTA file of DNA sequences of the loci of
interest and a tab-delimited file containing the coordinates of the CDS in each
sequence in the FASTA file. The input should represent a cluster of highly
similar (>= 70% sequence identity) proteins along with upstream and downstream
DNA context. Example input is the `examples` folder.

## Output

There are several output files. The main ones are
* `alignment3.fna` - Alignment of 3' end of the loci
* `alignment5.fna` - Alignment of 5' end of the loci
* `conservation.svg` - Visualization of the conservation
* `report.json` - JSON file of some of the results

## Installation/Usage

```
make build
make install
./predict.R [options] FASTA COORDINATES OUTPUT
```
