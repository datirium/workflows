Sequence reads are first cleaned from adapters and transformed into fully bisulfite-converted forward (C->T) and reverse
read (G->A conversion of the forward strand) versions, before they are aligned to similarly converted versions of the
genome (also C->T and G->A converted). Sequence reads that produce a unique best alignment from the four alignment processes
against the bisulfite genomes (which are running in parallel) are then compared to the normal genomic sequence and the
methylation state of all cytosine positions in the read is inferred. A read is considered to align uniquely if an alignment
has a unique best alignment score (as reported by the AS:i field). If a read produces several alignments with the same number
of mismatches or with the same alignment score (AS:i field), a read (or a read-pair) is discarded altogether.

On the next step we extract the methylation call for every single C analysed. The position of every single C will be written
out to a new output file, depending on its context (CpG, CHG or CHH), whereby methylated Cs will be labelled as forward
reads (+), non-methylated Cs as reverse reads (-). The output of the methylation extractor is then transformed into a bedGraph
and coverage file. The bedGraph counts output is then used to generate a genome-wide cytosine report which reports the number
on every single CpG (optionally every single cytosine) in the genome, irrespective of whether it was covered by any reads or not.
As this type of report is informative for cytosines on both strands the output may be fairly large (~46mn CpG positions or >1.2bn
total cytosine positions in the human genome).