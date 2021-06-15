What is MAnorm?
--------------

MAnorm is a robust model for quantitative comparison of ChIP-Seq data sets of TFs (transcription factors) or epigenetic modifications and you can use it for:

* Normalization of two ChIP-seq samples
* Quantitative comparison (differential analysis) of two ChIP-seq samples
* Evaluating the overlap enrichment of the protein binding sites(peaks)
* Elucidating underlying mechanisms of cell-type specific gene regulation

How MAnorm works?
----------------

MAnorm uses common peaks of two samples as a reference to build the rescaling model for normalization, which is based on the empirical assumption that if a chromatin-associated protein has a substantial number of peaks shared in two conditions, the binding at these common regions will tend to be determined by similar mechanisms, and thus should exhibit similar global binding intensities across samples. The observed differences on common peaks are presumed to reflect the scaling relationship of ChIP-Seq signals between two samples, which can be applied to all peaks.

What do the inputs mean?
----------------

### General

**Experiment short name/Alias**

* short name for you experiment to identify among the others

**ChIP-Seq PE sample 1**
* previously analyzed ChIP-Seq paired-end experiment to be used as Sample 1

**ChIP-Seq PE sample 2**

* previously analyzed ChIP-Seq paired-end experiment to be used as Sample 2

**Genome**

* Reference genome to be used for gene assigning

### Advanced

**Reads shift size for sample 1**

* This value is used to shift reads towards 3' direction to determine
the precise binding site. Set as half of the fragment length. Default 100

**Reads shift size for sample 2**

* This value is used to shift reads towards 5' direction to determine
the precise binding site. Set as half of the fragment length. Default 100

**M-value (log2-ratio) cutoff**

* Absolute M-value (log2-ratio) cutoff to define biased (differential binding)
peaks. Default: 1.0

**P-value cutoff**

* P-value cutoff to define biased peaks. Default: 0.01

**Window size**

* Window size to count reads and calculate read densities. 2000 is recommended for
sharp histone marks like H3K4me3 and H3K27ac, and 1000 for TFs or DNase-seq.
Default: 2000