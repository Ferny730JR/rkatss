# rKATS

<img src="https://github.com/user-attachments/assets/119df30c-61ab-43a5-9b66-62df67708578" alt="KATSS Logo of a blue cat with a RNA-shaped tail." width="300">

An R wrapper of KATSS for analyzing RNA sequences, identifying enriched
RNA-binding-protein motifs, and creating sequence logos.

## What rKATS provides

- k-mer counting and k-mer-presence analysis
- enrichment calculations with test/control or shuffled backgrounds
- bootstrap and iterative k-mer knockout workflows
- alignment and conversion of enriched k-mers into PWMs
- sequence-logo plotting and PWM comparison utilities
- support for FASTA, FASTQ, raw sequence, and gzip-compressed input files

## Installation

rKATS requires R >= 3.5 and CMake >= 3.12.0. Install the package from GitHub:

```r
install.packages("remotes")
remotes::install_github("Ferny730JR/rkatss")
```

Load the `rkats` package by using:

```r
library(rkats)
```

## FASTA-to-logo quick start

The basic workflow is:

```text
FASTA input -> count/enrichment -> PWM generation -> logo plot/export
```

For a single FASTA file, `algo = "shuffled"` creates a shuffled background for
the enrichment calculation. Here `bound.fa` is a FASTA file containing the
sequences of interest.

```r
library(rkats)

counts <- count_kmers("bound.fa", kmer = 5, sort = TRUE)

enrichments <- enrichments(
  "bound.fa",
  algo = "shuffled",
  kmer = 5,
  seed = 1,
  sort = TRUE
)

pwms <- get_pwms(enrichments, type = "RNA")
logo <- plot_logo(pwms, method = "bits", title = "Enriched RNA motif")

ggplot2::ggsave(
  "motif-logo.png",
  plot = logo,
  width = 6,
  height = 3,
  dpi = 300
)
```

Expected objects and files:

- `counts`: a data frame with `kmer` and `count` columns.
- `enrichments`: a data frame with `kmer` and `rval` columns. Bootstrap runs
  can add `stdev` and `pval` columns.
- `pwms`: a list containing `pwm` matrices and the number of contributing
  k-mers in `n`.
- `logo`: a `ggplot` object that can be printed or customized.
- `motif-logo.png`: the exported sequence-logo image.

`count_kmers()` is useful for inspecting k-mer abundance, but its `count`
column is not the input expected by `get_pwms()`. Use the `rval` output from
`enrichments()` or `ikke()` to generate PWMs.

## Input formats

All main sequence-analysis functions accept a file path or a character vector
or list of sequences. Files may contain FASTA, FASTQ, or raw sequences, and
gzip-compressed files are supported. A single sequence should be wrapped in a
list so it is not interpreted as a file path:

```r
count_kmers(list("ACGTACGT"), kmer = 3)
```

For a direct test/control comparison, pass the test and control FASTA files to
`enrichments()`:

```r
enrichment <- enrichments("bound.fa", "control.fa", kmer = 5)
```

## Documentation

- [Getting started](docs/getting-started.md): installation, inputs, and the
  first FASTA-to-logo workflow
- [Workflows](docs/workflows.md): test/control analysis, backgrounds,
  bootstrap, PWMs, plotting, and advanced analyses
- [Reference](docs/reference.md): exported functions, bundled data, and
  links to the package help pages

The bundled RBFOX2 examples can be loaded with:

```r
data(rbfox2_seqs)
data(rbfox2_enrichments)
data(rbfox2_pwms)
```

## License

rKATS is released under the [MIT License](LICENSE).
