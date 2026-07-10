# Getting started with rKATS

rKATS provides R functions for counting k-mers, measuring enrichment, turning
enriched k-mers into positional weight matrices (PWMs), and plotting sequence
logos.

## Installation

The package requires R >= 3.5 and CMake >= 3.12.0. Install it from GitHub with:

```r
install.packages("remotes")
remotes::install_github("Ferny730JR/rkatss")
library(rkats)
```

## Accepted sequence inputs

The sequence functions accept either a file path or sequences supplied in R.

| Input | Supported forms |
| --- | --- |
| File | FASTA, FASTQ, raw sequences, or gzip-compressed versions |
| Multiple sequences in R | Character vector or list |
| One sequence in R | A one-element list, such as `list("ACGTACGT")` |

When a single character string is supplied, rKATS treats it as a file path.
Wrap a single sequence in a list to avoid that interpretation.

## One-FASTA motif workflow

Assume `bound.fa` contains the RNA sequences of interest. The first call below
is optional and reports k-mer abundance. The enrichment result, rather than
the count result, is used to build the PWM.

```r
library(rkats)

# Stores the counts of the k-mers (optional)
counts <- count_kmers("bound.fa", kmer = 5, sort = TRUE)

# Compute the k-mer enrichment of the sequences
enrichments <- enrichments(
  "bound.fa",
  algo = "shuffled",
  kmer = 5,
  seed = 1,
  sort = TRUE
)

# Create sequence logo from enrichments
pwms <- get_pwms(enrichments, type = "RNA")
logo <- plot_logo(pwms, method = "bits", title = "Enriched RNA motif")

# Save the logo as a PNG
ggplot2::ggsave(
  "motif-logo.png",
  plot = logo,
  width = 6,
  height = 3,
  dpi = 300
)
```

The objects produced by this workflow have the following roles:

| Object | Expected result |
| --- | --- |
| `counts` | Data frame with `kmer` and `count` columns |
| `enrichments` | Data frame with `kmer` and `rval` columns |
| `pwms` | List with `pwm` matrices and contributing-k-mer counts in `n` |
| `logo` | A `ggplot` object |
| `motif-logo.png` | Exported sequence-logo image |

The exact enriched k-mers and numeric values depend on the input data. Setting
`seed = 1` makes the shuffled-background calculation reproducible.

## Using bundled example data

The package includes 1,000 bound and 1,000 control RBFOX2 sequences:

```r
data(rbfox2_seqs)

enrichment <- enrichments(
  rbfox2_seqs$bound,
  rbfox2_seqs$input,
  kmer = 5,
  sort = TRUE
)

pwms <- get_pwms(enrichment)
plot_logo(pwms, title = "RBFOX2 motif")
```

The same functions accept file paths, so vectors can be written to FASTA or
raw sequence files when building a file-based pipeline.

## Common input and output pitfalls

- `count_kmers()` returns `count`; `get_pwms()` expects `rval` from
  `enrichments()` or `ikke()`.
- Use `type = "RNA"` in `get_pwms()` when the logo should use `U`; the default
  is DNA and uses `T`.
- Bootstrap enrichment results include additional uncertainty columns. The
  standard non-bootstrap result contains `kmer` and `rval`.
- Use `ggplot2::ggsave()` to export the plot returned by `plot_logo()`.

Next, see [workflow options](workflows.md) or the [function reference](reference.md).
