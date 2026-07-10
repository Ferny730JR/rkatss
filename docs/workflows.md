# rKATS workflows

## Test/control enrichment

When a control set is available, pass the bound/test and control files to
`enrichments()` using the default `normal` algorithm:

```r
enrichment <- enrichments(
  "bound.fa",
  "control.fa",
  kmer = 5,
  sort = TRUE
)

pwms <- get_pwms(enrichment, type = "RNA")
plot_logo(pwms, title = "Motifs enriched in bound sequences")
```

The same call works with character vectors or lists of sequences. A control
file is optional when using a shuffled or probabilistic background.

## Counting versus enrichment

Use the function that matches the question being asked:

| Function | Measures | Main output |
| --- | --- | --- |
| `count_kmers()` | Total k-mer occurrences | `kmer`, `count` |
| `kmer_presence()` | Number of sequences containing each k-mer, counting each sequence at most once | `kmer`, `count` |
| `enrichments()` | Relative enrichment against a control or modeled background | `kmer`, `rval` |
| `ikke()` | Iterative k-mer knockout enrichments | `kmer`, `rval` |

Counts are useful for quality checks and abundance summaries. Use enrichment
results for `get_pwms()` because PWM weights are derived from the `rval` column.

## Background algorithms

`enrichments()` supports four algorithms:

| `algo` | Use case |
| --- | --- |
| `"normal"` | Compare test sequences with a supplied control set |
| `"shuffled"` | Estimate a background by shuffling the test sequences |
| `"probabilistic"` | Use the probabilistic background model |
| `"shuf+prob"` | Combine shuffled and probabilistic background calculations |

For shuffled analyses, set `seed` when a reproducible result is needed. The
`klet` argument controls the k-let length preserved during shuffling.

```r
enrichment <- enrichments(
  "bound.fa",
  algo = "shuffled",
  kmer = 5,
  klet = 2,
  seed = 1
)
```

## Bootstrap enrichment

Set `bootstrap_iters` to estimate variability. `sample` is the percentage of
sequences sampled per iteration:

```r
bootstrap <- enrichments(
  "bound.fa",
  "control.fa",
  kmer = 5,
  bootstrap_iters = 100,
  sample = 55.55,
  seed = 1
)

head(bootstrap)
```

Bootstrap results add `stdev` and, for enrichment calculations, `pval` columns.

## PWM generation and secondary motifs

`get_pwms()` aligns enriched k-mers to the most enriched k-mer and can produce
more than one PWM:

```r
pwms <- get_pwms(
  enrichment,
  num_pwms = 2,
  limit = 100,
  limit_per_logo = 20,
  type = "RNA"
)
```

- `limit` limits the total number of enrichment rows considered.
- `limit_per_logo` limits the number of rows contributing to each logo.
- `num_pwms` requests primary and secondary motif PWMs.
- `type = "RNA"` labels the fourth base as `U`; `type = "DNA"` uses `T`.

## Plotting and export

`plot_logo()` returns a `ggplot` object, so it can be printed, customized, or
exported with `ggsave()`:

```r
logo <- plot_logo(
  pwms,
  method = "probability",
  title = "RBFOX2 motifs",
  subtitle = "n"
)

print(logo)
ggplot2::ggsave("rbfox2-motifs.png", logo, width = 7, height = 4, dpi = 300)
```

Use `method = "bits"` for information-content logos or
`method = "probability"` for probability-based logos.

## Iterative knockout enrichment

`ikke()` repeatedly removes enriched k-mers to identify additional motifs:

```r
result <- ikke(
  "bound.fa",
  "control.fa",
  kmer = 5,
  iterations = 10,
  normalize = TRUE
)
```

For a probabilistic IKKE analysis, set `probabilistic = TRUE`; a control file
is not used in that mode.

## Reproducibility and performance

- Set `seed` for bootstrap and shuffled calculations.
- Keep the k-mer length as small as the analysis allows; the output space grows
  exponentially with k.
- Use `sort = TRUE` when the most enriched or abundant k-mers should appear
  first in the returned data frame.
- `threads` is available on the counting and enrichment functions, but its
  performance impact depends on the operation and platform.

For the complete argument lists and return values, see the [reference guide](reference.md)
and the linked Rd help pages.
