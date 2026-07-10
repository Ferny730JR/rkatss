# rKATS reference

The functions below are exported by the `rkats` package. Each entry links to
the corresponding package help page in `man/`.

## Counting

| Function | Purpose | Reference |
| --- | --- | --- |
| `count_kmers()` | Count all k-mer occurrences in sequence input | [`count_kmers`](../man/count_kmers.Rd) |
| `kmer_presence()` | Count how many sequences contain each k-mer | [`kmer_presence`](../man/kmer_presence.Rd) |

## Enrichment

| Function | Purpose | Reference |
| --- | --- | --- |
| `enrichments()` | Calculate k-mer enrichment against a control or modeled background | [`enrichments`](../man/enrichments.Rd) |
| `ikke()` | Perform iterative k-mer knockout enrichment analysis | [`ikke`](../man/ikke.Rd) |

## Motifs and logos

| Function | Purpose | Reference |
| --- | --- | --- |
| `align_kmers()` | Align enriched k-mers around the most enriched k-mer | [`align_kmers`](../man/align_kmers.Rd) |
| `get_pwms()` | Convert aligned enriched k-mers into one or more PWMs | [`get_pwms`](../man/get_pwms.Rd) |
| `plot_logo()` | Plot PWMs as bits or probability sequence logos | [`plot_logo`](../man/plot_logo.Rd) |
| `cor.pwm()` | Compare two PWMs by their best local correlation | [`cor.pwm`](../man/cor.pwm.Rd) |

## Sequence search

| Function | Purpose | Reference |
| --- | --- | --- |
| `seqseq()` | Find one sequence within another, with optional all-match output | [`seqseq`](../man/seqseq.Rd) |

## Bundled data

| Dataset | Contents | Reference |
| --- | --- | --- |
| `rbfox2_seqs` | 1,000 bound and 1,000 control RBFOX2 sequences | [`rbfox2_seqs`](../man/rbfox2_seqs.Rd) |
| `rbfox2_enrichments` | 5-mer enrichment data for RBFOX2 | [`rbfox2_enrichments`](../man/rbfox2_enrichments.Rd) |
| `rbfox2_pwms` | RBFOX2 5-mer PWM matrices | [`rbfox2_pwms`](../man/rbfox2_pwms.Rd) |

Load bundled data with `data()`:

```r
data(rbfox2_seqs)
data(rbfox2_enrichments)
data(rbfox2_pwms)
```

## Input and output conventions

The counting and enrichment functions accept FASTA, FASTQ, raw sequence, and
gzip-compressed files, as well as character vectors and lists. A single
sequence supplied in R should be wrapped in a list.

Typical outputs are:

- counting functions: a data frame with `kmer` and `count`
- enrichment functions: a data frame with `kmer` and `rval`
- bootstrap enrichment: additional `stdev` and `pval` columns
- `get_pwms()`: a list with `pwm` matrices and contributing-k-mer counts in
  `n`
- `plot_logo()`: a `ggplot` object

The main counting and enrichment functions support k-mers up to length 16;
`kmer_presence()` supports k-mers up to length 15. See the individual help
pages above for argument defaults and algorithm-specific parameters.
