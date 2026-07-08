#' K-mer Counting
#'
#' Count the k-mers in a FASTQ, FASTA, raw sequences file, or a vector/list
#' of sequences.
#'
#' @param input Name of the file which you want to count k-mers from, or a
#' character vector/list containing sequences. Files must contain raw sequences,
#' FASTA, or FASTQ data. Gzip-compressed files are supported. A single character
#' string is interpreted as a file path; vectors with multiple elements and
#' lists are interpreted as sequences.
#' @param kmer Length of the k-mer you want to count. Currently, only k-mers up
#' to length 16 are supported.
#' @param algo Whether to perform regular counts, or count shuffled sequences.
#' @param bootstrap_iters Number of iterations to bootstrap.
#' @param sample Percent to subsample during bootstrap (should be between 0-100%).
#' @param seed Specify the seed to be used by bootstrap. Since bootstrap
#' subsamples random sequences, seeding alters which random sequences will be
#' picked. This helps to ensure deterministic output which can be achieved by
#' using the same seed. To pick a random seed, set `seed=-1`.
#' @param klet Specify the k-let length to preserve during shuffling. This only
#' affects the output if `algo="shuffled"` is set. -1 chooses the default value.
#' @param sort Sort based on the counts from highest to lowest. Currently,
#' the output given is sorted based on kmers (AA... first, TT... last).
#' @param threads Number of threads to use. Currently not well optimized.
#'
#' @return Dataframe containing the counts for all k-mers.
#' @useDynLib rkats, .registration = TRUE
#' @export
#'
#' @examples
#' data(rbfox2_seqs)
#'
#' # Count di-nucleotides from a vector of sequences
#' count_kmers(rbfox2_seqs$bound, kmer = 2)
#'
#' # Count a single sequence by wrapping it in a list, otherwise it will treat
#' # it as a filename
#' count_kmers(list("ACGTACGT"), kmer = 2)
#'
#' # Create temporary file with sequences
#' tf <- tempfile()
#' writeLines(rbfox2_seqs$bound, tf)
#'
#' # Count di-nucleotides in file
#' count_kmers(tf, kmer = 2)
#'
#' # Count mono-nucleotides in file
#' count_kmers(tf, kmer = 1)
#'
#' # Count shuffled kmers
#' count_kmers(tf, kmer = 1, algo = "shuffled")
#'
#' # Specify k-let to preserve during shuffling
#' count_kmers(tf, kmer = 1, algo = "shuffled", klet = 2)
#'
#' # Count bootstrap kmers
#' result <- count_kmers(tf, bootstrap_iters = 100)
#' head(result)
#'
#' # Subsample 55.55% of the file per bootstrap iteration
#' result <- count_kmers(tf, bootstrap_iters = 100, sample = 55.55)
#' head(result)
#'
#' # Count bootstrap shuffled kmers
#' result <- count_kmers(tf, algo = "shuffled", bootstrap_iters = 100)
#' head(result)
#'
#' # Sort by count
#' result <- count_kmers(tf, kmer = 5, sort = TRUE)
#' head(result)
#'
#' # Cleanup file
#' unlink(tf)
count_kmers <- function(input, kmer = 3, algo=c("regular","shuffled"),
                        bootstrap_iters = 0, sample = 25, seed = -1, klet = -1,
                        sort = FALSE, threads = 1) {
  if(!is.character(input) && !is.list(input))
    stop("input must be a file path or a vector/list of sequences")
  if(!is.numeric(kmer) || kmer %% 1 != 0)
    stop("kmer must be an integer")
  if(!is.numeric(bootstrap_iters) || bootstrap_iters %% 1 != 0)
    stop("bootstrap_iters must be an integer")
  if(!is.numeric(sample) || sample <= 0 || 100 < sample)
    stop("sample must be a number between 0-100")
  if(!is.numeric(seed) || seed %% 1 != 0)
    stop("seed must be an integer")
  if(!is.numeric(klet) || klet %% 1 != 0)
    stop("klet must be an integer")
  if(!is.logical(sort))
    stop("sort must be logical")
  if(!is.numeric(threads) || threads %% 1 != 0)
    stop("threads must be an integer")

  if(is.list(input) || length(input) > 1) {
    input <- unlist(input, use.names = FALSE)
    if(!is.character(input) || anyNA(input))
      stop("sequences must be character strings")
    tmp <- tempfile()
    on.exit(unlink(tmp), add = TRUE)
    writeLines(input, tmp)
    input <- tmp
  } else {
    input <- path.expand(input)
  }

  sample <- as.integer((sample*1000) %% 100001)
  algo <- match.arg(algo)
  algo <- if(algo == "regular") 1 else 2

  .Call("count_kmers_R",
        input,
        as.integer(kmer),
        as.integer(klet),
        as.integer(sort),
        as.integer(bootstrap_iters),
        as.integer(sample),
        as.integer(algo),
        as.integer(seed),
        as.integer(threads))
}


#' Calculate k-mer enrichments
#'
#' @param test Test sequences. Either a file path or a character vector/list
#' containing sequences. Files must contain raw sequences, FASTA, or FASTQ data.
#' Gzip-compressed files are supported. A single character string is interpreted
#' as a file path; vectors with multiple elements and lists are interpreted as
#' sequences.
#' @param ctrl Control sequences (optional). Either `NULL`, a file path, or
#' a character vector/list containing sequences. The same formats and input
#' rules as `test` apply.
#' @param kmer Length of the k-mer to compute enrichments for. Currently, only
#' k-mers up to length 16 are supported.
#' @param algo The algorithm to use for computing enrichments.
#' @param bootstrap_iters Number of iterations to bootstrap.
#' @param sample Percent to subsample during bootstrap (should be between 0-100%).
#' @param seed Specify the seed to be used by bootstrap. Since bootstrap
#' subsamples random sequences, seeding alters which random sequences will be
#' picked. This helps to ensure deterministic output which can be achieved by
#' using the same seed. To pick a random seed, set `seed=-1`.
#' @param klet Specify the k-let length to preserve during shuffling. This only
#' affects the output if `algo="shuffled"` or `algo="shuf+prob"` is set. -1
#' chooses the default recommended value.
#' @param sort Sort data.frame based on the counts from highest to lowest.
#' Currently, the output given is sorted alphabetically based on k-mers.
#' @param threads Number of threads to use. Currently not well optimized/not
#' working.
#'
#' @return Data frame containing the k-mer enrichments.
#' @useDynLib rkats, .registration = TRUE
#' @export
#'
#' @examples
#' # Load data
#' data(rbfox2_seqs)
#'
#' # Use vectors of sequences directly
#' result <- enrichments(rbfox2_seqs$bound, rbfox2_seqs$input)
#' head(result)
#'
#' # Use a single sequence by wrapping it in a list
#' result <- enrichments(list("ACGTACGT"), algo = "shuffled")
#' head(result)
#'
#' # Create raw sequence files
#' test_file <- tempfile()
#' ctrl_file <- tempfile()
#' writeLines(rbfox2_seqs$bound, test_file)
#' writeLines(rbfox2_seqs$input, ctrl_file)
#'
#' ## Get enrichments when you have a test and control dataset
#' result <- enrichments(test_file, ctrl_file)
#' head(result)
#'
#' result <- enrichments(test_file, ctrl_file, kmer = 5)
#' head(result)
#'
#' ## Get enrichments without a control
#' result <- enrichments(test_file, algo = "shuffled", kmer = 5)
#' head(result)
#'
#' result <- enrichments(test_file, algo = "shuffled", kmer = 5, klet = 5)
#' head(result)
#'
#' result <- enrichments(test_file, algo = "probabilistic", kmer = 5)
#' head(result)
#'
#' ## Enabling bootstrap
#' result <- enrichments(test_file, ctrl_file, bootstrap_iters = 100)
#' head(result)
#'
#' result <- enrichments(test_file, ctrl_file, bootstrap_iters = 100,
#'                       sample = 55.55)
#' head(result)
#'
#' result <- enrichments(test_file, algo = "shuffled",
#'                       bootstrap_iters = 100, kmer = 5, klet = 5)
#' head(result)
#'
#' # Cleanup files
#' unlink(test_file)
#' unlink(ctrl_file)
enrichments <- function(test, ctrl = NULL, kmer = 3,
                        algo = c("normal", "shuffled", "probabilistic", "shuf+prob"),
                        bootstrap_iters = 0, sample = 25, seed = -1, klet = -1,
                        sort = TRUE, threads = 1)
{
  if(!is.character(test) && !is.list(test))
    stop("test must be a file path or a vector/list of sequences")
  if(!is.character(ctrl) && !is.list(ctrl) && !is.null(ctrl))
    stop("ctrl must be NULL, a file path, or a vector/list of sequences")
  if(!is.numeric(kmer) || kmer %% 1 != 0)
    stop("kmer must be an integer")
  if(!is.numeric(bootstrap_iters) || bootstrap_iters %% 1 != 0)
    stop("bootstrap_iters must be an integer")
  if(!is.numeric(sample) || sample <= 0 || 100 < sample)
    stop("sample must be a number between 0-100")
  if(!is.numeric(seed) || seed %% 1 != 0)
    stop("seed must be an integer")
  if(!is.numeric(klet) || klet %% 1 != 0)
    stop("klet must be an integer")
  if(!is.logical(sort))
    stop("sort must be logical")
  if(!is.numeric(threads) || threads %% 1 != 0)
    stop("threads must be an integer")
  if(16 >= kmer && kmer > 12) {
    menu_title <- paste(convert_bytes(4^kmer * 176),
                        "Are you sure you want to proceed?")
    if(utils::menu(c("Yes", "No! Fix your program!"),
                   title = menu_title) == 2)
      return(NULL)
  }

  if(is.list(test) || length(test) > 1) {
    test <- unlist(test, use.names = FALSE)
    if(!is.character(test) || anyNA(test))
      stop("test sequences must be character strings")
    test_tmp <- tempfile()
    on.exit(unlink(test_tmp), add = TRUE)
    writeLines(test, test_tmp)
    test <- test_tmp
  } else {
    test <- path.expand(test)
  }

  if(!is.null(ctrl)) {
    if(is.list(ctrl) || length(ctrl) > 1) {
      ctrl <- unlist(ctrl, use.names = FALSE)
      if(!is.character(ctrl) || anyNA(ctrl))
        stop("control sequences must be character strings")
      ctrl_tmp <- tempfile()
      on.exit(unlink(ctrl_tmp), add = TRUE)
      writeLines(ctrl, ctrl_tmp)
      ctrl <- ctrl_tmp
    } else {
      ctrl <- path.expand(ctrl)
    }
  }

  algo <- match.arg(algo)
  sample <- as.integer((sample*1000) %% 100001)
  if(algo == "normal") {
    algo <- 0
  } else if(algo == "shuffled") {
    algo <- 1
  } else if(algo == "probabilistic") {
    algo <- 2
  } else {
    algo <- 3
  }

  .Call("enrichments_R",
        test,
        ctrl,
        as.integer(kmer),
        as.integer(algo),
        as.integer(bootstrap_iters),
        as.integer(sample),
        as.integer(seed),
        as.integer(klet),
        as.integer(sort),
        as.integer(threads))
}


#' Iterative K-mer Knockout Enrichments
#'
#' @param testfile Test sequences. Either a file path or a character vector/list
#' containing sequences. Files must contain raw sequences, FASTA, or FASTQ data.
#' Gzip-compressed files are supported. A single character string is interpreted
#' as a file path; vectors with multiple elements and lists are interpreted as
#' sequences.
#' @param ctrlfile Control sequences (optional). Either `NULL`, a file path, or
#' a character vector/list containing sequences. The same formats and input
#' rules as `testfile` apply.
#' @param kmer Length of k-mer.
#' @param iterations Number of iterations to perform.
#' @param normalize Normalize enrichments to log2.
#' @param threads Number of threads to use. Specifying less than 1 thread sets
#' the number of threads as 1.
#' @param probabilistic Calculate probabilistic enrichments.
#'
#' @return Data frame containing the enrichments.
#' @useDynLib rkats, .registration = TRUE
#' @export
#'
#' @examples
#' # Load data
#' data(rbfox2_seqs)
#'
#' # Use a vector of sequences directly
#' result <- ikke(rbfox2_seqs$bound, probabilistic = TRUE)
#' head(result)
#'
#' # Use a single sequence by wrapping it in a list
#' result <- ikke(list("ACGTACGT"), kmer = 2, probabilistic = TRUE)
#' head(result)
#'
#' # Use test and control sequence vectors
#' result <- ikke(rbfox2_seqs$bound, rbfox2_seqs$input, kmer = 5)
#' head(result)
#'
#' # Specify the number of enrichments to obtain
#' result <- ikke(rbfox2_seqs$bound, rbfox2_seqs$input, kmer = 5,
#'                iterations = 2)
#' print(result)
#'
#' # Normalize enrichments to log2
#' result <- ikke(rbfox2_seqs$bound, rbfox2_seqs$input, kmer = 5,
#'                normalize = TRUE)
#' head(result)
#'
#' # File paths are also supported
#' test_file <- tempfile()
#' ctrl_file <- tempfile()
#' writeLines(rbfox2_seqs$bound, test_file)
#' writeLines(rbfox2_seqs$input, ctrl_file)
#'
#' result <- ikke(test_file, ctrl_file, kmer = 5)
#' head(result)
#'
#' unlink(test_file)
#' unlink(ctrl_file)
ikke <- function(testfile, ctrlfile = NULL, kmer = 3, iterations = 10,
                 probabilistic = FALSE, normalize = FALSE, threads = 1) {
  if(!is.character(testfile) && !is.list(testfile))
    stop("testfile must be a file path or a vector/list of sequences")
  if(!is.character(ctrlfile) && !is.list(ctrlfile) && !is.null(ctrlfile))
    stop("ctrlfile must be NULL, a file path, or a vector/list of sequences")
  if(!is.numeric(kmer) || kmer %% 1 != 0)
    stop("kmer must be an integer")
  if(!is.numeric(iterations) || iterations %% 1 != 0)
    stop("iterations must be an integer")
  if(!is.logical(probabilistic))
    stop("probabilistic must be logical")
  if(!is.logical(normalize))
    stop("normalize must be logical")
  if(!is.numeric(threads) || threads %% 1 != 0)
    stop("threads must be an integer")

  if(probabilistic && !is.null(ctrlfile)) {
    warning("Ignoring ctrlfile argument")
    ctrlfile <- NULL
  }
  if(!probabilistic && is.null(ctrlfile))
    stop("ctrlfile is required when using non-probabilistic ikke")

  if(is.list(testfile) || length(testfile) > 1) {
    testfile <- unlist(testfile, use.names = FALSE)
    if(!is.character(testfile) || anyNA(testfile))
      stop("test sequences must be character strings")
    test_tmp <- tempfile()
    on.exit(unlink(test_tmp), add = TRUE)
    writeLines(testfile, test_tmp)
    testfile <- test_tmp
  } else {
    testfile <- path.expand(testfile)
  }

  if(!is.null(ctrlfile)) {
    if(is.list(ctrlfile) || length(ctrlfile) > 1) {
      ctrlfile <- unlist(ctrlfile, use.names = FALSE)
      if(!is.character(ctrlfile) || anyNA(ctrlfile))
        stop("control sequences must be character strings")
      ctrl_tmp <- tempfile()
      on.exit(unlink(ctrl_tmp), add = TRUE)
      writeLines(ctrlfile, ctrl_tmp)
      ctrlfile <- ctrl_tmp
    } else {
      ctrlfile <- path.expand(ctrlfile)
    }
  }

  .Call("ikke_R",
        testfile,
        ctrlfile,
        as.integer(kmer),
        as.integer(iterations),
        probabilistic,
        normalize,
        as.integer(threads))
}


#' Count k-mer presence in a sequence
#'
#' Count the k-mer presence in a FASTQ, FASTA, raw sequences file, or a vector/list
#' of sequences. k-mer presence the number of sequences in which a k-mer appears.
#' That means that a k-mer is counted at most 1 time per sequnce (no repeats).
#' The final table then represents how many sequences a given k-mer is present in
#' given a dataset.
#'
#' @param input Name of the file which you want to count k-mers presence from, or a
#' character vector/list containing sequences. Files must contain raw sequences,
#' FASTA, or FASTQ data. Gzip-compressed files are supported. A single character
#' string is interpreted as a file path; vectors with multiple elements and
#' lists are interpreted as raw sequences.
#' @param kmer Length of the k-mer you want to count. Currently, only k-mers up
#' to length 15 are supported.
#' @param algo Whether to perform regular counts, or count shuffled sequences.
#' @param bootstrap_iters Number of iterations to bootstrap.
#' @param sample Percent to subsample during bootstrap (should be between 0-100%).
#' @param seed Specify the seed to be used by bootstrap. Since bootstrap
#' subsamples random sequences, seeding alters which random sequences will be
#' picked. This helps to ensure deterministic output which can be achieved by
#' using the same seed. To pick a random seed, set `seed=-1`.
#' @param klet Specify the k-let length to preserve during shuffling. This only
#' affects the output if `algo="shuffled"` is set. -1 chooses the default value.
#' @param sort Sort based on the counts from highest to lowest. By default,
#' the output given is sorted alphabetically on kmers (AA... first, TT... last).
#' @param threads Number of threads to use. Currently unsupported, changing the
#' value will make no difference in performance.
#'
#' @return Dataframe containing the counts for all k-mers.
#' @useDynLib rkats, .registration = TRUE
#' @export
#'
#' @examples
#' # Load data
#' data(rbfox2_seqs)
#'
#' # Count di-nucleotides from a vector of sequences
#' kmer_presence(rbfox2_seqs$bound, kmer = 2)
#'
#' # Count a single sequence by wrapping it in a list, otherwise it will treat
#' # it as a filename
#' kmer_presence(list("ACGTACGT"), kmer = 2)
#'
#' # Create temporary file with sequences
#' tf <- tempfile()
#' writeLines(rbfox2_seqs$bound, tf)
#'
#' # Count presence of di-nucleotides in file
#' kmer_presence(tf, kmer = 2)
#'
#' # Count presence of mono-nucleotides in file
#' kmer_presence(tf, kmer = 1)
#'
#' # Count kmer presence of shuffled sequences
#' kmer_presence(tf, kmer = 1, algo = "shuffled")
#'
#' # Specify k-let to preserve during shuffling
#' kmer_presence(tf, kmer = 1, algo = "shuffled", klet = 2)
#'
#' # Count kmer presence with bootstrap
#' result <- kmer_presence(tf, bootstrap_iters = 100)
#' head(result)
#'
#' # Subsample 55.55% of the file per bootstrap iteration
#' result <- kmer_presence(tf, bootstrap_iters = 100, sample = 55.55)
#' head(result)
#'
#' # Count presence of kmers in shuffled sequences with bootstrap
#' result <- kmer_presence(tf, algo = "shuffled", bootstrap_iters = 100)
#' head(result)
#'
#' # Sort by number of counts, descending
#' result <- kmer_presence(tf, kmer = 5, sort = TRUE)
#' head(result)
#'
#' # Cleanup file
#' unlink(tf)
kmer_presence <- function(input, kmer = 3, algo=c("regular","shuffled"),
                        bootstrap_iters = 0, sample = 25, seed = -1, klet = -1,
                        sort = FALSE, threads = 1) {
  if(!is.character(input) && !is.list(input))
    stop("input must be a file path or a vector/list of sequences")
  if(!is.numeric(kmer) || kmer %% 1 != 0)
    stop("kmer must be an integer")
  if(!is.numeric(bootstrap_iters) || bootstrap_iters %% 1 != 0)
    stop("bootstrap_iters must be an integer")
  if(!is.numeric(sample) || sample <= 0 || 100 < sample)
    stop("sample must be a number between 0-100")
  if(!is.numeric(seed) || seed %% 1 != 0)
    stop("seed must be an integer")
  if(!is.numeric(klet) || klet %% 1 != 0)
    stop("klet must be an integer")
  if(!is.logical(sort))
    stop("sort must be logical")
  if(!is.numeric(threads) || threads %% 1 != 0)
    stop("threads must be an integer")

  if(is.list(input) || length(input) > 1) {
    input <- unlist(input, use.names = FALSE)
    if(!is.character(input) || anyNA(input))
      stop("sequences must be character strings")
    tmp <- tempfile()
    on.exit(unlink(tmp), add = TRUE)
    writeLines(input, tmp)
    input <- tmp
  } else {
    input <- path.expand(input)
  }

  sample <- as.integer((sample*1000) %% 100001)
  algo <- match.arg(algo)
  algo <- if(algo == "regular") 1 else 2

  .Call("presence_R",
        input,
        as.integer(kmer),
        as.integer(klet),
        as.integer(sort),
        as.integer(bootstrap_iters),
        as.integer(sample),
        as.integer(algo),
        as.integer(seed),
        as.integer(threads))
}


#' Search for a sequence within a sequence
#'
#' This function returns the index in which a sub-sequence is found within a
#' larger sequence, or 0 if not found. Searching is case insensitive and 'T' and
#' 'U' are equivalent. The maximum length sequence to search for is currently
#' capped at 255.
#'
#' @param sequence    Big sequence
#' @param search      Small sequence to search for
#' @param all.matches Find all occurrences of search in sequence
#'
#' @return Index in which search was found, 0 if not found
#' @useDynLib rkats, .registration = TRUE
#' @export
#'
#' @examples
#' ## Find "AAA" in a random sequence
#' seq <-paste(sample(c("A","C","G","T"), 1000, TRUE), collapse="")
#' seqseq(seq, "AAA")
#'
#' ## Searching is case-insensitive
#' seq <- "accgtaagggtgccttac"
#' seqseq(seq, "GGGT")
#' seqseq(seq, "gCCtT")
#'
#' ## T and U are both interchangeable
#' seqseq(seq, "GCCUU")
#'
#' ## 0 is returned when search is not in sequence
#' seqseq(seq, "AAA")
seqseq <- function(sequence, search, all.matches = FALSE) {
  if(!is.character(sequence))
    stop("sequence must be a character string")
  if(!is.character(search))
    stop("search must be a character string")
  if(!is.logical(all.matches))
    stop("all.matches must be either TRUE or FALSE")

  return(.Call("seqseq_R", as.character(sequence), as.character(search), all.matches))
}


convert_bytes <- function(total_bytes) {
  # Conversion factors
  bytes_in_gb <- 1024^3
  bytes_in_mb <- 1024^2
  bytes_in_kb <- 1024

  # Calculating GB, MB, and KB
  gb <- floor(total_bytes / bytes_in_gb)
  remainder_after_gb <- total_bytes %% bytes_in_gb

  mb <- floor(remainder_after_gb / bytes_in_mb)
  remainder_after_mb <- remainder_after_gb %% bytes_in_mb

  kb <- ceiling(remainder_after_mb / bytes_in_kb)

  # Prepare output data frame
  result <- data.frame(gb = NA, mb = NA, kb = kb)
  if (gb > 0) {
    result$gb <- gb
  }
  if (mb > 0) {
    result$mb <- mb
  }

  # Creating the warning message with the sizes
  warning_msg <- paste("Dataframe being output will be ",
                       if (!is.na(result$gb)) paste(result$gb, "GB") else NULL,
                       if (!is.na(result$gb) && !is.na(result$mb)) " " else NULL,
                       if (!is.na(result$mb)) paste(result$mb, "MB") else NULL,
                       if ((!is.na(result$gb) || !is.na(result$mb)) && result$kb > 0) " " else NULL,
                       if (result$kb > 0) paste("",result$kb, "KB") else NULL,
                       ".", sep = "")

  return(warning_msg)
}

################################################################################
#                              KATSS SAMPLE DATA                               #
################################################################################

#' List containing the test and control sequences for RNA-binding protein RBFOX2
#'
#' @name rbfox2_seqs
#' @docType data
#' @keywords data
#' @usage data(rbfox2_seqs)
#' @format A list containing two lists: a list of control sequences, and a list
#' of test sequences.
#' \describe{
#'  \item{input}{1000 control sequences of 20 nt long}
#'  \item{bound}{1000 test sequences of 20 nt long}
#' }
"rbfox2_seqs"
