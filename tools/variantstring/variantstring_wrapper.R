#!/usr/bin/env Rscript

## Galaxy wrapper script for the variantstring R package.
## Reads arguments from the command line and dispatches to the
## appropriate variantstring function.

## Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

## args[1] = operation
## args[2] = input_file path
## args[3] = variant_column (integer)
## args[4] = has_header ("TRUE" or "FALSE")
## args[5] = output_file path
## args[6] = target_string (for compare_* operations, may be "")
## args[7] = position_string (for subset_position, may be "")
## args[8] = second_input_file (for overlay_variant, may be "")
## args[9] = second_variant_column (for overlay_variant, may be "")

op <- args[1]
input_file <- args[2]
variant_column <- as.integer(args[3])
has_header <- as.logical(args[4])
output_file <- args[5]
target_string <- if (length(args) >= 6 && nzchar(args[6])) args[6] else NULL
position_string <- if (length(args) >= 7 && nzchar(args[7])) args[7] else NULL
second_input_file <- if (length(args) >= 8 && nzchar(args[8])) args[8] else NULL
second_variant_column <- if (length(args) >= 9 && nzchar(args[9])) as.integer(args[9]) else 1L

suppressPackageStartupMessages(library(variantstring))

## Read input
df <- read.delim(input_file, header = has_header, stringsAsFactors = FALSE, quote = "")
variant_strings <- as.character(df[, variant_column])

## Dispatch on operation
if (op == "check_variant_string") {
    result <- sapply(variant_strings, function(s) {
        tryCatch(
            {
                check_variant_string(s)
                TRUE
            },
            error = function(e) FALSE
        )
    })
    out <- data.frame(variant_string = variant_strings, valid = result, stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "check_position_string") {
    result <- sapply(variant_strings, function(s) {
        tryCatch(
            {
                check_position_string(s)
                TRUE
            },
            error = function(e) FALSE
        )
    })
    out <- data.frame(position_string = variant_strings, valid = result, stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "variant_to_long") {
    long_list <- variant_to_long(variant_strings)
    out <- do.call(rbind, lapply(seq_along(long_list), function(i) {
        cbind(input_row = i, long_list[[i]], stringsAsFactors = FALSE)
    }))
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "position_from_variant_string") {
    result <- position_from_variant_string(variant_strings)
    out <- data.frame(variant_string = variant_strings, position_string = result, stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "subset_position") {
    result <- subset_position(position_string, variant_strings)
    out <- data.frame(variant_string = variant_strings, subset = result, stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "compare_variant_string") {
    result <- compare_variant_string(target_string, variant_strings)
    out <- cbind(variant_string = variant_strings, result, stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "compare_position_string") {
    result <- compare_position_string(target_string, variant_strings)
    out <- data.frame(variant_string = variant_strings, match = result, stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "get_component_variants") {
    result <- get_component_variants(variant_strings)
    out <- data.frame(
        variant_string = variant_strings,
        component_variants = sapply(result, function(x) paste(x, collapse = ";")),
        stringsAsFactors = FALSE
    )
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "extract_single_locus_variants") {
    result <- extract_single_locus_variants(variant_strings)
    out <- data.frame(
        variant_string = variant_strings,
        single_locus_variants = sapply(result, function(x) paste(x, collapse = ";")),
        stringsAsFactors = FALSE
    )
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "order_variant_string") {
    result <- order_variant_string(variant_strings)
    out <- data.frame(variant_string = variant_strings, ordered = result, stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "drop_read_counts") {
    result <- drop_read_counts(variant_strings)
    out <- data.frame(no_read_counts = result, stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "count_hets") {
    phased <- count_phased_hets(variant_strings)
    unphased <- count_unphased_hets(variant_strings)
    out <- data.frame(
        variant_string = variant_strings,
        phased_hets = phased, unphased_hets = unphased,
        stringsAsFactors = FALSE
    )
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else if (op == "overlay_variant") {
    df2 <- read.delim(second_input_file, header = has_header, stringsAsFactors = FALSE, quote = "")
    var2 <- as.character(df2[, second_variant_column])
    result <- overlay_variant(variant_strings, var2)
    out <- data.frame(var1 = variant_strings, var2 = var2, overlay = result, stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
} else {
    stop("Unknown operation: ", op)
}
