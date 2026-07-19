#!/usr/bin/env Rscript

## Wrapper script for the variantstring R package
## Exposes key functions as command-line operations on tabular input

library(variantstring)

args <- commandArgs(trailingOnly = TRUE)

op <- args[1]
input_file <- args[2]
output_file <- args[3]
variant_column <- as.integer(args[4])
has_header <- toupper(args[5]) == "TRUE"
target_string <- args[6]
position_string <- args[7]
second_input_file <- args[8]
second_variant_column <- as.integer(args[9])

## Read input
df <- read.delim(input_file, header = has_header, stringsAsFactors = FALSE, quote = "")
variant_strings <- as.character(df[, variant_column])

if (op == "check_variant_string") {
    result <- sapply(variant_strings, function(s) {
        tryCatch({
            check_variant_string(s)
            TRUE
        }, error = function(e) FALSE)
    })
    out <- data.frame(variant_string = variant_strings, valid = result,
                      stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else if (op == "check_position_string") {
    result <- sapply(variant_strings, function(s) {
        tryCatch({
            check_position_string(s)
            TRUE
        }, error = function(e) FALSE)
    })
    out <- data.frame(position_string = variant_strings, valid = result,
                      stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else if (op == "variant_to_long") {
    long_list <- variant_to_long(variant_strings)
    out <- do.call(rbind, lapply(seq_along(long_list), function(i) {
        cbind(input_row = i, long_list[[i]])
    }))
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else if (op == "position_from_variant_string") {
    result <- position_from_variant_string(variant_strings)
    out <- data.frame(variant_string = variant_strings, position_string = result,
                      stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else if (op == "subset_position") {
    result <- subset_position(position_string, variant_strings)
    out <- data.frame(variant_string = variant_strings, subset = result,
                      stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else if (op == "compare_variant_string") {
    result <- compare_variant_string(target_string, variant_strings)
    out <- cbind(variant_string = variant_strings, result)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else if (op == "compare_position_string") {
    result <- compare_position_string(target_string, variant_strings)
    out <- data.frame(variant_string = variant_strings, match = result,
                      stringsAsFactors = FALSE)
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
    out <- data.frame(variant_string = variant_strings, ordered = result,
                      stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else if (op == "drop_read_counts") {
    result <- drop_read_counts(variant_strings)
    out <- data.frame(variant_string = variant_strings, no_read_counts = result,
                      stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else if (op == "count_hets") {
    phased <- count_phased_hets(variant_strings)
    unphased <- count_unphased_hets(variant_strings)
    out <- data.frame(
        variant_string = variant_strings,
        phased_hets = phased,
        unphased_hets = unphased,
        stringsAsFactors = FALSE
    )
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else if (op == "overlay_variant") {
    df2 <- read.delim(second_input_file, header = has_header,
                      stringsAsFactors = FALSE, quote = "")
    var2 <- as.character(df2[, second_variant_column])
    result <- overlay_variant(variant_strings, var2)
    out <- data.frame(var1 = variant_strings, var2 = var2, overlay = result,
                      stringsAsFactors = FALSE)
    write.table(out, output_file, sep = "\t", row.names = FALSE, quote = FALSE)

} else {
    stop("Unknown operation: ", op)
}
