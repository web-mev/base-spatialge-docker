prep_stlist <- function(counts_file, coords_file, sample_name) {
    # STlist expects that the first column is the gene names- so we don't use row.names arg
    rnacounts <- read.table(counts_file, sep='\t', header=T, check.names=F)

    # Same as for the counts, the expectation is that the coordinates file does not
    # have row.names and instead has the barcodes in the first column. For now, however,
    # we set the row names and then later alter.
    spotcoords <- read.table(coords_file, sep='\t', row.names=1, header=T, check.names=T)

    # only take the first two columns for the (x,y) positions. Additional columns
    # can cause problems downstream
    spotcoords <- spotcoords[,c(1,2)]

    # the barcodes in coords dataframe can be a superset of the count matrix columns.
    # For example, if the matrix is filtered for QC, there may be poor quality spots
    # that were filtered out. 
    # The opposite is not the case since we cannot have a barcode without a position.
    barcodes_from_counts <- colnames(rnacounts)[2:dim(rnacounts)[2]] 
    diff_set <- setdiff(barcodes_from_counts, rownames(spotcoords))
    num_invalid_barcodes <- length(diff_set)
    if (num_invalid_barcodes > 0) {
        max_print <- 5
        if (num_invalid_barcodes < max_print) {
            invalid_barcodes <- paste(diff_set[1:num_invalid_barcodes], collapse=', ')
        } else {
            invalid_barcodes <- sprintf('%s, and %d others.', paste(diff_set[1:max_print], collapse=', '), num_invalid_barcodes-max_print)
        }
        message(sprintf('The set of barcodes in your count matrix must be a subset of those in your coordinates file. Problems include: %s', invalid_barcodes))
        quit(status=1)
    } else {
        # the STList constructor below will not accept coordinate files which are a superset of the 
        # count matrix barcodes. We handle that here:
        spotcoords <- spotcoords[barcodes_from_counts,]

        # we need the barcodes in the first col, not the row names. We used the rownames for indexing
        # convenience above, but need to change that now:
        spotcoords <- cbind(rownames(spotcoords), data.frame(spotcoords, row.names=NULL))
    }

    # Now, to avoid any unexpected issues downstream, we need to convert the column names, etc.
    # to preserve the barcodes/column names, we create a dataframe of the original and 'R mutated'
    # names. We then run through everything with the mutated names and finally map back.
    orig_col_names <- colnames(rnacounts)
    proper_names <- make.names(orig_col_names)
    colname_mapping = data.frame(
        orig_names = orig_col_names,
        row.names=proper_names,
        stringsAsFactors=F)
    colnames(rnacounts) <- proper_names


    spotcoords[,1] <- make.names(spotcoords[,1]) 

    # We will use a list of dataframes in the call to STlist
    rnacounts_list <- list()
    rnacounts_list[[sample_name]] <- rnacounts
    spotcoords_list <- list()
    spotcoords_list[[sample_name]] <- spotcoords

    # Import input data into spatialGE object
    spat <- STlist(
        rnacounts=rnacounts_list,
        spotcoords=spotcoords_list, 
        samples=c(sample_name)
    )
    return(spat)
}