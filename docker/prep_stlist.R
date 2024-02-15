library(tibble)
library(dplyr)
library(spatialGE)


prep_stlist <- function(counts_file, coords_file, sample_name, gene_mapping_df, source_gene_id, target_gene_id) {
    
    # This function returns a list containing the column-name mapping (so we can eventually map back to the
    # original barcodes) and an STList instance.
    
    # Note that the `gene_mapping_df`, `source_gene_id`, and `target_gene_id` are optional and 
    # are used to map from other identifier systems for situations when we need to compare
    # against a curated database like MSigDB, etc. If they are not supplied, no mapping
    # is performed.

    rnacounts <- read.table(counts_file, sep='\t', header=T, check.names=F, row.names=1)
    barcodes_from_counts <- colnames(rnacounts)[1:dim(rnacounts)[2]]
    # if the gene_mapping_df was passed:
    if (!missing(gene_mapping_df)) {

        rnacounts = merge(
                        rnacounts,
                        gene_mapping_df[,c(source_gene_id, target_gene_id)],
                        by.x=0,
                        by.y=source_gene_id
                    )
        # If no remaining rows, error out- the join didn't produce any results so it's most likely that the
        # gene identifier was not correct
        if(dim(rnacounts)[1] == 0){
            message('After mapping the gene identifiers, there were no remaining rows. Was the choice of gene identifier correct?')
            quit(status=1)
        }

        # remove the duplicate rows which can result:
        rnacounts <- rnacounts %>% distinct()

        # in the case of a trivial mapping (e.g. if source and target are BOTH "SYMBOL")
        # we get col names like 'Row.names' and 'SYMBOL.1'. This handles that edge case.
        if (target_gene_id == source_gene_id){
            gene_id_col <- 'Row.names'
        } else {
            gene_id_col <- target_gene_id
        }
        # reorder- we need the gene identifier as the first column.
        ordering <- c(gene_id_col, barcodes_from_counts)
    } else {
        # no remapping needed- simply take the row names as the gene ID.
        # STList expects the first column to be the gene, NOT the rownames.
        # This function is not explicit that the rownames become the first column
        # so we set the ordering and sort just to be certain.
        rnacounts <- tibble::rownames_to_column(rnacounts, var='__gene__')
        ordering <- c('__gene__', barcodes_from_counts)
    }

    # now sort the columns such that the first column is the gene ID:
    rnacounts <- rnacounts[,ordering]

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
    return(list(spat=spat, colname_mapping=colname_mapping))
}