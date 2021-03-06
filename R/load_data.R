

#' @title Read aparcstats or segstats table file created by FreeSurfer tool 'aparcstats2table' or 'asegstats2table'.
#'
#' @param filepath character string, the path to the aparcstats or asegstats file.
#'
#' @author T Schaefer
#'
#' @return data.frame with seg/parc stats
#'
#' @examples
#' \dontrun{
#' stats_file ="~/study1/lh.aparc_allsubjects.table";
#' df = read_segstats_table(stats_file);
#' }
#'
#' @export
read_segstats_table <- function(filepath) {
    return(read.table(filepath, sep='\t', header=TRUE, stringsAsFactors=FALSE));
}


#' @title Retrieve region-based morphometry data from a parcellation stats files created by 'aparcstats2table'.
#'
#' @author T Schaefer
#'
#' @param lh_segstatsfile character string, the path to the aparcstats file for the left hemisphere.
#'
#' @param rh_segstatsfile character string, the path to the aparcstats file for the right hemisphere.
#'
#' @return data.frame, see \code{\link[fsbrain]{group.agg.atlas.native}} for details as this one is designed to match that format.
#'
#' @note The parcellation stats files must be generated manually for your sample using 'aparcstats2table'.
#'
#' @examples
#' \dontrun{
#' lh_sfile ="~/study1/lh.aparc_allsubjects.table";
#' rh_sfile ="~/study1/rh.aparc_allsubjects.table";
#' df = region_data_from_segstat_tables(lh_sfile, rr_sfile);
#' }
#'
#' @export
region_data_from_segstat_tables <- function(lh_segstatsfile, rh_segstatsfile) {
    lh_df = convert.segstats.table.format(read_segstats_table(lh_segstatsfile));
    rh_df = convert.segstats.table.format(read_segstats_table(rh_segstatsfile));
    df_both_hemis = base::merge(lh_df, rh_df, by="subject");
    return(df_both_hemis);
}


#' @title Convert the parcstats data.frame to match the format used by fsbrain::group.agg.atlas.native.
#'
#' @author T Schaefer
#'
#' @param dt data.frame, as returned by \code{read.segstats.table}
#'
#' @return converted data.frame
#'
#' @keywords internal
convert.segstats.table.format <- function(dt) {

    # The column name of the subjects column contains a short table type
    # description, something like 'lh.aparc.area'.
    table_type = colnames(dt)[1];
    table_type_parts = strsplit(table_type, ".", fixed=T)[[1]];
    if(length(table_type_parts) == 3L) {
        hemi = table_type_parts[1];
        atlas = table_type_parts[2];
        measure = table_type_parts[3];
    } else {
        # no usable metadata in subjects column name.
        hemi = atlas = measure = NULL;
        stop("Could not parse measure from parcellation stats table");
    }

    # Rename first column to 'subject'.
    names(dt)[names(dt) == table_type] = 'subject';

    # We need to remove the mean column (named something like "lh_MeanThickness_thickness") if it exists.
    # However, it does not seem to exist in all cases (mabye it is descriptor-specific, and exists only for thickness)?
    dt$lh_MeanThickness_thickness = NULL; # this wont hurt even if the column does not exist.
    dt$rh_MeanThickness_thickness = NULL;

    # Strip measure from column names:
    # Change the column names from something like 'lh_lateraloccipital_thickness' to 'lh_lateraloccipital'.
    for(coln in colnames(dt)) {
        suffix = sprintf("_%s", measure);
        if(endsWith(coln, suffix)) {
            new_coln = substring(coln, 1L, nchar(coln)-nchar(suffix));
            names(dt)[names(dt) == coln] = new_coln;
        }
    }

    return(dt);
}

