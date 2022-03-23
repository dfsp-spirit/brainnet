

#' @title Plot all columns of a data.frame by group.
#'
#' @description Creates a plot with one subplot per data.frame column. Continuous data is plotted as histograms, factors as bar plots.
#'
#' @param df data.frame, the input data. Can contain factor columns and continuous columns.
#'
#' @param group_by character string, the column name that we should group by.
#'
#' @importFrom ggplot2 ggplot aes_string geom_boxplot
#'
#' @examples
#' \dontrun{
#' md = load_ABIDE_metadata_males();
#' dtf = md$demographics;
#' df_group_plot(dtf);
#' }
#'
#' @keywords internal
df_group_plot <- function(dtf, group_column="group") {

    if(requireNamespace("ggpubr", quietly = TRUE)) {

        plot_max_factor_levels = 20L; # Factors with more levels will not be plotted. Useful to skip columns like the subject ID, and also it makes little sense anyways.

        if(! is.data.frame(dtf)) {
            stop("Parameter 'df' must be a data.frame");
        }

        if(!(group_column %in% names(dtf))) {
            stop(sprintf("The data.frame df does not contain a column named '%s', please fix the 'group_column' parameter.", group_column));
        }

        num_dataframe_columns = length(names(dtf));
        num_data_columns = 0; # Gets incremented below
        #par(mfrow=c(num_plot_rows, num_plot_columns),mar=c(2,1,1,1));

        plot_list = list();

        for(i in 1:num_dataframe_columns) {
            current_col_name = names(dtf)[i];
            if(current_col_name != group_column) {
                if(is.character(dtf[,i])) {
                    dtf[,i] = as.factor(dtf[,i]);
                }
                if(is.factor(dtf[,i])) {
                    num_levels = length(levels(dtf[,i]));
                    if(num_levels < 2L | num_levels > plot_max_factor_levels) {
                        # TODO: Skipping at this time may lead to an empty row of subplots, we should adjust earlier.
                        if(num_levels < 2L) {
                            message(sprintf("Column #%d '%s' is a factor with %d levels. Need at least 2 for a meaningful bar plot. Skipping.\n", i, current_col_name, num_levels));
                        }
                        if(num_levels > plot_max_factor_levels) {
                            warning(sprintf("Column #%d '%s' is a factor with %d levels. Plotting makes little sense for more than %d. Skipping.\n", i, current_col_name, num_levels, plot_max_factor_levels));
                        }
                    } else {
                        #plot(dtf[,i], main = current_col_name);
                        num_data_columns = num_data_columns + 1L;
                        message(sprintf("Plotting categorical column #%d '%s'.\n", i, current_col_name));
                        e <- ggplot2::ggplot(dtf, ggplot2::aes_string(x=group_column, y=current_col_name, colour= group_column, group=group_column ));
                        e <- e + ggplot2::geom_boxplot();
                        plot_list[[current_col_name]] = e;
                    }

                } else{
                    if(is.numeric(dtf[,i])) {
                        message(sprintf("Plotting numerical column #%d '%s'\n", i, current_col_name));
                        #num_data_columns = num_data_columns + 1L;
                        #hist(dtf[,i], main = current_col_name);
                    } else {
                        # TODO: Skipping at this time may lead to an empty row of subplots, we should adjust earlier.
                        warning(sprintf("Column #%d '%s' is not numeric and not a factor. Skipping.\n", i, current_col_name));
                    }
                }
            } else {
                message(sprintf("Ingoring the group column #%d '%s'.\n", i, current_col_name));
            }
        }
        num_plot_columns = min(3L, num_data_columns); # We use 3 columns (or less, if there are less than 3 plots to make).
        num_plot_rows = ceiling(num_data_columns/num_plot_columns);

        cat(sprintf("Plotting %d subplots in %d rows and %d columns. There are %d plots in the plot_list.\n", num_data_columns, num_plot_rows, num_plot_columns, length(plot_list)));

        figure <- ggpubr::ggarrange(plot_list = unname(plot_list), labels = names(plot_list), ncol = num_plot_columns, nrow = num_plot_rows);
        return(list("fig"=figure, "plot_list"=plot_list));
    } else {
        stop("The optional dependency package 'ggpubr' is required to use this functionality. Please install the 'ggpubr' package.");
    }
}

