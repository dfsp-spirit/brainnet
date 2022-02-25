
#' @title Create effect size violin plots for all predictors at once.
#'
#' @param effects_mat n*m matrix in which the n rownames are the predictors (like age, sex, site, ...) and the m columns contain the effect sizes for the m vertices.
#'
#' @param plot_ylabel character string, the y label text for the plot.
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot coord_flip geom_violin geom_boxplot theme_light theme element_text ylab scale_x_discrete scale_fill_brewer scale_color_brewer
#'
#' @return a ggplot instance.
#'
#' @author Johanna Leyhausen
#'
#' @export
effect_size_violin_plots <- function(effects_mat, plot_ylabel="Effect Size (Cohen's f)") {
    tdata <- t(as.data.frame(effects_mat)); #transpose data
    mdata <- reshape2::melt(tdata, id=rownames(effects_mat)); # Reshape data to long format for usewith  ggplot.
    # Afterwards there are 3 columns in 'mdata': 'Var1' contains vertex name (like 'V1...Vn'), 'Var2' contains predictor name (Sex, age, ...) and 'value' contains the descriptor value for the respective vertex.

    #NOTE!: scale_x_discrete overwrites the axis labels with what you define in labs, make sure to have them in the order of your data
    #NOTE!: If you have more than 8 different things to plot, you have to use another color brewer palette (e.g. "Paired)
    labs <- rownames(effects_mat) #labels for axis
    p <- ggplot2::ggplot(mdata, ggplot2::aes(x= Var2, y= value, fill= Var2, alpha = 0.5, color = Var2)) +
        ggplot2::coord_flip() +
        ggplot2::geom_violin(width = 1.2) +
        ggplot2::geom_boxplot(width = 0.08) +
        ggplot2::theme_light() +
        ggplot2::theme(axis.title.y=ggplot2::element_blank(), legend.position = "none",
                       axis.text.y = ggplot2::element_text(color = "grey30", size = 24),
                       axis.text.x = ggplot2::element_text(color = "grey30", size = 24),
                       axis.title=ggplot2::element_text(size=26,face="bold")) +
        ggplot2::ylab(plot_ylabel) +
        ggplot2::scale_x_discrete(labels = labs) +
        ggplot2::scale_fill_brewer(palette = "Dark2") +
        ggplot2::scale_color_brewer(palette = "Dark2")
    return(p)
}
