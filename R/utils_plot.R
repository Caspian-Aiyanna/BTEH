# Plotting helpers shared across analyses

plot_kendall_heatmap <- function(cmat, cutoff, output_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not installed; skipping Kendall heatmap plot")
    return(invisible(NULL))
  }
  keep_vars <- colnames(cmat)
  kt_long <- as.data.frame(as.table(cmat))
  names(kt_long) <- c("Var1", "Var2", "tau")
  g <- ggplot2::ggplot(kt_long, ggplot2::aes(Var1, Var2, fill = tau)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Kendall τ",
                  title = sprintf("Kendall correlation (|τ| < %.2f)", cutoff)) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggplot2::ggsave(output_path, g, width = 8, height = 6, dpi = 300)
  invisible(output_path)
}

