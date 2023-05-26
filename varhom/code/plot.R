# Plot results from simulations 1 & 2

plot_sim_sd <- function(todo, plot_title, ...){
  df <- tibble(non_norm = as.vector(apply(todo, c(2,3), mean)),
               sampsize = rep(nvec, nv),
               SD = factor(rep(sdvec, each = nn)))
  
  g <- ggplot(data=df, aes(x=sampsize, y=non_norm*100, colour=SD)) + theme_gar +
    geom_abline(slope = 0, intercept = 0.05*100) +
    geom_line(linewidth=1.5) +
    geom_point(fill="white", shape=21, show.legend = FALSE) +
    scale_colour_viridis_d(option = "B", end = 0.9) + 
    coord_cartesian(ylim=c(0, 100)) +
    scale_x_continuous(breaks = nvec[c(TRUE, FALSE)]) +
    scale_y_continuous(breaks = c(0, 5, seq(25,100,25))) +
    labs(x = "Sample size",
         y = "% non-homogeneity decisions") +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          plot.title = element_text(size = 18),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(reverse=TRUE,override.aes = list(size = 2))) +
    ggtitle(plot_title)
  g
}

plot_sim_sd_comp <- function(vartest.bar, vartest.lev.m, vartest.lev.md, vartest.fk, vartest.pb, SD=1, plot_title, ...){
  df <- tibble(non_norm = c(as.vector(apply(vartest.bar[,,sdvec==SD], 2, mean)),
                            as.vector(apply(vartest.lev.m[,,sdvec==SD], 2, mean)),
                            as.vector(apply(vartest.lev.md[,,sdvec==SD], 2, mean)),
                            as.vector(apply(vartest.fk[,,sdvec==SD], 2, mean)),
                            as.vector(apply(vartest.pb[,,sdvec==SD], 2, mean))),
               sampsize = rep(nvec, 5),
               method = factor(rep(c("Bartlett", "Levene", "Brown-Forsythe", "Fligner-Killeen", "Perc. boot."), each = nn)))

  g <- ggplot(data=df, aes(x=sampsize, y=non_norm*100, colour=method)) + theme_gar +
    geom_abline(slope = 0, intercept = 0.05*100) +
    geom_line(linewidth=1.5) +
    geom_point(fill="white", shape=21, show.legend = FALSE) +
    scale_color_manual(values=c("#111111", "#56B4E9", "#E69F00", "#CC79A7", "#0072B2")) +
    # coord_cartesian(ylim=c(0, 100)) +
    scale_x_continuous(breaks = nvec[c(TRUE, FALSE)]) +
    scale_y_continuous(breaks = c(0, 5, seq(25,100,25))) +
    labs(x = "Sample size",
         y = "% non-homogeneity decisions") +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          plot.title = element_text(size = 18),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(reverse=FALSE,override.aes = list(linewidth = 3))) +
    ggtitle(plot_title)
  g
}

# Plot results from simulations 3 & 4
plot_sim_g <- function(todo, plot_title, ...){
  df <- tibble(non_norm = as.vector(apply(todo, c(2,3), mean)),
               sampsize = rep(nvec, nv),
               g = factor(rep(gvec, each = nn)))
  
g <- ggplot(data=df, aes(x=sampsize, y=non_norm*100, colour=g)) + theme_gar +
    geom_abline(slope = 0, intercept = 0.05*100) +
    geom_line(linewidth=1.5) +
    geom_point(fill="white", shape=21, show.legend = FALSE) +
    scale_colour_viridis_d(option = "B", end = 0.9) + 
    scale_x_continuous(breaks = nvec[c(TRUE, FALSE)]) +
    scale_y_continuous(breaks = c(0, 5, seq(25,100,25))) +
    labs(x = "Sample size") +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          plot.title = element_text(size = 18),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(reverse=TRUE,override.aes = list(linewidth = 3))) +
    ggtitle(plot_title)
  g
}

plot_sim_g_comp <- function(vartest.bar, vartest.lev.m, vartest.lev.md, vartest.fk, vartest.pb, G=0, plot_title, ...){
  df <- tibble(non_norm = c(as.vector(apply(vartest.bar[,,gvec==G], 2, mean)),
                            as.vector(apply(vartest.lev.m[,,gvec==G], 2, mean)),
                            as.vector(apply(vartest.lev.md[,,gvec==G], 2, mean)),
                            as.vector(apply(vartest.fk[,,gvec==G], 2, mean)),
                            as.vector(apply(vartest.pb[,,gvec==G], 2, mean))),
               sampsize = rep(nvec, 5),
               method = factor(rep(c("Bartlett", "Levene", "Brown-Forsythe", "Fligner-Killeen", "Perc. boot."), each = nn)))
  
  g <- ggplot(data=df, aes(x=sampsize, y=non_norm*100, colour=method)) + theme_gar +
    geom_abline(slope = 0, intercept = 0.05*100) +
    geom_line(linewidth=1.5) +
    geom_point(fill="white", shape=21, show.legend = FALSE) +
    scale_color_manual(values=c("#111111", "#56B4E9", "#E69F00", "#CC79A7", "#0072B2")) +
    scale_x_continuous(breaks = nvec[c(TRUE, FALSE)]) +
    scale_y_continuous(breaks = c(0, 5, seq(25,100,25))) +
    labs(x = "Sample size") +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          plot.title = element_text(size = 18),
          panel.grid.minor.y = element_blank(),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 12)) +
    guides(colour = guide_legend(reverse=FALSE,override.aes = list(linewidth = 3))) +
    ggtitle(plot_title)
  g
}

