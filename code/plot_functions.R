#########################################
# timeseries_plot, created by Scott McKinley
# Takes a matrix of path information and 
# plots them in a standardized style

plot_timeseries = function(t, paths, style = 0, title = "Stochastic Lab Plot") {
  if (style == 0) {
    ymax = max(abs(paths))
    ymin = -ymax
  }
  
  for (i in 1:length(paths[,1])) {
    if (i == 1) {
      plot(t,paths[1,],ylim = c(ymin,ymax),type="l",
           main = title,
           ylab = "x")
    } else
    {
      lines(t,paths[i,])
    }
  }
abline(h=0)
abline(v=0)
}