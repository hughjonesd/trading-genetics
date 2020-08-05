
# counterexample where assortative mating makes the correlation more negative

par(mar = rep(2.1, 4))
plot.new()
plot.window(c(-1.2, 1.2), c(-1.2, 1.2))
arrows(x0 = -1.1, x1 = 1.1, y0 = -1.1, length = .1)
arrows(y0 = -1.1, y1 = 1.1, x0 = -1.1, length = .1)
text(-1.1, 1.2, "S", cex = 0.8)
text(1.2, -1.1, "G", cex = 0.8)
polygon(x = c(-1, -.2, 1, .2), y = c(1, 1, -1, -1), col = rgb(.5, 1, .5, .4), border = NA)
polygon(x = c(-.7, -.5, .7, .5), y = c(1, 1, -1, -1), col = rgb(.5, 1, .5, .4), border = NA)

segments(y0 = seq(.9, -.9, -.3), 
         x0 = seq(-.9, by = 1.2/2*.3, length.out = 7),
         x1 = seq(-.2, by = 1.2/2 * .3, length.out = 7), 
         lty = 3, col = rgb(0, 0, 0))


