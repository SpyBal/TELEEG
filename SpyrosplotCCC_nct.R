SpyrosplotCCC.NCT = function (x, what = c("strength", "network", "edge", "centrality"), 
          ...) 
{
  what <- match.arg(what)
  if (what == "strength") {
    # hist(x$glstrinv.perm, main = paste("p =", x$glstrinv.pval), 
    #      xlab = "Difference in global strength", xlim = c(0, 
    #                                                       max(x$glstrinv.real, x$glstrinv.perm)))
    df= data.frame(stat=x$glstrinv.perm,obs = x$glstrinv.real)
    p=ggplot(df, aes(x=stat)) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      ylab("Density")+xlab("Difference in global strength")+
      geom_text(aes(label=paste("p-val =", x$glstrinv.pval), x=max(c(x$glstrinv.real, x$glstrinv.perm))- sd(c(x$glstrinv.real, x$glstrinv.perm)),y=max(density(df$stat)$y)),color="#00aedb", fontface="bold",size=7)+
      geom_density(alpha=.2, fill="#ffc425") + geom_point(aes(x=obs,y=0), size=10,shape=18,color="tomato")+
      xlim(c(0,max(x$glstrinv.real, x$glstrinv.perm)))+
      ggtitle("Network global strength comparison\n TRAILS - CCC2000")+
      theme_bw() +
      theme(axis.title = element_text(face="bold",size=15),
            plot.title = element_text(face="bold",size=20, hjust = 0.5),
            axis.text = element_text(face="italic",size=10))

    return(p)
    
    
    
  }
  if (what == "network") {
    # hist(x$nwinv.perm, main = paste("p =", x$nwinv.pval), 
    #      xlab = "Maximum of difference", xlim = c(0, max(x$nwinv.real, 
    #                                                      x$nwinv.perm)))
    df= data.frame(stat=x$nwinv.perm,obs = x$nwinv.real)
    p = ggplot(df, aes(x=stat)) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+ylab("Density")+xlab("Difference in global structure")+
      geom_text(aes(label=paste("p-val =", x$nwinv.pval), x=max(c(x$nwinv.real, x$nwinv.perm))- sd(c(x$nwinv.real, x$nwinv.perm)),y=max(density(df$stat)$y)),color="#00aedb", fontface="bold",size=7)+
      geom_density(alpha=.2, fill="#ffc425") + geom_point(aes(x=obs,y=0), size=10,shape=18,color="tomato")+
      xlim(c(0,max(x$nwinv.real, x$nwinv.perm)))+
      ggtitle("Network global structure comparison\n TRAILS - CCC2000")+
      theme_bw() +
      theme(axis.title = element_text(face="bold",size=15),
            plot.title = element_text(face="bold",size=20, hjust = 0.5),
            axis.text = element_text(face="italic",size=10))

    return(p)

    
  }
  if (what == "edge") {
    if (length(dim(x$einv.perm)) > 2) {
      extractLowTri <- function(einv.perm) {
        grid <- expand.grid(1:dim(einv.perm)[2], 1:dim(einv.perm)[2])
        grid <- grid[grid[, 1] != grid[, 2], ]
        grid <- grid[1:(nrow(grid)/2), ]
        out <- matrix(NA, nrow = dim(einv.perm)[3], 
                      ncol = nrow(grid))
        colnames(out) <- 1:nrow(grid)
        for (i in 1:nrow(grid)) {
          out[, i] <- einv.perm[grid[i, 1], grid[i, 
                                                 2], ]
          colnames(out)[i] <- paste(grid[i, 1], grid[i, 
                                                     2], sep = "-")
        }
        return(out)
      }
      x$einv.perm <- extractLowTri(x$einv.perm)
    }
    nedgetests <- ncol(x$einv.perm)
    for (i in 1:nedgetests) {
      hist(x$einv.perm[, i], main = paste("p =", x$einv.pval[i, 
                                                             3]), xlab = paste("Difference in edge strength:", 
                                                                               colnames(x$einv.perm)[i]), xlim = c(0, max(x$einv.real, 
                                                                                                                          x$einv.perm) * 1.1))
      points(x$einv.real[i], 0, col = "red", pch = 17)
    }
  }
  if (what == "centrality") {
    ncentests <- ncol(x$diffcen.perm)
    for (i in 1:ncentests) {
      hist(x$diffcen.perm[, i], main = paste("p =", reshape2::melt(x$diffcen.pval)$value[i]), 
           xlab = paste("Difference in '", reshape2::melt(x$diffcen.pval)$Var2[i], 
                        "' for node '", reshape2::melt(x$diffcen.pval)$Var1[i], 
                        "'", sep = ""), xlim = c(0, max(x$diffcen.perm[, 
                                                                       i], abs(reshape2::melt(x$diffcen.real)$value[i]))))
      points(abs(reshape2::melt(x$diffcen.real)$value[i]), 
             0, col = "red", pch = 17)
    }
  }
}
