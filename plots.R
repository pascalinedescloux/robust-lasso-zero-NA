### Producing plots of simulation results

## -- load pkgs
library(abind)
library(reshape2)
library(ggplot2)
library(dplyr)

simu.type <- "soracle" # "soracle" or "automatic"
rho <- 0 # fixed correlation parameter (single value)
s.values <- c(3, 10) # all sparsity indices to plot
pitot.values <- c(0.05, 0.2) # all pi.tot values (proportion of missing entries) to plot
a.values <- c(0, 5) # all a.values to plot

for (s in s.values) {
    for (pitot in pitot.values) {
        for(a  in a.values) {
            load(paste0(simu.type, ".s", s, ".pi", pitot, ".a", a, ".rho", rho, ".Rda"))
            assign(paste0(simu.type, ".s", s, ".pi", pitot, ".a", a, ".rho", rho),
                   get(paste0(simu.type, ".s", s, ".pi", pitot, ".a", a, ".rho", rho))$res)
        }
        assign(paste0("res.s", s, ".pi", pitot), 
               abind(mget(paste0(simu.type, ".s", s, ".pi", pitot, ".a", a.values, ".rho", rho)), along = 0))
    }
    assign(paste0("res.s", s),
           abind(mget(paste0("res.s", s, ".pi", pitot.values)), along = 0))
}
assign("results",
       abind(mget(paste0("res.s", s.values)), along = 0))
dimnames(results) <- list(sparsity = s.values,
                          pi.tot = pitot.values,
                          a = a.values,
                          estimator = dimnames(results)[[4]],
                          criterion = dimnames(results)[[5]],
                          rep = dimnames(results)[[6]])

## -- long format:
Lresults <- melt(results, value.name = "value")
# renaming facet labels:
#Lresults$sparsity <- factor(paste0("s=", Lresults$sparsity), 
#                            levels = c("s=3", "s=10"))
#Lresults$pi.tot <- factor(paste0(100* Lresults$pi.tot, "% NA"),
#                          levels = c("5% NA", "20% NA"))
# renaming criteria
levels(Lresults$criterion) <- c("proba of supp recovery", "TPR", "FDR",
                                "proba of sign recovery", "s-TPR", "s-FDR")

## -- plot for each criterion:
mean.narm <- function(x) {
    mean(x, na.rm = TRUE)
}
grouped <- group_by(Lresults, sparsity, pi.tot, a, estimator, criterion)
Lmeanres <- summarise(grouped, mean=mean.narm(value), sd=sd(value))

for(crit in levels(Lresults$criterion)) {
    g <- ggplot(data = Lmeanres[Lmeanres$criterion == crit, ],
                mapping = aes(x = a, y = mean, colour = estimator)) + 
        geom_line(size = 0.75) +
        geom_point(aes(shape = estimator), size = 2.5, stroke = 0.75) +
        facet_grid(pi.tot ~ sparsity) + 
        scale_shape_manual(values = c(1, 
                                      13, 
                                      9, 
                                      8)) +
        scale_colour_manual(values = c("dimgray", 
                                       "royalblue", 
                                       "chocolate", 
                                       "forestgreen")) +
        ylab(crit) +
        xlab("") +
        ylim(0, 1) +
        theme(plot.title = element_text(hjust = 0.5, size = 20, face="bold"),
              #axis.title.x = element_text(size=20),
              axis.text.x = element_text(hjust=c(0, 0, 0, 0, 0, 1), size = 15),
              axis.text.y = element_text(size = 12),
              axis.ticks = element_blank(),
              axis.title.y = element_text(size=20),
              legend.title = element_blank(),
              legend.text = element_text(size = 16),
              legend.position = "bottom",
              strip.text.y = element_text(size = 20),
              strip.text.x = element_text(size = 20))  +
        scale_x_continuous(labels=c("0" = "MCAR\n(a=0)", 
                                    "1" = "", "2" = "", "3" = "", "4" = "",
                                    "5" = "MNAR\n(a=5)"))
    print(g)
}



