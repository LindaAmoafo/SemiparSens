## Individual Curves
cbPalette_RP <- c("red", "#E69F00", "#009E73",  "#56B4E9", "#CC79A7" ,"grey")

Estimate.plot.data <- bind_rows(All.G_pava_EBRT_plot, All.G_pava_RP_plot, 
                                ref_dat) %>% 
  mutate(`Treatment` = factor(`Treatment`, levels= c("RP", "EBRT + AD")),
         Estimates = factor(Estimates, levels = c("KM Estimates","Gamma = 0", "Gamma = 0.5", "Gamma = 1", "Gamma = 1.5", "High Risk", "Gamma = -0.5", "Gamma = -1", "Gamma = -1.5", "Gamma =-2","Gamma =-2.5", "Low Risk")))

Estimate.plot.data_combine_split <- split(Estimate.plot.data, 
                                          f= Estimate.plot.data$`Treatment`)

Estimate.plot1 <- ggplot(Estimate.plot.data_combine_split$`RP` , 
                         aes(x= Time, y=`Survival Probability`, color=Estimates,
                             linetype = Estimates, size=Estimates)) +
  geom_step()+ 
  ylim(0, 1) +
  scale_x_continuous(breaks=seq(0,120, 24), limits = c(0, 130)) +
  theme_bw()+
  facet_grid(.~`Treatment`, scales = "fixed")+
  labs(x = "Time", y = "Survival Probability", color = "Estimate") +
  theme(legend.position = c(0.2, 0.25))+
  guides(colour = guide_legend(override.aes = list(linetype = c(1,1,1,1,1,5)))) + 
  scale_linetype_manual(values= c(rep("solid", 5), rep("longdash", 1)), guide = FALSE)+
  scale_size_manual(values = c(1,1,1,1,1 ,0.8), guide=FALSE)+
  scale_color_manual(values = cbPalette_RP, guide=FALSE, labels = c("KM Estimate", expression(paste(gamma[1], " = 0")), expression(paste(gamma[1], " = 0.5")),expression(paste(gamma[1], " = 1")),expression(paste(gamma[1], " = 1.5")), "High Risk Survival"))+
  theme(legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face= "bold"), #change legend title font size
        legend.text.align = 0,
        legend.text = element_text(size=12), #change legend text font size
        axis.text.x = element_text(size = 18, face="bold"),
        axis.title.x = element_text(size = 20, face="bold"),
        axis.text.y = element_text(size = 18, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.text.x = element_text(size = 18, face="bold"))


cbPalette_EBRT <- c("red", "#E69F00", "#009E73",  "#56B4E9", "#CC79A7","cyan", "purple" ,"grey")
plot2 <- Estimate.plot1  %+% Estimate.plot.data_combine_split$`EBRT + AD` + labs(y=NULL)+
  guides(colour = guide_legend(override.aes = list(linetype = c(1,1,1,1,1,1,1,5)))) + 
  scale_linetype_manual(values= c(rep("solid", 7), rep("longdash", 1)), guide = FALSE)+
  scale_size_manual(values = c(1,1,1,1,1,1,1,0.8), guide=FALSE)+
  scale_color_manual(values = cbPalette_EBRT, guide=FALSE, labels = c("KM Estimate",expression(paste(gamma[0], " = 0")), expression(paste(gamma[0], " = -0.5")),expression(paste(gamma[0], " = -1")),expression(paste(gamma[0], " = -1.5")),expression(paste(gamma[0], " = -2")),expression(paste(gamma[0], " = -2.5")), "Low Risk Survival"))



png("Survival Curve without Low Risk.png", width = 15, height = 10, units = 'in',res = 300)

grid.arrange(Estimate.plot1,plot2,ncol=2)

dev.off()


## Induced Counterfactual Plots
Counter.plot.data <- Counter.plot.dat %>% 
  mutate(`Treatment` = factor(`Treatment`, levels= c("RP", "EBRT + AD")),
         Estimates = factor(Estimates, levels = c("P[Y(t) \u2265 s | T=t]","Gamma = 0", "Gamma = 0.5", "Gamma = 1", "Gamma = 1.5","Gamma = -0.5", "Gamma = -1", "Gamma = -1.5", "Gamma = -2","Gamma = -2.5")))


Counter.plot.data_split <- split(Counter.plot.data, f= Counter.plot.data$Treatment)

Counter.plot1 <- ggplot(data=Counter.plot.data_split$`RP` , 
                        aes(x= Time, y=`Survival Probability`, color=Estimates,
                            linetype = Estimates, size=Estimates)) +
  geom_step() +
  ylim(0, 1) +
  #xlim(0,130) +
  scale_x_continuous(breaks=seq(0,120, 24), limits = c(0, 130)) +
  facet_wrap(~Treatment, scales = "fixed")+
  labs(x = "Time (in months)", y = "Survival Probability", color = "Estimate") +
  theme_bw()+
  theme(legend.position = c(0.35,0.32),
        legend.box = "horizontal") +
  guides(colour = guide_legend(override.aes = list(linetype = c(1, 2,2,2,2))))+ 
  scale_linetype_manual(values= c("solid", rep("dashed", 4)), guide = FALSE)+
  scale_size_manual(values = c(0.5,0.8,0.8,0.8,0.8), guide=FALSE)+
  scale_color_manual(values = c("red", "#045a8d","#2b8cbe", "#74a9cf", "#bdc9e1"), guide=FALSE, labels = c("P[Y(t) > s | T=t]",expression(paste(gamma[1], " = 0")), expression(paste(gamma[1], " = 0.5")),expression(paste(gamma[1], " = 1")),expression(paste(gamma[1], " = 1.5"))))+
  theme(legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=15, face= "bold"), #change legend title font size
        legend.text.align = 0,
        legend.text = element_text(size=12), #change legend text font size
        axis.text.x = element_text(size = 16, face="bold"),
        axis.title.x = element_text(size = 18, face="bold"),
        axis.text.y = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 18, face="bold"),
        strip.text.x = element_text(size = 16, face="bold"))


Counter.plot2 <- Counter.plot1  %+% Counter.plot.data_split$`EBRT + AD`+ labs(y=NULL) +
  guides(colour = guide_legend(override.aes = list(linetype = c(1,2,2,2,2,2,2)))) + 
  scale_linetype_manual(values= c(rep("solid", 1), rep("dashed", 6)), guide = FALSE)+
  scale_size_manual(values = c(0.5,0.8,0.8,0.8,0.8,0.8,0.8), guide=FALSE)+
  scale_color_manual(values = c("red", "#045a8d","#2b8cbe", "#74a9cf", "#bdc9e1", "#d0d1e6", "#f1eef6"), guide=FALSE, labels = c("P[Y(t) > s | T=t]",expression(paste(gamma[1], " = 0")), expression(paste(gamma[1], " = -0.5")),expression(paste(gamma[1], " = -1")),expression(paste(gamma[1], " = -1.5")), expression(paste(gamma[1], " = -2")), expression(paste(gamma[1], " = -2.5"))))


png("Counter Plot.png", width = 10, height = 7, units = 'in',res = 600)

grid.arrange(Counter.plot1,Counter.plot2, ncol=2)

dev.off()


## Contour plots
png("Diff Contour plot without Low Risk.png", width = 15, height = 10, units = 'in',res = 600)

contour(gamma.grid.Neg, gamma.grid, diff, nlevels = 12 , xlab = substitute(paste(bold(gamma[0]))), ylab = substitute(paste(bold(gamma[1]))), lty="dashed", lwd=1, cex=0.01, col="grey", cex.axis = 1.5, cex.lab = 1.5, axes = F)

axis(1, at = seq(0, 2.5,0.5), labels = -seq(0, 2.5,0.5))

# Add custom y-axis
axis(2)

contour(gamma.grid.Neg, gamma.grid, diff, levels = 0, lwd = 3, add = T, col = "red")

contour(gamma.grid.Neg, gamma.grid, Lower, lwd=3,label="", levels = 0, add = T, col = "blue")

contour(gamma.grid.Neg, gamma.grid, Upper, lwd=3,label="", levels = 0, add = T, col = "blue")

lines(x=c(0,2), y=c(1,1), lty = "dashed")
lines(x=c(2,2), y=c(0,1), lty = "dashed")

text(x=0.04*max(gamma.grid.Neg), y=0.05, cex=1, labels="naive estimate")
text(x=0.045*max(gamma.grid.Neg), y=0, cex=1,
     labels=paste(round(naive,2), "(",naive.lower,"-",naive.upper, ")"))
text(x=0.6, y=0.7, labels="Inconclusive", cex=2, col="grey")
text(x=0.5, y=0.15, labels="Favors RP", cex=2, col="grey")
text(x=0.8, y=1.35, labels="Favors EBRT + AD", cex=2, col="grey")
arrows(x0=0.75, y0=0.98, x1=0.8, y1=1.1, col="steelblue", lwd=2, length=0.15)
arrows(x0=0.6, y0=0.44, x1=0.55, y1=0.3, col="steelblue", lwd=2, length=0.15)

box()
dev.off()




