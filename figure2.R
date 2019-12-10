########## Performance w.r.t. adjusted Rand Index on Vincent Function
require(data.table)
require(ggplot2)
setwd("/home/npav/new_code/michalis/reprod/")
A <- fread("VincentZoom_de_rand1_k4_B2.txt")
A$method = factor(toupper(A$method), levels = c("GPBI","NBC2","SEED","BAPD"))
p <- ggplot(data=A, aes(x=iter,y=aRand, colour=method, shape=method, linetype=method)) + 
  geom_line(size=1.) + geom_point(size=2.5) +
  scale_x_continuous(breaks = c(1,seq(5,50,5)), limits = c(2,50)) + xlab("DE/rand/1 iteration") + 
  #ylab("Adjusted Rand Index") + 
  theme_bw()+ theme(text = element_text(size=16), panel.grid = element_blank()) +
  theme(text = element_text(size=28)) + 
  #theme(legend.title = element_blank(), legend.position = c(0.9,0.5))
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(axis.title.y = element_blank())
print(p)  

fname <- "/home/npav/new_code/michalis/latex/figs/VincentPerf."
ggsave(filename=paste0(fname,'png'), plot=p, width=12, height=8, units="in", device = "png")
dev.off()

print(p)  
ggsave(filename=paste0(fname,'pdf'), plot=p, width=12, height=8, units="in", device = "pdf")
dev.off()

# Summary information
A[method=="GPBI",.(iter,aRand, error=abs(numclusters-trueclusters))]

########## Performance w.r.t. adjusted Rand Index on Vincent Function
A[, Dclust := numclusters - trueclusters]
p <- ggplot(data=A, aes(x=iter,y=Dclust, colour=method, shape=method, linetype=method)) + 
  geom_line(size=1.) + geom_point(size=2.5) + geom_hline(yintercept=0, alpha=0.5, linetype="dashed") +
  scale_x_continuous(breaks = c(1,seq(5,50,5)), limits = c(2,50)) + xlab("DE/rand/1 iteration") + 
  #ylab("Adjusted Rand Index") +
  theme_bw()+ theme(text = element_text(size=16), panel.grid = element_blank()) +
  theme(text = element_text(size=28)) + 
  #theme(legend.title = element_blank(), legend.position = c(0.9,0.5))
  theme(legend.title = element_blank(), legend.position = "bottom") +
  theme(axis.title.y = element_blank())
print(p)  

fname <- "/home/npav/new_code/michalis/latex/figs/VincentClust."
ggsave(filename=paste0(fname,'png'), plot=p, width=12, height=8, units="in", device = "png")
dev.off()

print(p)  
ggsave(filename=paste0(fname,'pdf'), plot=p, width=12, height=8, units="in", device = "pdf")
dev.off()


# ARI at first iteration
th <- A[method=="GPBI"]$aRand[1]

# performance of NBC2 exceeds performance of GPBI on first iteration
A[method=="NBC2",]$aRand > th
which(cumsum(A[method=="NBC2",]$aRand > th) == 1)
which(cumsum(A[method=="NBC2",]$aRand > th) == 4)

A[method=="SEED",]$aRand > th
which(cumsum(A[method=="SEED",]$aRand > th) == 1)

# how many times GPBI has perfect performance
sum(A[method=="GPBI",.(aRand)] == 1)

# number of times GPBI estimates number of clusters correctly 
sum(A[method=="GPBI",.(Dclust)] ==0)
