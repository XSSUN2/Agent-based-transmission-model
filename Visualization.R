draw_aligned_legend <- function(labels, cols, ltys = rep(1, length(labels)),
                                x_line1 = 0.05, x_line2 = 0.25,
                                x_text = 0.32,
                                y_start = 0.80,
                                y_gap = 0.18,
                                lwd = 2,
                                cex = 1.5) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  y <- y_start - seq_along(labels) * 0 + 0
  y <- y_start - (seq_along(labels) - 1) * y_gap
  
  segments(
    x0 = x_line1,
    x1 = x_line2,
    y0 = y,
    y1 = y,
    col = cols,
    lty = ltys,
    lwd = lwd,
    xpd = NA
  )
  
  text(
    x = x_text,
    y = y,
    labels = labels,
    adj = 0,
    cex = cex,
    xpd = NA
  )
}


col_act <- "#4C78A8"
col_pq <- "#3FAE87"
bg_act <- "#C9D2EA33"
bg_pq <- "#C5E8DC33"

H_ACT <- read.csv("H_ACT_MDA.csv")
H_PQ <- read.csv("H_PQ_MDA.csv")
G_ACT <- read.csv("G_ACT_MDA.csv")
G_PQ <- read.csv("G_PQ_MDA.csv")
pm_ACT <- read.csv("pm_ACT_MDA.csv")
pm_PQ <- read.csv("pm_PQ_MDA.csv")

H_ACT_0.5 <- read.csv("H_ACT_0.5.csv")
H_PQ_0.5 <- read.csv("H_PQ_0.5.csv")
G_ACT_0.5 <- read.csv("G_ACT_0.5.csv")
G_PQ_0.5 <- read.csv("G_PQ_0.5.csv")
pm_ACT_0.5 <- read.csv("pm_ACT_0.5.csv")
pm_PQ_0.5 <- read.csv("pm_PQ_0.5.csv")
idx_0.5 <- which(H_ACT_0.5[,1825]!=0&H_PQ_0.5[,1825]!=0)

H_ACT_idpt <- read.csv("H_ACTidpt.csv")
H_PQ_idpt <- read.csv("H_PQidpt.csv")
G_ACT_idpt <- read.csv("G_ACTidpt.csv")
G_PQ_idpt <- read.csv("G_PQidpt.csv")
pm_ACT_idpt <- read.csv("pm_ACTidpt.csv")
pm_PQ_idpt <- read.csv("pm_PQidpt.csv")
idx_idpt <- which(H_ACT_idpt[,1825]!=0&H_PQ_idpt[,1825]!=0)

Difference <- H_ACT-H_PQ
median_diff <- vector()
upper_diff <- vector()
lower_diff <- vector()
for (i in 1:1825) {
  median_diff[i] <- median(Difference[,i])
  upper_diff[i] <- quantile(Difference[,i],0.975)
  lower_diff[i] <- quantile(Difference[,i],0.025)
}

Difference_0.5 <- H_ACT_0.5-H_PQ_0.5
median_diff_0.5 <- vector()
upper_diff_0.5 <- vector()
lower_diff_0.5 <- vector()
for (i in 1:1825) {
  median_diff_0.5[i] <- median(Difference_0.5[idx_0.5,i])
  upper_diff_0.5[i] <- quantile(Difference_0.5[idx_0.5,i],0.975)
  lower_diff_0.5[i] <- quantile(Difference_0.5[idx_0.5,i],0.025)
}

Difference_idpt <- H_ACT_idpt-H_PQ_idpt
median_diff_idpt <- vector()
upper_diff_idpt <- vector()
lower_diff_idpt <- vector()
for (i in 1:1825) {
  median_diff_idpt[i] <- median(Difference_idpt[idx_idpt,i])
  upper_diff_idpt[i] <- quantile(Difference_idpt[idx_idpt,i],0.975)
  lower_diff_idpt[i] <- quantile(Difference_idpt[idx_idpt,i],0.025)
}

##prob
H_ACT <- read.csv("H_ACTidpt.csv")
H_PQ <- read.csv("H_PQidpt.csv")
VI_ACT <- read.csv("VI_ACTidpt.csv")
VE_ACT <- read.csv("VE_ACTidpt.csv")
VI_PQ <- read.csv("VI_PQidpt.csv")
VE_PQ <- read.csv("VE_PQidpt.csv")

H_ACT_0.5 <- read.csv("H_ACT_0.5.csv")
H_PQ_0.5 <- read.csv("H_PQ_0.5.csv")
VI_ACT_0.5 <- read.csv("VI_ACT_0.5.csv")
VE_ACT_0.5 <- read.csv("VE_ACT_0.5.csv")
VI_PQ_0.5 <- read.csv("VI_PQ_0.5.csv")
VE_PQ_0.5 <- read.csv("VE_PQ_0.5.csv")

record_prob_ACT <- matrix(0,100,1825)
record_prob_PQ <- matrix(0,100,1825)
for (i in 1:100) {
  for (j in 1:1825) {
    record_prob_ACT[i,j] <- ifelse(H_ACT[i,j]==0&VE_ACT[i,j]==0&VI_ACT[i,j]==0,1,0)
    record_prob_PQ[i,j] <- ifelse(H_PQ[i,j]==0&VE_PQ[i,j]==0&VI_PQ[i,j]==0,1,0)
  }
  print(i)
}
record_prob_ACT_0.5 <- matrix(0,100,1825)
record_prob_PQ_0.5 <- matrix(0,100,1825)
for (i in 1:100) {
  for (j in 1:1825) {
    record_prob_ACT_0.5[i,j] <- ifelse(H_ACT_0.5[i,j]==0&VE_ACT_0.5[i,j]==0&VI_ACT_0.5[i,j]==0,1,0)
    record_prob_PQ_0.5[i,j] <- ifelse(H_PQ_0.5[i,j]==0&VE_PQ_0.5[i,j]==0&VI_PQ_0.5[i,j]==0,1,0)
  }
  print(i)
}


####Figure 3 prevalence of infections in MDA
pdf("Fig_MDA_prev.pdf", width = 18, height = 11)
layout(
  matrix(c(1, 2, 3, 10,
           4, 5, 6, 11,
           7, 8, 9, 12), nrow = 3, byrow = TRUE),
  widths = c(1, 1, 1, 0.45),
  heights = c(1, 1,1)
)
par(
  oma = c(1.2, 2.2, 2.0, 1.2),
  mgp = c(3, 1, 0),
  tcl = -0.3
)


par(mar = c(1.8, 7, 6.2, 0))
plot(1:(730)*12/365/2-12,H_ACT[1,1:730]/1000,type="n",ylab="Infection prevalence",xlab="",cex.axis=2,cex.lab=2.3,cex.main=2,xlim=c(-3,18),ylim=c(0,0.4),xaxt="n",yaxt="n",bty="l",xaxt="n")
mtext("a", side = 2, line = 0, at = 0.464, las = 1, font = 2, cex = 2)
axis(2,at=c(0,0.1,0.2,0.3,0.4),cex.axis=2.3)
title(main=expression(Full~adherence), line = 3,cex.main=2.5)


for (k in 1:100) {
  lines((366:730)*12/365/2-12,H_ACT[k,(366:730)]/1000,col=rgb(0.8, 0.8, 0.8, alpha = 0.2))
  lines((730:1825)*12/365/2-12,H_ACT[k,730:1825]/1000,col=bg_act)
  lines((730:1825)*12/365/2-12,H_PQ[k,730:1825]/1000,col=bg_pq)
}
segments(x0 = 0, y0 = 0.3, x1 = 0, y1 = 0, lty = 2)
segments(x0 = 1, y0 = 0.3, x1 = 1, y1 = 0, lty = 2)
segments(x0 = 2, y0 = 0.3, x1 = 2, y1 = 0, lty = 2)
text(0,0.33,"1",cex=2)
text(1,0.33,"2",cex=2)
text(2,0.33,"3",cex=2)
lines((366:730)*12/365/2-12,apply(H_ACT,2,median)[366:730]/1000)
lines((730:1825)*12/365/2-12,apply(H_ACT,2,median)[730:1825]/1000,col=col_act)
lines((730:1825)*12/365/2-12,apply(H_PQ,2,median)[730:1825]/1000,col=col_pq)




par(mar = c(1.8, 6, 6.2, 1.0))

plot(1:(730)*12/365/2-12,H_ACT_0.5[1,1:730]/1000,type="n",ylab="",xlab="",cex.axis=2,cex.lab=2.3,cex.main=2,xlim=c(-3,18),ylim=c(0,0.4),xaxt="n",yaxt="n",bty="l",xaxt="n")
mtext("b", side = 2, line = 0, at = 0.464, las = 1, font = 2, cex = 2)
title(main=expression(Half~adherence), line = 3,cex.main=2.5)

for (k in idx_0.5) {
  lines((366:730)*12/365/2-12,H_ACT_0.5[k,(366:730)]/1000,col=rgb(0.8, 0.8, 0.8, alpha = 0.2))
  lines((730:1825)*12/365/2-12,H_ACT_0.5[k,730:1825]/1000,col=bg_act)
  lines((730:1825)*12/365/2-12,H_PQ_0.5[k,730:1825]/1000,col=bg_pq)
}
segments(x0 = 0, y0 = 0.3, x1 = 0, y1 = 0, lty = 2)
segments(x0 = 1, y0 = 0.3, x1 = 1, y1 = 0, lty = 2)
segments(x0 = 2, y0 = 0.3, x1 = 2, y1 = 0, lty = 2)
text(0,0.33,"1",cex=2)
text(1,0.33,"2",cex=2)
text(2,0.33,"3",cex=2)
lines((366:730)*12/365/2-12,apply(H_ACT_0.5[idx_0.5,],2,median)[366:730]/1000)
lines((730:1825)*12/365/2-12,apply(H_ACT_0.5[idx_0.5,],2,median)[730:1825]/1000,col=col_act)
lines((730:1825)*12/365/2-12,apply(H_PQ_0.5[idx_0.5,],2,median)[730:1825]/1000,col=col_pq)


par(mar = c(1.8, 5, 6.2, 2.0))

plot(1:(730)*12/365/2-12,H_ACT_idpt[1,1:730]/1000,type="n",ylab="",xlab="",cex.axis=2,cex.lab=2.3,cex.main=2,xlim=c(-3,18),ylim=c(0,0.4),xaxt="n",yaxt="n",bty="l",xaxt="n")
mtext("c", side = 2, line = 0, at = 0.464, las = 1, font = 2, cex = 2)
title(main=expression(atop(Independent~participation,across~rounds)), line = 2,cex.main=2.5)
for (k in idx_idpt) {
  lines((366:730)*12/365/2-12,H_ACT_idpt[k,(366:730)]/1000,col=rgb(0.8, 0.8, 0.8, alpha = 0.2))
  lines((730:1825)*12/365/2-12,H_ACT_idpt[k,730:1825]/1000,col=bg_act)
  lines((730:1825)*12/365/2-12,H_PQ_idpt[k,730:1825]/1000,col=bg_pq)
}
segments(x0 = 0, y0 = 0.3, x1 = 0, y1 = 0, lty = 2)
segments(x0 = 1, y0 = 0.3, x1 = 1, y1 = 0, lty = 2)
segments(x0 = 2, y0 = 0.3, x1 = 2, y1 = 0, lty = 2)
text(0,0.33,"1",cex=2)
text(1,0.33,"2",cex=2)
text(2,0.33,"3",cex=2)
lines((366:730)*12/365/2-12,apply(H_ACT_idpt[idx_idpt,],2,median)[366:730]/1000)
lines((730:1825)*12/365/2-12,apply(H_ACT_idpt[idx_idpt,],2,median)[730:1825]/1000,col=col_act)
lines((730:1825)*12/365/2-12,apply(H_PQ_idpt[idx_idpt,],2,median)[730:1825]/1000,col=col_pq)


par(mar = c(3.8, 7, 4.2, 0))
plot(1:1825*12/365/2,Difference[1,]/1000,type="n",col="grey",main="",ylab="Paired difference in prevalence\n across simulations",xlab="",cex.axis=1.9,cex.lab=2.3,cex.main=2,xlim=c(-3,18),ylim=c(-0.3,0.3),xaxt="n",yaxt="n",bty="l")
mtext("d", side = 2, line = 0, at = 0.396, las = 1, font = 2, cex = 2)
axis(2,at=c(-0.3,0,0.3),cex.axis=2.3)

for (i in 1:100) {
  lines((366:1825)*12/365/2-12,Difference[i,(366:1825)]/1000,col=rgb(0.8, 0.8, 0.8, alpha = 0.2))
}
lines((366:1825)*12/365/2-12,median_diff[(366:1825)]/1000)
lines((366:1825)*12/365/2-12,upper_diff[(366:1825)]/1000,col="darkgrey")
lines((366:1825)*12/365/2-12,lower_diff[(366:1825)]/1000,col="darkgrey")
abline(h=0,lty=2)

par(mar = c(3.8, 6, 4.2, 1.0))
plot(1:1825*12/365/2,Difference_0.5[1,]/1000,type="n",col="grey",ylab="",xlab="",cex.axis=2,cex.lab=2.3,cex.main=2,xlim=c(-3,18),ylim=c(-0.3,0.3),xaxt="n",yaxt="n",bty="l")
mtext("e", side = 2, line = 0, at = 0.396, las = 1, font = 2, cex = 2)

for (i in idx_0.5) {
  lines((366:1825)*12/365/2-12,Difference_0.5[i,(366:1825)]/1000,col=rgb(0.8, 0.8, 0.8, alpha = 0.2))
}
lines((366:1825)*12/365/2-12,median_diff_0.5[(366:1825)]/1000)
lines((366:1825)*12/365/2-12,upper_diff_0.5[(366:1825)]/1000,col="darkgrey")
lines((366:1825)*12/365/2-12,lower_diff_0.5[(366:1825)]/1000,col="darkgrey")
abline(h=0,lty=2)


par(mar = c(3.8, 5, 4.2, 2.0))
plot(1:1825*12/365/2,Difference_idpt[1,]/1000,type="n",col="grey",ylab="",xlab="",cex.axis=2,cex.lab=2.3,cex.main=2,xlim=c(-3,18),ylim=c(-0.3,0.3),xaxt="n",yaxt="n",bty="l")
mtext("f", side = 2, line = 0, at = 0.396, las = 1, font = 2, cex = 2)


for (i in idx_idpt) {
  lines((366:1825)*12/365/2-12,Difference_idpt[i,(366:1825)]/1000,col=rgb(0.8, 0.8, 0.8, alpha = 0.2))
}
lines((366:1825)*12/365/2-12,median_diff_idpt[(366:1825)]/1000)
lines((366:1825)*12/365/2-12,upper_diff_idpt[(366:1825)]/1000,col="darkgrey")
lines((366:1825)*12/365/2-12,lower_diff_idpt[(366:1825)]/1000,col="darkgrey")
abline(h=0,lty=2)

par(mar = c(5.8, 7, 2.2, 0))
plot(1:1825*12/365/2-12,colSums(record_prob_PQ)/100,type="n",ylab="Proportion of simulations\n with extinction",xlab="Time post MDA (Months)",cex.axis=2,cex.lab=2.3,cex.main=2,xlim=c(-3,18),ylim=c(0,1),xaxt="n",yaxt="n",bty="l")
mtext("g", side = 2, line = 0, at = 1.2, las = 1, font = 2, cex = 2)
axis(2,at=c(0,0.5,1),cex.axis=2.3)
axis(1,at=c(-3,0,3,6,9,12,15,18),cex.axis=2.3)
lines(1:1825*12/365/2-12,rep(0,1825),col="#4C78A8",lwd=1.5)
lines(1:1825*12/365/2-12,rep(0,1825),col="#3FAE87B1",lwd=1.5)

par(mar = c(5.8, 6, 2.2, 1.0))
plot(1:1825*12/365/2-12,colSums(record_prob_PQ)/100,type="n",ylab="",xlab="Time post MDA (Months)",cex.axis=2,cex.lab=2.3,cex.main=2,xlim=c(-3,18),ylim=c(0,1),xaxt="n",yaxt="n",bty="l")
mtext("h", side = 2, line = 0, at = 1.2, las = 1, font = 2, cex = 2)
axis(1,at=c(-3,0,3,6,9,12,15,18),cex.axis=2.3)
lines(1:1825*12/365/2-12,colSums(record_prob_ACT_0.5)/100,col=col_act,lwd=1.5)
lines(1:1825*12/365/2-12,colSums(record_prob_PQ_0.5)/100,col=col_pq,lwd=1.5)

par(mar = c(5.8, 5, 2.2, 2.0))
plot(1:1825*12/365/2-12,colSums(record_prob_PQ)/100,type="n",ylab="",xlab="Time post MDA (Months)",cex.axis=2,cex.lab=2.3,cex.main=2,xlim=c(-3,18),ylim=c(0,1),xaxt="n",yaxt="n",bty="l")
mtext("i", side = 2, line = 0, at = 1.2, las = 1, font = 2, cex = 2)
axis(1,at=c(-3,0,3,6,9,12,15,18),cex.axis=2.3)
lines(1:1825*12/365/2-12,colSums(record_prob_ACT)/100,col=col_act,lwd=1.5)
lines(1:1825*12/365/2-12,colSums(record_prob_PQ)/100,col=col_pq,lwd=1.5)


par(mar = c(0, 0, 10.2, 0.1))
draw_aligned_legend(
  labels = c("Before MDA", "ACT alone", "ACT+PQ 0.25mg/kg"),
  cols = c("black", col_act, col_pq),
  ltys = c(1, 1, 1),
  x_line1 = 0.05,
  x_line2 = 0.25,
  x_text = 0.32,
  y_start = 0.80,
  y_gap = 0.1,
  lwd = 2,
  cex = 1.5
)

par(mar = c(5.2, 0, 5, 0.1))
draw_aligned_legend(
  labels = c("Trajectories", "Median", "95% PI"),
  cols = c(rgb(0.8, 0.8, 0.8, alpha = 0.6), "black", "darkgrey"),
  ltys = c(1, 1, 1),
  x_line1 = 0.05,
  x_line2 = 0.25,
  x_text = 0.32,
  y_start = 0.80,
  y_gap = 0.1,
  lwd = 2,
  cex = 1.5
)

par(mar = c(5.2, 0, 5, 0.1))
draw_aligned_legend(
  labels = c("ACT alone", "ACT+PQ 0.25mg/kg"),
  cols = c(col_act, col_pq),
  ltys = c( 1, 1),
  x_line1 = 0.05,
  x_line2 = 0.25,
  x_text = 0.32,
  y_start = 0.80,
  y_gap = 0.1,
  lwd = 2,
  cex = 1.5
)

dev.off()
system("pdfcrop Fig_MDA_prev.pdf Fig_MDA_prev.pdf")


###case management
H_PQ <- read.csv("H_PQ_TS.csv")
H_ACT <- read.csv("H_ACT_TS.csv")
G_PQ <- read.csv("G_PQ_TS.csv")
G_ACT <- read.csv("G_ACT_TS.csv")
pm_ACT <- read.csv("pm_ACT_TS.csv")
pm_PQ <- read.csv("pm_PQ_TS.csv")

col_act <- "#4C78A8"
col_pq <- "#3FAE87"
bg_act <- "#C9D2EA33"
bg_pq <- "#C5E8DC33"

####Figure 5 routine case management
pdf("Fig_RTS.pdf", width = 15, height = 4)
layout(
  matrix(c(1, 2, 3), nrow = 1, byrow = TRUE),
  widths = c(1, 1, 1),
  heights = c(1, 1)
)
par(
  oma = c(1.2, 2.2, 2.0, 1.2),
  mgp = c(3, 1, 0),
  tcl = -0.3
)

par(mar = c(1.8, 7, 6.2, 2))
plot(1:1825*12/365/2-12,H_ACT[1,]/1000,type="n",ylab="Infection prevalence",xlab="Time post treatment change (Months)",cex.axis=2,cex.lab=2,cex.main=2,xlim=c(-3,18),ylim=c(0.05,0.35),yaxt="n",xaxt="n",bty="l")
mtext("a", side = 2, line = 5, at = 0.398, las = 1, font = 2, cex = 1.8)
axis(2,at=c(0.05,0.2,0.35),cex.axis=2)
axis(1,at=c(-3,0,3,6,9,12,15,18),cex.axis=2)
for (k in 1:100) {
  lines((366:730)*12/365/2-12,H_ACT[k,(366:730)]/1000,col=bg_act)
  lines((730:1825)*12/365/2-12,H_ACT[k,730:1825]/1000,col=bg_act)
  lines((730:1825)*12/365/2-12,H_PQ[k,730:1825]/1000,col=bg_pq)
}
lines((366:730)*12/365/2-12,apply(H_ACT[,(366:730)],2,median)/1000,col=col_act)
lines((730:1825)*12/365/2-12,apply(H_ACT[,730:1825],2,median)/1000,col=col_act)
lines((730:1825)*12/365/2-12,apply(H_PQ[,730:1825],2,median)/1000,col=col_pq)
legend("topright",lty=c(1,1),lwd=2,col=c(col_act,col_pq),legend=c("ACT alone","ACT+PQ 0.25mg/kg"),cex=1.5,bty="n")


par(mar = c(1.8, 7, 6.2, 2.0))
plot(1:(1825)*12/365/2-12,G_ACT[1,]/4000,type="n",ylab="Population-averaged\n gametocytemia",xlab="",cex.axis=2,cex.lab=2,cex.main=2,xlim=c(-3,18),ylim=c(0,900),xaxt="n",yaxt="n",bty="l")
mtext("b", side = 2, line = 5, at = 1044, las = 1, font = 2, cex = 1.8)
axis(2,at=c(0,450,900),cex.axis=2)
axis(1,at=c(-3,0,3,6,9,12,15,18),cex.axis=2)
for (k in 1:100) {
  lines((366:730)*12/365/2-12,G_ACT[k,(366:730)],col=bg_act)
  lines((730:1825)*12/365/2-12,G_ACT[k,730:1825],col=bg_act)
  lines((730:1825)*12/365/2-12,G_PQ[k,730:1825],col=bg_pq)
}


lines((366:730)*12/365/2-12,apply(G_ACT[,(366:730)],2,median),col=col_act)
lines((730:1825)*12/365/2-12,apply(G_ACT[,730:1825],2,median),col=col_act)
lines((730:1825)*12/365/2-12,apply(G_PQ[,730:1825],2,median),col=col_pq)
legend("topright",lty=c(1,1),lwd=2,col=c(col_act,col_pq),legend=c("ACT alone","ACT+PQ 0.25mg/kg"),cex=1.5,bty="n")

par(mar = c(1.8, 7, 6.2, 2.0))
plot((366:730)*12/365/2-12,pm_ACT[1,1:365],type="n",ylab="Population-averaged\n infectiousness to mosquitoes",xlab="Time post treatment change (Months)",cex.axis=1.8,cex.lab=2,cex.main=2,xlim=c(-3,18),ylim=c(0,0.014),xaxt="n",yaxt="n",bty="l")
mtext("c", side = 2, line = 5, at = 0.01624, las = 1, font = 2, cex = 1.8)
axis(2,at=c(0,0.007,0.014),cex.axis=2)
axis(1,at=c(-3,0,3,6,9,12,15,18),cex.axis=2)
for (k in 1:100) {
  lines((366:730)*12/365/2-12,pm_ACT[k,(366:730)],col=bg_act)
  lines((730:1825)*12/365/2-12,pm_ACT[k,730:1825],col=bg_act)
  lines((730:1825)*12/365/2-12,pm_PQ[k,730:1825],col=bg_pq)
}


lines((366:730)*12/365/2-12,apply(pm_ACT[,(366:730)],2,median),col=col_act)
lines((730:1825)*12/365/2-12,apply(pm_ACT[,730:1825],2,median),col=col_act)
lines((730:1825)*12/365/2-12,apply(pm_PQ[,730:1825],2,median),col=col_pq)
legend("topright",lty=c(1,1),lwd=2,col=c(col_act,col_pq),legend=c("ACT alone","ACT+PQ 0.25mg/kg"),cex=1.5,bty="n")


dev.off()
system("pdfcrop Fig_RTS.pdf Fig_RTS.pdf")
