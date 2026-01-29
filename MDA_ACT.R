#MDA treatment

library(truncnorm)
library(deSolve)
odefunction <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    dD <- -kt*D
    dT1 <- kt*D-kt*T1
    dT2 <- kt*T1-kt*T2
    dC <- (kt*T2+q1*P1+q2*P2-q1*C-q2*C-qc*C)/Vc
    dP1 <- (q1*C-q1*P1)/V1
    dP2 <- (q2*C-q2*P2)/V2
    list(c(dD,dT1,dT2,dC,dP1,dP2))
  })
}

kd <- function(time,dr,kmax,ec50,gamma){
  dru <- cbind(time,dr*1000,0) 
  n <- dim(dru)[1]
  x <- dim(dru)[2]
  for (i in 2:n) {
    kt <- kmax * (dru[i,2])^gamma / ((dru[i,2])^gamma + ec50^gamma) #calculate k(t)
    kt_1 <- kmax * (dru[i-1,2])^gamma / ((dru[i-1,2])^gamma + ec50^gamma) #calculate k(t-1)
    kd <- (kt+kt_1)/2 #killing effect
    dru[i,x] <- kd #record
  }
  return(dru[,3])
}
####PPQ
dr_time <- 15*24
parameters <- c(kt = 1/2.79, qc = 51.4, Vc = 804, q1 = 2020, V1 = 3010, q2 = 149, V2 = 13300)
state <- c(D = 0, T1 = 0, T2 = 0, C = 0, P1 = 0, P2 = 0)
times <- seq(1,dr_time-1,by = 1)
out <- ode(state, times, odefunction, parameters)
drug <- data.frame(out)


state1 <- c(D = 960, T1 = 0, T2 = 0, C = 0, P1 = 0, P2 = 0)
times1 <- seq(dr_time,dr_time+23,by = 1)
out1 <- ode(state1, times1, odefunction, parameters)
drug1 <- data.frame(out1)
dim <- dim(drug1)[1]
state2 <- c(D=drug1[dim,2]+960, T1=drug1[dim,3],T2=drug1[dim,4],C=drug1[dim,5],P1=drug1[dim,6],P2=drug1[dim,7])
times2 <- seq(dr_time+24,dr_time+47,by=1)
out2 <- ode(state2, times2, odefunction, parameters)
drug2 <- data.frame(out2)
dim <- dim(drug2)[1]
state3 <- c(D=drug2[dim,2]+960, T1=drug2[dim,3],T2=drug2[dim,4],C=drug2[dim,5],P1=drug2[dim,6],P2=drug2[dim,7])
times3 <- seq(dr_time+48,dr_time+719,by=1)
out3 <- ode(state3, times3, odefunction, parameters)
drug3 <- data.frame(out3)

dim <- dim(drug3)[1]
state4 <- c(D=drug3[dim,2]+960, T1=drug3[dim,3],T2=drug3[dim,4],C=drug3[dim,5],P1=drug3[dim,6],P2=drug3[dim,7])
times4 <- seq(dr_time+720,dr_time+743,by=1)
out4 <- ode(state4, times4, odefunction, parameters)
drug4 <- data.frame(out4)
dim <- dim(drug4)[1]
state5 <- c(D=drug4[dim,2]+960, T1=drug4[dim,3],T2=drug4[dim,4],C=drug4[dim,5],P1=drug4[dim,6],P2=drug4[dim,7])
times5 <- seq(dr_time+744,dr_time+767,by=1)
out5 <- ode(state5, times5, odefunction, parameters)
drug5 <- data.frame(out5)
dim <- dim(drug5)[1]
state6 <- c(D=drug5[dim,2]+960, T1=drug5[dim,3],T2=drug5[dim,4],C=drug5[dim,5],P1=drug5[dim,6],P2=drug5[dim,7])
times6 <- seq(dr_time+768,dr_time+1439,by=1)
out6 <- ode(state6, times6, odefunction, parameters)
drug6 <- data.frame(out6)

dim <- dim(drug6)[1]
state7 <- c(D=drug6[dim,2]+960, T1=drug6[dim,3],T2=drug6[dim,4],C=drug6[dim,5],P1=drug6[dim,6],P2=drug6[dim,7])
times7 <- seq(dr_time+1440,dr_time+1463,by=1)
out7 <- ode(state7, times7, odefunction, parameters)
drug7 <- data.frame(out7)
dim <- dim(drug7)[1]
state8 <- c(D=drug7[dim,2]+960, T1=drug7[dim,3],T2=drug7[dim,4],C=drug7[dim,5],P1=drug7[dim,6],P2=drug7[dim,7])
times8 <- seq(dr_time+1464,dr_time+1487,by=1)
out8 <- ode(state8, times8, odefunction, parameters)
drug8 <- data.frame(out8)
dim <- dim(drug8)[1]
state9 <- c(D=drug8[dim,2]+960, T1=drug8[dim,3],T2=drug8[dim,4],C=drug8[dim,5],P1=drug8[dim,6],P2=drug8[dim,7])
times9 <- seq(dr_time+1488,3500,by=1)
out9 <- ode(state9, times9, odefunction, parameters)
drug9 <- data.frame(out9)

PQ<- rbind(drug,drug1,drug2,drug3,drug4,drug5,drug6,drug7,drug8,drug9)
t <- 1:3500
kd_PQ <- kd(t,PQ$C,0.2,9,12)
par(mar=c(5,7, 1.1,1.1))
plot((8761):(11900)/365/2-12,PQ$C[360:3499],ylab="PK curve of PPQ",type="l",xlab="Time post MDA (Months)",cex.axis=1.8,cex.lab=2,cex.main=2)
plot((8761):(11900)/365/2-12,kd_PQ[360:3499],ylab="PD curve of PPQ",type="l",xlab="Time post MDA (Months)",cex.axis=1.8,cex.lab=2,cex.main=2)


###DHA
ka <- 0.89
CL_F <- 48.8
V_F <- 44.4
D <- 120
Central <- D*ka/(V_F*ka-CL_F)*(exp(-(CL_F/V_F)*t)-exp(-ka*t))
DHA <- c(rep(0,dr_time),(Central+c(rep(0,24),Central[1:(length(Central)-24)])+c(rep(0,48),Central[1:(length(Central)-48)]))[1:720],(Central+c(rep(0,24),Central[1:(length(Central)-24)])+c(rep(0,48),Central[1:(length(Central)-48)]))[1:720],(Central+c(rep(0,24),Central[1:(length(Central)-24)])+c(rep(0,48),Central[1:(length(Central)-48)]))[1:3500])

plot((8761):(11900)/365/2-12,DHA[360:3499],ylab="PK curve of DHA",type="l",xlab="Time post MDA (Months)",cex.axis=1.8,cex.lab=2,cex.main=2)
kd_DHA <- kd(1:length(DHA),DHA,0.4,1.5,1.3)
plot((8761):(11900)/365/2-12,kd_DHA[360:3499],ylab="PD curve of DHA",type="l",xlab="Time post MDA (Months)",cex.axis=1.8,cex.lab=2,cex.main=2)

