#MDA treatment
library(truncnorm)
library(deSolve)
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

####Piperaquine
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

PPQ<- rbind(drug,drug1,drug2,drug3,drug4,drug5,drug6,drug7,drug8,drug9)
t <- 1:3500
kd_MDA_PPQ_three <- kd(t,PPQ$C,0.2,9,12)[360:3499]

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
state7 <- c(D=drug6[dim,2], T1=drug6[dim,3],T2=drug6[dim,4],C=drug6[dim,5],P1=drug6[dim,6],P2=drug6[dim,7])
times7 <- seq(dr_time+1440,dr_time+1463,by=1)
out7 <- ode(state7, times7, odefunction, parameters)
drug7 <- data.frame(out7)
dim <- dim(drug7)[1]
state8 <- c(D=drug7[dim,2], T1=drug7[dim,3],T2=drug7[dim,4],C=drug7[dim,5],P1=drug7[dim,6],P2=drug7[dim,7])
times8 <- seq(dr_time+1464,dr_time+1487,by=1)
out8 <- ode(state8, times8, odefunction, parameters)
drug8 <- data.frame(out8)
dim <- dim(drug8)[1]
state9 <- c(D=drug8[dim,2], T1=drug8[dim,3],T2=drug8[dim,4],C=drug8[dim,5],P1=drug8[dim,6],P2=drug8[dim,7])
times9 <- seq(dr_time+1488,3500,by=1)
out9 <- ode(state9, times9, odefunction, parameters)
drug9 <- data.frame(out9)

PPQ<- rbind(drug,drug1,drug2,drug3,drug4,drug5,drug6,drug7,drug8,drug9)
t <- 1:3500
kd_MDA_PPQ_two <- kd(t,PPQ$C,0.2,9,12)[360:3499]

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
times3 <- seq(dr_time+48,2400,by=1)
out3 <- ode(state3, times3, odefunction, parameters)
drug3 <- data.frame(out3)
PPQ<- rbind(drug,drug1,drug2,drug3)
t=1:2400
kd_PPQ <- kd(t,PPQ$C,0.2,9,12)[360:1799] ##length=1440


###DHA
ode_DHA <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    dD <- -kt*D
    dT1 <- kt*D-kt*T1
    dT2 <- kt*T1-kt*T2
    dT3 <- kt*T2-kt*T3
    dT4 <- kt*T3-kt*T4
    dA <- kt*T4-CLA_F*A/VA_F
    dDHA <- CLA_F*A/VA_F-CLD_F*DHA/VD_F-KE*DHA
    dC <- KE*DHA-KE*C
    
    list(c(dD,dT1,dT2,dT3,dT4,dA,dDHA,dC))
  })
}

dr_time <- 15*24
parameters <- c(kt = 5/0.869, CLA_F = 1890, VA_F = 778, CLD_F = 116, VD_F = 126, KE = 0.0717)
state <- c(D = 0, T1 = 0, T2 = 0, T3=0, T4=0, A = 0, DHA = 0, C = 0)
times <- seq(1,dr_time-1,by = 1)
out <- ode(state, times, ode_DHA, parameters)
drug <- data.frame(out)
state1 <- c(D = 120, T1 = 0, T2 = 0, T3=0, T4=0, A = 0, DHA = 0, C = 0)
times1 <- seq(dr_time,dr_time+23,by = 1)
out1 <- ode(state1, times1, ode_DHA, parameters)
drug1 <- data.frame(out1)
dim <- dim(drug1)[1]
state2 <- c(D=drug1[dim,2]+120, T1=drug1[dim,3],T2=drug1[dim,4],T3=drug1[dim,5],T4=drug1[dim,6],A=drug1[dim,7],DHA=drug1[dim,8],C=drug1[dim,9])
times2 <- seq(dr_time+24,dr_time+47,by=1)
out2 <- ode(state2, times2, ode_DHA, parameters)
drug2 <- data.frame(out2)
dim <- dim(drug2)[1]
state3 <- c(D=drug2[dim,2]+120, T1=drug2[dim,3],T2=drug2[dim,4],T3=drug2[dim,5],T4=drug2[dim,6],A=drug2[dim,7],DHA=drug2[dim,8],C=drug2[dim,9])
times3 <- seq(dr_time+48,2400,by=1)
out3 <- ode(state3, times3, ode_DHA, parameters)
drug3 <- data.frame(out3)
DHA<- rbind(drug,drug1,drug2,drug3)
kd_DHA <- kd(t,DHA$C/284.35/124*10^6,0.273,34.9,1)[360:799]## lenth=440


##primaquine
pk_PRQ <- function(t, state, parameters){
  with(as.list(c(state,parameters)),{
    dG <- -ka*G
    dT1 <- ka*G-ktr*T1
    dT2 <- ktr*T1-ktr*T2
    dT3 <- ktr*T2-ktr*T3
    dT4 <- ktr*T3-ktr*T4
    dT5 <- ktr*T4-ktr*T5
    dT6 <- ktr*T5-ktr*T6
    dP <- (1-Fm)*ktr*T6-CL_Fp*P/V_Fp-kdr*P
    dC <- Fm*ktr*T6+CL_Fp/V_Fp*P-CL_Fc/V_Fc*C
    dD1 <- kdr*P-kdr*D1
    dD2 <- kdr*D1-kdr*D2
    dD3 <- kdr*D2-kdr*D3
    list(c(dG,dT1,dT2,dT3,dT4,dT5,dT6,dP,dC,dD1,dD2,dD3))
  })
}
PRQ <- 0.25
parameters_PRQ <- c(ka = 2.27, ktr = 6/0.424, Fm = 0.359, kdr = 1.02, CL_Fp = 20.0, V_Fp =144 , CL_Fc = 1.06, V_Fc = 35.4)
state_PRQ <- c(G = 0, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times <- seq(1,dr_time-1,by = 1)
out_PRQ <- ode(state_PRQ, times, pk_PRQ, parameters_PRQ)
drug_PRQ <- data.frame(out_PRQ)
state_PRQ1 <- c(G = PRQ*60, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times1 <- seq(dr_time,2400,by = 1)
out_PRQ1 <- ode(state_PRQ1, times1, pk_PRQ, parameters_PRQ,method="radau",atol = 1e-9, rtol = 1e-9)
drug_PRQ1 <- data.frame(out_PRQ1)
PRQ_pk<- rbind(drug_PRQ,drug_PRQ1)

#length=440
kd_gf <- kd(2400,PRQ_pk$D3*1000/144/259.347,84,0.0001,1)[360:799]
kd_gm <- kd(2400,PRQ_pk$D3*1000/144/259.347,84,0.0001,1)[360:799]

kd_ACT <- kd_PPQ
kd_ACT[1:length(kd_DHA)] <- kd_ACT[1:length(kd_DHA)]+kd_DHA


