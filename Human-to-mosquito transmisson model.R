sto_model <- function(stage,N,mu,sd,rp,kd_ring,f,detp,detg,detgm,detgf,m,fg,kd_g,kd_gf,kd_gm,slop,dose,detgcl,al=42,as=25,t){
  rec <- matrix(0,t,(2*al+7))  
  if(N!=0){ #If inoculation occurs, initialise a parasite age distribution at the time of the mosquito bite
    pinit <- rtruncnorm(max(round(N,0),1),1,al,mu,sd) 
    pinit <- round(pinit,digits = 0)
    for (i in 1:al) { #record initial number of parasites
      rec[1,i] <- length(which(pinit==i)) 
    }
    rec[1,] <- rec[1,]+stage
  }else{ #If not, just use the intial setting of parasites in each age bins
    rec[1,] <- stage 
  }
  rec_f <- vector()
  rec_f[1] <- 0
  rec_fg <- vector()
  rec_fg[1] <- 0
  
  for (i in 2:t) { 
    for (j in 1:al) { #asexual parasites
      if (j==1){ #age = 1 hour
        # number of parasites in age bin a_L at t-1 * replication rate - number of death (detp--death rate ,kd--killing effect )
        p.a.t <- (rec[i-1,al]-rpois(1,(detp+kd_ring[i])*rec[i-1,al]))*rp[i]     # make sure it's not negative
        p.a.t <- max(0,round(p.a.t,0))
        rec[i,j] <- p.a.t
      }
      else{if (j == as+1){ #sexual conversion
        p.a.t <- rec[i-1,as]-rpois(1,(detp+kd_ring[i])*rec[i-1,as])
        p.a.t <- max(0,p.a.t)
        rec_f[i] <- rbinom(1,p.a.t,f)
        rec[i,j] <- p.a.t-rec_f[i]
      }
        else{
          p.a.t <- rec[i-1,j-1]-rpois(1,(detp+kd_ring[i])*rec[i-1,j-1])
          p.a.t <- max(0,p.a.t)
          rec[i,j] <- p.a.t
          
        }
      }
    }
    for (j in (al+1):(2*al-1)){ #sexual development
      #column al+1 is Pg(as+1), al+2 is Pg(as+2), ... 
      if (j==(2*al-as+1)){        # replication 
        p.a.t <- (rec[i-1,2*al-as]-rpois(1,(detp+kd_g[i])*rec[i-1,2*al-as]))*rp[i]
        p.a.t <- max(0,round(p.a.t,0))
        rec[i,j] <- p.a.t
      }else{if (j==(al+1)){  #sexual conversion
        p.a.t <- rec_f[i]-rpois(1,(detp+kd_g[i])*rec_f[i])
        p.a.t <- max(0,p.a.t)
        rec[i,j] <- p.a.t
      }
        else{
          p.a.t <- rec[i-1,j-1]-rpois(1,(detp+kd_g[i])*rec[i-1,j-1])
          p.a.t <- max(0,p.a.t)
          rec[i,j] <- p.a.t
        }
      }
    }
    # j in (2*al):(2*al+4) G1 to G5
    #random a number to be parasites going from G1 to G2, compare with number in G1
    #to make sure G1 has enough parasites to go to G2.
    g1tog2 <- min(rec[i-1,2*al],rpois(1,m*rec[i-1,2*al]))
    rec[i,2*al] <- max(0,rec[i-1,2*al]+rec[i-1,2*al-1]-rpois(1,(detg+kd_g[i]+detg)*rec[i-1,2*al])-rpois(1,detp*rec[i-1,2*al-1])-g1tog2)
    
    g2tog3 <- min(rec[i-1,2*al+1],rpois(1,m*rec[i-1,2*al+1]))
    rec[i,2*al+1] <- max(0,rec[i-1,2*al+1]+g1tog2-g2tog3-rpois(1,(detg+kd_g[i]+detg)*rec[i-1,2*al+1]))
    
    g3tog4 <- min(rec[i-1,2*al+2],rpois(1,m*rec[i-1,2*al+2]))
    rec[i,2*al+2] <- max(0,rec[i-1,2*al+2]+g2tog3-g3tog4-rpois(1,(detg+kd_g[i]+detg)*rec[i-1,2*al+2]))
    
    g4tog5 <- min(rec[i-1,2*al+3],rpois(1,m*rec[i-1,2*al+3]))
    rec[i,2*al+3] <- max(0,rec[i-1,2*al+3]+g3tog4-g4tog5-rpois(1,(detg+kd_g[i]+detg)*rec[i-1,2*al+3]))
    
    rec_fg[i] <- rbinom(1, g4tog5,fg)
    
    kill_f <- min(rpois(1,(detgf*(1+slop*dose*kd_gf[i]))*(rec[i-1,2*al+5])),rec[i-1,2*al+5])
    kill_m <- min(rpois(1,(detgm*(1+slop*dose*kd_gm[i]))*(rec[i-1,2*al+4])),rec[i-1,2*al+4])
    
    rec[i,2*al+4] <- max(0,rec[i-1,2*al+4]+rec_fg[i] -kill_m)
    rec[i,2*al+5] <- max(0,rec[i-1,2*al+5]+g4tog5-rec_fg[i]-kill_f)
    rec[i,2*al+6] <- max(0,kill_m+rec[i-1,2*al+6]-rpois(1,detgcl*rec[i-1,2*al+6]))
    rec[i,2*al+7] <- max(0,kill_f+rec[i-1,2*al+7]-rpois(1,detgcl*rec[i-1,2*al+7]))
    
    if (anyNA(rec[i,])==TRUE){
      rec=rec[-(i:t),]
      break
    }
  }
  
  rec <- data.frame(rec)
  return(rec)
}

###piperaquine
library(truncnorm)
library(deSolve)
#PK model of piperaquine
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
#PD model
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
####piperaquine
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

###PK model of DHA
ka <- 0.89
CL_F <- 48.8
V_F <- 44.4
t <- 1:2400
D <- 120
Central <- D*ka/(V_F*ka-CL_F)*(exp(-(CL_F/V_F)*t)-exp(-ka*t))
DHA <- c(rep(0,dr_time),(Central+c(rep(0,24),Central[1:(length(Central)-24)])+c(rep(0,48),Central[1:(length(Central)-48)]))[1:(length(t)-dr_time)])

##killing effect of piperaquine and DHA
kd_PPQ <- kd(t,PPQ$C,0.2,9,12)
kd_DHA <- kd(t,DHA,0.4,1.5,1.3)
kd_ring <- kd_PPQ+kd_DHA

##PK model of primaquine
pk_PQ <- function(t, state, parameters){
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

# 100 simulations for ACT only
benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0005,0.0022,0.02,0.2, 2.44,PQ,0.0247)
kd_gf = kd_gm = rep(0,2400)
record_ACTonly_ring <- matrix(0,100,720)
record_ACTonly_game <- matrix(0,100,720)
for (i in 1:100) {
  simulation_ACTonly <- try(sto_model(matrix(0,1,2*42+7),benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2400),kd_ring,benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],kd_DHA,kd_gf,kd_gm,benchmark[11],benchmark[12],benchmark[13],al=42,as=25,t=720),silent = TRUE)
  record_ACTonly_ring[i,] <- ifelse(rowSums(simulation_ACTonly[,1:25])==0,1,rowSums(simulation_ACTonly[,1:25]))
  record_ACTonly_game[i,] <- simulation_ACTonly[,88]+simulation_ACTonly[,89]
}
median_ACTonly <- apply(record_ACTonly_game,2,median) 
upper_ACTonly <- lower_ACTonly <- vector()

for (i in 1:720) {
  upper_ACTonly[i] <- quantile(record_ACTonly_game[,i],0.975) 
  lower_ACTonly[i] <- quantile(record_ACTonly_game[,i],0.025) 
}

##for ACT+0.125 mg/kg PQ
PQ <- 0.125
parameters_PQ <- c(ka = 2.27, ktr = 6/0.424, Fm = 0.359, kdr = 1.02, CL_Fp = 20.0, V_Fp =144 , CL_Fc = 1.06, V_Fc = 35.4)
state_PQ <- c(G = 0, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times <- seq(1,dr_time-1,by = 1)
out_PQ <- ode(state_PQ, times, pk_PQ, parameters_PQ)
drug_PQ <- data.frame(out_PQ)
state_PQ1 <- c(G = PQ*60/259.347*1000, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times1 <- seq(dr_time,2400,by = 1)
out_PQ1 <- ode(state_PQ1, times1, pk_PQ, parameters_PQ,method="radau",atol = 1e-9, rtol = 1e-9)
drug_PQ1 <- data.frame(out_PQ1)
PQ_pk<- rbind(drug_PQ,drug_PQ1)
kd_gf <- kd(2400,PQ_pk$D3,84,0.0001,1)
kd_gm <- kd(2400,PQ_pk$D3,84,0.0001,1)
benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0005,0.0022,0.02,0.2, 2.44,PQ,0.0247)

record_PQ0.125_ring <- matrix(0,100,720)
record_PQ0.125_game <- matrix(0,100,720)
for (i in 1:100) {
  simulation_PQ0.125 <- try(sto_model(matrix(0,1,2*42+7),benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2400),kd_ring,benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],kd_DHA,kd_gf,kd_gm,benchmark[11],benchmark[12],benchmark[13],al=42,as=25,t=720),silent = TRUE)
  record_PQ0.125_ring[i,] <- length(ifelse(rowSums(simulation_PQ0.125[,1:25])==0,1,rowSums(simulation_PQ0.125[,1:25])))
  record_PQ0.125_game[i,] <- simulation_PQ0.125[,88]+simulation_PQ0.125[,89]
}
median_PQ0.125 <- apply(record_PQ0.125_game,2,median) 
upper_PQ0.125 <- lower_PQ0.125 <- vector()

for (i in 1:720) {
  upper_PQ0.125[i] <- quantile(record_PQ0.125_game[,i],0.975) 
  lower_PQ0.125[i] <- quantile(record_PQ0.125_game[,i],0.025) 
}

##ACT + 0.25mg/kg PQ
PQ <- 0.25
parameters_PQ <- c(ka = 2.27, ktr = 6/0.424, Fm = 0.359, kdr = 1.02, CL_Fp = 20.0, V_Fp =144 , CL_Fc = 1.06, V_Fc = 35.4)
state_PQ <- c(G = 0, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times <- seq(1,dr_time-1,by = 1)
out_PQ <- ode(state_PQ, times, pk_PQ, parameters_PQ)
drug_PQ <- data.frame(out_PQ)
state_PQ1 <- c(G = PQ*60/259.347*1000, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times1 <- seq(dr_time,2400,by = 1)
out_PQ1 <- ode(state_PQ1, times1, pk_PQ, parameters_PQ,method="radau",atol = 1e-9, rtol = 1e-9)
drug_PQ1 <- data.frame(out_PQ1)
PQ_pk<- rbind(drug_PQ,drug_PQ1)
kd_gf <- kd(2400,PQ_pk$D3,84,0.0001,1)
kd_gm <- kd(2400,PQ_pk$D3,84,0.0001,1)
benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0005,0.0022,0.02,0.2, 2.44,PQ,0.0247)

record_PQ0.25_ring <- matrix(0,100,720)
record_PQ0.25_game <- matrix(0,100,720)
for (i in 1:100) {
  simulation_PQ0.25 <- try(sto_model(matrix(0,1,2*42+7),benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2400),kd_ring,benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],kd_DHA,kd_gf,kd_gm,benchmark[11],benchmark[12],benchmark[13],al=42,as=25,t=720),silent = TRUE)
  record_PQ0.25_ring[i,] <- length(ifelse(rowSums(simulation_PQ0.25[,1:25])==0,1,rowSums(simulation_PQ0.25[,1:25])))
  record_PQ0.25_game[i,] <- simulation_PQ0.25[,88]+simulation_PQ0.25[,89]
}
median_PQ0.25 <- apply(record_PQ0.25_game,2,median) 
upper_PQ0.25 <- lower_PQ0.25 <- vector()

for (i in 1:720) {
  upper_PQ0.25[i] <- quantile(record_PQ0.25_game[,i],0.975) 
  lower_PQ0.25[i] <- quantile(record_PQ0.25_game[,i],0.025) 
}


##ACT + 0.5mg/kg PQ
PQ <- 0.5
parameters_PQ <- c(ka = 2.27, ktr = 6/0.424, Fm = 0.359, kdr = 1.02, CL_Fp = 20.0, V_Fp =144 , CL_Fc = 1.06, V_Fc = 35.4)
state_PQ <- c(G = 0, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times <- seq(1,dr_time-1,by = 1)
out_PQ <- ode(state_PQ, times, pk_PQ, parameters_PQ)
drug_PQ <- data.frame(out_PQ)
state_PQ1 <- c(G = PQ*60/259.347*1000, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times1 <- seq(dr_time,2400,by = 1)
out_PQ1 <- ode(state_PQ1, times1, pk_PQ, parameters_PQ,method="radau",atol = 1e-9, rtol = 1e-9)
drug_PQ1 <- data.frame(out_PQ1)
PQ_pk<- rbind(drug_PQ,drug_PQ1)
kd_gf <- kd(2400,PQ_pk$D3,84,0.0001,1)
kd_gm <- kd(2400,PQ_pk$D3,84,0.0001,1)
benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0005,0.0022,0.02,0.2, 2.44,PQ,0.0247)

record_PQ0.5_ring <- matrix(0,100,720)
record_PQ0.5_game <- matrix(0,100,720)
for (i in 1:100) {
  simulation_PQ0.5 <- try(sto_model(matrix(0,1,2*42+7),benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2400),kd_ring,benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],kd_DHA,kd_gf,kd_gm,benchmark[11],benchmark[12],benchmark[13],al=42,as=25,t=720),silent = TRUE)
  record_PQ0.5_ring[i,] <- length(ifelse(rowSums(simulation_PQ0.5[,1:25])==0,1,rowSums(simulation_PQ0.5[,1:25])))
  record_PQ0.5_game[i,] <- simulation_PQ0.5[,88]+simulation_PQ0.5[,89]
}
median_PQ0.5 <- apply(record_PQ0.5_game,2,median) 
upper_PQ0.5 <- lower_PQ0.5 <- vector()

for (i in 1:720) {
  upper_PQ0.5[i] <- quantile(record_PQ0.5_game[,i],0.975) 
  lower_PQ0.5[i] <- quantile(record_PQ0.5_game[,i],0.025) 
}

##ACT + 0.75mg/kg PQ
PQ <- 0.75
parameters_PQ <- c(ka = 2.27, ktr = 6/0.424, Fm = 0.359, kdr = 1.02, CL_Fp = 20.0, V_Fp =144 , CL_Fc = 1.06, V_Fc = 35.4)
state_PQ <- c(G = 0, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times <- seq(1,dr_time-1,by = 1)
out_PQ <- ode(state_PQ, times, pk_PQ, parameters_PQ)
drug_PQ <- data.frame(out_PQ)
state_PQ1 <- c(G = PQ*60/259.347*1000, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times1 <- seq(dr_time,2400,by = 1)
out_PQ1 <- ode(state_PQ1, times1, pk_PQ, parameters_PQ,method="radau",atol = 1e-9, rtol = 1e-9)
drug_PQ1 <- data.frame(out_PQ1)
PQ_pk<- rbind(drug_PQ,drug_PQ1)
kd_gf <- kd(2400,PQ_pk$D3,84,0.0001,1)
kd_gm <- kd(2400,PQ_pk$D3,84,0.0001,1)
benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0005,0.0022,0.02,0.2, 2.44,PQ,0.0247)

record_PQ0.75_ring <- matrix(0,100,720)
record_PQ0.75_game <- matrix(0,100,720)
for (i in 1:100) {
  simulation_PQ0.75 <- try(sto_model(matrix(0,1,2*42+7),benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2400),kd_ring,benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],kd_DHA,kd_gf,kd_gm,benchmark[11],benchmark[12],benchmark[13],al=42,as=25,t=720),silent = TRUE)
  record_PQ0.75_ring[i,] <- length(ifelse(rowSums(simulation_PQ0.75[,1:25])==0,1,rowSums(simulation_PQ0.75[,1:25])))
  record_PQ0.75_game[i,] <- simulation_PQ0.75[,88]+simulation_PQ0.75[,89]
}
median_PQ0.75 <- apply(record_PQ0.75_game,2,median) 
upper_PQ0.75 <- lower_PQ0.75 <- vector()

for (i in 1:720) {
  upper_PQ0.75[i] <- quantile(record_PQ0.75_game[,i],0.975) 
  lower_PQ0.75[i] <- quantile(record_PQ0.75_game[,i],0.025) 
}

median_ring <- apply(record_ACTonly_ring,2,median) 
upper_ring <- lower_ring <- vector()
for (i in 1:720) {
  upper_ring[i] <- quantile(record_ACTonly_ring[,i],0.95) 
  lower_ring[i] <- quantile(record_ACTonly_ring[,i],0.05) 
}

par(mar=c(5,7, 2,1.1))
plot(1:720/24,median_ring/4000,type="n",log = 'y',xlab = "Time post infection (Days)",ylab = "Parasitemia /ml",ylim=c(1e1,1e+11),xlim=c(0,30),cex.axis=2,cex.lab=2,yaxt='n',bty="l")
axis(2,at=c(100,10000,10^6,10^8,10^10),labels = c(expression(10^2),expression(10^4),expression(10^6),expression(10^8),expression(10^10)),cex.axis=1.8)
polygon(polygon(c(1:720/24, rev(1:720)/24), c(upper_ring/4000, rev(lower_ring)/4000),
                col = rgb(255, 100, 100, maxColorValue = 255, alpha = 150), border = NA)  )
lines(1:720/24,median_ring/4000,col="red")
arrows(15,10^9.8,15,10^9.2,length=0.1)
text(15,10^10.5,"Treatment\n on Day 15",cex=1)

plot(1:720/24,median_ACTonly/4000,type="n",log = 'y',xlab = "Time post infection (Days)",ylab = "Gametocytemia /ml",ylim=c(1e1,1e+5),xlim=c(0,30),cex.axis=2,cex.lab=2,yaxt='n',bty="l")
axis(2,at=c(10,100,1000,10^4,10^5),labels = c(expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5)),cex.axis=1.8)

polygon(polygon(c(1:720/24, rev(1:720)/24), c(upper_PQ0.75/4000, rev(lower_PQ0.75)/4000),
                col = rgb(173, 216, 230, maxColorValue = 255, alpha = 150), border = NA)  )
lines(1:720/24,median_PQ0.75/4000,col=rgb(173, 216, 230, maxColorValue = 255))

polygon(polygon(c(1:720/24, rev(1:720)/24), c(upper_PQ0.125/4000, rev(lower_PQ0.125)/4000),
                col = rgb(0, 0, 205, maxColorValue = 255, alpha = 150), border = NA)  )
lines(1:720/24,median_PQ0.125/4000,col=rgb(0, 0, 205, maxColorValue = 255))

polygon(polygon(c(1:720/24, rev(1:720)/24), c(upper_PQ0.25/4000, rev(lower_PQ0.25)/4000),
                col = rgb(30, 144, 255, maxColorValue = 255, alpha = 150), border = NA)  )
lines(1:720/24,median_PQ0.25/4000,col=rgb(30, 144, 255, maxColorValue = 255))

polygon(polygon(c(1:720/24, rev(1:720)/24), c(upper_PQ0.5/4000, rev(lower_PQ0.5)/4000),
                col = rgb(135, 206, 235, maxColorValue = 255, alpha = 150), border = NA)  )
lines(1:720/24,median_PQ0.5/4000,col=rgb(135, 206, 235, maxColorValue = 255))

polygon(polygon(c(1:720/24, rev(1:720)/24), c(upper_ACTonly/4000, rev(lower_ACTonly)/4000),
                col = rgb(0, 0, 139, maxColorValue = 255, alpha = 150)
                , border = NA)  )
lines(1:720/24,median_ACTonly/4000,col=rgb(0, 0, 139, maxColorValue = 255))

legend("topleft",lty=1,lwd=3,bty="n",cex=1.2,col=c(rgb(0, 0, 139, maxColorValue = 255),col=rgb(0, 0, 205, maxColorValue = 255),
                                                   rgb(30, 144, 255, maxColorValue = 255),rgb(135, 206, 235, maxColorValue = 255),rgb(173, 216, 230, maxColorValue = 255)),legend=c("ACT alone","ACT+PQ 0.125mg/kg","ACT+PQ 0.25mg/kg", "ACT+PQ 0.5mg/kg","ACT+PQ 0.75mg/kg"))
arrows(15,10^4.5,15,10^4.1,length=0.1)
text(15,10^4.9,"Treatment\n on Day 15",cex=1)




###sensitivity dose
transmission_model <- function(N, Gm, Gf, V, rho, r){
  M <- rbinom(N,round(Gm,0),V/1000)
  F <- rbinom(N,round(Gf,0),V/1000)
  n <- vector()
  for (i in 1:N) {
    n[i] <- ifelse(min(round(rho*M[i],0),F[i])==0,0,rbinom(1,min(round(rho*M[i],0),F[i]),r))
  }
  record <- cbind(M,F,n)
  return(record)
}

N=100
simulation_record <-list()
ring_record <- list()
game_record <- list()
ring <- matrix(0,100,720)
game <- matrix(0,100,720)

benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0005,0.0022,0.02,0.2, 2.44,0,0.0247)
record <- matrix(0,100,720)
for (k in 1:100) {
  set.seed(k)
  kd_gf <- rep(0,720)
  kd_gm <- rep(0,720)
  simulation <- try(sto_model(matrix(0,1,2*42+7),benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2400),kd_ring,benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],kd_DHA,kd_gf,kd_gm,benchmark[11],benchmark[12],benchmark[13],al=42,as=25,720),silent = TRUE)
  ring[k,] <- ifelse(rowSums(simulation[,1:25])==0,1,rowSums(simulation[,1:25]))
  game[k,] <- simulation[,88]+simulation[,89]
  Phm <- rep(0,720)
  for (j in c(dr_time+12,24+dr_time)) {
    a <- transmission_model(N,simulation[j,88]/4000,simulation[j,89]/4000,0.87,2.4,0.085)
    Phm[j] <- sum(a[,3]>0)/N
  }
  record[k,] <- Phm
  print(k)
}
ring_record[[1]] <- ring
game_record[[1]] <- game
simulation_record[[1]] <- record
PQ_record <- c(0.125,0.25,0.5,0.75)

for (i in 1:4) {
  PQ <- PQ_record[i]
  parameters_PQ <- c(ka = 2.27, ktr = 6/0.424, Fm = 0.359, kdr = 1.02, CL_Fp = 20.0, V_Fp =144 , CL_Fc = 1.06, V_Fc = 35.4)
  state_PQ <- c(G = 0, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
  times <- seq(1,dr_time-1,by = 1)
  out_PQ <- ode(state_PQ, times, pk_PQ, parameters_PQ)
  drug_PQ <- data.frame(out_PQ)
  state_PQ1 <- c(G = PQ*60/259.347*1000, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
  times1 <- seq(dr_time,2400,by = 1)
  out_PQ1 <- ode(state_PQ1, times1, pk_PQ, parameters_PQ,method="radau",atol = 1e-9, rtol = 1e-9)
  drug_PQ1 <- data.frame(out_PQ1)
  PQ_pk<- rbind(drug_PQ,drug_PQ1)
  kd_gf <- kd(t,PQ_pk$D3,84,0.0001,1)
  kd_gm <- kd(t,PQ_pk$D3,84,0.0001,1)
  benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0005,0.0022,0.02,0.2, 2.44,PQ,0.0247)
  record <- matrix(0,100,720)
  for (k in 1:100) {
    set.seed(k)
    simulation <- try(sto_model(matrix(0,1,2*42+7),benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2400),kd_ring,benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],kd_DHA,kd_gf,kd_gm,benchmark[11],benchmark[12],benchmark[13],al=42,as=25,720),silent = TRUE)
    ring[k,] <- ifelse(rowSums(simulation[,1:25])==0,1,rowSums(simulation[,1:25]))
    game[k,] <- simulation[,88]+simulation[,89]
    Phm <- rep(0,720)
    for (j in c(dr_time+12,24+dr_time)) {
      a <- transmission_model(N,simulation[j,88]/4000,simulation[j,89]/4000,0.87,2.4,0.085)
      Phm[j] <- sum(a[,3]>0)/N
    }
    record[k,] <- Phm
    print(k)
  }
  ring_record[[i+1]] <- ring
  game_record[[i+1]] <- game
  simulation_record[[i+1]] <- record
}




result_dose <- matrix(0,1000,3)
result_dose[,1] <- c(rep(0,200),rep(0.125,200),rep(0.25,200),rep(0.5,200),rep(0.75,200))
result_dose[,2] <- c(rep(c(rep(0.5,100),rep(1,100)),5))

for (i in 1:5) {
  result_dose[((i-1)*200+1):((i-1)*200+100),3] <- simulation_record[[i]][,dr_time+24*0.5]
  result_dose[((i-1)*200+101):((i-1)*200+200),3] <- simulation_record[[i]][,dr_time+24*1]
}
result_dose <- as.data.frame(result_dose)
names(result_dose) <- c("Dose","Day","Infectivity")
result_dose$Day <- factor(result_dose$Day)


upper_game <- lower_game <- median_game<- matrix(0,5,720)

for (j in 1:5) {
  for (i in 1:720) {
    median_game[j,i] <- median(game_record[[j]][,i])
    upper_game[j,i] <- quantile(game_record[[j]][,i],0.975) 
    lower_game[j,i] <- quantile(game_record[[j]][,i],0.025) 
  }
}

library(ggplot2)
ggplot(result_dose, aes(x=Dose, y=Infectivity, fill=as.factor(Day), group=interaction(Dose, Day))) + 
  geom_boxplot() + 
  xlab("Dose (mg/kg)") + 
  ylab(expression(p[M])) + 
  labs(fill = "Days post treatment") + 
  theme(
    axis.title.x = element_text(size = 24),
    axis.text.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 24),
    panel.background = element_rect(colour="white", fill="white", linetype="solid"),
    panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey"), 
    panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey"),
    axis.line = element_line(size = 0.5, linetype = "solid", colour = "black")
  ) + 
  scale_fill_manual(
    values = c( "#1E90FF", "#87CEEB")  # light → dark blue
  ) + 
  scale_x_continuous(breaks=seq(0, 1, by=0.25)) + 
  scale_y_continuous(breaks=seq(0, 1, by=0.2), limits=c(0, 0.4))



