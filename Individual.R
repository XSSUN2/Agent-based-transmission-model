sto_model <- function(stage,N,mu,sd,rp,kd_ring,f,detp,detg,detgm,detgf,m,fg,kd_g,kd_gf,kd_gm,slop,dose,detgcl,al=42,as=25,t){
  rec <- matrix(0,t,(2*al+7))
  if(N!=0){
    pinit <- rtruncnorm(max(round(N,0),1),1,al,mu,sd) 
    pinit <- round(pinit,digits = 0)
    for (i in 1:al) { #record initial number of parasites
      rec[1,i] <- length(which(pinit==i)) 
    }
    rec[1,] <- rec[1,]+stage
  }else{
    rec[1,] <- stage
  }
  rec_f <- vector()
  rec_f[1] <- 0
  rec_fg <- vector()
  rec_fg[1] <- 0
  
  for (i in 2:t) { 
    for (j in 1:al) { #first compute asexual parasites and in this part drug is effective
      if (j==1){ #age = 1 hour
        # number of parasites in state 48 last time * replication rate - number of death (detp--death rate ,kd--killing effect )
        # Here kd will be explained in the PQP file.
        p.a.t <- (rec[i-1,al]-rpois(1,(detp+kd_ring[i])*rec[i-1,al]))*rp[i]
        # make sure it's not negative
        p.a.t <- max(0,round(p.a.t,0))
        rec[i,j] <- p.a.t
      }
      else{if (j == as+1){
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
    for (j in (al+1):(2*al-1)){ # 47 states for sexual development
      #column al+1 is Pg(as+1), al+2 is Pg(as+2), ... 
      if (j==(2*al-as+1)){ # age 1 of the sexually committed parasites
        # replication 
        p.a.t <- (rec[i-1,2*al-as]-rpois(1,(detp+kd_g[i])*rec[i-1,2*al-as]))*rp[i]
        p.a.t <- max(0,round(p.a.t,0))
        rec[i,j] <- p.a.t
      }else{if (j==(al+1)){ 
        #fraction
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
    #random a number to be parasites going from g1 to g2, compare with number in g1
    #to make sure g1 has enough parasites to go to g2.
    #the maturation rate is m
    g1tog2 <- min(rec[i-1,2*al],rpois(1,m*rec[i-1,2*al]))
    #all the parasites in (Pg(as-1) - nature death (death rate is detp)) will go to state g1 
    #+ parasites in g1 last time - nature death (death rate is detg)
    #- parasites goes to g2
    #= number in g1 state
    rec[i,2*al] <- max(0,rec[i-1,2*al]+rec[i-1,2*al-1]-rpois(1,(detg+kd_g[i]+detg)*rec[i-1,2*al])-rpois(1,detp*rec[i-1,2*al-1])-g1tog2)
    
    #next g2, g3, g4 
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
times <- seq(1,dr_time-1,by=1)
out <- ode(state, times, odefunction, parameters)
drug <- data.frame(out)
state1 <- c(D = 960, T1 = 0, T2 = 0, C = 0, P1 = 0, P2 = 0)
times1 <- seq(dr_time,dr_time+23,by=1)
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
PQ<- rbind(drug,drug1,drug2,drug3)

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
times <- seq(1,dr_time-1,by=1)
out <- ode(state, times, ode_DHA, parameters)
drug <- data.frame(out)
state1 <- c(D = 120, T1 = 0, T2 = 0, T3=0, T4=0, A = 0, DHA = 0, C = 0)
times1 <- seq(dr_time,dr_time+23,by=1)
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

##kd_ring
t=2400
kd_PQ <- kd(t,PQ$C,0.2,9,12)

kd_DHA <- kd(t,DHA$C/284.35/124*10^6,0.273,34.9,1)

kd_ring <- kd_PQ+kd_DHA

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

###ACT only
PRQ=0
kd_gf = kd_gm = rep(0,2400)
benchmark <- c(0.08,7,8,71,0.006,0.03,0.0005,0.0022,0.02,0.2, 2.44,PRQ,0.0247)
set.seed(5)
simulation_ACTonly <- try(sto_model(matrix(0,1,2*42+7),benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2400),kd_ring,benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],kd_DHA,kd_gf,kd_gm,benchmark[11],benchmark[12],benchmark[13],al=42,as=25,t=1200),silent = TRUE)

ring <- ifelse(rowSums(simulation_ACTonly[,1:24])==0,1,rowSums(simulation_ACTonly[,1:24]))
game <- simulation_ACTonly[,88]+simulation_ACTonly[,89]


###ACT+0.25 primaquine (PRQ is primaquine in the code)


PRQ <- 0.25
parameters_PRQ <- c(ka = 2.27, ktr = 6/0.424, Fm = 0.359, kdr = 1.02, CL_Fp = 20.0, V_Fp =144 , CL_Fc = 1.06, V_Fc = 35.4)
state_PRQ <- c(G = 0, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times <- seq(1,dr_time-1,by=1)
out_PRQ <- ode(state_PRQ, times, pk_PRQ, parameters_PRQ)
drug_PRQ <- data.frame(out_PRQ)
state_PRQ1 <- c(G = PRQ*60, T1 = 0, T2 = 0,T3 = 0, T4 = 0, T5 = 0, T6 = 0, P = 0, C = 0, D1 = 0, D2 = 0, D3 = 0)
times1 <- seq(dr_time,2400,by=1)
out_PRQ1 <- ode(state_PRQ1, times1, pk_PRQ, parameters_PRQ,method="radau",atol = 1e-9, rtol = 1e-9)
drug_PRQ1 <- data.frame(out_PRQ1)
PRQ_pk<- rbind(drug_PRQ,drug_PRQ1)

kd_gf <- kd(2400,PRQ_pk$D3*1000/144/259.347,84,0.0001,1)[-(1:360)]
kd_gm <- kd(2400,PRQ_pk$D3*1000/144/259.347,84,0.0001,1)[-(1:360)]
kd_DHA_PQ <- kd_DHA[-(1:360)]
kd_ring_PQ <- kd_ring[-(1:360)]



benchmark <- c(0.08,7,5.5,71,0.006,0.03,0.0005,0.0022,0.02,0.2, 2.44,PRQ,0.0247)
set.seed(5)

simulation_PQ <- try(sto_model(as.numeric(simulation_ACTonly[360,]),0,benchmark[2],benchmark[3],rep(benchmark[4],2400),kd_ring_PQ,benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],kd_DHA_PQ,kd_gf,kd_gm,benchmark[11],benchmark[12],benchmark[13],al=42,as=25,t=840),silent = TRUE)
ring_PQ <- ifelse(rowSums(simulation_PQ[,1:24])==0,1,rowSums(simulation_PQ[,1:24]))
game_PQ <- simulation_PQ[,88]+simulation_PQ[,89]

transmission_model <- function(N, Gm, Gf, V, rho, r){
  M <- rbinom(N,round(Gm,0),V/1000)
  F <- rbinom(N,round(Gf,0),V/1000)
  n <- vector()
  for (i in 1:N) {
    n[i] <- ifelse(min(round(rho*M[i],0),F[i])==0,0,rbinom(1,min(round(rho*M[i],0),F[i]),r))
  }
  record <- sum(n>0)/N
  return(record)
}

trans_prob <- matrix(0,2,960)

game_male_ACTonly <- simulation_ACTonly[241:1200,88]/4000
game_female_ACTonly <- simulation_ACTonly[241:1200,89]/4000

for (i in 1:960) {
  trans_prob[1,i] <- transmission_model(100000,game_male_ACTonly[i],game_female_ACTonly[i],3.4,0.8,0.03)
  print(i)
}
game_male_PQ <- simulation_PQ[,88]/4000
game_female_PQ <- simulation_PQ[,89]/4000
for (i in 1:840) {
  trans_prob[2,i+120] <- transmission_model(100000,game_male_PQ[i],game_female_PQ[i],3.4,0.8,0.03)
  print(i)
}
trans_prob[2,1:120] <- trans_prob[1,1:120]





####Figure 2
col_act <- "#4C78A8"
col_pq <- "#3FAE87"

pdf("Fig_htom.pdf", width = 14, height = 4)
layout(
  matrix(c(1, 2 , 3), nrow = 1, byrow = TRUE),
  widths = c(1, 1,1 )
)
par(
  oma = c(1.2, 2.2, 1.0, 1.2),
  mgp = c(3, 1, 0),
  tcl = -0.3
)


par(mar = c(4.2, 7, 2.5, 1.0))
plot(1:2400/24,DHA$DHA/124*10^3,type="l",ylab="Drug concentration (ng/mL)",
     xlab="Time post infection (Days)",cex.axis=2,cex.lab=2,cex.main=2,xlim=c(0,30),
     lwd=1.4,col="#A65628",lty=1,bty="l",yaxt="n",ylim=c(0,500))
mtext("a", side = 2, line = 5, at = 580, las = 1, font = 2, cex = 1.8)
axis(2,at=c(0,250,500),cex.axis=2)
lines(1:2400/24,PQ$C*10^6/804,col="#FF7F00",lwd=1.5)
lines(1:2400/24,PRQ_pk$P*1000/144,col="#8B36A4",lwd=1.5)
legend("topright",bty="n",lwd=1.5,cex=1.2,col=c("#A65628","#FF7F00", "#8B36A4"),legend=c("DHA","PPQ","PQ"))
arrows(15,440.3,15,410,length=0.1)
text(15,480,"Administration\n on Day 15",cex=1.4)

par(mar = c(4.2, 7, 2.5, 1.0))
plot(1:1200/24,ring/4000,type="n",log = 'y',xlab="Time post infection (Days)",ylab = "Parasitemia /ml",ylim=c(1e0,10^16),xlim=c(0,30),cex.axis=2,cex.lab=2,yaxt='n',bty="l")
mtext("b", side = 2, line = 5, at = 10^(16*1.16), las = 1, font = 2, cex = 1.8)
axis(2,at=c(1,10^4,10^8,10^12,10^16),labels = c(expression(10^0),expression(10^4),expression(10^8),expression(10^12),expression(10^16)),cex.axis=2)
lines(1:1200/24,ring/4000,lwd=1.2)

lines(1:1200/24,game/4000,col=col_act,lwd=1.5)
lines(1:1200/24,c(game[1:360]/4000,game_PQ/4000),col="#3FAE87C8",lwd=1.5)
legend("topright",lty=1,lwd=1.5,bty="n",cex=1.2,col=c(1,col_act,"#3FAE87C8"),legend=c("Ring-stage parasitemia","Gametocytemia\n ACT alone","Gametocytemia\n ACT+0.25mg/kg PQ"),y.intersp = 1.5)



par(mar = c(4.2, 7, 2.5, 1.0))
plot((1:150)*8/24,c(rep(0,30),trans_prob[1,(1:120)*8]),type="l",ylim=c(0,0.3),ylab=expression(Infectiousness~(P[M])),xlab="Time post infection (Days)",cex.axis=1.8,cex.lab=2,xlim=c(0,30),
     lwd=1.4,col=col_act,lty=1,yaxt="n",bty="l")
mtext("c", side = 2, line = 5, at =0.3*1.16, las = 1, font = 2, cex = 1.8)

box(bty="l",lwd=1.2)
lines((1:150)*8/24,c(rep(0,30),trans_prob[2,(1:120)*8]),col="#3FAE87C8",lwd=1.5)
axis(2,at=c(0,0.1,0.2,0.3),cex.axis=2)
legend("topright",bty="n",lwd=1.5,cex=1.2,col=c(col_act,"#3FAE87C8"),legend=c("ACT alone","ACT+0.25 mg/kg PQ"))

par(
  fig = c(0.09, 0.2, 0.25, 0.7),
  new = TRUE,
  mar = c(2.0, 0.3, 0.3, 2.0),
  oma = c(0, 0, 0, 0)
)

plot(1:2400/24,DHA$DHA/124*10^3, type="l", lwd=1.5, col="#A65628",
     xlim=c(14.9,15.6), ylim=c(0,55),
     xlab="", ylab="", axes=FALSE)
lines(1:2400/24,PQ$C*10^6/804, lwd=1.5, col="#FF7F00")
lines(1:2400/24,PRQ_pk$P*1000/144, lwd=1.5, col="#8B36A4")
box()
axis(1, cex.axis=1.5,at=c(15,15.5))
axis(4, cex.axis=1.5,las=1,at=c(0,50))


dev.off()

system("pdfcrop Fig_htom.pdf Fig_htom.pdf")


