library(truncnorm)
##initial parameter settings (increasing infectious mosquitoes (Iv) can speed up convergence to steady state)
H=1000
Sv=9900
Ev=0
Iv=100 ##set more infectious mosquitoes can speed up reaching steady state
theta=0.2 #proportion of asymptomatic infections
Ah <- H*theta

total_bench <- log(10^9)  #threshold of parasitemia for giving treatment
At <- 194*24  #duration of asymptomatic infections (Days)

b=0.4/24 #biting rate
rho = 0.8 #male gametocyte-to-gamete ratio
pf = 0.03 #probability of successful fertilization 
V=3.4 #mosquito blood volume
mu=1/12/24 #death rate of mosquitoes
Lambda =2000*mu #birth rate of mosquitoes
kv=1/12/24 #latent period of mosquitoes
P_mh=0.55 #probability of mosquito-to-human transmission
tmax <- 2.5*365*24 #total time simulated
At <- 194*24
## simulate 100 trajectories before MDA/ changing treatment strategy in case management

source("MDA_ACT.R")

for (l in 1:100) {
  set.seed(500*l)
  ##record matrix
  record_vec <- matrix(0,tmax,3)  #number of S, E and I mosquitoes
  record_vec[1,] <- c(Sv,Ev,Iv)   
  bite_mh <- matrix(0,H,tmax)     #the time since blood stage infection (after the liver stage following the infectious bite)
  treat_ACT <- matrix(0,H,tmax)   #killing rate of ACT (for asexual parsites)
  treat_DHA <- matrix(0,H,tmax)   #killing rate of DHA (for sexual-committed parasites and G1-G4)
  treat_PRQ_m <- matrix(0,H,tmax) #killing rate of primaquine (for G5 male and female gametocytes)
  treat_PRQ_f <- matrix(0,H,tmax)
  record_Ah <- matrix(0,Ah,tmax)  
  infected <- matrix(0,H,tmax)   
  
  treat_ACT[1,1:length(kd_ACT)] <- kd_ACT[1:length(kd_ACT)]
  treat_DHA[1,1:length(kd_DHA)] <- kd_DHA[1:length(kd_DHA)]
  
  record_H_all <- list()
  for (i in 1:H){
    record_H_all[[i]] <- matrix(0,tmax,91)  #record parasitemia in each age bins for each human
  }
  record_H_all[[1]][1,] <- as.numeric(simulation_ACTonly[dr_time,] ) #start with one infected person
  for (t in 2:(24*365)){ 
    ##human to mosquito
    times_hm <- rpois(1,b*record_vec[t-1,1]) #number of susceptible mosquito bites
    idx_hm <- sample.int(H, times_hm, replace = TRUE) #decide who receive these bites
    n_hm <- 0
    for (i in idx_hm) { #For each bite, if at least one developed zygote developed, the mosquito is infected 
      M <- rbinom(1,round(record_H_all[[i]][t-1,88]/4000,0),V/1000)
      F <- rbinom(1,round(record_H_all[[i]][t-1,89]/4000,0),V/1000)
      success <- rbinom(1,round(min(rho*M,F),0),pf)
      n_hm <- n_hm+success 
    }
    ##updates of mosquito population
    record_vec[t,1] <- max(record_vec[t-1,1]+rpois(1,Lambda)-rpois(1,mu*record_vec[t-1,1])-n_hm,0)
    latent_infect <- rpois(1,kv*record_vec[t-1,2])
    record_vec[t,2] <- max(record_vec[t-1,2]+n_hm-rpois(1,mu*record_vec[t-1,2])-latent_infect,0)
    record_vec[t,3] <- max(record_vec[t-1,3]+latent_infect-rpois(1,mu*record_vec[t-1,3]),0)
    
    ##mosquito to human 
    times_mh <- rbinom(1,rpois(1,0.4/24*record_vec[t-1,3]),P_mh) #the number of successful infectious mosquito bites
    if (times_mh>0){
      idx_mh <- sample.int(H, times_mh, replace = TRUE) #decide who receive these bites
      latent_time <- 7*24 #a latent time for liver stage
      for (i in 1:length(idx_mh)) {
        if (idx_mh[i]>(H-Ah)){ #not considered if the person is already in asymptomatic infection, 
          if (latent_time+t+25*24<=tmax & sum(record_Ah[idx_mh[i]-(H-Ah),t:(min(t+240,tmax))])==0){
            record_Ah[idx_mh[i]-(H-Ah),(latent_time+t):min(tmax,latent_time+t+At-1)]<-1:min(At,length((latent_time+t):tmax))
          }
        }else{
          if (latent_time+t<=tmax&idx_mh[i]<=(H-Ah)){ 
            bite_mh[idx_mh[i],latent_time+t] <- bite_mh[idx_mh[i],latent_time+t]+1
          }
        }
      }
    }
    for (i in 1:(H-Ah)) { #update the parasitemia in each bins for each person (symptomatic)
      total <-log(sum(record_H_all[[i]][t-1,c(1:24,61:83,88,89)])/4000) #if the parasitemia>10^9, give treatment
      if (total > total_bench & treat_ACT[i,min(t+2,tmax)]<0.1){
        treat_ACT[i,t:min(tmax,length(kd_ACT)+t-1)] <- treat_ACT[i,t:min(tmax,length(kd_ACT)+t-1)]+kd_ACT[1:min(tmax-t+1,length(kd_ACT))]
        treat_DHA[i,t:min(tmax,length(kd_DHA)+t-1)] <- treat_DHA[i,t:min(tmax,length(kd_DHA)+t-1)]+kd_DHA[1:min(tmax-t+1,length(kd_DHA))]
        }
      simulation <- try(sto_model(record_H_all[[i]][t-1,],bite_mh[i,t]*benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2),treat_ACT[i,(t-1):t],benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],treat_DHA[i,(t-1):t],treat_PRQ_f[i,(t-1):t],treat_PRQ_m[i,(t-1):t],benchmark[11],benchmark[12],benchmark[13],al=42,as=25,2),silent = TRUE)
      record_H_all[[i]][t,] <- as.numeric(simulation[2,])
      infected[i,t] <- ifelse(sum(as.numeric(simulation[2,1:42]))>320,1,0)
    }
    for (i in (H-Ah+1):H) {#asymptomatic
      if (record_Ah[i-H+Ah,t]>0){
        simulation <- try(sto_model(record_H_all[[i]][t-1,],ifelse(record_Ah[i-H+Ah,t]==1,1,0)*benchmark[1]*4000,benchmark[2],benchmark[3],rep(ifelse(record_Ah[i-H+Ah,t]<asym_bench,benchmark[4],benchmark[4]-67.38),2),treat_ACT[i,(t-1):t],benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],treat_DHA[i,(t-1):t],treat_PRQ_f[i,(t-1):t],treat_PRQ_m[i,(t-1):t],benchmark[11],benchmark[12],benchmark[13],al=42,as=25,2),silent = TRUE)
        record_H_all[[i]][t,] <- as.numeric(simulation[2,])
        infected[i,t] <- ifelse(sum(as.numeric(simulation[2,1:42]))>320,1,0)
      }else{
        record_H_all[[i]][t,] <- 0 #If the duration of asymptomatic infection is exceeded, parasites are cleared
      }
    }
    print(t)
  }
  file <- paste0("transmissionMDA",theta,"base",l,".RData")
  save.image(file)
}



##MDA
library(truncnorm)
rm(list = ls())
for (l in 1:100){
  base <- paste0("transmissionMDA",theta,"base",l,".RData")
  load(base)
  p_c <- 0.90 #MDA coverage 
  ad <- 0.5 #adherence
  
  coverage_mat <- matrix(0, H, 3)
  set.seed(l*1000)
  always_n <- round(p_c * H * ad)
  
  if (always_n > 0) {
    always_people <- sample.int(H, always_n, replace = FALSE)
    coverage_mat[always_people, ] <- 1
  } else {
    always_people <- integer(0)
  }
  
  remaining_people <- setdiff(seq_len(H), always_people)
  
  each_n <- round(p_c * H * (1 - ad))
  
  if (each_n > 0) {
    for (j in 1:3) {
      sampled <- sample(remaining_people, each_n, replace = FALSE)
      coverage_mat[sampled, j] <- 1
    }
  }
  
  for (i in 1:H) {
    cov_i <- coverage_mat[i,]
    idx_cov <- sum(coverage_mat[i,])
    
    t1 <- t
    t2 <- t + 720
    t3 <- t + 1440
    
    if (idx_cov==3){
      treat_ACT[i, t1:(t1 + length(kd_MDA_PPQ_three) - 1)] <- kd_MDA_PPQ_three
      
      treat_ACT[i, t1:(t1 + 439)] <- kd_DHA + treat_ACT[i, t1:(t1 + 439)]
      treat_ACT[i, t2:(t2 + 439)] <- kd_DHA + treat_ACT[i, t2:(t2 + 439)]
      treat_ACT[i, t3:(t3 + 439)] <- kd_DHA + treat_ACT[i, t3:(t3 + 439)]
      
      treat_DHA[i, t1:(t1 + 439)] <- kd_DHA
      treat_DHA[i, t2:(t2 + 439)] <- kd_DHA
      treat_DHA[i, t3:(t3 + 439)] <- kd_DHA
      
      treat_PRQ_m[i, t1:(t1 + 439)] <- kd_gm
      treat_PRQ_m[i, t2:(t2 + 439)] <- kd_gm
      treat_PRQ_m[i, t3:(t3 + 439)] <- kd_gm
      
      treat_PRQ_f[i, t1:(t1 + 439)] <- kd_gf
      treat_PRQ_f[i, t2:(t2 + 439)] <- kd_gf
      treat_PRQ_f[i, t3:(t3 + 439)] <- kd_gf

    }else{
      if (idx_cov==1){
        idx_i <- which(cov_i == 1)
        t_i <- t + (idx_i - 1) * 720
        
        treat_ACT[i, t_i:(t_i + length(kd_PPQ) - 1)] <-  kd_PPQ
        treat_ACT[i, t_i:(t_i + 439)] <- kd_DHA + treat_ACT[i, t_i:(t_i + 439)]
        
        treat_DHA[i, t_i:(t_i + 439)] <- kd_DHA
        treat_PRQ_m[i, t_i:(t_i + 439)] <- kd_gm
        treat_PRQ_f[i, t_i:(t_i + 439)] <- kd_gf
        
      }else{
        if (idx_cov==2){
          idx_i <- which(cov_i == 1)
          if (all(diff(idx_i) == 1)) {
            t_i <- t + (idx_i[1] - 1) * 720
            treat_ACT[i, t_i:(t_i + length(kd_MDA_PPQ_two) - 1)] <- kd_MDA_PPQ_two
            for (r in idx_i) {
              t_r <- t + (r - 1) * 720
              treat_ACT[i, t_r:(t_r + 439)] <- kd_DHA + treat_ACT[i, t_r:(t_r + 439)]
              treat_DHA[i, t_r:(t_r + 439)] <- kd_DHA
              treat_PRQ_m[i, t_r:(t_r + 439)] <- kd_gm
              treat_PRQ_f[i, t_r:(t_r + 439)] <- kd_gf
            }
            
          }else{
            for (r in idx_i) {
              t_r <- t + (r - 1) * 720
              treat_ACT[i, t_r:(t_r + length(kd_PPQ) - 1)] <- kd_PPQ
              treat_ACT[i, t_r:(t_r + 439)] <- kd_DHA + treat_ACT[i, t_r:(t_r + 439)]
              treat_DHA[i, t_r:(t_r + 439)] <- kd_DHA
              treat_PRQ_m[i, t_r:(t_r + 439)] <- kd_gm
              treat_PRQ_f[i, t_r:(t_r + 439)] <- kd_gf
            }
            
          }
        }
      }
    }
  }
  benchmark[12] <- 0 #if ACT+0.25 primaquine in MDA, set it to be 0.25
  
  
  for (t in 8761:(2.5*24*365)){
    ##human to mosquito
    times_hm <- rpois(1,0.4/24*record_vec[t-1,1])
    idx_hm <- sample.int(H, times_hm, replace = TRUE)
    n_hm <- 0
    for (i in idx_hm) {
      M <- rbinom(1,round(record_H_all[[i]][t-1,88]/4000,0),V/1000)
      F <- rbinom(1,round(record_H_all[[i]][t-1,89]/4000,0),V/1000)
      success <- rbinom(1,round(min(rho*M,F),0),pf)
      n_hm <- n_hm+success
    }
    ##mosquito population
    record_vec[t,1] <- max(record_vec[t-1,1]+rpois(1,Lambda)-rpois(1,mu*record_vec[t-1,1])-n_hm,0)
    latent_infect <- rpois(1,kv*record_vec[t-1,2])
    record_vec[t,2] <- max(record_vec[t-1,2]+n_hm-rpois(1,mu*record_vec[t-1,2])-latent_infect,0)
    record_vec[t,3] <- max(record_vec[t-1,3]+latent_infect-rpois(1,mu*record_vec[t-1,3]),0)
    
    ##mosquito to human 
    times_mh <- rbinom(1,rpois(1,b*record_vec[t-1,3]),P_mh)
    if (times_mh>0){
      idx_mh <- sample.int(H, times_mh, replace = TRUE)
      latent_time <- 7*24
      for (i in 1:length(idx_mh)) {
        if (idx_mh[i]>(H-Ah)){
          if (latent_time+t+25*24<=tmax & sum(record_Ah[idx_mh[i]-(H-Ah),t:(min(t+240,tmax))])==0){
            record_Ah[idx_mh[i]-(H-Ah),(latent_time+t):min(tmax,latent_time+t+At-1)]<-1:min(At,length((latent_time+t):tmax))
          }
        }else{
          if (latent_time+t<=tmax&idx_mh[i]<=(H-Ah)){
            bite_mh[idx_mh[i],latent_time+t] <- bite_mh[idx_mh[i],latent_time+t]+1
          }
        }
      }
    }
    
    for (i in 1:(H-Ah)) {
      total <-log(sum(record_H_all[[i]][t-1,c(1:24,61:83,88,89)])/4000)
      #total <-sum(record_H_all[[i]][t-1,c(88,89)])/4000
      if (total > total_bench & treat_ACT[i,min(t+26,tmax)]<0.1){
        treat_ACT[i,t:min(tmax,length(kd_ACT)+t-1)] <- treat_ACT[i,t:min(tmax,length(kd_ACT)+t-1)]+kd_ACT[1:min(tmax-t+1,length(kd_ACT))]
        treat_DHA[i,t:min(tmax,length(kd_DHA)+t-1)] <- treat_DHA[i,t:min(tmax,length(kd_DHA)+t-1)]+kd_DHA[1:min(tmax-t+1,length(kd_DHA))]
        }
      simulation <- try(sto_model(record_H_all[[i]][t-1,],bite_mh[i,t]*benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2),treat_ACT[i,(t-1):t],benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],treat_DHA[i,(t-1):t],treat_PRQ_f[i,(t-1):t],treat_PRQ_m[i,(t-1):t],benchmark[11],benchmark[12],benchmark[13],al=42,as=25,2),silent = TRUE)
      record_H_all[[i]][t,] <- as.numeric(simulation[2,])
      infected[i,t] <- ifelse(sum(as.numeric(simulation[2,1:42]))>320,1,0)
    }
    for (i in (H-Ah+1):H) {
      if (record_Ah[i-H+Ah,t-1]==At){
        record_H_all[[i]][t,] <- 0
      }else{
        simulation <- try(sto_model(record_H_all[[i]][t-1,],ifelse(record_Ah[i-H+Ah,t]==1,1,0)*benchmark[1]*4000,benchmark[2],benchmark[3],rep(ifelse(record_Ah[i-H+Ah,t]<asym_bench,benchmark[4],benchmark[4]-67.38),2),treat_ACT[i,(t-1):t],benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],treat_DHA[i,(t-1):t],treat_PRQ_f[i,(t-1):t],treat_PRQ_m[i,(t-1):t],benchmark[11],benchmark[12],benchmark[13],al=42,as=25,2),silent = TRUE)
        record_H_all[[i]][t,] <- as.numeric(simulation[2,])
        infected[i,t] <- ifelse(sum(as.numeric(simulation[2,1:42]))>320,1,0)
        if (sum(as.numeric(simulation[2,1:42]))==0&record_Ah[i-H+Ah,t]>0){
          record_Ah[i-H+Ah,(t+1):(t+At-record_Ah[i-H+Ah,t])] <- 0
        }
      }
    }
    print(t)
    if (sum(infected[,t])==0&sum(record_vec[t,2:3])==0){
      break
    }
  }
  file <- paste0("transmissionMDA,",theta,"c0.9ACT_",l,".RData")
  save.image(file)
  file_HI <- paste0("infectedMDA,",theta,"c0.9ACT_",l,".csv")
  file_IV <- paste0("infectiousMosMDA,",theta,",c0.9ACT_",l,".csv")
  write.csv(infected[,1:1825*12],file_HI,row.names = FALSE)
  write.csv(record_vec[1:1825*12,c(2,3)],file_IV,row.names=FALSE)
  game_record_m <- matrix(0,H,21900)
  game_record_f <- matrix(0,H,21900)
  for (z in 1:H) {
    game_record_m[z,] <- record_H_all[[z]][,88]
    game_record_f[z,] <- record_H_all[[z]][,89]
  }
  file <- paste0("gameMDA_ACT_m",l,".csv")
  write.csv(game_record_m[,1:1825*12],file,row.names = FALSE)
  file <- paste0("gameMDA_ACT_f",l,".csv")
  write.csv(game_record_f[,1:1825*12],file,row.names = FALSE)
}

H_ACT <- matrix(0,100,1825)
VI_ACT <- matrix(0,100,1825)
VE_ACT <- matrix(0,100,1825)
for (l in 1:100) {
  name_H <- paste0("infectedMDA",theta,"c0.9ACT_",l,".csv")
  name_V <- paste0("infectiousMosMDA",theta,"c0.9ACT_",l,".csv")
  file_H <- read.csv(name_H)
  file_V <- read.csv(name_V)
  H_ACT[l,] <- colSums(file_H)
  VI_ACT[l,] <- file_V[,2]
  VE_ACT[l,] <- file_V[,1]
  print(l)
}
write.csv(H_ACT,"H_ACT.csv",row.names = FALSE)
write.csv(VI_ACT,"VI_ACT.csv",row.names = FALSE)
write.csv(VE_ACT,"VE_ACT.csv",row.names = FALSE)
H_PQ <- matrix(0,100,1825)
VI_PQ <- matrix(0,100,1825)
VE_PQ <- matrix(0,100,1825)
for (l in 1:100) {
  name_H <- paste0("infectedMDA",theta,"c0.9PQ_",l,".csv")
  name_V <- paste0("infectiousMosMDA",theta,"c0.9PQ_",l,".csv")
  file_H <- read.csv(name_H)
  file_V <- read.csv(name_V)
  H_PQ[l,] <- colSums(file_H)
  VI_PQ[l,] <- file_V[,2]
  VE_PQ[l,] <- file_V[,1]
  print(l)
}
write.csv(H_PQ,"H_PQ.csv",row.names = FALSE)
write.csv(VI_PQ,"VI_PQ.csv",row.names = FALSE)
write.csv(VE_PQ,"VE_PQ.csv",row.names = FALSE)




##routine case management
for (l in 1:100) {
  set.seed(500*l)
  base <- paste0("transmissionMDA",theta,"base",l,".RData")
  load(base)
  benchmark[12] <- 0.25 #if change to ACT+0.25PQ, set it to be 0.25
  for (t in (24*365+1):(2.5*365*24)){
    ##human to mosquito
    times_hm <- rpois(1,b*record_vec[t-1,1])
    idx_hm <- sample.int(H, times_hm, replace = TRUE)
    n_hm <- 0
    for (i in idx_hm) {
      M <- rbinom(1,round(record_H_all[[i]][t-1,88]/4000,0),V/1000)
      F <- rbinom(1,round(record_H_all[[i]][t-1,89]/4000,0),V/1000)
      success <- rbinom(1,round(min(rho*M,F),0),pf)
      n_hm <- n_hm+success
    }
    ##mosquito population
    record_vec[t,1] <- max(record_vec[t-1,1]+rpois(1,Lambda)-rpois(1,mu*record_vec[t-1,1])-n_hm,0)
    latent_infect <- rpois(1,kv*record_vec[t-1,2])
    record_vec[t,2] <- max(record_vec[t-1,2]+n_hm-rpois(1,mu*record_vec[t-1,2])-latent_infect,0)
    record_vec[t,3] <- max(record_vec[t-1,3]+latent_infect-rpois(1,mu*record_vec[t-1,3]),0)
    
    ##mosquito to human 
    times_mh <- rbinom(1,rpois(1,0.4/24*record_vec[t-1,3]),P_mh)
    if (times_mh>0){
      idx_mh <- sample.int(H, times_mh, replace = TRUE)
      latent_time <- 7*24
      for (i in 1:length(idx_mh)) {
        if (idx_mh[i]>(H-Ah)){
          if (latent_time+t+25*24<=tmax & sum(record_Ah[idx_mh[i]-(H-Ah),t:(min(t+240,tmax))])==0){
            record_Ah[idx_mh[i]-(H-Ah),(latent_time+t):min(tmax,latent_time+t+At-1)]<-1:min(At,length((latent_time+t):tmax))
          }
        }else{
          if (latent_time+t<=tmax&idx_mh[i]<=(H-Ah)){
            bite_mh[idx_mh[i],latent_time+t] <- bite_mh[idx_mh[i],latent_time+t]+1
          }
        }
      }
    }
    for (i in 1:(H-Ah)) {
      total <-log(sum(record_H_all[[i]][t-1,c(1:24,61:83,88,89)])/4000)
      if (total > total_bench & treat_ACT[i,min(t+2,tmax)]<0.1){
        treat_ACT[i,t:min(tmax,length(kd_ACT)+t-1)] <- treat_ACT[i,t:min(tmax,length(kd_ACT)+t-1)]+kd_ACT[1:min(tmax-t+1,length(kd_ACT))]
        treat_DHA[i,t:min(tmax,length(kd_DHA)+t-1)] <- treat_DHA[i,t:min(tmax,length(kd_DHA)+t-1)]+kd_DHA[1:min(tmax-t+1,length(kd_DHA))]
        treat_PRQ_m[i,t:min(tmax,length(kd_gm)+t-1)] <- treat_PRQ_m[i,t:min(tmax,length(kd_gm)+t-1)]+kd_gm[1:min(tmax-t+1,length(kd_gm))]
        treat_PRQ_f[i,t:min(tmax,length(kd_gf)+t-1)] <- treat_PRQ_f[i,t:min(tmax,length(kd_gf)+t-1)]+kd_gf[1:min(tmax-t+1,length(kd_gf))]
      }
      simulation <- try(sto_model(record_H_all[[i]][t-1,],bite_mh[i,t]*benchmark[1]*4000,benchmark[2],benchmark[3],rep(benchmark[4],2),treat_ACT[i,(t-1):t],benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],treat_DHA[i,(t-1):t],treat_PRQ_f[i,(t-1):t],treat_PRQ_m[i,(t-1):t],benchmark[11],benchmark[12],benchmark[13],al=42,as=25,2),silent = TRUE)
      record_H_all[[i]][t,] <- as.numeric(simulation[2,])
      infected[i,t] <- ifelse(sum(as.numeric(simulation[2,1:42]))>320,1,0)
    }
    for (i in (H-Ah+1):H) {
      if (record_Ah[i-H+Ah,t]>0){
        simulation <- try(sto_model(record_H_all[[i]][t-1,],ifelse(record_Ah[i-H+Ah,t]==1,1,0)*benchmark[1]*4000,benchmark[2],benchmark[3],rep(ifelse(record_Ah[i-H+Ah,t]<asym_bench,benchmark[4],benchmark[4]-67.38),2),treat_ACT[i,(t-1):t],benchmark[5],benchmark[6],benchmark[7],benchmark[8],benchmark[8],benchmark[9],benchmark[10],treat_DHA[i,(t-1):t],treat_PRQ_f[i,(t-1):t],treat_PRQ_m[i,(t-1):t],benchmark[11],benchmark[12],benchmark[13],al=42,as=25,2),silent = TRUE)
        record_H_all[[i]][t,] <- as.numeric(simulation[2,])
        infected[i,t] <- ifelse(sum(as.numeric(simulation[2,1:42]))>320,1,0)
      }else{
        record_H_all[[i]][t,] <- 0
      }
    }
    print(t)
  }
  file <- paste0("transmissionTS",theta,"PQ_",l,".RData")
  save.image(file)
  file_HI <- paste0("infectedTS,",theta,"PQ_",l,".csv")
  file_IV <- paste0("infectiousMosTS,",theta,",PQ_",l,".csv")
  write.csv(infected[,1:1825*12],file_HI,row.names = FALSE)
  write.csv(record_vec[1:1825*12,c(2,3)],file_IV,row.names=FALSE)
  game_record_m <- matrix(0,H,21900)
  game_record_f <- matrix(0,H,21900)
  for (z in 1:H) {
    game_record_m[z,] <- record_H_all[[z]][,88]
    game_record_f[z,] <- record_H_all[[z]][,89]
  }
  file <- paste0("gameMDA_PQ_m",l,".csv")
  write.csv(game_record_m[,1:1825*12],file,row.names = FALSE)
  file <- paste0("gameMDA_PQ_f",l,".csv")
  write.csv(game_record_f[,1:1825*12],file,row.names = FALSE)
}

average_game <- matrix(0,100,1825)
for (l in 1:100) {
  name_m <- paste0("gameMDA_PQ_m",l,".csv")
  name_f <- paste0("gameMDA_PQ_f",l,".csv")
  file_m <- read.csv(name_m)
  file_f <- read.csv(name_f)
  total_game <- file_m+file_f
  average_game[l,] <- colSums(total_game)/H/4000
  print(l)
}
write.csv(average_game,"G_PQ.csv",row.names = FALSE)


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


library(foreach)
library(doParallel)
cl <- parallel::makeCluster(8)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)


for (l in 1:100) {
  name_m <- paste0("gameMDA_PQ_m",l,".csv")
  name_f <- paste0("gameMDA_PQ_f",l,".csv")
  file_m <- read.csv(name_m)
  file_f <- read.csv(name_f)
  for (j in 1:1000) {
    result <- foreach (i=1:1825) %dopar% {
      if (file_m[j,i]==0){
        simulation=0
      }else{
        simulation <- transmission_model(100000,file_m[j,i]/4000,file_f[j,i]/4000,3.4,0.8,0.03)
      }
      simulation
    }
    for (i in 1:1825) {
      pm[j,i] <- result[[i]]
    }
    print(j)
  }
  file_pm <- paste0("pm_PQ",l,".csv")
  write.csv(pm,file_pm,row.names = FALSE)
}
