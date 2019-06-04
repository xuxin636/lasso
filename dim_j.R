library(mvtnorm)
cond <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
E <- as.matrix(read.csv("/rigel/home/xx2319/lasso/0-1categoryofdata.csv"))
ww <- 40
w <- E[,2:(ww+1)]
J = ncol(w)
N = nrow(w)
K = 3
response <- w;
###my code###
Q <- matrix(1,J,K);Q[J,2:3] <- 0;Q[(J-1),3] <- 0;
##initial value###
A_initial <- matrix(0,J,K);A_initial[,1] <- runif(J,1,2);A_initial[,2] <- runif(J,1,2);A_initial[,3] <- runif(J,1,2);
A_initial <- A_initial*Q;
d_initial <- rnorm(J,0,1);D_initial <- cbind(d_initial,A_initial);
KK <- 20;theta_min <- -4;theta_max <- 4;mm1 <- seq(theta_min,theta_max,(theta_max-theta_min)/KK);mm <- mm1[-1]
THETA_tuta <- matrix(0,nrow=KK*KK*KK,ncol=3);THETA_tuta[,3] <- rep(mm,KK*KK);
THETA_tuta[,2] <-rep(c(rep(1,KK)%*%t(mm)),KK);THETA_tuta[,1] <-c(rep(1,KK*KK)%*%t(mm))#针对K <- 3的theta分块,获取theta的分块
THETA_tuta <- cbind(rep(1,nrow(THETA_tuta)),THETA_tuta)
theta_tmp <- rowSums(THETA_tuta[,2:4]*THETA_tuta[,2:4])/2
xx <- seq(0.001,0.03,0.001);xx1 <- matrix(0,nrow = length(xx)*length(xx),ncol=2);xx1[,2] <- rep(xx,length(xx));xx1[,1] <- c(rep(1,length(xx))%*%t(xx))
lammda <- c(rep(xx1[cond,1],J/2),rep(xx1[cond,2],J/2))*N;
soft <- function(a,b){
  if(a>0&a>b){return(a-b)}
  else{return(0)}
}
response <- t(response);
ll <- mm[2]-mm[1]
A_0 <- t(D_initial)
cc <- exp(THETA_tuta%*%A_0%*%response-theta_tmp)/apply(1+exp((THETA_tuta%*%A_0)),1,prod);
theta_post <- sweep(cc, 2, colSums(cc), "/") 
ss_1 <- colSums(exp(log(exp(THETA_tuta%*%A_0)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%A_0)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%A_0)))%*%(1-response)
)) 
log_lik_1 <-sum(log(ss_1))
likelihood_0 <- log_lik_1-sum(sweep(A_0,2,lammda,"*"));
for(j in 1:J){
  ss_0 <- ss_1
  A_grad <- sum(sweep(THETA_tuta[,1]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,1]*theta_post* c(1/(1+exp(-THETA_tuta%*%A_0[,j]))))
  A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%A_0[,j])/(1+exp(THETA_tuta%*%A_0[,j]))/(1+exp(THETA_tuta%*%A_0[,j]))))*(THETA_tuta[,1]*THETA_tuta[,1]))
  d_tuta <- A_0[1,j]-A_grad/A_grad_2
  aa <- A_0;aa[1,j] <- d_tuta;
  
  ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
  ))
  eps <- sum(log(ss_1/ss_0))
  while(eps>5){
    ss_0 <- ss_1
    A_grad <- sum(sweep(THETA_tuta[,1]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,1]*theta_post* c(1/(1+exp(-THETA_tuta%*%aa[,j]))))
    A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%aa[,j])/(1+exp(THETA_tuta%*%aa[,j]))/(1+exp(THETA_tuta%*%aa[,j]))))*(THETA_tuta[,1]*THETA_tuta[,1]))
    d_tuta <- d_tuta-A_grad/A_grad_2
    aa <- A_0;aa[1,j] <- d_tuta;
    ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
    ))
    eps <- sum(log(ss_1/ss_0))
  }
  
  A_0[1,j] <- d_tuta;
  A_grad <- sum(sweep(THETA_tuta[,2]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,2]*theta_post* c(1/(1+exp(-THETA_tuta%*%A_0[,j]))))  
  A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%A_0[,j])/(1+exp(THETA_tuta%*%A_0[,j]))/(1+exp(THETA_tuta%*%A_0[,j]))))*(THETA_tuta[,2]*THETA_tuta[,2]))
  A1_tuta <-A_0[2,j]-A_grad/A_grad_2
  aa <- A_0;aa[2,j] <- A1_tuta;
  ss_0 <- ss_1 
  ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
  ))
  eps <- sum(log(ss_1/ss_0))
  
  while(eps>5){
    ss_0 <- ss_1
    A_grad <- sum(sweep(THETA_tuta[,2]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,2]*theta_post* c(1/(1+exp(-THETA_tuta%*%aa[,j]))))
    A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%aa[,j])/(1+exp(THETA_tuta%*%aa[,j]))/(1+exp(THETA_tuta%*%aa[,j]))))*(THETA_tuta[,2]*THETA_tuta[,2]))
    A1_tuta <-A1_tuta-A_grad/A_grad_2
    aa <- A_0;aa[2,j] <- A1_tuta;
    ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
    ))
    eps <- sum(log(ss_1/ss_0))
  }
  
  A_0[2,j]<-soft(A1_tuta,-lammda[j]/A_grad_2)
  A_grad <- sum(sweep(THETA_tuta[,3]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,3]*theta_post* c(1/(1+exp(-THETA_tuta%*%A_0[,j]))))  
  A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%A_0[,j])/(1+exp(THETA_tuta%*%A_0[,j]))/(1+exp(THETA_tuta%*%A_0[,j]))))*(THETA_tuta[,3]*THETA_tuta[,3]))
  A2_tuta <-A_0[3,j]-A_grad/A_grad_2
  aa <- A_0;aa[3,j] <- A2_tuta;
  ss_0 <- ss_1
  ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
  ))
  eps <- sum(log(ss_1/ss_0))
  
  while(eps>5){
    ss_0 <- ss_1
    A_grad <- sum(sweep(THETA_tuta[,3]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,3]*theta_post* c(1/(1+exp(-THETA_tuta%*%aa[,j]))))
    A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%aa[,j])/(1+exp(THETA_tuta%*%aa[,j]))/(1+exp(THETA_tuta%*%aa[,j]))))*(THETA_tuta[,3]*THETA_tuta[,3]))
    A2_tuta <-A2_tuta-A_grad/A_grad_2
    aa <- A_0;aa[3,j] <- A2_tuta;
    ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
    ))
    eps <- sum(log(ss_1/ss_0))
  }
  
  A_0[3,j] <- soft(A2_tuta,-lammda[j]/A_grad_2)
  A_grad <- sum(sweep(THETA_tuta[,4]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,4]*theta_post* c(1/(1+exp(-THETA_tuta%*%A_0[,j]))))  
  A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%A_0[,j])/(1+exp(THETA_tuta%*%A_0[,j]))/(1+exp(THETA_tuta%*%A_0[,j]))))*(THETA_tuta[,4]*THETA_tuta[,4]))
  A3_tuta <-A_0[4,j]-A_grad/A_grad_2
  aa <- A_0;aa[4,j] <- A3_tuta;
  ss_0 <- ss_1 
  ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
  ))
  eps <- sum(log(ss_1/ss_0))
  
  while(eps>5){
    ss_0 <- ss_1
    A_grad <- sum(sweep(THETA_tuta[,4]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,4]*theta_post* c(1/(1+exp(-THETA_tuta%*%aa[,j]))))
    A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%aa[,j])/(1+exp(THETA_tuta%*%aa[,j]))/(1+exp(THETA_tuta%*%aa[,j]))))*(THETA_tuta[,4]*THETA_tuta[,4]))
    A3_tuta <-A3_tuta-A_grad/A_grad_2
    aa <- A_0;aa[4,j] <- A3_tuta;
    ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
    ))
    eps <- sum(log(ss_1/ss_0))
  }
  
  A_0[4,j] <- soft(A3_tuta,-lammda[j]/A_grad_2)
  
}


A_0 <- A_0*t(cbind(rep(1,J),Q))
ss_1 <- colSums(exp(log(exp(THETA_tuta%*%A_0)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%A_0)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%A_0)))%*%(1-response)
)) 
log_lik_1 <-sum(log(ss_1))
likelihood_1 <- log_lik_1-sum(lammda*A_0[3:4,]);
eps <- likelihood_1-likelihood_0

while(eps>5){
  likelihood_0 <- likelihood_1;
  cc <- exp(THETA_tuta%*%A_0%*%response-theta_tmp)/apply(1+exp((THETA_tuta%*%A_0)),1,prod);
  theta_post <- sweep(cc, 2, colSums(cc), "/") 
  for(j in 1:J){
    ss_0 <- ss_1
    A_grad <- sum(sweep(THETA_tuta[,1]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,1]*theta_post* c(1/(1+exp(-THETA_tuta%*%A_0[,j]))))
    A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%A_0[,j])/(1+exp(THETA_tuta%*%A_0[,j]))/(1+exp(THETA_tuta%*%A_0[,j]))))*(THETA_tuta[,1]*THETA_tuta[,1]))
    d_tuta <- A_0[1,j]-A_grad/A_grad_2
    aa <- A_0;aa[1,j] <- d_tuta;
    
    ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
    ))
    eps <- sum(log(ss_1/ss_0))
    while(eps>5){
      ss_0 <- ss_1
      A_grad <- sum(sweep(THETA_tuta[,1]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,1]*theta_post* c(1/(1+exp(-THETA_tuta%*%aa[,j]))))
      A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%aa[,j])/(1+exp(THETA_tuta%*%aa[,j]))/(1+exp(THETA_tuta%*%aa[,j]))))*(THETA_tuta[,1]*THETA_tuta[,1]))
      d_tuta <- d_tuta-A_grad/A_grad_2
      aa <- A_0;aa[1,j] <- d_tuta;
      ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
      ))
      eps <- sum(log(ss_1/ss_0))
    }
    
    A_0[1,j] <- d_tuta;
    A_grad <- sum(sweep(THETA_tuta[,2]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,2]*theta_post* c(1/(1+exp(-THETA_tuta%*%A_0[,j]))))  
    A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%A_0[,j])/(1+exp(THETA_tuta%*%A_0[,j]))/(1+exp(THETA_tuta%*%A_0[,j]))))*(THETA_tuta[,2]*THETA_tuta[,2]))
    A1_tuta <-A_0[2,j]-A_grad/A_grad_2
    aa <- A_0;aa[2,j] <- A1_tuta;
    ss_0 <- ss_1 
    ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
    ))
    eps <- sum(log(ss_1/ss_0))
    
    while(eps>5){
      ss_0 <- ss_1
      A_grad <- sum(sweep(THETA_tuta[,2]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,2]*theta_post* c(1/(1+exp(-THETA_tuta%*%aa[,j]))))
      A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%aa[,j])/(1+exp(THETA_tuta%*%aa[,j]))/(1+exp(THETA_tuta%*%aa[,j]))))*(THETA_tuta[,2]*THETA_tuta[,2]))
      A1_tuta <-A1_tuta-A_grad/A_grad_2
      aa <- A_0;aa[2,j] <- A1_tuta;
      ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
      ))
      eps <- sum(log(ss_1/ss_0))
    }
    
    A_0[2,j]<-soft(A1_tuta,-lammda[j]/A_grad_2)
    A_grad <- sum(sweep(THETA_tuta[,3]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,3]*theta_post* c(1/(1+exp(-THETA_tuta%*%A_0[,j]))))  
    A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%A_0[,j])/(1+exp(THETA_tuta%*%A_0[,j]))/(1+exp(THETA_tuta%*%A_0[,j]))))*(THETA_tuta[,3]*THETA_tuta[,3]))
    A2_tuta <-A_0[3,j]-A_grad/A_grad_2
    aa <- A_0;aa[3,j] <- A2_tuta;
    ss_0 <- ss_1
    ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
    ))
    eps <- sum(log(ss_1/ss_0))
    
    while(eps>5){
      ss_0 <- ss_1
      A_grad <- sum(sweep(THETA_tuta[,3]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,3]*theta_post* c(1/(1+exp(-THETA_tuta%*%aa[,j]))))
      A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%aa[,j])/(1+exp(THETA_tuta%*%aa[,j]))/(1+exp(THETA_tuta%*%aa[,j]))))*(THETA_tuta[,3]*THETA_tuta[,3]))
      A2_tuta <-A2_tuta-A_grad/A_grad_2
      aa <- A_0;aa[3,j] <- A2_tuta;
      ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
      ))
      eps <- sum(log(ss_1/ss_0))
    }
    
    A_0[3,j] <- soft(A2_tuta,-lammda[j]/A_grad_2)
    A_grad <- sum(sweep(THETA_tuta[,4]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,4]*theta_post* c(1/(1+exp(-THETA_tuta%*%A_0[,j]))))  
    A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%A_0[,j])/(1+exp(THETA_tuta%*%A_0[,j]))/(1+exp(THETA_tuta%*%A_0[,j]))))*(THETA_tuta[,4]*THETA_tuta[,4]))
    A3_tuta <-A_0[4,j]-A_grad/A_grad_2
    aa <- A_0;aa[4,j] <- A3_tuta;
    ss_0 <- ss_1 
    ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
    ))
    eps <- sum(log(ss_1/ss_0))
    
    while(eps>5){
      ss_0 <- ss_1
      A_grad <- sum(sweep(THETA_tuta[,4]*theta_post,2,response[j,],"*"))-sum(THETA_tuta[,4]*theta_post* c(1/(1+exp(-THETA_tuta%*%aa[,j]))))
      A_grad_2 <- -sum(rowSums(theta_post*c(exp(THETA_tuta%*%aa[,j])/(1+exp(THETA_tuta%*%aa[,j]))/(1+exp(THETA_tuta%*%aa[,j]))))*(THETA_tuta[,4]*THETA_tuta[,4]))
      A3_tuta <-A3_tuta-A_grad/A_grad_2
      aa <- A_0;aa[4,j] <- A3_tuta;
      ss_1 <- colSums(exp(log(exp(THETA_tuta%*%aa)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%aa)))%*%(1-response)
      ))
      eps <- sum(log(ss_1/ss_0))
    }
    
    A_0[4,j] <- soft(A3_tuta,-lammda[j]/A_grad_2)
    
  }
  
  
  A_0 <- A_0*t(cbind(rep(1,J),Q))
  ss_1 <- colSums(exp(log(exp(THETA_tuta%*%A_0)*dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%A_0)))%*%response+log(dmvnorm(THETA_tuta[,2:4],rep(0,K),diag(1,K))*ll*ll*ll/(1+exp(THETA_tuta%*%A_0)))%*%(1-response)
  )) 
  log_lik_1 <-sum(log(ss_1))
  likelihood_1 <- log_lik_1-sum(lammda*A_0[3:4,]);
  eps <- likelihood_1-likelihood_0
}

bic <- -2*log_lik_1+log(N)*(J*K)
RESULT <- rbind(c(bic,0,0,0),t(A_0))
write.csv(RESULT, file =paste0('dim_j',cond,'.csv'))
