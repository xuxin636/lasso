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
theta_square <- THETA_tuta[,2:4]*THETA_tuta[,2:4]
theta_tmp <- rowSums(theta_square)/2
xx <- seq(0.001,0.03,0.001);xx1 <- matrix(0,nrow = length(xx)*length(xx),ncol=2);xx1[,2] <- rep(xx,length(xx));xx1[,1] <- c(rep(1,length(xx))%*%t(xx))
lammda <- c(rep(xx1[cond,1],J/2),rep(xx1[cond,2],J/2))*N;
soft <- function(a,b){
  if(a>0&a>b){return(a-b)}
  else{return(0)}
}
response <- t(response);
A_0 <- t(D_initial)
temp_0 <- like_temp_0 <- THETA_tuta%*%A_0
cc1 <- exp(like_temp_0%*%response-theta_tmp-rowSums(log(1+exp(like_temp_0))))
theta_post <- sweep(cc1, 2, colSums(cc1), "/") 
penalty<- -sum(sweep(A_0,2,lammda,"*"))
th0 <- rowSums(theta_post);
th1 <- THETA_tuta[,2]*th0;th2 <- THETA_tuta[,3]*th0;th3 <- THETA_tuta[,4]*th0;
uu <- theta_post%*%t(response);uu0 <-colSums(uu);uu1 <-colSums(THETA_tuta[,2]*uu);uu2 <-colSums(THETA_tuta[,3]*uu);uu3 <-colSums(THETA_tuta[,4]*uu);
frac_00 <- log(colSums(exp(temp_0%*%response-rowSums(log(1+exp(temp_0)))-theta_tmp)))
timstart <- Sys.time()
for(j in 1:J){
  xx <- c(1/(1+exp(-temp_0[,j])))
  A_grad <- uu0[j]-sum(th0* xx)
  A_grad_2 <- -sum(th0*xx*(1-xx))
  d_tuta <- A_0[1,j]-A_grad/A_grad_2
  temp_1 <- temp_0;temp_1[,j] <- THETA_tuta%*%c(d_tuta,A_0[-1,j])
  frac_0 <- log(colSums(exp(temp_0%*%response-rowSums(log(1+exp(temp_0)))-theta_tmp)));
  frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
  eps <- sum(frac_1-frac_0)
  while(eps>5){
    frac_0 <-frac_1
    temp_0 <- temp_1;
    A_0[1,j] <- d_tuta
    xx <- c(1/(1+exp(-temp_0[,j])))
    A_grad <- uu0[j]-sum(th0* xx)
    A_grad_2 <- -sum(th0*xx*(1-xx))
    d_tuta <- d_tuta-A_grad/A_grad_2
    temp_1[,j] <- THETA_tuta%*%c(d_tuta,A_0[-1,j])
    frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
    eps <-  sum(frac_1-frac_0)
  }
  
  frac_0 <-frac_1
  A_grad <- uu1[j]-sum(th1* xx)  
  A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,1])
  A1_tuta <-A_0[2,j]-A_grad/A_grad_2
  temp_1[,j] <- THETA_tuta%*%c(A_0[1,j],A1_tuta,A_0[3:4,j])
  frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
  eps <- sum(frac_1-frac_0) 
  while(eps>5){
    frac_0 <-frac_1
    temp_0 <- temp_1;
    xx <- c(1/(1+exp(-temp_0[,j])))
    A_0[2,j] <- A1_tuta
    A_grad <- uu1[j]-sum(th1* xx)  
    A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,1])
    A1_tuta <-A_0[2,j]-A_grad/A_grad_2
    temp_1[,j] <- THETA_tuta%*%c(A_0[1,j],A1_tuta,A_0[3:4,j])
    frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
    eps <- sum(frac_1-frac_0) 
  }
  A_0[2,j]<-soft(A_0[2,j],-lammda[j]/A_grad_2)
  
  
  frac_0 <-frac_1
  A_grad <- uu2[j]-sum(th2* xx)  
  A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,2])
  A2_tuta <-A_0[3,j]-A_grad/A_grad_2
  temp_1[,j] <- THETA_tuta%*%c(A_0[1:2,j],A2_tuta,A_0[4,j])
  frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
  eps <- sum(frac_1-frac_0) 
  while(eps>5){
    frac_0 <-frac_1
    temp_0 <- temp_1;
    xx <- c(1/(1+exp(-temp_0[,j])))
    A_0[3,j] <- A2_tuta
    A_grad <- uu2[j]-sum(th2* xx)  
    A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,2])
    A2_tuta <-A_0[3,j]-A_grad/A_grad_2
    temp_1[,j] <- THETA_tuta%*%c(A_0[1:2,j],A2_tuta,A_0[4,j])
    frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
    eps <- sum(frac_1-frac_0) 
  }
  A_0[3,j] <- soft(A_0[3,j],-lammda[j]/A_grad_2)
  
  frac_0 <-frac_1
  A_grad <- uu3[j]-sum(th3* xx)
  A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,3])
  A3_tuta <-A_0[4,j]-A_grad/A_grad_2
  temp_1[,j] <- THETA_tuta%*%c(A_0[1:3,j],A3_tuta)
  frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
  eps <- sum(frac_1-frac_0)  
  
  while(eps>5){
    frac_0 <-frac_1
    temp_0 <- temp_1;
    xx <- c(1/(1+exp(-temp_0[,j])))
    A_0[4,j] <- A3_tuta
    A_grad <- uu3[j]-sum(th3* xx)
    A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,3])
    A3_tuta <-A3_tuta-A_grad/A_grad_2
    temp_1[,j] <- THETA_tuta%*%c(A_0[1:3,j],A3_tuta)
    frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
    eps <- sum(frac_1-frac_0) 
  }
  A_0[4,j] <- soft(A_0[4,j],-lammda[j]/A_grad_2)
}

A_0 <- A_0*t(cbind(rep(1,J),Q))
eps <-sum(frac_0-frac_00) -sum(lammda*A_0[3:4,])-penalty
while(eps>5){
  cc1 <- exp(temp_0%*%response-theta_tmp-rowSums(log(1+exp(temp_0))))
  theta_post <- sweep(cc1, 2, colSums(cc1), "/") 
  penalty<- -sum(sweep(A_0,2,lammda,"*"))
  frac_00 <- frac_0
  for(j in 1:J){
    xx <- c(1/(1+exp(-temp_0[,j])))
    A_grad <- uu0[j]-sum(th0* xx)
    A_grad_2 <- -sum(th0*xx*(1-xx))
    d_tuta <- A_0[1,j]-A_grad/A_grad_2
    temp_1 <- temp_0;temp_1[,j] <- THETA_tuta%*%c(d_tuta,A_0[-1,j])
    frac_0 <- log(colSums(exp(temp_0%*%response-rowSums(log(1+exp(temp_0)))-theta_tmp)));
    frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
    eps <- sum(frac_1-frac_0)
    while(eps>5){
      frac_0 <-frac_1
      temp_0 <- temp_1;
      A_0[1,j] <- d_tuta
      xx <- c(1/(1+exp(-temp_0[,j])))
      A_grad <- uu0[j]-sum(th0* xx)
      A_grad_2 <- -sum(th0*xx*(1-xx))
      d_tuta <- d_tuta-A_grad/A_grad_2
      temp_1[,j] <- THETA_tuta%*%c(d_tuta,A_0[-1,j])
      frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
      eps <-  sum(frac_1-frac_0)
    }
    
    frac_0 <-frac_1
    A_grad <- uu1[j]-sum(th1* xx)  
    A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,1])
    A1_tuta <-A_0[2,j]-A_grad/A_grad_2
    temp_1[,j] <- THETA_tuta%*%c(A_0[1,j],A1_tuta,A_0[3:4,j])
    frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
    eps <- sum(frac_1-frac_0) 
    while(eps>5){
      frac_0 <-frac_1
      temp_0 <- temp_1;
      xx <- c(1/(1+exp(-temp_0[,j])))
      A_0[2,j] <- A1_tuta
      A_grad <- uu1[j]-sum(th1* xx)  
      A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,1])
      A1_tuta <-A_0[2,j]-A_grad/A_grad_2
      temp_1[,j] <- THETA_tuta%*%c(A_0[1,j],A1_tuta,A_0[3:4,j])
      frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
      eps <- sum(frac_1-frac_0) 
    }
    A_0[2,j]<-soft(A_0[2,j],-lammda[j]/A_grad_2)
    
    
    frac_0 <-frac_1
    A_grad <- uu2[j]-sum(th2* xx)  
    A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,2])
    A2_tuta <-A_0[3,j]-A_grad/A_grad_2
    temp_1[,j] <- THETA_tuta%*%c(A_0[1:2,j],A2_tuta,A_0[4,j])
    frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
    eps <- sum(frac_1-frac_0) 
    while(eps>5){
      frac_0 <-frac_1
      temp_0 <- temp_1;
      xx <- c(1/(1+exp(-temp_0[,j])))
      A_0[3,j] <- A2_tuta
      A_grad <- uu2[j]-sum(th2* xx)  
      A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,2])
      A2_tuta <-A_0[3,j]-A_grad/A_grad_2
      temp_1[,j] <- THETA_tuta%*%c(A_0[1:2,j],A2_tuta,A_0[4,j])
      frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
      eps <- sum(frac_1-frac_0) 
    }
    A_0[3,j] <- soft(A_0[3,j],-lammda[j]/A_grad_2)
    
    frac_0 <-frac_1
    A_grad <- uu3[j]-sum(th3* xx)
    A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,3])
    A3_tuta <-A_0[4,j]-A_grad/A_grad_2
    temp_1[,j] <- THETA_tuta%*%c(A_0[1:3,j],A3_tuta)
    frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
    eps <- sum(frac_1-frac_0)  
    
    while(eps>5){
      frac_0 <-frac_1
      temp_0 <- temp_1;
      xx <- c(1/(1+exp(-temp_0[,j])))
      A_0[4,j] <- A3_tuta
      A_grad <- uu3[j]-sum(th3* xx)
      A_grad_2 <- -sum(th0*xx*(1-xx)*theta_square[,3])
      A3_tuta <-A3_tuta-A_grad/A_grad_2
      temp_1[,j] <- THETA_tuta%*%c(A_0[1:3,j],A3_tuta)
      frac_1 <- log(colSums(exp(temp_1%*%response-rowSums(log(1+exp(temp_1)))-theta_tmp)))
      eps <- sum(frac_1-frac_0) 
    }
    A_0[4,j] <- soft(A_0[4,j],-lammda[j]/A_grad_2)
  }
  A_0 <- A_0*t(cbind(rep(1,J),Q))
  eps <-sum(frac_0-frac_00) -sum(sweep(A_0,2,lammda,"*"))-penalty
}

timend <- Sys.time()
tt <- timend-timstart
tt
bic <- -2*sum(log(colSums(exp(temp_0%*%response-rowSums(log(1+exp(temp_0)))-theta_tmp))))+log(N)*(J*K)
RESULT <- rbind(c(bic,0,0,0),t(A_0))
write.csv(RESULT, file =paste0('dim_j',cond,'.csv'))
