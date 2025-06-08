library(pracma)
library(expm)
library(TruncatedDistributions)
library(DiscreteWeibull)
clear()

#Discrete Weibull Parameters
q <- 0.9
b <- 1.5
p <- 0.8 #Percentage of MR
trunc <- 30 #Maximum Observation Time – Discrete

#########################
##MODEL Brown-Proschan###
#########################

E1 <- c()

for (i in (0:trunc)){ #Falhas
  for (j in (0:trunc)){ #Tempo
    for(k in (0:trunc)){ #quando ?ltima falha
      
      
      if (i == 0 & j >= 0 & k == 0){
        E <- c(i,j,k)
        E1 <- rbind(E1,E)
      }
      
      if (i == j & i < trunc & j < trunc & i > 0 & k <= i){
        E <- c(i,j,k)
        E1 <- rbind(E1,E)
      }
      
      if (i < j & i < trunc & j < trunc & k < trunc & i > 0 & k <= j){
        E <- c(i,j,k)
        E1 <- rbind(E1,E)
      }
      
      if (i < j & j == trunc & i > 0 & k <= j){
        E <- c(i,j,k)
        E1 <- rbind(E1,E)
      }
      
      if (i == trunc & j == trunc){
        E <- c(i,j,k)
        E1 <- rbind(E1,E)
      }
      
    }
  }
}


E2 <- E1

num_row = nrow(E1)
mat_prob  <- zeros(num_row)


for (i in 1:num_row){     
  num_falha1 <- E1[i,1]
  tempo1  <- E1[i,2]
  ult_falha1 <- E1[i,3]
  
  for (j in 1:num_row){   
    num_falha2 <- E2[j,1]
    tempo2 <- E2[j,2]
    ult_falha2 <- E2[j,3]
    
    if ((num_falha2 == num_falha1 + 1) & (tempo2 == tempo1 + 1) & num_falha1 < trunc & tempo1 < trunc & 
        ult_falha2 == tempo2){
      sobr <- tempo1 - ult_falha1
      temp <- tempo2 - ult_falha1
      mat_prob[i,j] <- (q^((temp-1)^b)-(q^(temp^b)))/(q^(sobr^b))*(1-p)
    }
    
    if ((num_falha2 == num_falha1 + 1) & (tempo2 == tempo1 + 1) & num_falha1 < trunc & tempo1 < trunc & 
        ult_falha2 == ult_falha1){
      sobr <- tempo1 - ult_falha1
      temp <- tempo2 - ult_falha1
      mat_prob[i,j] <- (q^((temp-1)^b)-(q^(temp^b)))/(q^(sobr^b))*p
    }
    
    if ((num_falha2 == num_falha1) & (tempo2 == tempo1 + 1) & num_falha1 < trunc & tempo1 < trunc 
        & ult_falha2 == ult_falha1){
      sobr <- tempo1 - ult_falha1
      temp <- tempo2 - ult_falha1
      mat_prob[i,j] <- 1 - (q^((temp-1)^b)-(q^(temp^b)))/(q^(sobr^b))
    }
    
    if ((num_falha2 == num_falha1) & (tempo2 == tempo1) & tempo1 == trunc & 
        ult_falha1 == ult_falha2){
      
      mat_prob[i,j] <- 1 
    }
  }
}

P2 <- mat_prob%^%(trunc)
O <- rep(0,num_row)
O[1] <- 1
final <- (O%*%P2)
media <- sum(E1[,1]*final)

cat("Mean failures=",media,"\n")
