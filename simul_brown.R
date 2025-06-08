library(pracma)
clear()
tic()
corridas <- 200000
p <- 0.8#percentual de MR

q <- 0.9
b <- 1.5
trunc <- 30 #Truncamento
n <- 40


falhas <- matrix(0,corridas,1)
cont_mr <- matrix(0,corridas,1)
saida <- c()

for (j in 1:corridas){
  tempo_entre <- c()
  tempo_global <- c()
  trunc_falha <- 0
  saidaA <- c()
  renov <- c()
  for (i in 1:n){
    if (i == 1){
      
      r <- runif(1, trunc_falha, 1)
      y <- (log(1-r)/log(q))^(1/b)
      y <- ceiling(y)
      trunc_falha <- 1 - q^(y^(b))
      tempo_entre[i] <- y
      tempo_global[i] <- y
      
    }else{
      z <- rbinom(1, 1, p)
      renov[i-1] <- z
      if (z == 1){
        r <- runif(1, trunc_falha, 1)
        y <- (log(1-r)/log(q))^(1/b)
        y <- ceiling(y)
        if (length(which(renov == 0)) == 0){
          ult_renov <- 0
        }else{
          ult_renov <- tempo_global[max(which(renov == 0))]
        }
        tempo_global[i] <- y + ult_renov
        trunc_falha <- 1 - q^(y^(b))
        tempo_entre[i] <- tempo_global[i]+tempo_global[i-1]
        
      }else{
        r1 <- runif(1, 0, 1)
        y <- (log(1-r1)/log(q))^(1/b)
        y <- ceiling(y)
        tempo_entre[i] <- y
        trunc_falha <- 1 - q^(y^(b))
        tempo_global[i] <- tempo_global[i-1]+tempo_entre[i]
      }
    }
  }
  falhas[j,1] <- sum(tempo_global <= trunc)
  if (tempo_global[1] > trunc){
    cont_mr[j,1] <- 0
  }else{
    cont_mr[j,1] <- sum(renov[1:sum(tempo_global <= trunc)])
    
  }
  
}

#tempo_global
#falhas
#cont_mr

#media_naofalhas <- length(which(falhas[,1]==0))/corridas


print(mean(falhas[,1]))     
print(mean(cont_mr[,1]))
mean(cont_mr[,1])/(mean(falhas[,1]))#-media_naofalhas)
  
toc()