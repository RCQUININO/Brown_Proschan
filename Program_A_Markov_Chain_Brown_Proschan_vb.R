library(pracma)
clear()
tic()

#-----------------------------
# Parameters
#-----------------------------
q <- 0.9
b <- 1.5
p <- 0.8 #MR
trunc <- 30

log_q <- log(q)

#-----------------------------
# State (i,a)
#-----------------------------
O <- matrix(0, nrow = trunc+1, ncol = trunc+1)
O[1,1] <- 1

#-----------------------------
# Temporal evolution
#-----------------------------
for (t in 1:trunc){
  
  O_new <- matrix(0, nrow = trunc+1, ncol = trunc+1)
  
  for (i in 0:t){
    for (a in 0:t){
      
      prob <- O[i+1, a+1]
      
      
      if (!is.finite(prob) || prob == 0) next
      
      sobr <- a
      temp <- a + 1
      
      #-----------------------------
      # Steady-state probability (LOG)
      #-----------------------------
      a1 <- (temp-1)^b
      c1 <- temp^b
      d1 <- sobr^b
      
      log_num1 <- a1 * log_q
      log_num2 <- c1 * log_q
      log_den  <- d1 * log_q
      
      prob_fail <- exp(log_num1 - log_den) * (-expm1(log_num2 - log_num1))
      
      if (!is.finite(prob_fail)) prob_fail <- 0
      prob_fail <- max(min(prob_fail,1),0)
      
      #-------------------
      # 1. Failure and MR
      #-------------------
      if (i+1 <= trunc && a+1 <= t){
        O_new[i+2, a+2] <- O_new[i+2, a+2] + prob * prob_fail * p
      }
      
      #-------------------
      # 2. Failure and PR
      #-------------------
      if (i+1 <= trunc){
        O_new[i+2, 1] <- O_new[i+2, 1] + prob * prob_fail * (1 - p)
      }
      
      #-------------------
      # 3. No failure
      #-------------------
      if (a+1 <= t){
        O_new[i+1, a+2] <- O_new[i+1, a+2] + prob * (1 - prob_fail)
      }
      
    }
  }
  
  O <- O_new
}

#-----------------------------
# Result
#-----------------------------
media <- sum((0:trunc) %*% O)

cat("Mean failures =", media, "\n")
toc()