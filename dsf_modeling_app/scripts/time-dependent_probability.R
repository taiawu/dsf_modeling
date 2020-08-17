# kl for light, kd for dark
#
#               Unfolded:DyeBound  Irrev.:DyeBound
#                   ^ ^              ^  ^
#                 kl1 |             kl2 |
#                     | kd1             | kd2
#           < kf      v v               v  v
#    Folded <---> Unfolded ---> Irrev. unfolded
#             ku >        kir > 
#


#declaring rates (time unit^-1)
kf <- 100
ku <- 100
kl1 <- 1000
kl2 <- kl1
kd1 <- 1000
kd2 <- kd1
kir <- 0.1

#rate matrix: rows represent states, columns represent paths out of states
rates <- matrix(c(
  0,ku,0,0,0,
  kf,0,kl1,kir,0,
  0,kd1,0,0,0,
  0,0,0,0,kl2,
  0,0,0,kd2,0), 
  nrow = 5, ncol = 5, byrow = TRUE
  )

#transition boundary matrix
#boundary fraction likelihood of taking path out of state
# 0|xxxxxxxway1xxxxxxxx|yyyway2yyy|1
bounds <- matrix(c(
  0,1,0,0,0,
  kf/sum(kf,kl1,kir),0,sum(kf,kl1)/sum(kf,kl1,kir),1,0,
  0,1,0,0,0,
  0,0,0,0,1,
  0,0,0,1,0), 
  nrow = 5, ncol = 5, byrow = TRUE
)


#Summary of rates, for my own brain-functioning
# dfold/dt = kf(rev) - ku(fold)
# drev/dt = ku(fold) + kd1(revstar) - kl1(rev) - kf(rev) - kir(rev)
# drevstar/dt = kl1(rev) - kd1(revstar)
# dirrev/dt = kir(rev) + kd2(irrevstar) - kl2(irrev)
# dirrevstar/dt = kl2(irrev) - kd2(irrevstar)

#foldS = exp(-ku * t)
#revS = exp(-(kl1 + kf + kir) * t)
#revstarS = exp(-kd1 * t)
#irrevS = exp(-kl2 * t)
#irrevstarS = exp(-kd2 * t)

#set up algorithm
first_irrev <- c()
cycles <- 0

while(cycles < 11){
  explen <- 500
  
  t <- 0
  state <- matrix(c(
    1,0,0,0,0), 
    nrow = 1, ncol = 5)
  
  #declare "time in each state" record
  state_times <- c(0,0,0,0,0)
  x <- 0
  
  
  while(t < explen) {
    #random number generation
    u1 <- runif(1)
    u2 <- runif(1)
    
    #set rates/bounds for that state
    step_param <- state%*%rates
    step_bound <- state%*%bounds
    
    #set next state to "unknown"
    state <- state*0
    
    #step time
    tstep <- -1/sum(step_param) * log(u1)
    
    #determine which state to move to, as determined by bounds
    s <- 1
    for(i in step_bound){
      if(u2 < i){
        state[1,s] <- 1
        break
      }
      s <- s + 1
    }
    #add time spent in that state
    state_times[s] <- state_times[s] + tstep
    
    #at what t does first transition to irreversibly unfolded occur?
    if(s == 4 && x == 0){
      print(t)
      first_irrev <- c(first_irrev, t)
      x <- 1
    }
    t <- t + tstep
  }
  
  cycles <- cycles + 1
}



  







