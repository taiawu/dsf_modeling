#           < kf
#    Folded <---> Unfolded ---> Irrev. unfolded
#             ku >        kir > 

#declaring rates (time unit^-1)

kf <- 100
ku <- 100
kir <- 0.01


#rate matrix: rows represent states, columns represent paths out of states
rates <- matrix(c(
  0,ku,0,
  kf,0,kir,
  0,0,0), 
  nrow = 3, ncol = 3, byrow = TRUE
)

#boundary fraction likelihood of taking path out of state
# 0|xxxxxxxway1xxxxxxxx|yyyway2yyy|1
bounds <- matrix(c(
  0,1,0,
  kf/sum(kf,kir),0,1,
  0,0,0), 
  nrow = 3, ncol = 3, byrow = TRUE
)

#set up algorithm
first_irrev <- c()
cycles <- 0
explen <- 500

#for (cycles < n) number of cycles, 
#set molecule in "folded" state
#run gillespie: 
#while total time elapsed is less than experiment length
#wait random timestep (based on u1 and rates out of state) in current state
#move state based on random u2 (once leaving state, likelihood of entering
#another if multiple options)


while(cycles < 1001){
  
  t <- 0
  #state as matrix allows using matrix multiplication to set rates/bounds
  state <- matrix(c(
    1,0,0), 
    nrow = 1, ncol = 3)
  
  #state_times <- c(0,0,0)
  
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
    #state_times[s] <- state_times[s] + tstep
    
    #condition: end experiment if "irreversible unfolded" state entered
    if(s == 3){
      break
    }
    #advance time
    t <- t + tstep
  }
  
  first_irrev <- c(first_irrev, t)
  cycles <- cycles + 1
}

#viz results
hist(first_irrev, breaks = 20)
