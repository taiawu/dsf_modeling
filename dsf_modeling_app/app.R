library(SciViews) # for ln function used in LE models
library(signal) # for sgolayfilt
library(tidyverse)
library(shinyBS) # drop-down panels
library(shinyalert) # for pop-up messages
library(shinycssloaders) # spinning plot loading icon
library(shiny) # for shiny web-apps
source("scripts/data_modeling.R")


### time dependent probability function
calc_time_dep_prob <- function(kf = 10,
                               ku = 10,
                               kl1 = 100,
                               kl2 = 100,
                               kd1 = 100,
                               kd2 = 100,
                               kir = 0.1) {
    
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
    
    list(state_times = state_times,
         first_irrev = first_irrev,
         state = s)
}

#### tdp3 functions
calc_tdp3 <- function(kf, ku, kir) {
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
}



server <- function(session, input, output) { 

# #### time dependent model -----
    time_dep_list <- eventReactive( input$simulate_time_dep, {
                                            calc_time_dep_prob(kf = input$time_dep_kf,
                                                               ku = input$time_dep_ku,
                                                               kl1 = input$time_dep_kl1,
                                                               kl2 = input$time_dep_kl2,
                                                               kd1 = input$time_dep_kd1,
                                                               kd2 = input$time_dep_kd2,
                                                               kir = input$time_dep_kir)
                                        })
    
    

    output$state_time_plot <- renderPlot({
        tibble(y  = time_dep_list()$state_times, cycle = c(1:length(time_dep_list()$state_times))) %>%
            ggplot() + 
            geom_point(aes(cycle, y), size = 5) +
            labs(title = "Time Dependent Probability, State Times", y = "State Times", x = "Cycles")

                                        })

    output$first_irrev_plot <- renderPlot({
        
        tibble(y = time_dep_list()$first_irrev, cycle = c(1:length(time_dep_list()$first_irrev))) %>%
            ggplot() + 
            geom_point(aes(cycle, y ), size = 5) +
            labs(title = "Time Dependent Probability, First Irreversible", y = "Time to first irreversible", x = "Cycles")
                                            })
    
    
########### tdp_3 model  ---------- 
    tdp3_result <- eventReactive(input$simulate_tdp3, {
        calc_tdp3(input$tdp_kf, input$tdp_ku, input$tdp_kir)
    })
    
output$tdp3_plot <- renderPlot({ tdp3_result() })
    
########### LE models ----------
df_model <-  reactive({
    #make_model_df(start_T = (273+25), end_T = (273+95), dHu_ = 270, T_half_ = (55+273), dCp_ = 8, Ea_ = 100, T_star_ = 85+273, v_ = 1, nat_dye = 0, unf_dye = 1, fin_dye = 1, decay_rate = 0.8) -> df
    make_model_df(start_T = (273+25), end_T = (273+95), 
                  dHu_ = input$dHu_, 
                  T_half_ = input$T_half_ + 273, #55+273, 
                  dCp_ = input$dCp_, 
                  Ea_ = input$Ea_, 
                  T_star_ = input$T_star_ + 273, 
                  v_ = input$v_, 
                  nat_dye =input$nat_dye_, 
                  unf_dye = input$unf_dye_, 
                  fin_dye = input$fin_dye_, 
                  decay_rate = input$decay_rate_
    )
    ########### end LE models ----------
})

p_model <- reactive(make_model_plot(df_model()))
output$plot_model <- renderPlot(p_model())

}



ui <- navbarPage(useShinyalert(),

                 tabPanel("TDP3",
                          sidebarLayout(
                              sidebarPanel(
                                  #input$tdp_kf, input$tdp_ku, input$tdp_kir
                                sliderInput("tdp_kf", "K_f", min = 1, max = 1000, value = 100, step = 1),  
                                sliderInput("tdp_ku", "K_u", min = 1, max = 1000, value = 100, step = 1),
                                sliderInput("tdp_kir", "K_ir", min = 0.01, max = 100, value = 1, step = 0.01),
                                actionButton("simulate_tdp3", "Simulate!")
                                  
                              ), # end sidebarPanel
                              mainPanel(
                                  plotOutput("tdp3_plot") %>% withSpinner(color="#525252")
                              )
                          ) # end sidebar layout
                          ),
                 tabPanel("Time Dependent",
                          sidebarLayout(
                              sidebarPanel(

                                  sliderInput("time_dep_kf", "K_f", min = 1, max = 100, value = 10, step = 1),
                                  sliderInput("time_dep_ku", "K_u", min = 1, max = 100, value = 10, step = 1),
                                  sliderInput("time_dep_kl1", "K_l1", min = 0.01, max = 1000, value = 10, step = 1),
                                  sliderInput("time_dep_kl2", "Kl2", min = 0.01, max = 1000, value = 10, step = 1),
                                  sliderInput("time_dep_kd1", "K_d1", min = 0.01, max = 1000, value = 10, step = 1),
                                  sliderInput("time_dep_kd2", "K_d2", min = 0.01, max = 1000, value = 10, step = 1),
                                  sliderInput("time_dep_kir", "K_ir", min = 0.01, max = 100, value = 1, step = 0.01),
                 actionButton("simulate_time_dep", "Simulate!")
                
                              ), # end sidebarPanel
                              mainPanel(
                                  plotOutput("state_time_plot") %>% withSpinner(color="#525252"),
                                  plotOutput("first_irrev_plot") %>% withSpinner(color="#525252")
                              )
                          ) # end sidebar layout
                 ),
                 tabPanel("Lumry-Eyring",
                          tags$head( # set the slider aesthetic
                              tags$style(
                                  ".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: grey; border-color: transparent;}",
                                  ".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: grey; border-color: transparent;}",
                                  ".js-irs-2 .irs-single, .js-irs-2 .irs-bar-edge, .js-irs-2 .irs-bar {background: grey; border-color: transparent;}",
                                  
                                  ".js-irs-3 .irs-single, .js-irs-3 .irs-bar-edge, .js-irs-3 .irs-bar {background: grey; border-color: transparent;}",
                                  ".js-irs-4 .irs-single, .js-irs-4 .irs-bar-edge, .js-irs-4 .irs-bar {background: grey; border-color: transparent;}",
                                  ".js-irs-5 .irs-single, .js-irs-5 .irs-bar-edge, .js-irs-5 .irs-bar {background: grey; border-color: transparent;}",
                                  
                                  ".js-irs-6 .irs-single, .js-irs-6 .irs-bar-edge, .js-irs-6 .irs-bar {background: grey; border-color: transparent;}",
                                  ".js-irs-7 .irs-single, .js-irs-7 .irs-bar-edge, .js-irs-7 .irs-bar {background: grey; border-color: transparent;}",
                                  ".js-irs-8 .irs-single, .js-irs-8 .irs-bar-edge, .js-irs-8 .irs-bar {background: grey; border-color: transparent;}",
                                  ".js-irs-9 .irs-single, .js-irs-9 .irs-bar-edge, .js-irs-9 .irs-bar {background: grey; border-color: transparent;}",
                                  ".js-irs-10 .irs-single, .js-irs-9 .irs-bar-edge, .js-irs-10 .irs-bar {background: grey; border-color: transparent;}"
                                  
                              )
                          ),
                          sidebarLayout(
                              sidebarPanel(
                                  p("Tune model parameters below", style = "font-family: 'Avenir Next'; font-size: 20px; color: black",align = "center"),
                                  # # "Tune model parameters",
                                  bsCollapse(id = "thermo_pars", open = "Panel 1",
                                             bsCollapsePanel(p("Thermodynamic parameters", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                             sliderInput("T_half_", "Thermodynamic melting temperature (C)", min = 25, max = 95, value = 55, step = 1),
                                                             bsTooltip("T_half_", "At this temperature, the equilibrium ratio of folded to reversibly unfolded states is 1:1.",
                                                                       "right", options = list(container = "body")),
                                                             
                                                             sliderInput("dHu_", "Enthalpy of unfolding (kJ/mol)", min = 1, max = 600, value = 250, step = 10),
                                                             bsTooltip("dHu_", "Change in enthalpy between the folded and reversibly unfolded state. This model assumes the enthalpy of reversibly and irreversibly unfolded states are equal, a common simplification.",
                                                                       "right", options = list(container = "body", color = "white")),
                                                             
                                                             sliderInput("dCp_", "Change in heat capacity with unfolding (kJ/mol)", min = 0.1, max = 50, value = 8, step = 0.5),
                                                             bsTooltip("dCp_", "When a protein unfolds, the heat capacity of the solution changes. The magnitude of this change influences the Tm of the protein.",
                                                                       "right", options = list(container = "body", color = "white"))
                                                             , style = "default")
                                             
                                             
                                  ),
                                  
                                  bsCollapse(id = "kinetic_pars", open = "Panel 2",
                                             bsCollapsePanel(p("Kinetic parameters", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                             sliderInput("Ea_", "Activation energy of unfolding (kJ/mol)", min = 10, max = 1000, value = 200, step = 10),
                                                             bsTooltip("Ea_", "The energy barrier between the reversibly and irreversibly unfolded states. In this model, this value changes with temperature according to the Arrhenius Equation.",
                                                                       "right", options = list(container = "body")),
                                                             
                                                             sliderInput("T_star_", "Temperature at which irreversible unfolding becomes significant, T* (C)", min = 25, max = 95, value = 55, step = 1),
                                                             bsTooltip("T_star_", "Specifically, the temperature at which the rate constant of irreversible unfolding is 1/min.",
                                                                       "right", options = list(container = "body")),
                                                             
                                                             sliderInput("v_", "DSF experiment heating rate (C/min)", min = 0.1, max = 15, value = 1, step = 0.1),
                                                             bsTooltip("v_", "The thermocycling ramp rate used in the simulated DSF experiment.",
                                                                       "right", options = list(container = "body")),
                                                             style = "default")
                                  ),
                                  
                                  bsCollapse(id = "dye_pars", open = "Panel 3",
                                             bsCollapsePanel(p("Dye parameters", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"), 
                                                             # the heading above causes this warning: Warning in if (getAttribs(panels[[i]])$value %in% open) { : the condition has length > 1 and only the first element will be used
                                                             sliderInput("nat_dye_", "Folded state", min = 0, max = 1, value = 0, step = 0.1),
                                                             bsTooltip("nat_dye_", "The degree of dye binding and activation observed to the folded state of the protein. Detection of the folded state underlies many high background issues with DSF for hydrophobic proteins.",
                                                                       "right", options = list(container = "body")),
                                                             
                                                             sliderInput("unf_dye_", "Reversibly unfolded state", min = 0, max = 1, value = 1, step = 0.1),
                                                             bsTooltip("unf_dye_", "The degree of dye binding and activation observed to the reverisibly unfolded state of the protein. Low detection of this state may underlie the invisibility of some proteins in DSF.",
                                                                       "right", options = list(container = "body")),
                                                             
                                                             sliderInput("fin_dye_", "Irreversibly unfolded state", min = 0, max = 1, value = 1, step = 0.1),
                                                             bsTooltip("fin_dye_", "The degree of dye binding and activation observed to the irreversibly unfolded state of the protein. This parameter is non-zero for most proteins.",
                                                                       "right", options = list(container = "body")),
                                                             
                                                             sliderInput("decay_rate_", "Temperature sensitivity of dye activation", min = 0, max = 1, value = 0.85, step = 0.1),
                                                             bsTooltip("decay_rate_", "If all dye binding sites remained constant over the experiment, this is the rate at which fluorescence would decrease with temperature. Empirically this is close to 0.8, and is related to the general temperature-sensitivity of fluorophore quantum yields and the strength of hydrophobic effect.",
                                                                       "right", options = list(container = "body"))
                                                             
                                                             ,style = "default")
                                  )
                                  
                              ),
                              mainPanel(
                                  plotOutput("plot_model", width = "100%", height = "400px")
                              ))
                          
                 )

                 )


shinyApp(ui, server)