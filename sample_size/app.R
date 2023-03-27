#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(shinyhelper)
library(dplyr)
source("../helper.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
    withMathJax(),

    # Application title
    titlePanel("Version 1: Power calculations to detect treatment 
               effect heterogeneity by a single binary effect modifier in a cluster randomized trial with binary outcomes"),
    h5("This application faciliates the design of a CRT randomized to two arms with a constant \\(n\\) individuals per cluster. Let \\(Y_{ki}\\) be a binary outcome for individual i in cluster k,
       \\(X_{ki}\\) be a univariate binary effect modifier for individual i in cluster k with prevalence \\(\\pi_1\\) and covariate ICC \\(\\rho_{x}\\), 
       \\(W_k\\) be a treatment indicator for cluster k. Assuming the data-generating process 
       \\(logit (P(Y_{ki}=1 | W_k, X_{ki}, \\alpha_k))=\\beta_1 + \\beta_2 W_k + \\beta_3 X_{ki} + \\beta_4 X_{ki} W_k + \\alpha_k\\), where \\(\\alpha_k\\) is Normally distributed
       with mean 0 and variance \\(\\sigma_\\alpha^2\\), this application faciliates the design of a trial testing
       the null hypothesis \\(\\beta_4 = 0\\) using a Wald test."),
    sidebarLayout(
        sidebarPanel(
          h2("Trial parameters"),
          sliderInput("type_1_error",
                      "Type I Error (%)",
                      min = 0.01,
                      max = 20,
                      value = 5),
                numericInput("inter_eff",
                "Interaction odds ratio (OR)",
                1.179, min = 0, max = 10) %>% 
            helper(type = "inline",
                   title = "Interaction odds ratio (OR)",
                   content = withMathJax(helpText("Defined as \\(exp(\\beta_4)\\) above."))),
            numericInput("n",
                       "Cluster size",
                       1584, min = 2, max = 5000) %>% 
            helper(type = "inline",
                   title = "Cluster size",
                   content = withMathJax(helpText("Defined as \\(n\\) above."))),
            sliderInput("bar_W",
                      "Treatment allocation",
                      min = 0.01,
                      max = 0.99,
                      value = 0.5)  %>% 
            helper(type = "inline",
                   title = "Treatment allocation",
                   content = withMathJax(helpText("Proportion of clusters allocated to the treatment condition."))),
            sliderInput("p_x",
                      "Overall prevalence of binary effect modifier",
                      min = 0.01,
                      max = 0.99,
                      value = 0.4418) %>% 
            helper(type = "inline",
                   title = "Overall prevalence of binary effect modifier",
                   content = withMathJax(helpText("Defined as \\(\\pi_1\\) above."))),
            sliderInput("icc_x",
                      "Covariate Intraclass Correlation Coefficient (ICC)",
                      min = 0,
                      max = 1,
                      value = 0.007) %>% 
            helper(type = "inline",
                   title = "Covariate ICC",
                   content = withMathJax(helpText("ICC of binary effect modifier, or \\(\\rho_{x}\\) above."))),
            sliderInput("icc_yx",
                      "Outcome ICC",
                      min = 0,
                      max = 0.2,
                      value = 0.09)%>% 
            helper(type = "inline",
                   title = "Outcome ICC",
                   content = withMathJax(helpText("Defined as \\(\\rho_{y|x}\\), assuming the approximation \\(\\rho_{y|x}=\\frac{\\sigma_\\alpha^2}{\\sigma_\\alpha^2 + \\pi^3/2}\\)"))),
            sliderInput("baseline_prev",
                        "Baseline prevalence",
                        min = 0,
                        max = 1,
                        value = 0.139)  %>% 
            helper(type = "inline",
                   title = "Baseline prevalence",
                   content = withMathJax(helpText("Defined as the prevalence in the control group 
                                                  when the effect modifier is equal to 0, or \\(expit(\\beta_1)\\)."))),
            numericInput("trt_eff",
                        "Treatment OR",
                        1.39, min = 0, max = 10) %>%
            helper(type = "inline",
                   title = "Treatment OR",
                   content = withMathJax(helpText("Defined as \\(exp(\\beta_2)\\) above"))),
            numericInput("main_cov_eff",
                        "Main binary effect modifier OR",
                       1, min = 0, max = 10)  %>%
            helper(type = "inline",
                   title = "Main binary effect modifier OR",
                   content = withMathJax(helpText("Defined as \\(exp(\\beta_3)\\) above"))),
            h2("Plot attributes"),
            selectInput("to_vary", label="To vary in plot",
                      choices=c("Interaction OR",
                                "Cluster size",
                                "Covariate ICC",
                                "Outcome ICC"
                                ))%>% 
            helper(type = "inline",
                   title = "To vary in plot",
                   content = helpText("Please select the attribute you would like on the x-axis in the 
                                      \"Plot of parameter vs. number of clusters\" on the right hand side")),
          numericInput("minimum",
                       "Minimum value for",
                       1.1, min = 0, max = 4)%>% 
            helper(type = "inline",
                   title = "To vary in plot",
                   content = helpText("Please select the minimum value you would like on the x-axis in the 
                                      \"Plot of parameter vs. number of clusters\" on the right hand side")),
           numericInput("maximum",
                       "Maximum value for",
                       1.5, min = 0, max = 4) %>%
            helper(type = "inline",
                   title = "To vary in plot",
                   content = helpText("Please select the maximum value you would like on the x-axis in the 
                                      \"Plot of parameter vs. number of clusters\" on the right hand side")),
          width = 4
          ),
        
       
        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
              tabPanel("Number of clusters", fluid = TRUE,
                  sidebarLayout(
                      sliderInput("power",
                                  "Desired power (%)",
                                  min = 1,
                                  max = 99,
                                  value = 80),
                      mainPanel(
                          h3("Results"),
                          textOutput("text"),
                          h3("Plot of parameter vs. number of clusters"),
                          plotOutput("plotCluster", dblclick = "plot_reset"),
                          downloadButton('downloadDataCluster', 'Download Data'),
                          downloadButton('downloadPlotCluster', 'Download Plot')
                      )
                  )
              ),
              tabPanel("Power", fluid = TRUE,
                  sidebarLayout(
                      numericInput("K",
                                    "Number of clusters",
                                    value = 22, min = 0, max = 100),
                      mainPanel(h3("Results"),
                                textOutput("text_power"),
                                h3("Plot of parameter vs. power"),
                                plotOutput("plotPower", dblclick = "plot_reset"),
                                downloadButton('downloadDataPower', 'Download Data'),
                                downloadButton('downloadPlotPower', 'Download Plot')
                                
                    )
                  ),
              )
           )
        )
    )
)

logit <- function(x) {
  return(log(x/(1-x)))
}

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

logit_var_icc_1 <- function(beta_1,
                            beta_2,
                            beta_3,
                            beta_4,
                            n,
                            p_x,
                            bar_W,
                            sigma_gamma_sq) {
  sigma <- function(W, X) {
    mu <- expit(beta_1+beta_2*W+beta_3*X+beta_4*X*W)
    return(1/(mu*(1-mu)))
  }
  
  kappa <- function(W,X,n) {
    sigma.inv = 1/sigma(W,X)
    d = -1*sigma_gamma_sq/(1+n*sigma_gamma_sq*sigma.inv)
    return(sigma.inv+sigma.inv^2*n*d)
  }
  num1 = (1-p_x)*kappa(0,0,n) + (p_x)*kappa(0,1,n)
  denom1 = (1-bar_W)*p_x*(1-p_x)*kappa(0,1,n)*kappa(0,0,n)
  num2 = (1-p_x)*kappa(1,0,n) + (p_x)*kappa(1,1,n)
  denom2 = (bar_W)*p_x*(1-p_x)*kappa(1,1,n)*kappa(1,0,n)
  
  returned = (1/n * (num1/denom1 + num2/denom2))
  return(returned)
}

logit_var <- function(beta_1,
                      beta_2,
                      beta_3,
                      beta_4,
                      n,
                      p_x,
                      icc_x,
                      bar_W,
                      sigma_gamma_sq,
                      use_exact = TRUE) {
  
  sigma <- function(W, X) {
    mu <- expit(beta_1+beta_2*W+beta_3*X+beta_4*X*W)
    return(1/(mu*(1-mu)))
  }
  
  if (icc_x == 1 & use_exact) {
    return(logit_var_icc_1(beta_1,
                           beta_2,
                           beta_3,
                           beta_4,
                           n,
                           p_x,
                           bar_W,
                           sigma_gamma_sq))
  }
  eta_2 = (1/n)*((1+(n-1)*icc_x)*p_x + (n-1)*((1-icc_x)*p_x^2))
  P_1 <- function(s) {
    num1 = sigma_gamma_sq*((1-p_x)/sigma(s,0 )+p_x/sigma(s,1))^2
    denom1 <- 1+sigma_gamma_sq*n*((1-p_x)/sigma(s,0)+p_x/sigma(s,1))
    num2 = sigma_gamma_sq*(1/sigma(s,0)-1/sigma(s,1))^2*(eta_2-p_x^2)
    denom2 = (sigma_gamma_sq*(p_x-1)*n/sigma(s,0)-
                sigma_gamma_sq*p_x*n/sigma(s,1) - 1)^3
    return((-1*num1/denom1+num2/denom2))
  }
  
  P_2 <- function(s) {
    num1 = sigma_gamma_sq*p_x*((1-p_x)/sigma(s,0)+p_x/sigma(s,1))/sigma(s,1)
    denom1 <- 1+sigma_gamma_sq*n*((1-p_x)/sigma(s,0)+p_x/sigma(s,1))
    num2 = sigma_gamma_sq*(1/sigma(s,0)-1/sigma(s,1))*
      (1/sigma(s,1))*((sigma_gamma_sq * n / sigma(s,0)) + 1)*(eta_2-p_x^2)
    denom2 = (sigma_gamma_sq*(p_x-1)*n/sigma(s,0)-
                sigma_gamma_sq*p_x*n/sigma(s,1) - 1)^3
    return((-1*num1/denom1+num2/denom2))
  }
  
  P_3 <- function(s) {
    num1 = sigma_gamma_sq*eta_2*(1/sigma(s,1))^2
    denom1 <- 1+sigma_gamma_sq*n*((1-p_x)/sigma(s,0)+p_x/sigma(s,1))
    
    return((-1*num1/denom1))
  }
  
  a_1 = bar_W*p_x/sigma(1,1) + bar_W*(1-p_x)/sigma(1,0)  + 
    (1-bar_W)*p_x/sigma(0,1)  + 
    (1-bar_W)*(1-p_x)/sigma(0,0)  + (n*(bar_W*P_1(1) + (1-bar_W)*P_1(0)))
  
  a_2 = bar_W*(p_x/sigma(1,1) + (1-p_x)/sigma(1,0)) +
    n*bar_W*P_1(1) 
  b_1 = p_x*((1-bar_W)/sigma(0,1)  + bar_W/sigma(1,1)) +
    (n*(bar_W*P_2(1) + (1-bar_W)*P_2(0)))
  b_2 = bar_W*p_x/sigma(1,1) + n*bar_W*P_2(1)
  
  d_1 = p_x*((1-bar_W)/sigma(0,1)  + bar_W/sigma(1,1)) + 
    (n*(bar_W*P_3(1) + (1-bar_W)*P_3(0)))
  d_2 = bar_W*p_x/sigma(1,1) + n*bar_W*P_3(1)
  
  U = n*matrix(c(a_1,a_2,b_1,b_2,
                 a_2,a_2,b_2,b_2,
                 b_1,b_2,d_1,d_2,
                 b_2,b_2,d_2,d_2),byrow=TRUE,nrow=4)
  
  num1 = a_1-a_2
  denom1 = (a_1-a_2)*(d_1-d_2)-(b_1-b_2)^2
  
  num2 = a_2
  denom2 = a_2*d_2-b_2^2
  returned = (num1/denom1 + num2/denom2) * (1/(n))
  return(returned)
}

get_K_logit = function(type_1_error,
                       power, 
                       K,
                       beta_1,
                       beta_2,
                       beta_3,
                       beta_4,
                       n,
                       p_x,
                       icc_x,
                       bar_W,
                       sigma_gamma_sq) {
  if (is.na(K)) {
    sigma_4_sq = logit_var(beta_1,
                           beta_2,
                           beta_3,
                           beta_4,
                           n,
                           p_x,
                           icc_x,
                           bar_W,
                           sigma_gamma_sq)
    K_raw = sigma_4_sq*(qnorm(1-type_1_error/200)+qnorm(power/100))^2/beta_4^2
    return(ceiling(K_raw))
  }
  
  if (is.na(power)) {
    sigma_4_sq = logit_var(beta_1,
                           beta_2,
                           beta_3,
                           beta_4,
                           n,
                           p_x,
                           icc_x,
                           bar_W,
                           sigma_gamma_sq)
    return((pnorm(sqrt(K/sigma_4_sq)*abs(beta_4) - qnorm(1-type_1_error/200)) + 
              pnorm(-1*sqrt(K/sigma_4_sq)*abs(beta_4) - qnorm(1-type_1_error/200)))*100 )
  }
 
}

server <- function(input, output,session) {
    observe_helpers()
    observe({
        input$to_vary
        if (input$to_vary == "Interaction OR") {
            updateNumericInput(session, 
                               "minimum", 
                               label = paste("Minimum plot value for", input$to_vary),
                              1.1, min = 0, max = 4)
          
          updateNumericInput(session, 
                             "maximum", 
                             label = paste("Maximum plot value for", input$to_vary),
                             1.5, min = 0, max = 4)
        }
      
      if (input$to_vary == "Cluster size") {

        updateNumericInput(session, 
                           "minimum", 
                           label = paste("Minimum plot value for", input$to_vary),
                           max(n()-20,2), min = 0, max = 5000)
        updateNumericInput(session, 
                          "maximum", 
                          label = paste("Maximum plot value for", input$to_vary),
                          n()+20, min = 0, max = 5000)
      }
      
        if (input$to_vary == "Outcome ICC") {
          updateNumericInput(session, 
                           "minimum", 
                           label = paste("Minimum plot value for", input$to_vary),
                           0, min = 0, max = 1)
          updateNumericInput(session, 
                            "maximum", 
                            label = paste("Maximum plot value for", input$to_vary),
                            0.2, min = 0, max = 1)
        }
      
      if (input$to_vary == "Covariate ICC") {
        updateNumericInput(session, 
                           "minimum", 
                           label = paste("Minimum plot value for", input$to_vary),
                           0, min = 0, max = 1)
        updateNumericInput(session, 
                          "maximum", 
                          label = paste("Maximum plot value for", input$to_vary),
                          0.2, min = 0, max = 1)
      }
    })
  
    inter_eff <- reactive({ 
      validate(
        need(input$inter_eff > 0 & input$inter_eff < 10, 
             "Error: Interaction OR must be between 0 and 10"),
      )
      input$inter_eff
    })
    
    trt_eff <- reactive({ 
      validate(
        need(input$trt_eff > 0 & input$trt_eff < 10, 
             "Error: Treatment OR must be between 0 and 10"),
      )
      input$trt_eff
    })
    
    main_cov_eff <- reactive({ 
      validate(
        need(input$main_cov_eff > 0 & input$main_cov_eff < 10, 
             "Error: Covariate OR must be between 0 and 10"),
      )
      input$main_cov_eff
    })
    
    n <- reactive({ 
      validate(
        need(input$n > 2 & input$main_cov_eff < 5000, 
             "Error: Cluster size must be between 2 and 5,000"),
      )
      input$n
    })
    
    minimum <- reactive({ 
      validate(
        need(input$minimum < input$maximum, "Error: Plot minimum must be less than plot maximum"),
        need(if(input$to_vary == "Interaction OR") {
          input$minimum > 0} else {TRUE}, "Error: Minimum OR must be positive"),
        need(if(input$to_vary == "Cluster size") {
          input$minimum > 2} else {TRUE}, "Error: Minimum cluster size must be greater than 2"),
        need(if(input$to_vary == "Covariate ICC") {
          input$minimum >= 0} else {TRUE}, "Error: Minimum covariate ICC must be positive"),
        need(if(input$to_vary == "Outcome ICC") {
          input$minimum >= 0} else {TRUE}, "Error: Minimum outcome ICC must be positive")
      )
      input$minimum
    })
    
    maximum <- reactive({ 
      validate(
        need(input$minimum < input$maximum, "Minimum must be less than maximum"),
        need(if(input$to_vary == "Interaction OR") {
          input$maximum < 10} else {TRUE}, "Error: Maximum OR must be less than 10"),
        need(if(input$to_vary == "Cluster size") {
          input$minimum < 5000} else {TRUE}, "Error: Maximum cluster size must be less than 5000"),
        need(if(input$to_vary == "Covariate ICC") {
          input$minimum <= 1} else {TRUE}, "Error: Maximum covariate ICC cannot be greater than 1"),
        need(if(input$to_vary == "Outcome ICC") {
          input$minimum < 0.2} else {TRUE}, "Error: Maximum outcome ICC must be less than 0.2")
        
      )
      input$maximum
    })
    
    output$text <- renderText({
      beta_1 <- logit(input$baseline_prev)
      beta_2 <- log(trt_eff())
      beta_3 <- log(main_cov_eff())
      beta_4 <- log(inter_eff())
      sigma_gamma_sq = (pi^2/3*input$icc_yx)/(1-input$icc_yx)
      paste("This CRT requires",
                  get_K(beta_1,
                        beta_2,
                        beta_3,
                        beta_4,
                        n(),
                        input$p_x,
                        input$icc_x,
                        input$bar_W,
                        sigma_gamma_sq,
                        sigma_sq_logit,
                        input$power/100,
                        input$type_1_error/100), "clusters to detect an odds ratio of",
            inter_eff(), 
            "with", paste(input$power, "%",sep=""), "power.")
      
     
    })
    
    output$text_power <- renderText({
      beta_1 <- logit(input$baseline_prev)
      beta_2 <- log(trt_eff())
      beta_3 <- log(main_cov_eff())
      beta_4 <- log(inter_eff())
      sigma_gamma_sq = (pi^2/3*input$icc_yx)/(1-input$icc_yx)
      paste("This CRT with ", input$K, " clusters has ",
            round(get_power(
                        input$K,
                        beta_1,
                        beta_2,
                        beta_3,
                        beta_4,
                        n(),
                        input$p_x,
                        input$icc_x,
                        input$bar_W,
                        sigma_gamma_sq,
                        sigma_sq_logit,
                        input$type_1_error/100),2), "% power detect an odds ratio of ",
            inter_eff(), ".",sep="")
      
      
    })
    
    dataInputCluster <- reactive({
      min = minimum()
      max = maximum()
      beta_1 <- logit(input$baseline_prev)
      beta_2 <- log(trt_eff())
      beta_3 <- log(main_cov_eff())
      beta_4 <- log(inter_eff())
      sigma_gamma_sq = (pi^2/3*input$icc_yx)/(1-input$icc_yx)
      if (input$to_vary == "Interaction OR") {
        if (min < 0.01) {
          min = 0.01
        } 
        if (max > 10) {
          max = 10
        }
        data = data.frame(x=seq(min, max,by=0.01))
        data$y = sapply(data$x, function(x) {
          return(get_K(
            beta_1,
            beta_2,
            beta_3,
            log(x),
            n(),
            input$p_x,
            input$icc_x,
            input$bar_W,
            sigma_gamma_sq,
            sigma_sq_logit,
            input$power/100,
            input$type_1_error/100))})
      }
      if (input$to_vary == "Cluster size") {
        if (min < 2) {
          min = 2
        } 
        data = data.frame(x=seq(min, max,by=1))
        data$y = sapply(data$x, function(x) {
          return(get_K_logit(
            input$type_1_error,
            input$power,
            NA,
            beta_1,
            beta_2,
            beta_3,
            beta_4,
            x,
            input$p_x,
            input$icc_x,
            input$bar_W,
            sigma_gamma_sq))})
      }
      
      if (input$to_vary == "Covariate ICC") {
        if (min < 0) {
          min = 0
        } 
        if (max > 1) {
          max = 1
        }
        data = data.frame(x=seq(min, max,by=0.01))
        data$y = sapply(data$x, function(x) {
          return(get_K_logit(input$type_1_error,
                             input$power,
                             NA,
                             beta_1,
                             beta_2,
                             beta_3,
                             beta_4,
                             n(),
                             input$p_x,
                             x,
                             input$bar_W,
                             sigma_gamma_sq))})
      }
      
      if (input$to_vary == "Outcome ICC") {
        if (min < 0) {
          min = 0
        } 
        if (max > 0.2) {
          max = 0.2
        }
        data = data.frame(x=seq(min, max,by=0.01))
        data$y = sapply(data$x, function(x) {
          return(get_K_logit(input$type_1_error,
                             input$power,
                             NA,
                             beta_1,
                             beta_2,
                             beta_3,
                             beta_4,
                             n(),
                             input$p_x,
                             input$icc_x,
                             input$bar_W,
                             (pi^2/3*x)/(1-x)))})
      }
      names(data)[2] <- "cluster"
      
      return(data)
    })
    
    dataInputPower <- reactive({
      min = minimum()
      max = maximum()
      beta_1 <- logit(input$baseline_prev)
      beta_2 <- log(trt_eff())
      beta_3 <- log(main_cov_eff())
      beta_4 <- log(inter_eff())
      sigma_gamma_sq = (pi^2/3*input$icc_yx)/(1-input$icc_yx)
      
      if (input$to_vary == "Interaction OR") {
        if (min < 0.01) {
          min = 0.01
        } 
        if (max > 10) {
          max = 10
        }
        data = data.frame(x=seq(min, max,by=0.01))
        data$y = sapply(data$x, function(x) {
          return(get_power(
            K,
            beta_1,
            beta_2,
            beta_3,
            log(x),
            n(),
            input$p_x,
            input$icc_x,
            input$bar_W,
            sigma_gamma_sq,
            sigma_sq_logit,
            input$type_1_error/100)) })
      }
      if (input$to_vary == "Cluster size") {
        if (min < 2) {
          min = 2
        } 
        data = data.frame(x=seq(min, max,by=1))
        data$y = sapply(data$x, function(x) {
          return(get_K_logit(
            input$type_1_error,
            NA,
            input$K,
            beta_1,
            beta_2,
            beta_3,
            beta_4,
            x,
            input$p_x,
            input$icc_x,
            input$bar_W,
            sigma_gamma_sq))})
      }
      
      if (input$to_vary == "Covariate ICC") {
        if (min < 0) {
          min = 0
        } 
        if (max > 1) {
          max = 1
        }
        data = data.frame(x=seq(min, max,by=0.01))
        data$y = sapply(data$x, function(x) {
          return(get_K_logit(
                            input$type_1_error,
                             NA,
                             input$K,
                             beta_1,
                             beta_2,
                             beta_3,
                             beta_4,
                             n(),
                             input$p_x,
                             x,
                             input$bar_W,
                             sigma_gamma_sq))})
      }
      
      if (input$to_vary == "Outcome ICC") {
        if (min < 0) {
          min = 0
        } 
        if (max > 0.2) {
          max = 0.2
        }
        data = data.frame(x=seq(min, max,by=0.01))
        data$y = sapply(data$x, function(x) {
          return(get_K_logit(input$type_1_error,
                             NA,
                             input$K,
                             beta_1,
                             beta_2,
                             beta_3,
                             beta_4,
                             n(),
                             input$p_x,
                             input$icc_x,
                             input$bar_W,
                             (pi^2/3*x)/(1-x)))})
      }
      names(data)[2] <- "power"
      return(data)
    })
    
    plotInputCluster <- reactive({
      df <- dataInputCluster()
      p <- ggplot(df,aes(x=x,y=cluster))+
        geom_point()+
        labs(x=input$to_vary, y="Number of clusters") 
    })
    
    output$plotCluster <- renderPlot({
      print(plotInputCluster())
    }, res = 96)
    
    output$downloadDataCluster <- downloadHandler(
      filename = function() { paste(input$to_vary,'cluster', '.csv', sep='') },
      content = function(file) {
        df <- dataInputCluster()
        names(df) <- c(input$to_vary,"cluster")
        write.csv(df, file)
      }
    )
    output$downloadPlotCluster <- downloadHandler(
      filename = function() { paste(input$to_vary,'cluster', '.png', sep='') },
      content = function(file) {
        ggsave(file,plotInputCluster())
      }
    )
    
    plotInputPower <- reactive({
      df <- dataInputPower()
      p <- ggplot(df,aes(x=x,y=power))+
        geom_point()+
        labs(x=input$to_vary, y="Power") 
    })
    
    output$plotPower <- renderPlot({
      print(plotInputPower())
    }, res = 96)
    
    output$downloadDataPower <- downloadHandler(
      filename = function() { paste(input$to_vary,'power', '.csv', sep='') },
      content = function(file) {
        df <- dataInputPower()
        names(df) <- c(input$to_vary,"power")
        write.csv(df, file)
      }
    )
    output$downloadPlotPower <- downloadHandler(
      filename = function() { paste(input$to_vary,'power', '.png', sep='') },
      content = function(file) {
        ggsave(file,plotInputPower())
      }
    )
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
