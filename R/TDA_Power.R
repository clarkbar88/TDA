#' Complete TDA power analyis
#'
#' @param value vector of values concentration, BIBI, etc.
#' @param year vector of years - if multiple values for same year, values are averaged arithmetically
#' @param site vector of sites
#' @param trend.magnitude.threshold absolute value of ecologically meaningful trend magnitude. if not provided, defaults to 10% of average value
#' @param power 1-B - defaults to 0.8
#' @param power.simulation how many permutations to run to calculate power - default to 1000
#' @return A tibble of trends and power.
#' @examples
#' TDA_Power(1, 1)

#TDA Power Function
#value - vector of values (concentration, BIBI, etc.)
#year - vector of years - if multiple values for same year, values are averaged (arithmetically)
#site - vector of sites
#trend.magnitude.threshold - absolute value of ecologically meaningful trend magnitude,
#if not provided, defaults to 10% of average value
#power - 1-B - defaults to 0.8
#power.simulations.n - how many permutations to run to calculate power - default to 1000
#requireNamespace(c('tidyverse','lubridate','EnvStats'))

TDA_Power<-function(value,year,site,trend.magnitude.threshold,power=0.8,power.simulations.n=1000){
  if(missing(trend.magnitude.threshold)) trend.magnitude.threshold=0.1*mean(value) #default to 10% of average
  trend.magnitude.threshold<-abs(trend.magnitude.threshold)
  set.seed(2)
  
  tibble(value=value,year=year,site=site) %>%
    group_by(site,year) %>%
    summarise(value=mean(value)) %>% ##average duplicates
    group_by(site) %>%
    nest() %>%
    mutate(n=map(.x=data,.f=~length(.x$year)),
           Trend=map(.x=data,.f=~TDA_kendall(.x$value,.x$year)$sen.slope),
           ProbabilityofNegativeSlope=map(.x=data,.f=~TDA_kendall(.x$value,.x$year)$`Probability of Negative Slope`),
           ProbabilityofPostiveSlope=map(.x=data,.f=~TDA_kendall(.x$value,.x$year)$`Probability of Positive Slope`),
           Power_to_Detect_Trend=map(.x=data,.f=~{
             X=.x$value
             n.years<-length(X)
             years<-.x$year
             sample.matrix<-matrix(sample(X,n.years*power.simulations.n,replace=T)+
                                     abs(trend.magnitude.threshold)*(years-years[1]),
                                   power.simulations.n,
                                   n.years,
                                   byrow=T)
             
             tda.out<-plyr::adply(sample.matrix,1,function(samp.row) TDA_kendall(samp.row,years)$`Probability of Positive Slope`)$V1
             #use probability (2/3) threshold to determine if a trend is found
             length(which(tda.out>=(2/3)))/power.simulations.n
           })) %>%
    select(site,n,Trend,ProbabilityofNegativeSlope,ProbabilityofPostiveSlope,Power_to_Detect_Trend) %>%
    unnest() %>%
    mutate(Result=ifelse(Trend>=trend.magnitude.threshold&ProbabilityofPostiveSlope>=(2/3),'Increasing',
                         ifelse(Trend<=(-1*trend.magnitude.threshold)&ProbabilityofNegativeSlope>=(2/3),'Decreasing',
                                ifelse(Power_to_Detect_Trend>=power,'Maintaining',
                                       ifelse(ProbabilityofPostiveSlope>=(2/3)|ProbabilityofNegativeSlope>=(2/3),'Maintaining',
                                       'Inadequate Data')))))
  
}
