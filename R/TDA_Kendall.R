requireNamespace(c('tidyverse','lubridate','EnvStats'))

TDA_Kendall<-function(conc,time,kendall.conf.level=.95){
  # if(!(trend.direction %in% c('both','positive','negative'))){
  #   stop("trend.direction must be 'positive', 'negative', or 'both'")
  # }
  #require(EnvStats)
  if(missing(time)) time<-1:length(conc)
  trend <- kendallTrendTest(y=conc, x=time, conf.level= kendall.conf.level, correct = FALSE)
  n=length(time)
  V=n*(n-1)*(2*n+5)/18 # this is the same as Var.S
  nC2<-factorial(n)/(factorial(2)*factorial(n-2))
  
  
  likeli_func<-function(x){
    ifelse(x>=.99,'Virtually certain',ifelse(x>=.95,'Extremely likely',ifelse(x>=.9,'Very likely',
                                                                              ifelse(x>=2/3,'Likely',ifelse(x>=1/3,'About as likely as not',ifelse(x>=.1,'Unlikely',
                                                                                                                                                   ifelse(x>=.05,'Very unlikely',ifelse(x>=.01,'Extremely unlikely','Exceptionally unlikely'))))))))
  }
  
  rank_slope_zero<-approx(y=1:nC2,x=trend$slopes,xout=0)$y
  z1a_zero<-(rank_slope_zero*2-nC2)/-sqrt(V)
  neg.prob<-round(1-pnorm(z1a_zero),3)  #negative slope
  pos.prob<-round(pnorm(z1a_zero),3) #positive slope
  
  list(sen.slope=trend$estimate[2],
       kendall.p=as.numeric(round(trend$p.value,4)),
       'Probability of Negative Slope'=neg.prob,
       'Likelihood of Negative Slope'=likeli_func(neg.prob),
       'Probability of Positive Slope'=pos.prob,
       'Likelihood of Positive Slope'=likeli_func(pos.prob)
  )
}