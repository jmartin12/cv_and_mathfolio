#change this to determine how many simulations we run
iter = 100000

# change file names for different stocks
data<-read.csv("C:\\Users\\ocorn\\Documents\\MTH-553\\capstone\\TNET.csv")

#change this for each strike price for the stock
strike_price = 80
s0 = 76.74

# used to store all results of simulations
all_prices_in_all_simulations = list()
all_closing_prices = list()

mu <- avg_log_daily_ret(data) 
vol <- variance_of_stock(data) 

for (i in 1:iter) {
  # annualized
  up <- (exp(-mu/252)+ exp((mu+vol^2)/252)+sqrt((exp(-mu/252)+exp((mu+vol^2)/252))^2-4))/2
  dn <- 1/up
  p <- (up*exp(mu/252)-1)/(up^2-1)
  s <- sample(c(-1,1),size=30,replace=TRUE,prob=c(1-p,p))
  simbin<-s0*up^cumsum(c(0,s))
  
  # add individual simulation to result set
  all_prices_in_all_simulations[[i]] <- simbin
  
  # add final closing price to result set
  all_closing_prices[i] <- simbin[31]
}

predicted_close = mean(as.numeric(unlist(all_closing_prices)))
profit = predicted_close - strike_price


all_call_values = list()
for (i in 1:length(all_closing_prices)) {
  if (all_closing_prices[[i]] <= strike_price) {
    all_call_values[i] <- 0
  }
  else {
    all_call_values[i] <- all_closing_prices[[i]] - strike_price
  }
  
}

average_call_value <- mean(as.numeric(unlist(all_call_values)))

 for our graphs
 max_boundary = max(as.numeric(unlist(all_prices_in_all_simulations)))
 min_boundary = min(as.numeric(unlist(all_prices_in_all_simulations)))
 
 for (i in 1:iter) {
   if (i == 1) {
     plot(
       unlist(all_prices_in_all_simulations[i]),
       ylim=c(min_boundary, max_boundary),
       main = "Binomial model - SNAP",
       type="l",
       ylab="Price",
       xlab= "April 2023 -- Days of month",
       col="blue")
   }
   else if (i == 2){
     lines(unlist(all_prices_in_all_simulations[i]), col = "red")
   }
   else {
     lines(unlist(all_prices_in_all_simulations[i]), col = "green")
   }
 }
