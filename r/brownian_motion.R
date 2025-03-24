#change this to determine how many simulations we run
iter = 1000000
# change file names for different stocks
data<-read.csv("C:\\Users\\ocorn\\Documents\\MTH-553\\capstone\\SNAP.csv") 

#change this for each strike price for the stock
strike_price = 14
s0 = 11.74

# used to store all results of simulations
all_prices_in_all_simulations = list()
all_closing_prices = list()

# divide by 100 because the function definition converts the daily returns to a percentage basis
mean <- avg_log_daily_ret(data)/100
stdev <- sqrt(variance_of_stock(data))/100

for (i in 1:iter) {
  # annualized
  simreturns <- rnorm(30, mean=mean/252, sd=stdev)
  simcumreturns <- cumsum(c(0, simreturns))
  simstock <- s0* exp(simcumreturns)
  
  all_prices_in_all_simulations[[i]] <- simstock
  all_closing_prices[i] <- simstock[31]
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
 
  max_boundary = max(as.numeric(unlist(all_prices_in_all_simulations)))
  min_boundary = min(as.numeric(unlist(all_prices_in_all_simulations)))
  
  for (i in 1:iter) {
    if (i == 1) {
      plot(
        unlist(all_prices_in_all_simulations[i]),
        ylim=c(min_boundary, max_boundary),
        main = "Geometric brownian motion - SNAP",
        type="l",
        ylab="Price",
        xlab= "April 2023 - Day of month",
        col="blue")
    }
    else if (i == 2){
      lines(unlist(all_prices_in_all_simulations[i]), col = "red")
    }
    else {
      lines(unlist(all_prices_in_all_simulations[i]), col = "green")
    }
  }
