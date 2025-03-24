determine_percent_days_stock_went_up <- function(stock_data) {
  closing_prices <- stock_data[,5] # 5 is the index for close_price, for all files used in this project
  total_days <- length(closing_prices)
  days_positive <- 0
  
  for (i in 1:(total_days - 1)) {
    if (closing_prices[i + 1] > closing_prices[i]) {
      days_positive <- days_positive + 1
    }
  }
  
  percentage_of_days_stock_moved_up <- days_positive / total_days
}


avg_daily_abs_change <- function(stock_data) {
  closing_prices <- stock_data[,5]
  changes = rep(NA, times=length(closing_prices))
  
  for (i in 1:(length(closing_prices) - 1)) {
    changes[i] <- abs(closing_prices[i + 1] - closing_prices[i])
  }
  
  mean(changes, na.rm = TRUE)
}


avg_daily_stock_movement_on_positive_days <- function(stock_data) {
  closing_prices <- stock_data[,5]
  changes = rep(NA, times = length(closing_prices))
  
  for (i in 1:(length(closing_prices) - 1)) {
    if (closing_prices[i + 1] > closing_prices[i]) {
      changes[i] <- closing_prices[i + 1] - closing_prices[i]
    }
  }
  
  mean(changes, na.rm = TRUE)
}

avg_daily_stock_movement_on_negative_days <- function(stock_data) {
  closing_prices <- stock_data[,5]
  changes = rep(NA, times = length(closing_prices))
  
  for (i in 1:(length(closing_prices) - 1)) {
    if (closing_prices[i + 1] < closing_prices[i]) {
      changes[i] <- closing_prices[i + 1] - closing_prices[i]
    }
  }
  
  mean(changes, na.rm = TRUE)
}

variance_of_stock <- function(stock_data) {
  var(diff(log(stock_data[,5]))) * 100 
}