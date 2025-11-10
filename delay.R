rates_time$time[2]

t <- round(rates_time$time, 3)
d <- diff(t)
d==d[10]


max(df_avg$time-rates_time$time[1:986])
max(df_avg$time-neff_repl$time[1:986])
