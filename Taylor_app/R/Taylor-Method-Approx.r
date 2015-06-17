file.create("/Users/kylemcgrath/Documents/SDE-REU-2015/SUNY-SDE-2015/REU15/Taylor_app/R/data.csv");
data <- file("/Users/kylemcgrath/Documents/SDE-REU-2015/SUNY-SDE-2015/REU15/Taylor_app/R/data.csv");
write(c("Time", "Approximation", "Theoretical", "Experimental"), file = data, ncolumns = 4, sep = ", ");