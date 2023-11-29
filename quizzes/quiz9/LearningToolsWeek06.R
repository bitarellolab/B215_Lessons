######## Analysis of Biological Data Labs -- Learning the tools --  Tutorial 6 ######
# Frequency data
# This file contains all of the commands in the Learning the Tools section of tutorial 6. 

#  Confidence intervals for a proportion

# install.packages("binom", dependencies = TRUE)
# Remove the # from the previous line if you need to install the package "binom" for the first time.

library(binom)

binom.confint(x = 30, n = 87, method = "ac")



#  Binomial test

binom.test(x = 14, n = 18, p = 0.5)


# Goodness of fit

# table
MMlist <- read.csv("DataForLabs/MandMlist.csv")
MMlist$color

MMtable <- table(MMlist$color)
MMtable

#  chisq.test()

expected_proportions <- c(0.24, 0.14, 0.16, 0.20, 0.13, 0.13)

sum(MMtable)

55 * expected_proportions

chisq.test(MMtable, p = expected_proportions)


#  Goodness of fit to a Poisson distribution

extinctData <- read.csv("DataForLabs/MassExtinctions.csv", stringsAsFactors = TRUE)
head(extinctData)

number_of_extinctions <- extinctData$numberOfExtinctions

table(number_of_extinctions)

#Generating expected frequencies
mean_extinctions <- mean(number_of_extinctions)
mean_extinctions

dpois(x = 3, lambda = 4.21)

0:20

expected_probability = dpois(x = 0:20, lambda = 4.21)
expected_probability 

length(number_of_extinctions)

76 * expected_probability 

expected_combined <- c(5.878568, 9.999264, 14.032300, 14.768996, 12.435494,  8.725572, 5.247808, 4.911998)
observed_combined <- c(13, 15, 16, 7, 10, 4, 2, 9)

# chi-squared statistic
chisq.test(observed_combined, p = expected_combined, rescale.p = TRUE)$statistic

#calculating the p-value
pchisq(q = 23.939, df = 6, lower.tail = FALSE)



####  Making frequency distribution plot for Activity

freqs <- dbinom(0:5, size = 5, p = 1/6)
num3s <- 0:5
diceFreqTable <-data.frame(num3s,freqs)

ggplot(data = diceFreqTable, aes(x=num3s, y=freqs) )+ 
  geom_bar(colour="black", fill="orange", width=1, stat="identity") + 
  xlab("Number of 3's") + ylab("Predicted frequency") +
  theme_classic() + 
  scale_x_continuous(breaks=seq(0, 5, 1)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"))