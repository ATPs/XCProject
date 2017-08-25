getwd()
setwd('E:\\Study\\coursera\\Udacity\\ud651DataAnalysisWithR\\')

statesInfo = read.csv('stateData.csv')

subset(statesInfo, state.region == 1)

#dataSet[ROWS,COLUMNS]

stateSubsetBracket = statesInfo[statesInfo$state.region ==1,]
head(stateSubsetBracket,2)
dim(stateSubsetBracket)


#Factor Variables 20160511

reddit = read.csv('E:\\Study\\coursera\\Udacity\\ud651DataAnalysisWithR\\reddit.csv')

table(reddit$employment.status)

summary(reddit)

str(reddit)

levels(reddit$age.range)

library(ggplot2)
qplot(data = reddit, x = age.range)

qplot(data = reddit, x = income.range)

#to order
qplot(data = reddit, x = factor(age.range, levels = c("Under 18", "18-24", "25-34", "35-44", "45-54", "55-64", "65 or Above")))

#or
reddit$age.range = ordered(reddit$age.range, levels = c("Under 18", "18-24", "25-34", "35-44", "45-54", "55-64", "65 or Above"))
#or
reddit$age.range = factor(reddit$age.range, levels = c("Under 18", "18-24", "25-34", "35-44", "45-54", "55-64", "65 or Above"), odered = T)

reddit$income.range = ordered(reddit$income.range, levels = c("Under $20,000", "$20,000 - $29,999", "$30,000 - $39,999","$40,000 - $49,999","$50,000 - $69,999","$70,000 - $99,999", "$150,000 or more"))
qplot(data = reddit, x = income.range)

#20160613
#pseudo Facebook User Data
