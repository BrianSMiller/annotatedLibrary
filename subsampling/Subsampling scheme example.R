#Subsampling scheme - example
#written by Danielle Harris, 13 Apr 2017

#Assuming 1 year of data, in hourly files and we want to analyse ~900 hours of data

#######################################################################################
##Step 1: Work out total number of hours in the whole dataset
#NB: in ISOdatatime GMT as a timezone is treated as UTC so no daylight savings are observed.  The default is the "current" (I assume the computer's) timezone.
start.date.time<-ISOdatetime(2014, 1, 1, 0, 0, 0, tz="GMT")
end.date.time<-ISOdatetime(2014, 12, 31, 23, 0, 0, tz="GMT") 
total.days<-353
total.hours<-total.days*24

########################################################################################
#Step 2: assume that we will use a 1 hour file as a sampling unit - work out the spacing required to generate 900 hours in a sample.  

no.samples.required<-200
spacing.required<-total.hours/(no.samples.required) 
#note that the seq function in R uses no.samples.required-1, but this is without considering that a random start will be used between 1 and the spacing number.  Therefore, the spacing.required is used as a non-integer number of hours but then will be rounded to the nearest hour across the subsample.

#######################################################################################
#Step 3: ask whether this spacing could create a potentially meaningful temporal pattern in the sampling scheme i.e., aliasing and adjust if necessary.

#9.73 hours is not divisible by 24, so we would not expect any aliasing to occur.  However, we can check this at the end. 

###################################################################################
#Step 4: create a sequence of hours sampled using the above spacing and with a random start point for the first number

#select a random number between 1 and the required spacing. 
random.hour.start<-runif(1,1,spacing.required)

#convert this decimal number of hours to a proper date time object using the start date
random.start.date.time<-start.date.time + as.difftime(random.hour.start, unit="hours")

#need to turn the required spacing into seconds for the sequence function
spacing.required.secs<-spacing.required*60*60

#so this sequence starts with the random start hour on the first day, then increases using the spacing required as a step lengh, and generates the required number of samples.  Finally, the whole sequence is rounded to the nearest integer.
random.sample<-round(seq(random.start.date.time,by = spacing.required.secs,length.out= no.samples.required),"hours")
#note that the seq splits the data by the second but then we round to the nearest hour. This means that sometimes the spacing is 9 hours, sometimes 10...but this helps to break any temporal repetitive patterning.

#extract the hours to check whether there is any temporal aliasing
extract.hours<-format(random.sample, "%H")
extract.months<-format(random.sample, "%m")
extract.days<-format(random.sample, "%d")

summary.hours<-table(extract.hours)
plot(summary.hours,xlab="Sampled hours",ylab="Frequency")

summary.months<-table(extract.months)
plot(summary.months,xlab="Sampled months",ylab="Frequency")

summary.days<-table(extract.days)
plot(summary.days,xlab="Sampled days",ylab="Frequency")

#check the relationship between chosen days, hours, and months.  The most important plot here is the first one of hours against day of the month.  This helps to flag if there is any temporal repetetive patterning to be aware of.
plot(extract.days,extract.hours)
plot(extract.months,extract.days)
plot(extract.months,extract.hours)
