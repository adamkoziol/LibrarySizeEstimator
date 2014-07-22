#!/usr/bin/env Rscript

# http://www.cbs.dtu.dk/courses/27626/Exercises/Alignment_exercise.php
# samtools view HG00418_A.bam | cut -f9 > OLC-795_insertsizes.txt

library("ggplot2", lib.loc="/usr/local/lib/R/site-library")

cmdArgs <- commandArgs(trailingOnly=TRUE)
name <- cmdArgs[2]

csvFile <- paste(cmdArgs[1], "/", name, "_insertsizes.csv", sep='')

a = read.csv(csvFile)

a.v = a[a[,1]>0,1]

b = data.frame(x = a[a[,1]>0,1])
c = data.frame(x = b[b[,1]<1000,1]) 
# hospitalName <- append(hospitalName, sortedResults[ranking, 1], after = length(hospitalName))

# Lets get rid of outliers and use the 5% and 95% intervals to calculate mean and standard deviation:
mn = quantile(a.v, seq(0,1,0.05))[2]
mx = quantile(a.v, seq(0,1,0.05))[20]
tmp1 <- mean(a.v[a.v >= mn & a.v <= mx])    # Mean
tmp2 <- sd(a.v[a.v >= mn & a.v <= mx])      # SD

annotateMean <- paste("Mean: " , round(tmp1, digits=1), " bp")
annotateSD <- paste("Standard Deviation" , round(tmp2, digits=1), " bp")
title <- paste("Distribution of library insert sizes for sequencing run of strain", cmdArgs[2])

p <-ggplot(c, aes(x=c[,1])) 
p + geom_histogram(binwidth=20, fill="white", color="black") + 
  annotate("text", x=750, y=Inf, vjust = 6, label = annotateMean) + 
  annotate("text", x=750, y=Inf, vjust = 7.5, label = annotateSD) +
  scale_x_continuous(name="Insert size (bp)") + 
  scale_y_continuous(name="Frequency") +
  ggtitle(title)

filename <- paste(name, "_insert_sizes.pdf", sep='')

# saves the chart as insert_sizes.pdf, in landscape formate with letter paper dimensions
ggsave(filename, width=11, height=8.5)

textOutput <- paste(name, "_insert_sizes.txt", sep='')
sink(textOutput) + cat(name) + cat("\t") + cat(tmp1) + cat("\t") + cat(tmp2)  + sink()
