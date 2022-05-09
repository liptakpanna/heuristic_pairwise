library(tidyverse)
library(DBI)
library(MonetDB.R)

con <- dbConnect(MonetDB.R(), 
                 host="localhost", 
                 dbname="thesis", 
                 port="50000",
                 user="monetdb", 
                 password="monetdb")

getlookup_simple <- function(id1, id2, score_match, score_mismatch, score_gap){
  start <- proc.time()
  score <- dbGetQuery(con, paste(
    "select use_lookup_simple( ", id1, "," , id2 ," ,", 
    score_match,',', score_mismatch,',',score_gap,
    ");", sep='')
  )
  end <- proc.time()
  
  score <- str_split(score, ',', simplify=TRUE)
  
  #There was no match from the lookup table, it would fallback on the needleman
  if(length(score) < 3){
    return(list(used=FALSE))
  }
  
  result <- list(score=strtoi(trimws(score[1])), time=(end-start)[["elapsed"]], used=TRUE)
  
  return(result)
}

getneedle <- function(id1, id2, score_match, score_mismatch, score_gap){
  start <- proc.time()
  score <- dbGetQuery(con, paste(
    "select * from needleman((select id, seq,", 
    score_match,',', score_mismatch,',',score_gap,
    " from seq_data where id in ( ", id1, "," , id2 ,")));", sep='')
  )
  end <- proc.time()
  
  #print(score$align1)
  #print(score$align2)
  
  result <- list(score=score$score, time=(end-start)[["elapsed"]])
  
  return(result)
}

#getneedle(3,210)
getlookup_simple(3,20, 1, -1, -2)

also <- 1
db <- 100
y <- c(also:(also+db-1))
x <- combn(y, 2)

measure_used_score <- function(felsohatar,score_match, score_mismatch, score_gap){
  ossz <- 0
  identical <- 0
  ossztime_lookup <- 0
  ossztime_needle <- 0
  max_score <- 0
  used <- 0
  for (i in 1:felsohatar){
    t <- x[,i]
    id1 <- t[1]
    id2 <- t[2]
    #print(paste(id1, ',', id2))
    
    res1 <- getlookup_simple(id1,id2,score_match, score_mismatch, score_gap)
    
    
    if(res1$used == TRUE){
      res2 <- getneedle(id1,id2,score_match, score_mismatch, score_gap)
      
      score1 <- res1$score
      score2 <- res2$score
      
      if(score2 > max_score){
        max_score = score2
      }
      
      ossztime_lookup = ossztime_lookup + res1$time
      ossztime_needle = ossztime_needle + res2$time
      
      #print(paste("RES1: ", score1, " RES2: ", score2, " elteres: ", (score2-score1), " negyzetes: ", (score2-score1)^2))
      
      if(score1 == score2){
        identical = identical + 1
      }
      
      ossz = ossz + (score2-score1)^2
      
      #ossz = ossz + score2-score1
      
      used = used + 1
    }
  }
  
  avg = ossz / used
  print(paste("HANYSZOR HASZNÁLTA A LOOKUPOT: ", used))
  print(paste("ATLAG négyzetes eltérés: ", avg))
  print(paste("EGYEZESEK szama", identical))
  print(paste("OSSZ LOOKUP TIME: ", ossztime_lookup))
  print(paste("OSSZ NEEDLE TIME: ", ossztime_needle))
  print(paste("MAXSCORE: ", max_score))
  
  result <- list(used=used, avg=avg, identical=identical, time_lookup = ossztime_lookup, time_needle = ossztime_needle)
  return(result)
}

library(dplyr)

measure_overall_used <- function(count, limit, k, sample_size,  upperlimit, 
                                 score_match, score_mismatch, score_gap, 
                                 lookup=1, isLog=FALSE, div=100) {
  avgscore <- 0
  identical <- 0
  sumtime_lookup <- 0
  sumtime_needle <- 0
  sumused <- 0
  for(i in 1:count) {
    print(paste("FUTTATÁS ", i))

    if (lookup == 1){
      dbSendQuery(con,paste("call update_lookup(",
                            score_match,',', score_mismatch, ',',score_gap, ','
                            , k, ',', sample_size,',', upperlimit, ");"))
    }
    else if (lookup == 2){
      dbSendQuery(con,paste("call update_lookup2(1,-1,-2,", k, ',', sample_size,',', upperlimit, ");"))
    }
    
    result <- measure_used_score(limit, score_match, score_mismatch, score_gap)
    avgscore = avgscore + result$avg
    identical = identical + result$identical
    sumtime_lookup = sumtime_lookup + result$time_lookup
    sumtime_needle = sumtime_needle + result$time_needle
    sumused = sumused + result$used
    
    if(isLog){
      df <- data.frame(size=c(sample_size),
                       k=c(k),
                       diversity=c(div),
                       used=c(result$used),
                       dist_sq=c(result$avg),
                       dist=c(sqrt(result$avg)),
                       match=c(result$identical),
                       t1=c(result$time_lookup),
                       l_time=c(result$time_lookup/result$used),
                       t2=c(result$time_needle),
                       n_time=c(result$time_needle/result$used)
                      )
      df <- df %>% 
        mutate(across(where(is.numeric), round, 3))
      write.table( df,  
                   file="/home/panna/Projects/thesis/pairwise/meresek_50dist.csv", 
                   append = T, 
                   sep=',', 
                   row.names=F, 
                   col.names=F )
    }
  }
  print("===================================================")
  print(paste("ATLAG négyzetes eltérés: ", avgscore/count))
  print(paste("ATLAG EGYEZESEK szama", identical/count))
  print(paste("ATLAG OSSZ LOOKUP TIME: ", sumtime_lookup/count))
  print(paste("ATLAG OSSZ NEEDLE TIME: ", sumtime_needle/count))
  print(paste("ATLAG USED LOOKUP: ", sumused/count))
}

also <- 1
db <- 100
y <- c(also:(also+db-1))
x <- combn(y, 2)

score_match <- 1
score_mismatch <- -1
score_gap <- -2

samplesize <- 10
k <- 7
lookup <- 1

ks <- c(5,7,10,15)
samplesizes <- c(10,15,20)

for(s in samplesizes) {
  for(k1 in ks){
    measure_overall_used(5, length(x[1,]), k1, s, db, score_match, score_mismatch, score_gap, 1, TRUE)
  }
}

measure_overall_used(1, length(x[1,]), k, samplesize, db,score_match, score_mismatch, score_gap, 1, TRUE)

measure_overall_used(1, length(x[1,]), 40, 2, db,score_match, score_mismatch, score_gap, 2, TRUE)

measure_overall_used(5,4950,7)
#change sample size to 20
measure_overall_used(5,4950,5)
measure_overall_used(5,4950,7)

measure_overall_used(5,4950,10, TRUE, 30, 100)

measure_overall_used(5,4950,5, TRUE, 30, 100)
measure_overall_used(5,4950,7, TRUE, 30, 100)


measure_overall_used(5,44850,10,FALSE)

measure_used_score(100)
