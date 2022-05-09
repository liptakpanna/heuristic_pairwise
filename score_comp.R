library(tidyverse)
library(DBI)
library(MonetDB.R)

con <- dbConnect(MonetDB.R(), 
                 host="localhost", 
                 dbname="thesis", 
                 port="50000",
                 user="monetdb", 
                 password="monetdb")

getlookup <- function(id1, id2){
  start <- proc.time()
  score <- dbGetQuery(con, paste(
    "select use_lookup( ", id1, "," , id2 ," ,1,-2,-1);", sep='')
  )
  end <- proc.time()
  
  score <- str_split(score, ',', simplify=TRUE)
  result <- list(score=strtoi(score[3]), time=(end-start)[["elapsed"]])
  
  return(result)
}

getneedle <- function(id1, id2){
  start <- proc.time()
  score <- dbGetQuery(con, paste(
    "select * from needleman((select id, seq, 1, -1, -2 from seq_data where id in ( ", id1, "," , id2 ,")));", sep='')
  )
  end <- proc.time()
  
  result <- list(score=score$score, time=(end-start)[["elapsed"]])
  
  return(result)
}

also <- 51
db <- 100
y <- c(also:(also+db-1))
x <- combn(y, 2)

measure_score <- function(felsohatar){
  ossz <- 0
  identical <- 0
  ossztime_lookup <- 0
  ossztime_needle <- 0
  max_score <- 0
  for (i in 1:felsohatar){
    t <- x[,i]
    id1 <- t[1]
    id2 <- t[2]
    #print(paste(id1, ',', id2))
    
    res1 <- getlookup(id1,id2)
    res2 <- getneedle(id1,id2)
    
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
  }
  
  avg = ossz / felsohatar
  print(paste("ATLAG négyzetes eltérés: ", avg))
  print(paste("EGYEZESEK szama", identical))
  print(paste("OSSZ LOOKUP TIME: ", ossztime_lookup))
  print(paste("OSSZ NEEDLE TIME: ", ossztime_needle))
  print(paste("MAXSCORE: ", max_score))
  
  result <- list(avg=avg, identical=identical, time_lookup = ossztime_lookup, time_needle = ossztime_needle)
  return(result)
}

measure_overall <- function(count, limit, k) {
  avgscore <- 0
  identical <- 0
  sumtime_lookup <- 0
  sumtime_needle <- 0
  
  for(i in 1:count) {
    dbSendQuery(con,paste("call update_lookup(1,-2,-1,", k, ");"))
    
    result <- measure_score(limit)
    avgscore = avgscore + result$avg
    identical = identical + result$identical
    sumtime_lookup = sumtime_lookup + result$time_lookup
    sumtime_needle = sumtime_needle + result$time_needle
  }
  print("===================================================")
  print(paste("ATLAG négyzetes eltérés: ", avgscore/count))
  print(paste("ATLAG EGYEZESEK szama", identical/count))
  print(paste("ATLAG OSSZ LOOKUP TIME: ", sumtime_lookup/count))
  print(paste("ATLAG OSSZ NEEDLE TIME: ", sumtime_needle/count))
}


measure_overall(10,4950,5)


