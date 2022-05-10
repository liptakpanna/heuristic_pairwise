library(tidyverse)
library(DBI)
library(MonetDB.R)

con <- dbConnect(MonetDB.R(), 
                 host="localhost", 
                 dbname="thesis", 
                 port="50000",
                 user="monetdb", 
                 password="monetdb")

getlookup_simple_affine <- function(id1, id2, score_match, score_mismatch, score_gap_open, score_gap_extend){
  start <- proc.time()
  score <- dbGetQuery(con, paste(
    "select use_lookup_simple_affine( ", id1, "," , id2 ," ,", 
    score_match,',', score_mismatch,',',score_gap_open,',',score_gap_extend,
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

get_needle_affine <- function(id1, id2, score_match, score_mismatch, score_gap_open, score_gap_extend){
  res <- dbGetQuery(con, paste("select score, time from alignments_help_affine where (id1 =", id1, " and id2=", id2, ") or (id2 =",id1, "and id1 =",id2, ");"));
  if(nrow(res) > 0){
      return(list(score=res$score, time=res$time))
  }
  
  start <- proc.time()
  score <- dbGetQuery(con, paste(
    "select * from needleman_affine((select id, seq, ", 
    score_match,',', score_mismatch,',',score_gap_open,',',score_gap_extend,"
    from seq_data where id in ( ", id1, "," , id2 ,")));", sep='')
  )
  end <- proc.time()
  
  t <- (end-start)[["elapsed"]]

  dbSendQuery(con, paste("insert into alignments_help_affine values (",id1, ",", id2, ",", score$score, ",", t,");"))
  
  return(list(score=score$score, time=t))
}

measure_used_score_affine <- function(count, score_match, score_mismatch, score_gap_open, score_gap_extend){
  scorediff <- 0
  identical <- 0
  time_lookup <- 0
  time_needle <- 0
  max_score <- 0
  used <- 0
  for (i in 1:count){
    t <- x[,i]
    id1 <- t[1]
    id2 <- t[2]
    
    res1 <- getlookup_simple_affine(id1,id2, score_match, score_mismatch, score_gap_open, score_gap_extend)
    
    
    if(res1$used == TRUE){
      res2 <- get_needle_affine(id1,id2, score_match, score_mismatch, score_gap_open, score_gap_extend)
      
      score1 <- res1$score
      score2 <- res2$score
      
      if(score2 > max_score){
        max_score = score2
      }
      
      time_lookup = time_lookup + res1$time
      time_needle = time_needle + res2$time
            
      if(score1 == score2){
        identical = identical + 1
      }
      
      scorediff = scorediff + (score2-score1)^2
      
      used = used + 1
    }
  }
  
  avg = scorediff / used
  print(paste("USED LOOKUP: ", used))
  print(paste("AVG SQUARE SUM SCORE DIFF: ", avg))
  print(paste("IDENTICAL", identical))
  print(paste("SUM LOOKUP TIME: ", time_lookup))
  print(paste("SUM NEEDLE TIME: ", time_needle))
  print(paste("MAXSCORE: ", max_score))
  
  result <- list(used=used, avg=avg, identical=identical, time_lookup = time_lookup, time_needle = time_needle)
  return(result)
}

measure_overall_used_affine <- function(count, limit, k, sample_size,  upperlimit, 
                                        score_match, score_mismatch, score_gap_open, score_gap_extend,
                                        isLog=FALSE, div=100) {
  avgscore <- 0
  identical <- 0
  sumtime_lookup <- 0
  sumtime_needle <- 0
  sumused <- 0
  t <- 0
  for(i in 1:count) {
    print(paste("LOOP ", i, ", ", Sys.time()))
    
    start <- proc.time()

    dbSendQuery(con,paste("call update_lookup_affine(", score_match,',', score_mismatch,',',score_gap_open,',',score_gap_extend,','
                          , k, ',', sample_size,',', upperlimit, ");"))

    end <- proc.time()
    t = t + (end-start)[["elapsed"]]

    result <- measure_used_score_affine(limit,score_match, score_mismatch, score_gap_open, score_gap_extend)
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
                   file=<filepath>, 
                   append = T, 
                   sep=',', 
                   row.names=F, 
                   col.names=F )
    }
  }
  print("===================================================")
  print(paste("AVG SQUARE DIFF: ", avgscore/count))
  print(paste("AVG IDENTICAL", identical/count))
  print(paste("AVG SUM LOOKUP TIME: ", sumtime_lookup/count))
  print(paste("AVG SUM NEEDLE TIME: ", sumtime_needle/count))
  print(paste("AVG USED LOOKUP: ", sumused/count))
  print(paste("AVG UPDATE LOOKUP TIME: ", t/count))
}

minid <- 1
maxid <- 100
y <- c(minid:(minid+maxid-1))
x <- combn(y, 2)
score_match <- 1
score_mismatch <- -1
score_gap_open <- -3
score_gap_extend <- -1

for(s in samplesizes) {
  for(k in ks){
    measure_overall_used_affine(5, length(x[1,]), k, s, db, score_match, score_mismatch, score_gap_open, score_gap_extend, TRUE)
  }
}