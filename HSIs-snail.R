# Clear all variables before new model run #
rm(list=ls(all=TRUE))

## Read in depth info
d.400 <- as.matrix(read.csv("Depths-400.csv", header = T))
d.800 <- as.matrix(read.csv("Depths-800.csv", header = T))
d.1600 <- as.matrix(read.csv("Depths-1600.csv", header = T))
d.3200 <- as.matrix(read.csv("Depths-3200.csv", header = T))
d.6400 <- as.matrix(read.csv("Depths-6400.csv", header = T))

## Create empty natrices for HSI values
s.400 <- matrix(nrow = 4015, ncol = 3502)
s.800 <- matrix(nrow = 4015, ncol = 909)
s.1600 <- matrix(nrow = 4015, ncol = 239)
s.3200 <- matrix(nrow = 4015, ncol = 65)
s.6400 <- matrix(nrow = 4015, ncol = 20)

### APPLE SNAILS ###

# .sn denotes snail suitabilities
s.sn.400 <- s.400
s.sn.800 <- s.800
s.sn.1600 <- s.1600
s.sn.3200 <- s.3200
s.sn.6400 <- s.6400

# Daily totals are summed for each year for each grid
t.sn.400 <- matrix(nrow = 11, ncol = ncol(s.400))
t.sn.800 <- matrix(nrow = 11, ncol = ncol(s.800))
t.sn.1600 <- matrix(nrow = 11, ncol = ncol(s.1600))
t.sn.3200 <- matrix(nrow = 11, ncol = ncol(s.3200))
t.sn.6400 <- matrix(nrow = 11, ncol = ncol(s.6400))

## Apple Snail parameters

# Egg-laying season
stt.day <- 61
end.day <- 306

# Egg production factors
ep.1 <- 0.00176 # for days 61-91
ep.2 <- 0.01053 # for days 92-121
ep.3 <- 0.00702 # for days 122-152
ep.4 <- 0.00527 # for days 153-182
ep.5 <- 0.00351 # for days 183-214
ep.6 <- 0.00246 # for days 215-245
ep.7 <- 0.00140 # for days 246-275
ep.8 <- 0.00070 # for days 276-306

tick <- Sys.time()
print(Sys.time())
## Begin day-by-day analysis of Apple Snails at 6400m
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.6400)) {
      
      if (d.6400[((y-1)*365 + d), g] >= 0.1) { # if water depth is 10cm+, conditions are considered suitable
        
        s.sn.6400[((y-1)*365 + d), g] <- 1
        
      } else {
        
        s.sn.6400[((y-1)*365 + d), g] <- 0
        
      }
      
      if (d.6400[((y-1)*365 + d), g] > (0.2 + max(d.6400[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]))) { # if today's depth is more than 20cm above the high for the previous 30 days, young from that period are flooded out
        
        s.sn.6400[((y-1)*365 + d - 31):((y-1)*365 + d - 1), g] <- 0
        
      }
      
      if (max(d.6400[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 30 days, all recruits less than a month old at the start of the drought die
        
        s.sn.6400[((y-1)*365 + d - 60):((y-1)*365 + d - 31), g] <- 0
        
      }
      
      if (max(d.6400[((y-1)*365 + d - 61):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 60 days, all recruits less than two months old at the start of the drought die
        
        s.sn.6400[((y-1)*365 + d - 120):((y-1)*365 + d - 91), g] <- 0
        
      }
      
    }
    
  }
  
}

## Multiply daily 0 or 1 values by egg production factors
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.6400)) {
      
      if (d >= 61 && d <= 91) {
        
        s.sn.6400[((y-1)*365 + d), g] <- s.sn.6400[((y-1)*365 + d), g]*ep.1
        
      }
      
      if (d >= 92 && d <= 121) {
        
        s.sn.6400[((y-1)*365 + d), g] <- s.sn.6400[((y-1)*365 + d), g]*ep.2
        
      }
      
      if (d >= 122 && d <= 152) {
        
        s.sn.6400[((y-1)*365 + d), g] <- s.sn.6400[((y-1)*365 + d), g]*ep.3
        
      }
      
      if (d >= 153 && d <= 182) {
        
        s.sn.6400[((y-1)*365 + d), g] <- s.sn.6400[((y-1)*365 + d), g]*ep.4
        
      }
      
      if (d >= 183 && d <= 214) {
        
        s.sn.6400[((y-1)*365 + d), g] <- s.sn.6400[((y-1)*365 + d), g]*ep.5
        
      }
      
      if (d >= 215 && d <= 245) {
        
        s.sn.6400[((y-1)*365 + d), g] <- s.sn.6400[((y-1)*365 + d), g]*ep.6
        
      }
      
      if (d >= 246 && d <= 275) {
        
        s.sn.6400[((y-1)*365 + d), g] <- s.sn.6400[((y-1)*365 + d), g]*ep.7
        
      }
      
      if (d >= 276 && d <= 306) {
        
        s.sn.6400[((y-1)*365 + d), g] <- s.sn.6400[((y-1)*365 + d), g]*ep.8
        
      }
      
    }
    
  }
  
}

## Sum yearly totals for each grid
for (y in 1:11) {
  
  for (g in 1:ncol(s.6400)) {
    
      t.sn.6400[y, g] <- round(sum(s.sn.6400[((y-1)*365 + stt.day):((y-1)*365 + end.day), g]), 2)
    
  }
  
}
r.avg <- matrix(nrow = 11, ncol = 1); c.avg <- matrix(nrow = 1, ncol = ncol(s.6400))

for (r in 1:11) {
  
  r.avg[r,1] <- round(mean(t.sn.6400[r,]), 2)
  
}

for (c in 1:ncol(s.6400)) {
  
  c.avg[1,c] <- round(mean(t.sn.6400[,c]), 2)
  
}


write.table(t.sn.6400, file = "Apple_Snail-6400.csv", sep = ",", 
            row.names = c(2005:2015), col.names = c(1:ncol(s.6400)))

## Begin day-by-day analysis of Apple Snails at 3200m
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.3200)) {
      
      if (d.3200[((y-1)*365 + d), g] >= 0.1) { # if water depth is 10cm+, conditions are considered suitable
        
        s.sn.3200[((y-1)*365 + d), g] <- 1
        
      } else {
        
        s.sn.3200[((y-1)*365 + d), g] <- 0
        
      }
      
      if (d.3200[((y-1)*365 + d), g] > (0.2 + max(d.3200[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]))) { # if today's depth is more than 20cm above the high for the previous 30 days, young from that period are flooded out
        
        s.sn.3200[((y-1)*365 + d - 31):((y-1)*365 + d - 1), g] <- 0
        
      }
      
      if (max(d.3200[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 30 days, all recruits less than a month old at the start of the drought die
        
        s.sn.3200[((y-1)*365 + d - 60):((y-1)*365 + d - 31), g] <- 0
        
      }
      
      if (max(d.3200[((y-1)*365 + d - 61):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 60 days, all recruits less than two months old at the start of the drought die
        
        s.sn.3200[((y-1)*365 + d - 120):((y-1)*365 + d - 91), g] <- 0
        
      }
      
    }
    
  }
  
}

## Multiply daily 0 or 1 values by egg production factors
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.3200)) {
      
      if (d >= 61 && d <= 91) {
        
        s.sn.3200[((y-1)*365 + d), g] <- s.sn.3200[((y-1)*365 + d), g]*ep.1
        
      }
      
      if (d >= 92 && d <= 121) {
        
        s.sn.3200[((y-1)*365 + d), g] <- s.sn.3200[((y-1)*365 + d), g]*ep.2
        
      }
      
      if (d >= 122 && d <= 152) {
        
        s.sn.3200[((y-1)*365 + d), g] <- s.sn.3200[((y-1)*365 + d), g]*ep.3
        
      }
      
      if (d >= 153 && d <= 182) {
        
        s.sn.3200[((y-1)*365 + d), g] <- s.sn.3200[((y-1)*365 + d), g]*ep.4
        
      }
      
      if (d >= 183 && d <= 214) {
        
        s.sn.3200[((y-1)*365 + d), g] <- s.sn.3200[((y-1)*365 + d), g]*ep.5
        
      }
      
      if (d >= 215 && d <= 245) {
        
        s.sn.3200[((y-1)*365 + d), g] <- s.sn.3200[((y-1)*365 + d), g]*ep.6
        
      }
      
      if (d >= 246 && d <= 275) {
        
        s.sn.3200[((y-1)*365 + d), g] <- s.sn.3200[((y-1)*365 + d), g]*ep.7
        
      }
      
      if (d >= 276 && d <= 306) {
        
        s.sn.3200[((y-1)*365 + d), g] <- s.sn.3200[((y-1)*365 + d), g]*ep.8
        
      }
      
    }
    
  }
  
}

## Sum yearly totals for each grid
for (y in 1:11) {
  
  for (g in 1:ncol(s.3200)) {
    
    t.sn.3200[y, g] <- round(sum(s.sn.3200[((y-1)*365 + stt.day):((y-1)*365 + end.day), g]), 2)
    
  }
  
}
r.avg <- matrix(nrow = 11, ncol = 1); c.avg <- matrix(nrow = 1, ncol = ncol(s.3200))

for (r in 1:11) {
  
  r.avg[r,1] <- round(mean(t.sn.3200[r,]), 2)
  
}

for (c in 1:ncol(s.3200)) {
  
  c.avg[1,c] <- round(mean(t.sn.3200[,c]), 2)
  
}


write.table(t.sn.3200, file = "Apple_Snail-3200.csv", sep = ",", 
            row.names = c(2005:2015), col.names = c(1:ncol(s.3200)))

## Begin day-by-day analysis of Apple Snails at 1600m
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.1600)) {
      
      if (d.1600[((y-1)*365 + d), g] >= 0.1) { # if water depth is 10cm+, conditions are considered suitable
        
        s.sn.1600[((y-1)*365 + d), g] <- 1
        
      } else {
        
        s.sn.1600[((y-1)*365 + d), g] <- 0
        
      }
      
      if (d.1600[((y-1)*365 + d), g] > (0.2 + max(d.1600[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]))) { # if today's depth is more than 20cm above the high for the previous 30 days, young from that period are flooded out
        
        s.sn.1600[((y-1)*365 + d - 31):((y-1)*365 + d - 1), g] <- 0
        
      }
      
      if (max(d.1600[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 30 days, all recruits less than a month old at the start of the drought die
        
        s.sn.1600[((y-1)*365 + d - 60):((y-1)*365 + d - 31), g] <- 0
        
      }
      
      if (max(d.1600[((y-1)*365 + d - 61):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 60 days, all recruits less than two months old at the start of the drought die
        
        s.sn.1600[((y-1)*365 + d - 120):((y-1)*365 + d - 91), g] <- 0
        
      }
      
    }
    
  }
  
}

## Multiply daily 0 or 1 values by egg production factors
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.1600)) {
      
      if (d >= 61 && d <= 91) {
        
        s.sn.1600[((y-1)*365 + d), g] <- s.sn.1600[((y-1)*365 + d), g]*ep.1
        
      }
      
      if (d >= 92 && d <= 121) {
        
        s.sn.1600[((y-1)*365 + d), g] <- s.sn.1600[((y-1)*365 + d), g]*ep.2
        
      }
      
      if (d >= 122 && d <= 152) {
        
        s.sn.1600[((y-1)*365 + d), g] <- s.sn.1600[((y-1)*365 + d), g]*ep.3
        
      }
      
      if (d >= 153 && d <= 182) {
        
        s.sn.1600[((y-1)*365 + d), g] <- s.sn.1600[((y-1)*365 + d), g]*ep.4
        
      }
      
      if (d >= 183 && d <= 214) {
        
        s.sn.1600[((y-1)*365 + d), g] <- s.sn.1600[((y-1)*365 + d), g]*ep.5
        
      }
      
      if (d >= 215 && d <= 245) {
        
        s.sn.1600[((y-1)*365 + d), g] <- s.sn.1600[((y-1)*365 + d), g]*ep.6
        
      }
      
      if (d >= 246 && d <= 275) {
        
        s.sn.1600[((y-1)*365 + d), g] <- s.sn.1600[((y-1)*365 + d), g]*ep.7
        
      }
      
      if (d >= 276 && d <= 306) {
        
        s.sn.1600[((y-1)*365 + d), g] <- s.sn.1600[((y-1)*365 + d), g]*ep.8
        
      }
      
    }
    
  }
  
}

## Sum yearly totals for each grid
for (y in 1:11) {
  
  for (g in 1:ncol(s.1600)) {
    
    t.sn.1600[y, g] <- round(sum(s.sn.1600[((y-1)*365 + stt.day):((y-1)*365 + end.day), g]), 2)
    
  }
  
}
r.avg <- matrix(nrow = 11, ncol = 1); c.avg <- matrix(nrow = 1, ncol = ncol(s.1600))

for (r in 1:11) {
  
  r.avg[r,1] <- round(mean(t.sn.1600[r,]), 2)
  
}

for (c in 1:ncol(s.1600)) {
  
  c.avg[1,c] <- round(mean(t.sn.1600[,c]), 2)
  
}


write.table(t.sn.1600, file = "Apple_Snail-1600.csv", sep = ",", 
            row.names = c(2005:2015), col.names = c(1:ncol(s.1600)))

## Begin day-by-day analysis of Apple Snails at 800m
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.800)) {
      
      if (d.800[((y-1)*365 + d), g] >= 0.1) { # if water depth is 10cm+, conditions are considered suitable
        
        s.sn.800[((y-1)*365 + d), g] <- 1
        
      } else {
        
        s.sn.800[((y-1)*365 + d), g] <- 0
        
      }
      
      if (d.800[((y-1)*365 + d), g] > (0.2 + max(d.800[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]))) { # if today's depth is more than 20cm above the high for the previous 30 days, young from that period are flooded out
        
        s.sn.800[((y-1)*365 + d - 31):((y-1)*365 + d - 1), g] <- 0
        
      }
      
      if (max(d.800[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 30 days, all recruits less than a month old at the start of the drought die
        
        s.sn.800[((y-1)*365 + d - 60):((y-1)*365 + d - 31), g] <- 0
        
      }
      
      if (max(d.800[((y-1)*365 + d - 61):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 60 days, all recruits less than two months old at the start of the drought die
        
        s.sn.800[((y-1)*365 + d - 120):((y-1)*365 + d - 91), g] <- 0
        
      }
      
    }
    
  }
  
}

## Multiply daily 0 or 1 values by egg production factors
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.800)) {
      
      if (d >= 61 && d <= 91) {
        
        s.sn.800[((y-1)*365 + d), g] <- s.sn.800[((y-1)*365 + d), g]*ep.1
        
      }
      
      if (d >= 92 && d <= 121) {
        
        s.sn.800[((y-1)*365 + d), g] <- s.sn.800[((y-1)*365 + d), g]*ep.2
        
      }
      
      if (d >= 122 && d <= 152) {
        
        s.sn.800[((y-1)*365 + d), g] <- s.sn.800[((y-1)*365 + d), g]*ep.3
        
      }
      
      if (d >= 153 && d <= 182) {
        
        s.sn.800[((y-1)*365 + d), g] <- s.sn.800[((y-1)*365 + d), g]*ep.4
        
      }
      
      if (d >= 183 && d <= 214) {
        
        s.sn.800[((y-1)*365 + d), g] <- s.sn.800[((y-1)*365 + d), g]*ep.5
        
      }
      
      if (d >= 215 && d <= 245) {
        
        s.sn.800[((y-1)*365 + d), g] <- s.sn.800[((y-1)*365 + d), g]*ep.6
        
      }
      
      if (d >= 246 && d <= 275) {
        
        s.sn.800[((y-1)*365 + d), g] <- s.sn.800[((y-1)*365 + d), g]*ep.7
        
      }
      
      if (d >= 276 && d <= 306) {
        
        s.sn.800[((y-1)*365 + d), g] <- s.sn.800[((y-1)*365 + d), g]*ep.8
        
      }
      
    }
    
  }
  
}

## Sum yearly totals for each grid
for (y in 1:11) {
  
  for (g in 1:ncol(s.800)) {
    
    t.sn.800[y, g] <- round(sum(s.sn.800[((y-1)*365 + stt.day):((y-1)*365 + end.day), g]), 2)
    
  }
  
}
r.avg <- matrix(nrow = 11, ncol = 1); c.avg <- matrix(nrow = 1, ncol = ncol(s.800))

for (r in 1:11) {
  
  r.avg[r,1] <- round(mean(t.sn.800[r,]), 2)
  
}

for (c in 1:ncol(s.800)) {
  
  c.avg[1,c] <- round(mean(t.sn.800[,c]), 2)
  
}


write.table(t.sn.800, file = "Apple_Snail-800.csv", sep = ",", 
            row.names = c(2005:2015), col.names = c(1:ncol(s.800)))

## Begin day-by-day analysis of Apple Snails at 400m
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.400)) {
      
      if (d.400[((y-1)*365 + d), g] >= 0.1) { # if water depth is 10cm+, conditions are considered suitable
        
        s.sn.400[((y-1)*365 + d), g] <- 1
        
      } else {
        
        s.sn.400[((y-1)*365 + d), g] <- 0
        
      }
      
      if (d.400[((y-1)*365 + d), g] > (0.2 + max(d.400[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]))) { # if today's depth is more than 20cm above the high for the previous 30 days, young from that period are flooded out
        
        s.sn.400[((y-1)*365 + d - 31):((y-1)*365 + d - 1), g] <- 0
        
      }
      
      if (max(d.400[((y-1)*365 + d - 31):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 30 days, all recruits less than a month old at the start of the drought die
        
        s.sn.400[((y-1)*365 + d - 60):((y-1)*365 + d - 31), g] <- 0
        
      }
      
      if (max(d.400[((y-1)*365 + d - 61):((y-1)*365 + d  - 1), g]) <= 0) { # if the ground has been dry for at least 60 days, all recruits less than two months old at the start of the drought die
        
        s.sn.400[((y-1)*365 + d - 120):((y-1)*365 + d - 91), g] <- 0
        
      }
      
    }
    
  }
  
}

## Multiply daily 0 or 1 values by egg production factors
for (y in 1:11) {
  
  for (d in stt.day:end.day) {
    
    for (g in 1:ncol(s.400)) {
      
      if (d >= 61 && d <= 91) {
        
        s.sn.400[((y-1)*365 + d), g] <- s.sn.400[((y-1)*365 + d), g]*ep.1
        
      }
      
      if (d >= 92 && d <= 121) {
        
        s.sn.400[((y-1)*365 + d), g] <- s.sn.400[((y-1)*365 + d), g]*ep.2
        
      }
      
      if (d >= 122 && d <= 152) {
        
        s.sn.400[((y-1)*365 + d), g] <- s.sn.400[((y-1)*365 + d), g]*ep.3
        
      }
      
      if (d >= 153 && d <= 182) {
        
        s.sn.400[((y-1)*365 + d), g] <- s.sn.400[((y-1)*365 + d), g]*ep.4
        
      }
      
      if (d >= 183 && d <= 214) {
        
        s.sn.400[((y-1)*365 + d), g] <- s.sn.400[((y-1)*365 + d), g]*ep.5
        
      }
      
      if (d >= 215 && d <= 245) {
        
        s.sn.400[((y-1)*365 + d), g] <- s.sn.400[((y-1)*365 + d), g]*ep.6
        
      }
      
      if (d >= 246 && d <= 275) {
        
        s.sn.400[((y-1)*365 + d), g] <- s.sn.400[((y-1)*365 + d), g]*ep.7
        
      }
      
      if (d >= 276 && d <= 306) {
        
        s.sn.400[((y-1)*365 + d), g] <- s.sn.400[((y-1)*365 + d), g]*ep.8
        
      }
      
    }
    
  }
  
}

## Sum yearly totals for each grid
for (y in 1:11) {
  
  for (g in 1:ncol(s.400)) {
    
    t.sn.400[y, g] <- round(sum(s.sn.400[((y-1)*365 + stt.day):((y-1)*365 + end.day), g]), 2)
    
  }
  
}

write.table(t.sn.400, file = "Apple_Snail-400.csv", sep = ",", 
            row.names = c(2005:2015), col.names = c(1:ncol(s.400)))

