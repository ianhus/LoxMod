# Clear all variables before new model run #
rm(list=ls(all=TRUE))
#tick <- Sys.time()
#print(Sys.time())

## Set the working directory
setwd("E:/#pubz/4. Eco")

## Read in depth info
d.400 <- as.matrix(read.csv("Depths-400.csv", header = T))
d.800 <- as.matrix(read.csv("Depths-800.csv", header = T))
d.1600 <- as.matrix(read.csv("Depths-1600.csv", header = T))
d.3200 <- as.matrix(read.csv("Depths-3200.csv", header = T))
d.6400 <- as.matrix(read.csv("Depths-6400.csv", header = T))

## Create empty natrices for HSI values
s.400 <- matrix(nrow = 17*11, ncol = 3502)
s.800 <- matrix(nrow = 17*11, ncol = 909)
s.1600 <- matrix(nrow = 17*11, ncol = 239)
s.3200 <- matrix(nrow = 17*11, ncol = 65)
s.6400 <- matrix(nrow = 17*11, ncol = 20)

### WADING BIRDS ###

# .wd denotes water depth suitability
s.wd.400 <- s.400
s.wd.800 <- s.800
s.wd.1600 <- s.1600
s.wd.3200 <- s.3200
s.wd.6400 <- s.6400

# .dh denotes recession rate suitabilities
s.dh.400 <- s.400[1:(16*11),]
s.dh.800 <- s.800[1:(16*11),]
s.dh.1600 <- s.1600[1:(16*11),]
s.dh.3200 <- s.3200[1:(16*11),]
s.dh.6400 <- s.6400[1:(16*11),]

# Minimum of wd and dh is used
s.wb.400 <- matrix(nrow = nrow(s.dh.400), ncol = ncol(s.dh.400))
s.wb.800 <- matrix(nrow = nrow(s.dh.800), ncol = ncol(s.dh.800))
s.wb.1600 <- matrix(nrow = nrow(s.dh.1600), ncol = ncol(s.dh.1600))
s.wb.3200 <- matrix(nrow = nrow(s.dh.3200), ncol = ncol(s.dh.3200))
s.wb.6400 <- matrix(nrow = nrow(s.dh.6400), ncol = ncol(s.dh.6400)) 

# The top 23% of grids for each week is considered the landscape suitability
wb.ls.400 <- matrix(nrow = nrow(s.dh.400), ncol = 1)
wb.ls.800 <- matrix(nrow = nrow(s.dh.800), ncol = 1)
wb.ls.1600 <- matrix(nrow = nrow(s.dh.1600), ncol = 1)
wb.ls.3200 <- matrix(nrow = nrow(s.dh.3200), ncol = 1)
wb.ls.6400 <- matrix(nrow = nrow(s.dh.6400), ncol = 1)

# Wood Stork suitability gets summed to one value a year
t.ws.400 <- matrix(nrow = 11, ncol = 1)
t.ws.800 <- matrix(nrow = 11, ncol = 1)
t.ws.1600 <- matrix(nrow = 11, ncol = 1)
t.ws.3200 <- matrix(nrow = 11, ncol = 1)
t.ws.6400 <- matrix(nrow = 11, ncol = 1)

# White Ibis and small heron suitability gets summed to one value a year
t.ih.400 <- matrix(nrow = 11, ncol = 1)
t.ih.800 <- matrix(nrow = 11, ncol = 1)
t.ih.1600 <- matrix(nrow = 11, ncol = 1)
t.ih.3200 <- matrix(nrow = 11, ncol = 1)
t.ih.6400 <- matrix(nrow = 11, ncol = 1)

tick <- Sys.time()
print(Sys.time())

## Begin week-by-week analysis of Wading Birds at 6400m
for (y in 1:11) {
  
  for (w in 1:17) { # HSIs for wading birds go from Jan-Apr, approximately 17 weeks
    
    for (g in 1:ncol(s.6400)) {
      
      s.wd.6400[(y-1)*17 + w, g] <- mean(d.6400[((y-1)*365 + (w - 1)*7 + 1):((y-1)*365 + w*7), g]) # find average weekly water depth
      
      if (w > 1) {
      
        s.dh.6400[(y-1)*16 + w - 1, g] <- s.wd.6400[(y-1)*17 + w, g] - s.wd.6400[(y-1)*17 + w - 1, g] # find change in average water depth (since this requires a previous week, I'm only doing it for weeks 2-17)
          
      }
      
    }
    
  }
  
}

## Determine HSI_depth for 6400 based on formula
for (w in 1:nrow(s.wd.6400)) {
  
  for (g in 1:ncol(s.wd.6400)) {
    
    if (s.wd.6400[w, g] <= -0.3/3.28) {
      
      s.wd.6400[w, g] <- 0; next()
      
    }
    
    if (s.wd.6400[w, g] > -0.3/3.28 && s.wd.6400[w, g] <= 0) {
      
      s.wd.6400[w, g] <- s.wd.6400[w, g]/(-0.3/3.28) + 1; next()
      
    }
    
    if (s.wd.6400[w, g] > 0 && s.wd.6400[w, g] <= 0.5/3.28) {
      
      s.wd.6400[w, g] <- 1; next()
      
    }
    
    if (s.wd.6400[w, g] > 0.5/3.28 && s.wd.6400[w, g] <= 0.8/3.28) {
      
      s.wd.6400[w, g] <- (0.8/3.28 - s.wd.6400[w, g])/(0.3/3.28); next()
      
    }
    
    if (s.wd.6400[w, g] > 0.8/3.28) {
      
      s.wd.6400[w, g] <- 0; next()
      
    }
    
  }
  
}

## Determine HSI_recession for 6400 based on formula
for (w in 1:nrow(s.dh.6400)) {
  
  for (g in 1:ncol(s.dh.6400)) {
    
    if (s.dh.6400[w, g] <= -0.6/3.28) {
      
      s.dh.6400[w, g] <- 0; next()
      
    }
    
    if (s.dh.6400[w, g] > -0.6/3.28 && s.dh.6400[w, g] <= -0.16/3.28) {
      
      s.dh.6400[w, g] <- (s.dh.6400[w, g] + 0.6/3.28)/(0.44/3.28); next()
      
    }
    
    if (s.dh.6400[w, g] > -0.16/3.28 && s.dh.6400[w, g] <= -0.05/3.28) {
      
      s.dh.6400[w, g] <- 1; next()
      
    }
    
    if (s.dh.6400[w, g] > -0.05/3.28 && s.dh.6400[w, g] <= 0.05/3.28) {
      
      s.dh.6400[w, g] <- (0.5 - 10*s.dh.6400[w, g]*3.28); next()
      
    }
    
    if (s.dh.6400[w, g] > 0.05/3.28) {
      
      s.dh.6400[w, g] <- 0; next()
      
    }
    
  }
  
}

## Wading bird suitability is the minimum of the wd and dh for that grid for that week
s.wd.6400 <- s.wd.6400[-c(171,154,137,120,103,86,69,52,35,18,1),]

for (w in 1:nrow(s.dh.6400)) {
  
  for (g in 1:ncol(s.dh.6400)) {
    
    s.wb.6400[w, g] <- min(s.wd.6400[w, g], s.dh.6400[w, g])
    
  }
  
}

# write.table(s.wb.6400, file = "Wading_Birds_wklygrd-6400.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = c(1:ncol(s.wb.6400)))

## Landscape suitability is based on the top 23% of grids each week
for (w in 1:nrow(s.dh.6400)) {
  
  wb.ls.6400[w, 1] <- mean(head(sort(s.wb.6400[w,], decreasing = T),5))
  
}
  
# write.table(wb.ls.6400, file = "Wading_Birds_wklyldscp-6400.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = F)
  
## Wood Stork yearly suitability is based on SI_land for Jan-Mar
for (y in 1:11) {
  
  t.ws.6400[y, 1] <- mean(wb.ls.6400[((y - 1)*16 + 1):((y - 1)*16 + 12), 1])
  
}
  
write.table(t.ws.6400, file = "Wood_Stork-6400.csv", sep = ",", row.names = c(2005:2015),
            col.names = F)

## Ibis and small heron yearly suitability is based on SI_land for Mar-Apr
for (y in 1:11) {
  
  t.ih.6400[y, 1] <- max(0, (1 -  length(which(wb.ls.6400[((y - 1)*16 + 8):(y*16), 1] <= 0.5))/6))
  
}

write.table(t.ih.6400, file = "Ibis_Heron-6400.csv", sep = ",", row.names = c(2005:2015),
            col.names = F)

tock <- Sys.time()
runtime <- tock - tick
print(runtime)
print(Sys.time())
tick <- Sys.time()

## Begin week-by-week analysis of Wading Birds at 3200m
for (y in 1:11) {
  
  for (w in 1:17) { # HSIs for wading birds go from Jan-Apr, approximately 17 weeks
    
    for (g in 1:ncol(s.3200)) {
      
      s.wd.3200[(y-1)*17 + w, g] <- mean(d.3200[((y-1)*365 + (w - 1)*7 + 1):((y-1)*365 + w*7), g]) # find average weekly water depth
      
      if (w > 1) {
        
        s.dh.3200[(y-1)*16 + w - 1, g] <- s.wd.3200[(y-1)*17 + w, g] - s.wd.3200[(y-1)*17 + w - 1, g] # find change in average water depth (since this requires a previous week, I'm only doing it for weeks 2-17)
        
      }
      
    }
    
  }
  
}

## Determine HSI_depth for 3200 based on formula
for (w in 1:nrow(s.wd.3200)) {
  
  for (g in 1:ncol(s.wd.3200)) {
    
    if (s.wd.3200[w, g] <= -0.3/3.28) {
      
      s.wd.3200[w, g] <- 0; next()
      
    }
    
    if (s.wd.3200[w, g] > -0.3/3.28 && s.wd.3200[w, g] <= 0) {
      
      s.wd.3200[w, g] <- s.wd.3200[w, g]/(-0.3/3.28) + 1; next()
      
    }
    
    if (s.wd.3200[w, g] > 0 && s.wd.3200[w, g] <= 0.5/3.28) {
      
      s.wd.3200[w, g] <- 1; next()
      
    }
    
    if (s.wd.3200[w, g] > 0.5/3.28 && s.wd.3200[w, g] <= 0.8/3.28) {
      
      s.wd.3200[w, g] <- (0.8/3.28 - s.wd.3200[w, g])/(0.3/3.28); next()
      
    }
    
    if (s.wd.3200[w, g] > 0.8/3.28) {
      
      s.wd.3200[w, g] <- 0; next()
      
    }
    
  }
  
}

## Determine HSI_recession for 3200 based on formula
for (w in 1:nrow(s.dh.3200)) {
  
  for (g in 1:ncol(s.dh.3200)) {
    
    if (s.dh.3200[w, g] <= -0.6/3.28) {
      
      s.dh.3200[w, g] <- 0; next()
      
    }
    
    if (s.dh.3200[w, g] > -0.6/3.28 && s.dh.3200[w, g] <= -0.16/3.28) {
      
      s.dh.3200[w, g] <- (s.dh.3200[w, g] + 0.6/3.28)/(0.44/3.28); next()
      
    }
    
    if (s.dh.3200[w, g] > -0.16/3.28 && s.dh.3200[w, g] <= -0.05/3.28) {
      
      s.dh.3200[w, g] <- 1; next()
      
    }
    
    if (s.dh.3200[w, g] > -0.05/3.28 && s.dh.3200[w, g] <= 0.05/3.28) {
      
      s.dh.3200[w, g] <- (0.5 - 10*s.dh.3200[w, g]*3.28); next()
      
    }
    
    if (s.dh.3200[w, g] > 0.05/3.28) {
      
      s.dh.3200[w, g] <- 0; next()
      
    }
    
  }
  
}

## Wading bird suitability is the minimum of the wd and dh for that grid for that week
s.wd.3200 <- s.wd.3200[-c(171,154,137,120,103,86,69,52,35,18,1),]

for (w in 1:nrow(s.dh.3200)) {
  
  for (g in 1:ncol(s.dh.3200)) {
    
    s.wb.3200[w, g] <- min(s.wd.3200[w, g], s.dh.3200[w, g])
    
  }
  
}

# write.table(s.wb.3200, file = "Wading_Birds_wklygrd-3200.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = c(1:ncol(s.wb.3200)))

## Landscape suitability is based on the top 23% of grids each week
for (w in 1:nrow(s.dh.3200)) {
  
  wb.ls.3200[w, 1] <- mean(head(sort(s.wb.3200[w,], decreasing = T),15))
  
}

# write.table(wb.ls.3200, file = "Wading_Birds_wklyldscp-3200.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = F)

## Wood Stork yearly suitability is based on SI_land for Jan-Mar
for (y in 1:11) {
  
  t.ws.3200[y, 1] <- mean(wb.ls.3200[((y - 1)*16 + 1):((y - 1)*16 + 12), 1])
  
}

write.table(t.ws.3200, file = "Wood_Stork-3200.csv", sep = ",", row.names = c(1:11),
            col.names = F)

## Ibis and small heron yearly suitability is based on SI_land for Mar-Apr
for (y in 1:11) {
  
  t.ih.3200[y, 1] <- max(0, (1 -  length(which(wb.ls.3200[((y - 1)*16 + 8):(y*16), 1] <= 0.5))/6))
  
}

write.table(t.ih.3200, file = "Ibis_Heron-3200.csv", sep = ",", row.names = c(1:11),
            col.names = F)

tock <- Sys.time()
runtime <- tock - tick
print(runtime)
print(Sys.time())
tick <- Sys.time()

## Begin week-by-week analysis of Wading Birds at 1600m
for (y in 1:11) {
  
  for (w in 1:17) { # HSIs for wading birds go from Jan-Apr, approximately 17 weeks
    
    for (g in 1:ncol(s.1600)) {
      
      s.wd.1600[(y-1)*17 + w, g] <- mean(d.1600[((y-1)*365 + (w - 1)*7 + 1):((y-1)*365 + w*7), g]) # find average weekly water depth
      
      if (w > 1) {
        
        s.dh.1600[(y-1)*16 + w - 1, g] <- s.wd.1600[(y-1)*17 + w, g] - s.wd.1600[(y-1)*17 + w - 1, g] # find change in average water depth (since this requires a previous week, I'm only doing it for weeks 2-17)
        
      }
      
    }
    
  }
  
}

## Determine HSI_depth for 1600 based on formula
for (w in 1:nrow(s.wd.1600)) {
  
  for (g in 1:ncol(s.wd.1600)) {
    
    if (s.wd.1600[w, g] <= -0.3/3.28) {
      
      s.wd.1600[w, g] <- 0; next()
      
    }
    
    if (s.wd.1600[w, g] > -0.3/3.28 && s.wd.1600[w, g] <= 0) {
      
      s.wd.1600[w, g] <- s.wd.1600[w, g]/(-0.3/3.28) + 1; next()
      
    }
    
    if (s.wd.1600[w, g] > 0 && s.wd.1600[w, g] <= 0.5/3.28) {
      
      s.wd.1600[w, g] <- 1; next()
      
    }
    
    if (s.wd.1600[w, g] > 0.5/3.28 && s.wd.1600[w, g] <= 0.8/3.28) {
      
      s.wd.1600[w, g] <- (0.8/3.28 - s.wd.1600[w, g])/(0.3/3.28); next()
      
    }
    
    if (s.wd.1600[w, g] > 0.8/3.28) {
      
      s.wd.1600[w, g] <- 0; next()
      
    }
    
  }
  
}

## Determine HSI_recession for 1600 based on formula
for (w in 1:nrow(s.dh.1600)) {
  
  for (g in 1:ncol(s.dh.1600)) {
    
    if (s.dh.1600[w, g] <= -0.6/3.28) {
      
      s.dh.1600[w, g] <- 0; next()
      
    }
    
    if (s.dh.1600[w, g] > -0.6/3.28 && s.dh.1600[w, g] <= -0.16/3.28) {
      
      s.dh.1600[w, g] <- (s.dh.1600[w, g] + 0.6/3.28)/(0.44/3.28); next()
      
    }
    
    if (s.dh.1600[w, g] > -0.16/3.28 && s.dh.1600[w, g] <= -0.05/3.28) {
      
      s.dh.1600[w, g] <- 1; next()
      
    }
    
    if (s.dh.1600[w, g] > -0.05/3.28 && s.dh.1600[w, g] <= 0.05/3.28) {
      
      s.dh.1600[w, g] <- (0.5 - 10*s.dh.1600[w, g]*3.28); next()
      
    }
    
    if (s.dh.1600[w, g] > 0.05/3.28) {
      
      s.dh.1600[w, g] <- 0; next()
      
    }
    
  }
  
}

## Wading bird suitability is the minimum of the wd and dh for that grid for that week
s.wd.1600 <- s.wd.1600[-c(171,154,137,120,103,86,69,52,35,18,1),]

for (w in 1:nrow(s.dh.1600)) {
  
  for (g in 1:ncol(s.dh.1600)) {
    
    s.wb.1600[w, g] <- min(s.wd.1600[w, g], s.dh.1600[w, g])
    
  }
  
}

# write.table(s.wb.1600, file = "Wading_Birds_wklygrd-1600.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = c(1:ncol(s.wb.1600)))

## Landscape suitability is based on the top 23% of grids each week
for (w in 1:nrow(s.dh.1600)) {
  
  wb.ls.1600[w, 1] <- mean(head(sort(s.wb.1600[w,], decreasing = T),55))
  
}

# write.table(wb.ls.1600, file = "Wading_Birds_wklyldscp-1600.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = F)

## Wood Stork yearly suitability is based on SI_land for Jan-Mar
for (y in 1:11) {
  
  t.ws.1600[y, 1] <- mean(wb.ls.1600[((y - 1)*16 + 1):((y - 1)*16 + 12), 1])
  
}

write.table(t.ws.1600, file = "Wood_Stork-1600.csv", sep = ",", row.names = c(1:11),
            col.names = F)

## Ibis and small heron yearly suitability is based on SI_land for Mar-Apr
for (y in 1:11) {
  
  t.ih.1600[y, 1] <- max(0, (1 -  length(which(wb.ls.1600[((y - 1)*16 + 8):(y*16), 1] <= 0.5))/6))
  
}

write.table(t.ih.1600, file = "Ibis_Heron-1600.csv", sep = ",", row.names = c(1:11),
            col.names = F)

tock <- Sys.time()
runtime <- tock - tick
print(runtime)
print(Sys.time())
tick <- Sys.time()

## Begin week-by-week analysis of Wading Birds at 800m
for (y in 1:11) {
  
  for (w in 1:17) { # HSIs for wading birds go from Jan-Apr, approximately 17 weeks
    
    for (g in 1:ncol(s.800)) {
      
      s.wd.800[(y-1)*17 + w, g] <- mean(d.800[((y-1)*365 + (w - 1)*7 + 1):((y-1)*365 + w*7), g]) # find average weekly water depth
      
      if (w > 1) {
        
        s.dh.800[(y-1)*16 + w - 1, g] <- s.wd.800[(y-1)*17 + w, g] - s.wd.800[(y-1)*17 + w - 1, g] # find change in average water depth (since this requires a previous week, I'm only doing it for weeks 2-17)
        
      }
      
    }
    
  }
  
}

## Determine HSI_depth for 800 based on formula
for (w in 1:nrow(s.wd.800)) {
  
  for (g in 1:ncol(s.wd.800)) {
    
    if (s.wd.800[w, g] <= -0.3/3.28) {
      
      s.wd.800[w, g] <- 0; next()
      
    }
    
    if (s.wd.800[w, g] > -0.3/3.28 && s.wd.800[w, g] <= 0) {
      
      s.wd.800[w, g] <- s.wd.800[w, g]/(-0.3/3.28) + 1; next()
      
    }
    
    if (s.wd.800[w, g] > 0 && s.wd.800[w, g] <= 0.5/3.28) {
      
      s.wd.800[w, g] <- 1; next()
      
    }
    
    if (s.wd.800[w, g] > 0.5/3.28 && s.wd.800[w, g] <= 0.8/3.28) {
      
      s.wd.800[w, g] <- (0.8/3.28 - s.wd.800[w, g])/(0.3/3.28); next()
      
    }
    
    if (s.wd.800[w, g] > 0.8/3.28) {
      
      s.wd.800[w, g] <- 0; next()
      
    }
    
  }
  
}

## Determine HSI_recession for 800 based on formula
for (w in 1:nrow(s.dh.800)) {
  
  for (g in 1:ncol(s.dh.800)) {
    
    if (s.dh.800[w, g] <= -0.6/3.28) {
      
      s.dh.800[w, g] <- 0; next()
      
    }
    
    if (s.dh.800[w, g] > -0.6/3.28 && s.dh.800[w, g] <= -0.16/3.28) {
      
      s.dh.800[w, g] <- (s.dh.800[w, g] + 0.6/3.28)/(0.44/3.28); next()
      
    }
    
    if (s.dh.800[w, g] > -0.16/3.28 && s.dh.800[w, g] <= -0.05/3.28) {
      
      s.dh.800[w, g] <- 1; next()
      
    }
    
    if (s.dh.800[w, g] > -0.05/3.28 && s.dh.800[w, g] <= 0.05/3.28) {
      
      s.dh.800[w, g] <- (0.5 - 10*s.dh.800[w, g]*3.28); next()
      
    }
    
    if (s.dh.800[w, g] > 0.05/3.28) {
      
      s.dh.800[w, g] <- 0; next()
      
    }
    
  }
  
}

## Wading bird suitability is the minimum of the wd and dh for that grid for that week
s.wd.800 <- s.wd.800[-c(171,154,137,120,103,86,69,52,35,18,1),]

for (w in 1:nrow(s.dh.800)) {
  
  for (g in 1:ncol(s.dh.800)) {
    
    s.wb.800[w, g] <- min(s.wd.800[w, g], s.dh.800[w, g])
    
  }
  
}

# write.table(s.wb.800, file = "Wading_Birds_wklygrd-800.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = c(1:ncol(s.wb.800)))

## Landscape suitability is based on the top 23% of grids each week
for (w in 1:nrow(s.dh.800)) {
  
  wb.ls.800[w, 1] <- mean(head(sort(s.wb.800[w,], decreasing = T),210))
  
}

# write.table(wb.ls.800, file = "Wading_Birds_wklyldscp-800.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = F)

## Wood Stork yearly suitability is based on SI_land for Jan-Mar
for (y in 1:11) {
  
  t.ws.800[y, 1] <- mean(wb.ls.800[((y - 1)*16 + 1):((y - 1)*16 + 12), 1])
  
}

write.table(t.ws.800, file = "Wood_Stork-800.csv", sep = ",", row.names = c(1:11),
            col.names = F)

## Ibis and small heron yearly suitability is based on SI_land for Mar-Apr
for (y in 1:11) {
  
  t.ih.800[y, 1] <- max(0, (1 -  length(which(wb.ls.800[((y - 1)*16 + 8):(y*16), 1] <= 0.5))/6))
  
}

write.table(t.ih.800, file = "Ibis_Heron-800.csv", sep = ",", row.names = c(1:11),
            col.names = F)

tock <- Sys.time()
runtime <- tock - tick
print(runtime)
print(Sys.time())
tick <- Sys.time()

## Begin week-by-week analysis of Wading Birds at 400m
for (y in 1:11) {
  
  for (w in 1:17) { # HSIs for wading birds go from Jan-Apr, approximately 17 weeks
    
    for (g in 1:ncol(s.400)) {
      
      s.wd.400[(y-1)*17 + w, g] <- mean(d.400[((y-1)*365 + (w - 1)*7 + 1):((y-1)*365 + w*7), g]) # find average weekly water depth
      
      if (w > 1) {
        
        s.dh.400[(y-1)*16 + w - 1, g] <- s.wd.400[(y-1)*17 + w, g] - s.wd.400[(y-1)*17 + w - 1, g] # find change in average water depth (since this requires a previous week, I'm only doing it for weeks 2-17)
        
      }
      
    }
    
  }
  
}

## Determine HSI_depth for 400 based on formula
for (w in 1:nrow(s.wd.400)) {
  
  for (g in 1:ncol(s.wd.400)) {
    
    if (s.wd.400[w, g] <= -0.3/3.28) {
      
      s.wd.400[w, g] <- 0; next()
      
    }
    
    if (s.wd.400[w, g] > -0.3/3.28 && s.wd.400[w, g] <= 0) {
      
      s.wd.400[w, g] <- s.wd.400[w, g]/(-0.3/3.28) + 1; next()
      
    }
    
    if (s.wd.400[w, g] > 0 && s.wd.400[w, g] <= 0.5/3.28) {
      
      s.wd.400[w, g] <- 1; next()
      
    }
    
    if (s.wd.400[w, g] > 0.5/3.28 && s.wd.400[w, g] <= 0.8/3.28) {
      
      s.wd.400[w, g] <- (0.8/3.28 - s.wd.400[w, g])/(0.3/3.28); next()
      
    }
    
    if (s.wd.400[w, g] > 0.8/3.28) {
      
      s.wd.400[w, g] <- 0; next()
      
    }
    
  }
  
}

## Determine HSI_recession for 400 based on formula
for (w in 1:nrow(s.dh.400)) {
  
  for (g in 1:ncol(s.dh.400)) {
    
    if (s.dh.400[w, g] <= -0.6/3.28) {
      
      s.dh.400[w, g] <- 0; next()
      
    }
    
    if (s.dh.400[w, g] > -0.6/3.28 && s.dh.400[w, g] <= -0.16/3.28) {
      
      s.dh.400[w, g] <- (s.dh.400[w, g] + 0.6/3.28)/(0.44/3.28); next()
      
    }
    
    if (s.dh.400[w, g] > -0.16/3.28 && s.dh.400[w, g] <= -0.05/3.28) {
      
      s.dh.400[w, g] <- 1; next()
      
    }
    
    if (s.dh.400[w, g] > -0.05/3.28 && s.dh.400[w, g] <= 0.05/3.28) {
      
      s.dh.400[w, g] <- (0.5 - 10*s.dh.400[w, g]*3.28); next()
      
    }
    
    if (s.dh.400[w, g] > 0.05/3.28) {
      
      s.dh.400[w, g] <- 0; next()
      
    }
    
  }
  
}

## Wading bird suitability is the minimum of the wd and dh for that grid for that week
s.wd.400 <- s.wd.400[-c(171,154,137,120,103,86,69,52,35,18,1),]

for (w in 1:nrow(s.dh.400)) {
  
  for (g in 1:ncol(s.dh.400)) {
    
    s.wb.400[w, g] <- min(s.wd.400[w, g], s.dh.400[w, g])
    
  }
  
}

# write.table(s.wb.400, file = "Wading_Birds_wklygrd-400.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = c(1:ncol(s.wb.400)))

## Landscape suitability is based on the top 23% of grids each week
for (w in 1:nrow(s.dh.400)) {
  
  wb.ls.400[w, 1] <- mean(head(sort(s.wb.400[w,], decreasing = T),805))
  
}

# write.table(wb.ls.400, file = "Wading_Birds_wklyldscp-400.csv", sep = ",", row.names = c(rep(1:16, 11)),
#             col.names = F)

## Wood Stork yearly suitability is based on SI_land for Jan-Mar
for (y in 1:11) {
  
  t.ws.400[y, 1] <- mean(wb.ls.400[((y - 1)*16 + 1):((y - 1)*16 + 12), 1])
  
}

write.table(t.ws.400, file = "Wood_Stork-400.csv", sep = ",", row.names = c(1:11),
            col.names = F)

## Ibis and small heron yearly suitability is based on SI_land for Mar-Apr
for (y in 1:11) {
  
  t.ih.400[y, 1] <- max(0, (1 -  length(which(wb.ls.400[((y - 1)*16 + 8):(y*16), 1] <= 0.5))/6))
  
}

write.table(t.ih.400, file = "Ibis_Heron-400.csv", sep = ",", row.names = c(1:11),
            col.names = F)

tock <- Sys.time()
runtime <- tock - tick
print(runtime)
print(Sys.time())

