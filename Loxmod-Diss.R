###############################################
##############  Ian Hahus, 2018  ##############
####     Water Balance Model for WCA 1     ####
###############################################

# Clear all variables before new model run #
rm(list=ls(all=TRUE))

## Read inputs from file ##
ins <- read.csv("input.csv", header = F)

## Assign parameter values based on model inputs ##
L <- ins[1,1] 
per.day <- ins[2,1]
a <- ins[3,1]
b <- ins[4,1]
GM <- ins[5,1]
ds.max <- ins[6,1]
d.dry <- ins[7,1]
ndsmax <- ins[8,1]

StgB <- 3.5 
d.det <- 0.03

## Read in other data ##
Q <- read.csv("Flows05-.csv")/35.3  # [m3s-1]  
Q[is.na(Q)] <- 0                  # replaces NA with 0   
Q$G301 <- -Q$G301   
P <- read.csv("RainT05-.csv", header = T)*0.0254 # [m]
Tpoly <- read.csv(paste("Thiessen_", L, "m.csv", sep = ""), header = F) # Thiessen polygon assignments
ET <- read.csv("LOXWS_ET05-.csv", header = T)/1000 # convert mm to m    
ET[is.na(ET)] <- 0                  # replaces NA with 0 
grid.data <- read.csv(paste("DATA_MATRIX_", L, "m.csv", sep = ""))   # the attributes of each grid cell
grid.data$Area <- grid.data$Area*((143238*4047)/sum(grid.data$Area)) # adjusting total area to match real value per Lox CCP
gnsadj <- read.csv(paste("Ga_and_Str_adj_", L, "m.csv", sep = ""), header = F) # the adjacency file for stage gauges and flow strucs
f2 <- gnsadj[1:6,1]                 # lists the grid nearest each marsh stage gauge
g2 <- gnsadj[,2]                    # lists the grid nearest each flow structure
stage <- read.csv("Stage05-.csv")/3.28       #measured daily mean stage at each marsh gauge [m]
adj <- read.csv(paste("ADJ_", L, "m.csv", sep = ""), header = F)  #the adjacency matrix created in Grid_Resize.R !!!!!!!!!!!!!!!!!
elevs <- grid.data$elev_ft_NG/3.28           #cell elevations [m] 

## Define constants and parameters ##
days <- nrow(Q)
ncell <- nrow(grid.data) 

## Create an empty stage vector for each cell 
stage.day <- matrix(nrow = (days), ncol = ncell) # a column for each cell, a row for each day
stage.day[1,] <- grid.data$StgJan05/3.28 # set the first stage equal to the initial value (1/1/05) 
stage.day.temp <- grid.data$StgJan05/3.28 # initialize temporary, sub-daily stage vector


ord <- seq(1, ncell, 1) # the order in which cells are to be evaluated

for (d in 1:days) {   
  
  for (subd in 1:per.day) {
    
    d.stage <- mat.or.vec(ncell, 1)
    
    for (c in ord) {
    
      prow <- P[d,] # the values for each rain gauge on this day
      prow2 <- prow[c(Tpoly[c,1], Tpoly[c,2], Tpoly[c,3], Tpoly[c,4], Tpoly[c,5], Tpoly[c,6])] # ordering the values based on how close they are the current cell
      prow3 <- prow2[!is.na(prow2)] # removing NAs from the ordered list (needed if the nearest gauge(s) do not have values for this day)
      P.d <- prow3[1]
      
      d.stage[c] <- (P.d - GM*(stage.day.temp[c] - StgB) - ET2[d,1])/per.day #calculate incremental change in stage
      
      if (c %in% g2) {                 #if current cell is nearest a flow structure(s), give/take its flow
        matches <- which(match(g2, c) == 1) # identifies which structure(s) is nearest (needed is n>1)
        d.stage[c] <- d.stage[c] + sum(Q[d,c(matches)])*86400/grid.data$Area[c]/per.day   #find appropriate flow rate [cmd] based on day and grid cell, convert to meters to find stage change
      }
      
      if (abs(d.stage[c]) <= ds.max) { # if ds.max is not exceeded in the cell nearest a flow structure
        
        stage.day.temp[c] <- stage.day.temp[c] + d.stage[c] #update stage for current cell in current day
        
        if ((elevs[c] - d.dry) > stage.day.temp[c]) { # if current stage is below allowable elevation
          
          fill <- (elevs[c] - d.dry) - stage.day.temp[c]
          stage.day.temp[c] <- stage.day.temp[c] + fill
          hole <- -fill*grid.data$Area[c]
          
          nbr <- 8 - sum(as.integer(is.na(adj[,c]))) # find out how many neighbors the current cell has
          
          nbs <- matrix(nr = 1, nc = nbr) # create vector to hold neighboring grid numbers
          j <- 1
          
          for (i in 1:8) {  # store the adjacent cells in a vector w/o NAs
            if (is.na(adj[i,c]) == F) {
              nbs[j] <- adj[i,c]
              j <- j + 1
            }
          }
          
          for (i in 1:length(nbs)) { # spread the excess water to neighboring cells
            
            if (hole < 0) { # it will be intitially, but might "run out" before all neighbors are subtracted from
              
              d.stage.nbs <- min(max(ndsmax, (stage.day.temp[c] - stage.day.temp[nbs[i]]),
                                     -(stage.day.temp[nbs[i]] - (elevs[nbs[i]] - d.dry)), hole/grid.data$Area[nbs[i]]), 0)
              stage.day.temp[nbs[i]] <- stage.day.temp[nbs[i]] + d.stage.nbs
              hole <- hole - d.stage.nbs*grid.data$Area[nbs[i]]
              
            }
            
          }
          
          if (hole < 0) { # if the neighbors of the current cell are not enough to spread out excess
            
            front <- matrix(nr = 1, nc = 8*length(nbs)) # initialize vector for neighbors of neighbors
            nbs.old <- nbs # save old values to check for overlap
            
            for (i in 1:length(nbs)) { # populate neighbor neighbor vector, including NAs
              front[(8*(i-1)+1):(8*i)] <- adj[,nbs[i]]
            }
            
            nbr <- length(front) - sum(as.integer(is.na(front))) # find out how many non-NA neighbors there are
            k <- 1
            nbs <- mat.or.vec(1,nbr) # initialize vector of unique neighboring grids
            
            for (i in 1:length(front)) {  # store the adjacent cells in a vector w/o NAs
              if (is.na(front[i]) == F) {
                nbs[k] <- front[i]
                k <- k + 1
              }
            }
            
            nbs <- unique(as.vector(nbs)) # remove duplicate entries
            nbs.old[(length(nbs.old)+1)] <- c # add the original cell to list of potential overlaps
            
            for (i in 1:length(nbs.old)) {
              
              if (nbs.old[i] %in% nbs == T) { # remove the original cell if included as a neighbor
                nbs <- nbs[-match(nbs.old[i],nbs)]
              }
              
            }
            
            for (i in 1:length(nbs)) { # allocate extra flow to neighbor neighbors up to ds.max
              
              if (hole < 0) { # it will be intitially, but might "run out" before all neighbors are added to
                
                d.stage.nbs <- min(max(ndsmax, (stage.day.temp[c] - stage.day.temp[nbs[i]]),
                                       -(stage.day.temp[nbs[i]] - (elevs[nbs[i]] - d.dry)), hole/grid.data$Area[nbs[i]]), 0)
                stage.day.temp[nbs[i]] <- stage.day.temp[nbs[i]] + d.stage.nbs
                hole <- hole - d.stage.nbs*grid.data$Area[nbs[i]]
                
              }
              
              if ((elevs[nbs[i]] - d.dry) > stage.day.temp[nbs[i]]) {
                
                fill <- (elevs[nbs[i]] - d.dry) - stage.day.temp[nbs[i]]
                stage.day.temp[nbs[i]] <- stage.day.temp[nbs[i]] + fill
                hole <- hole - fill*grid.data$Area[nbs[i]]
                
              }
              
            } # for each neighbor neighboring cell
            
          } # if the neighbors of the current cell are not enough to spread out excess
          
        } # if current stage is below allowable elevation
        
      } # if ds.max is not exceeded in the cell nearest a flow structure
      ##############################################################################################################################
      ######################################################################################
      if (d.stage[c] > ds.max) {   #0.2m set as max daily stage change for cell
        
        sgn <- d.stage[c]/abs(d.stage[c])
        d..stage <- sgn*ds.max # scaling d.stage down to ds.max while retaining sign
        stage.day.temp[c] <- stage.day.temp[c] + d..stage   # updating stage for the current cell
        
        pile <- (d.stage[c] - ds.max)*grid.data$Area[c] # the amount of excess from cell in question, used in logical below 
        
        nbr <- 8 - sum(as.integer(is.na(adj[,c]))) # find out how many neighbors the current cell has
        
        nbs <- matrix(nr = 1, nc = nbr) # create vector to hold neighboring grid numbers
        j <- 1
        
        for (i in 1:8) {  # store the adjacent cells in a vector w/o NAs
          if (is.na(adj[i,c]) == F) {
            nbs[j] <- adj[i,c]
            j <- j + 1
          }
        }
        
        for (i in 1:length(nbs)) { # spread the excess water to neighboring cells
          
          if (pile > 0) { # it will be intitially, but might "run out" before all neighbors are added to
            
            d.stage.nbs <- max(min(ds.max, (stage.day.temp[c] - stage.day.temp[nbs[i]]), pile/grid.data$Area[nbs[i]]), 0)
            stage.day.temp[nbs[i]] <- stage.day.temp[nbs[i]] + d.stage.nbs
            pile <- pile - d.stage.nbs*grid.data$Area[nbs[i]]
            
          } 
          
        }
        
        if (pile > 0) { # if the neighbors of the current cell are not enough to spread out excess   
          
          front <- matrix(nr = 1, nc = 8*length(nbs)) # initialize vector for neighbors of neighbors
          nbs.old <- nbs # save old values to check for overlap
          
          for (i in 1:length(nbs)) { # populate neighbor neighbor vector, including NAs
            front[(8*(i-1)+1):(8*i)] <- adj[,nbs[i]]
          }
          
          nbr <- length(front) - sum(as.integer(is.na(front))) # find out how many non-NA neighbors there are
          k <- 1
          nbs <- mat.or.vec(1,nbr) # initialize vector of unique neighboring grids
          
          for (i in 1:length(front)) {  # store the adjacent cells in a vector w/o NAs
            if (is.na(front[i]) == F) {
              nbs[k] <- front[i]
              k <- k + 1
            }
          }
          
          nbs <- unique(as.vector(nbs)) # remove duplicate entries
          nbs.old[(length(nbs.old)+1)] <- c # add the original cell to list of potential overlaps
          
          for (i in 1:length(nbs.old)) {
            
            if (nbs.old[i] %in% nbs == T) { # remove the original cell if included as a neighbor
              nbs <- nbs[-match(nbs.old[i],nbs)] 
            }
            
          }
          
          for (i in 1:length(nbs)) { # allocate extra flow to neighbor neighbors up to ds.max
            
            if (pile > 0) { # it will be intitially, but might "run out" before all neighbors are added to
              
              d.stage.nbs <- max(min(ds.max, (stage.day.temp[c] - stage.day.temp[nbs[i]]), pile/grid.data$Area[nbs[i]]), 0)
              stage.day.temp[nbs[i]] <- stage.day.temp[nbs[i]] + d.stage.nbs
              pile <- pile - d.stage.nbs*grid.data$Area[nbs[i]]
              
            } 
            
          } # for each neighbor neighboring cell
          
        } # if the first round of adjacent cells isn't enough to spread it out
        
      } # if d.stage is greater than ds.max  
      
      ########################################################################################################
      ########################################################################################################
      ########################################################################################################
      
      if (d.stage[c] < -ds.max) { 
        
        
        sgn <- d.stage[c]/abs(d.stage[c])
        d..stage <- sgn*ds.max # scaling d.stage down to ds.max while retaining sign
        stage.day.temp[c] <- stage.day.temp[c] + max(d..stage, -(stage.day.temp[c] - (elevs[c] - d.dry)))    # updating stage for the current cell
        
        hole <- (d.stage[c] + min(ds.max, (stage.day.temp[c] - (elevs[c] - d.dry))))*grid.data$Area[c] # the amount of debt from cell in question, used in logical below 
        
        if ((elevs[c] - d.dry) > stage.day.temp[c]) {
          
          fill <- (elevs[c] - d.dry) - stage.day.temp[c]
          stage.day.temp[c] <- stage.day.temp[c] + fill
          hole <- hole - fill*grid.data$Area[c]
          
        }
        
        nbr <- 8 - sum(as.integer(is.na(adj[,c]))) # find out how many neighbors the current cell has
        
        nbs <- matrix(nr = 1, nc = nbr) # create vector to hold neighboring grid numbers
        j <- 1
        
        for (i in 1:8) {  # store the adjacent cells in a vector w/o NAs
          if (is.na(adj[i,c]) == F) {
            nbs[j] <- adj[i,c]
            j <- j + 1
          }
        }
        
        for (i in 1:length(nbs)) { # spread the excess water to neighboring cells
          
          if (hole < 0) { # it will be intitially, but might "run out" before all neighbors are subtracted from
            
            d.stage.nbs <- min(max(ndsmax, (stage.day.temp[c] - stage.day.temp[nbs[i]]), 
                                   -(stage.day.temp[nbs[i]] - (elevs[nbs[i]] - d.dry)), hole/grid.data$Area[nbs[i]]), 0)
            stage.day.temp[nbs[i]] <- stage.day.temp[nbs[i]] + d.stage.nbs
            hole <- hole - d.stage.nbs*grid.data$Area[nbs[i]]
            
          }
          
          if ((elevs[nbs[i]] - d.dry) > stage.day.temp[nbs[i]]) {
            
            fill <- (elevs[nbs[i]] - d.dry) - stage.day.temp[nbs[i]]
            stage.day.temp[nbs[i]] <- stage.day.temp[nbs[i]] + fill
            hole <- hole - fill*grid.data$Area[nbs[i]]
            
          }
          
        } 
        
        if (hole < 0) { # if the neighbors of the current cell are not enough to spread out excess   
          
          front <- matrix(nr = 1, nc = 8*length(nbs)) # initialize vector for neighbors of neighbors
          nbs.old <- nbs # save old values to check for overlap
          
          for (i in 1:length(nbs)) { # populate neighbor neighbor vector, including NAs
            front[(8*(i-1)+1):(8*i)] <- adj[,nbs[i]]
          }
          
          nbr <- length(front) - sum(as.integer(is.na(front))) # find out how many non-NA neighbors there are
          k <- 1
          nbs <- mat.or.vec(1,nbr) # initialize vector of unique neighboring grids
          
          for (i in 1:length(front)) {  # store the adjacent cells in a vector w/o NAs
            if (is.na(front[i]) == F) {
              nbs[k] <- front[i]
              k <- k + 1
            }
          }
          
          nbs <- unique(as.vector(nbs)) # remove duplicate entries
          nbs.old[(length(nbs.old)+1)] <- c # add the original cell to list of potential overlaps
          
          for (i in 1:length(nbs.old)) {
            
            if (nbs.old[i] %in% nbs == T) { # remove the original cell if included as a neighbor
              nbs <- nbs[-match(nbs.old[i],nbs)] 
            }
            
          }
          
          for (i in 1:length(nbs)) { # allocate extra flow to neighbor neighbors up to ds.max
            
            if (hole < 0) { # it will be intitially, but might "run out" before all neighbors are added to
              
              d.stage.nbs <- min(max(ndsmax, (stage.day.temp[c] - stage.day.temp[nbs[i]]), 
                                     -(stage.day.temp[nbs[i]] - (elevs[nbs[i]] - d.dry)), hole/grid.data$Area[nbs[i]]), 0)
              stage.day.temp[nbs[i]] <- stage.day.temp[nbs[i]] + d.stage.nbs
              hole <- hole - d.stage.nbs*grid.data$Area[nbs[i]]
              
            }
            
            if ((elevs[nbs[i]] - d.dry) > stage.day.temp[nbs[i]]) {
              
              fill <- (elevs[nbs[i]] - d.dry) - stage.day.temp[nbs[i]]
              stage.day.temp[nbs[i]] <- stage.day.temp[nbs[i]] + fill
              hole <- hole - fill*grid.data$Area[nbs[i]]
              
            }
            
          } # for each neighbor neighboring cell
          
        } # if the first round of adjacent cells isn't enough to spread it out
        
      } # if d.stage is less than negative ds.max
      ########################################################################################################
      ########################################################################################################
      ########################################################################################################  
      
      zone <- c(adj[,c], c) # define zone as current grid and all adjacent grids
      zone.stage <- stage.day.temp[zone]
      c.new <- zone[which.max(zone.stage)] # find which grid has highest stage and make this the "exporting" grid
      adj.new <- zone[zone!=c.new] # designate new adjacency matrix for this grid
      
      if (stage.day.temp[c.new] - elevs[c.new] > d.det) { # if water depth is above detention (Min et al.)
        
        ## Find water surface elevation difference between adjacent cells ##
        dh <- mat.or.vec(nr = 8, nc = 1)
        for (m in 1:8) {
          dh[m] <- stage.day.temp[c.new] - stage.day.temp[(adj.new[m])] #find drop to stage or elev, whichever is higher
        }
        dh[is.na(dh)] <- -9999       #re-assigning NAs to work with logical statements
        
        ## Find relative slope between adjacent cells ##
        if (max(dh, na.rm = T) > 0.005) {        #If there is positive slope from current cell, distribute water from cell
          
          dh[dh<0.005] <- NA
          dh[dh==0] <- NA
          
          slope <- mat.or.vec(nr = 8, nc = 1)
          
          for (n in 1:length(dh)) {        #for adjacent "corners", the distance is greater
            slope[n] <- dh[n]/sqrt((grid.data$Yutm[c.new] - grid.data$Yutm[adj.new[n]])^2 + 
                                     (grid.data$Xutm[c.new] - grid.data$Xutm[adj.new[n]])^2) 
          }
          for (s in 1:length(slope)) {   #re-re-assign NAs to find proportion of slope (and flow) to allot to each downhill neighbor
            if (is.na(slope[s]) == T) {
              slope[s] <- 0  
            }
            if (slope[s] < 0) {
              slope[s] <- 0  #maybe change this...have to see how I do "pre-calculation"...
            }
          }
          prop <- matrix(nrow = 8, ncol = 1)
          for (p in 1:length(slope)) {   #find the portion of "total slope" for each downhill cell
            prop[p] <- slope[p]/sum(slope)
          }
          
          #My idea is to determine total outflow based on the largest head difference, then allocate
          #that flow proportionally based on the head gradients with the downhill neighbors
          Q.coeff <- cbind(prop,adj.new)     #link proportion of flow to corresponding neighbors
          Q.tot <- min(grid.data$Area[c.new]*(ds.max), grid.data$Area[c.new]*(stage.day.temp[c.new] - elevs[c.new]),  
                       ((stage.day.temp[c.new] - elevs[c.new])^(1.67 - b)*max(slope)^(0.5))*L*86400/per.day/a) #d-dependent Mannings
          
          Q.nbr <- length(which(Q.coeff[,1] > 0))  # number of downhill neighbors
          Q.ind <- which(Q.coeff[,1] > 0)                   # indices of the downhill neighbors
          
          if (Q.nbr == 1) { # if there is only one downhill neighbor, give either predicted flow or the flow that would equilibrate their stages, whichever is least
            
            Q.tot.1 <- min(Q.tot, grid.data$Area[Q.coeff[Q.ind,2]]*(ds.max), 
                           ((grid.data$Area[c.new]*grid.data$Area[Q.coeff[Q.ind,2]]*dh[Q.ind])/(grid.data$Area[c.new]+grid.data$Area[Q.coeff[Q.ind,2]]))) # determine which of the possible flow( volume)s is the least
            stage.day.temp[c.new] <- stage.day.temp[c.new] - Q.tot.1/grid.data$Area[c.new] # decrease original cell stage by appropraite amount
            stage.day.temp[Q.coeff[Q.ind,2]] <- stage.day.temp[Q.coeff[Q.ind,2]] - Q.tot/grid.data$Area[Q.coeff[Q.ind,2]] # increase destination cell stage by appropraite amount
            
          } else { # if there is more than one downhill neighbor to receive leftover flow
            
            Q.vol <- matrix(nrow = Q.nbr, ncol = 4)           # create matrix to store data
            Q.vol[,1] <- Q.ind                                # populate first col with cell numbers
            Q.vol[,3] <- c(dh[c(Q.ind)])
            
            for (o in 1:Q.nbr) {
              Q.vol[o,2] <- Q.coeff[Q.ind[o],2] # grid number corresponding to each downhill neighbor
              Q.vol[o,4] <- (grid.data$Area[Q.vol[o,2]]*grid.data$Area[c.new]*Q.vol[o,3])/
                (grid.data$Area[Q.vol[o,2]] + grid.data$Area[c.new]*Q.coeff[Q.vol[o,1],1]) # total flow volume required to equilibrate given neighbor with original cell
            }
            
            Q.eq <- min(Q.vol[,4])
            Q.rat <- Q.tot/Q.eq
            
            if (Q.rat <= 1) { # if calculated flow is less than least equlibrating volume
              
              stage.day.temp[c.new] <- stage.day.temp[c.new] - Q.tot/grid.data$Area[c.new] 
              
              for (q in Q.ind) {
                
                stage.day.temp[Q.coeff[q,2]] <- stage.day.temp[Q.coeff[q,2]] + Q.tot*Q.coeff[q,1]/grid.data$Area[Q.coeff[q,2]]
                
              } # if calculated flow is less than least equlibrating volume
              
            } 
            
            if (Q.rat > 1) { # if calculated flow is greater than the least equilibrating volume
              
              stage.day.temp[c.new] <- stage.day.temp[c.new] - Q.eq/grid.data$Area[c.new] 
              
              for (q in Q.ind) {
                
                stage.day.temp[Q.coeff[q,2]] <- stage.day.temp[Q.coeff[q,2]] + Q.eq*Q.coeff[q,1]/grid.data$Area[Q.coeff[q,2]] # give each downhill neighbor its proportion of reduced total flow
                
              }
              
              Q.2 <- Q.eq
              lev <- matrix(Q.vol[which.min(Q.vol[,4]),], ncol = 4, byrow = T)
              Q.vol.2 <- matrix(Q.vol[-which.min(Q.vol[,4]),], ncol = 4, byrow = F)
              new.area <- grid.data$Area[lev[1,2]] + grid.data$Area[c.new]
              
              if (Q.2 < Q.tot && nrow(Q.vol.2) > 0) {
                
                dh.2 <- mat.or.vec(nr = nrow(Q.vol.2), nc = 1)
                for (m in 1:nrow(Q.vol.2)) {
                  dh.2[m] <- stage.day.temp[c.new] - max(stage.day.temp[Q.vol.2[m,2]], elevs[Q.vol.2[m,2]]) #find drop to stage or elev, whichever is higher
                }
                
                ## Find relative slope between adjacent cells ##
                
                slope.2 <- mat.or.vec(nr = nrow(Q.vol.2), nc = 1)
                
                for (n in 1:nrow(Q.vol.2)) {        
                  
                  if (Q.vol.2[n,1] %in% c(1,3,5,7)) { # for adjacent "corners", the distance is greater
                    slope.2[n] <- dh.2[n]/sqrt(2)
                  } else{
                    slope.2[n] <- dh.2[n]
                  }
                }
                
                prop.2 <- matrix(nrow = nrow(Q.vol.2), ncol = 1)
                for (p in 1:length(slope.2)) {   #find the portion of "total slope" for each downhill cell
                  prop.2[p] <- slope.2[p]/sum(slope.2)
                  if (sum(slope.2) == 0) {
                    prop.2[p] <- 0 # in case all grids are at same level as grid [c]
                  }
                }
                
                Q.coeff.2 <- cbind(prop.2, Q.vol.2[,2])
                
                for (r in 1:nrow(Q.vol.2)) {
                  
                  Q.vol.2[r,3] <- dh.2[r]
                  
                  Q.vol.2[r,4] <- (grid.data$Area[Q.vol.2[r,2]]*(new.area)*Q.vol.2[r,3])/
                    (grid.data$Area[Q.vol.2[r,2]] + (new.area)*Q.coeff.2[r]) # total flow volume required to equilibrate given neighbor with original cell
                  
                }
                
                Q.eq.2 <- min(Q.vol.2[,4])
                Q.rat.2 <- (Q.tot - Q.2)/Q.eq.2
                
                if (Q.eq.2 == 0) { # if equilibrating volume is 0, skip allocation steps
                  Q.rat.2 <- 0
                }
                
                if (Q.rat.2 <= 1 && Q.eq.2 > 0) { # if calculated flow is less than least equlibrating volume
                  
                  stage.day.temp[c.new] <- stage.day.temp[c.new] - (Q.tot - Q.2)/new.area 
                  stage.day.temp[c(lev[,2])] <- stage.day.temp[c(lev[,2])] - (Q.tot - Q.2)/new.area
                  
                  for (q in 1:nrow(Q.vol.2)) {
                    
                    stage.day.temp[Q.vol.2[q,2]] <- stage.day.temp[Q.vol.2[q,2]] + (Q.tot - Q.2)*prop.2[q]/grid.data$Area[Q.vol.2[q,2]]
                    
                  } # if calculated flow is less than least equlibrating volume
                  
                } 
                
                if (Q.rat.2 > 1) { # if calculated flow is greater than the least equilibrating volume
                  
                  stage.day.temp[c.new] <- stage.day.temp[c.new] - Q.eq.2/new.area
                  stage.day.temp[c(lev[,2])] <- stage.day.temp[c(lev[,2])] - (Q.tot - Q.2)/new.area
                  
                  for (q in 1:nrow(Q.vol.2)) {
                    
                    stage.day.temp[Q.vol.2[q,2]] <- stage.day.temp[Q.vol.2[q,2]] + Q.eq.2*prop.2[q]/grid.data$Area[Q.vol.2[q,2]] # give each downhill neighbor its proportion of reduced total flow
                    
                  }
                  
                  Q.2 <- Q.2 + Q.eq.2
                  
                }
                
                lev <- matrix(c(lev, Q.vol.2[which.min(Q.vol.2[,4]),]), ncol = 4, byrow = T)
                Q.vol.2 <- matrix(Q.vol.2[-which.min(Q.vol.2[,4]),], ncol = 4, byrow = F)
                new.area <- grid.data$Area[c(lev[,2])] + grid.data$Area[c.new]
                
                
              } # if allocated flow is less than original calculated flow AND there are still neighbors to receive it
              
            } # if there is more than one downhill neighbor to receive leftover flow
            
          } # if calculated flow is greater than the least equilibrating volume
          
        } # if there is at least one positive drop
        
      } # if there is standing water in current cell
      
    } # for every cell
    
  } # for each subdaily timestep
  
  if (d < days) {
    stage.day[d+1,] <- stage.day.temp
  }
  
}     # for every day

#####################################################################
#####################################################################
#####################################################################
#####################################################################

## Create .csv with daily stage values for each grid
stage.day <- round(stage.day, 4)
write.table(stage.day, file = paste("Stages-", L, ".csv", sep = ""), sep = ",", row.names = F, 
            col.names = c(1:nrow(grid.data)))

## Create .csv with daily water depth values for each grid
depth.day <- matrix(nrow = (days), ncol = ncell) # a column for each cell, a row for each day

for (i in 1:length(elevs)) {
  
  depth.day[, i] <- stage.day[, i] - elevs[i]
  
}

depth.day <- round(depth.day, 4)
write.table(depth.day, file = paste("Depths-", L, ".csv", sep = ""), sep = ",", row.names = F, 
            col.names = c(1:nrow(grid.data)))
