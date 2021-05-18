# Convert FIBR census data into demography matrix for simulations
lib <- c("data.table")
lapply(lib,library,character.only=T)

# Fetch census csv
census <- read.csv("data/genomicscensusdata.csv")
streams <- unique(census$stream)

# What generation time are we using in months?
gen_time <- 8

# Build a demography for each river
river_demogs <- lapply(streams,function(stream){
  
  # Subset and chunk into gen_time month windows...
  census_sub <- census[census$stream == stream,]
  chunk_seq <- seq(1,nrow(census_sub),gen_time)
  chunk_seq2 <- chunk_seq + gen_time - 1
  chunk_seq2[length(chunk_seq2)] <- nrow(census_sub)
  
  # For each chunk, we want to take the maximum size of the population
  chunk_demos <- data.frame(rbindlist(lapply(1:length(chunk_seq),function(i){
    #print(i)
    
    # Take the largest size estimate during the chunk
    largest <- which(rowSums(census_sub[chunk_seq[i]:chunk_seq2[i],c("F.N.hat","M.N.hat")]) == max(rowSums(census_sub[chunk_seq[i]:chunk_seq2[i],c("F.N.hat","M.N.hat")]),na.rm = T))
    
    # Approximate standard deviation
    male_sd <- (census_sub[largest,c("M.N.ucl")] - census_sub[largest,c("M.N.lcl")]) / 3.92
    female_sd <- (census_sub[largest,c("F.N.ucl")] - census_sub[largest,c("F.N.lcl")]) / 3.92
    
    # Use this as our estimate
    out <- data.frame(generation = i,
                      female_mean = census_sub[largest,c("F.N.hat")],
                      female_sd = female_sd,
                      male_mean = census_sub[largest,c("M.N.hat")],
                      male_sd = male_sd)
  })))
  
  # Save
  write.table(chunk_demos,
              paste0("data/",stream,"_simulation_demography.txt"),
              row.names = F,quote = F,sep = "\t")
})
