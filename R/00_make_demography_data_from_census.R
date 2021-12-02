# Convert FIBR census data into demography matrix for simulations
lib <- c("data.table","ggplot2")
lapply(lib,library,character.only=T)

# Fetch census csv
census <- read.csv("data/genomicscensusdata.csv")
streams <- unique(census$stream)

# Plot all together first
ggplot(census[census$year <= 2013,],aes(x=sampling,y=F.N.hat+M.N.hat,colour=stream))+
  geom_line()

# What generation time are we using in months?
gen_time <- 8

# First trim for pre 2013, and then work out how many whole generations we can run
census <- census[census$year <= 2013,]

# Build a demography for each river
river_demogs <- lapply(streams,function(stream){
  
  # Subset and chunk into gen_time month windows...
  census_sub <- census[census$stream == stream,]
  chunk_seq <- seq(1,nrow(census_sub),gen_time)
  chunk_seq2 <- chunk_seq + gen_time - 1
  chunk_seq2[length(chunk_seq2)] <- nrow(census_sub)
  
  # Remove any sets where we don't have a full generation...
  chunk_seq <- chunk_seq[which((chunk_seq2 - chunk_seq) == gen_time-1)]
  chunk_seq2 <- chunk_seq2[1:length(chunk_seq)]
  
  # For each chunk, we want to take the maximum size of the population
  chunk_demos <- data.frame(rbindlist(lapply(1:length(chunk_seq),function(i){
    #print(i)
    
    # # Take the largest size estimate during the chunk
    # largest <- (chunk_seq[i]:chunk_seq2[i])[which(rowSums(census_sub[chunk_seq[i]:chunk_seq2[i],c("F.N.hat","M.N.hat")]) == max(rowSums(census_sub[chunk_seq[i]:chunk_seq2[i],c("F.N.hat","M.N.hat")]),na.rm = T))]
    # 
    # # Approximate standard deviation
    # male_sd <- (census_sub[largest,c("M.N.ucl")] - census_sub[largest,c("M.N.lcl")]) / 3.92
    # female_sd <- (census_sub[largest,c("F.N.ucl")] - census_sub[largest,c("F.N.lcl")]) / 3.92
    # 
    # # Use this as our estimate
    # out <- data.frame(generation = i,
    #                   female_mean = census_sub[largest,c("F.N.hat")],
    #                   female_sd = female_sd,
    #                   male_mean = census_sub[largest,c("M.N.hat")],
    #                   male_sd = male_sd)
    
    # Approximate standard deviation
    male_sd <- (census_sub[chunk_seq[i]:chunk_seq2[i],c("M.N.ucl")] - census_sub[chunk_seq[i]:chunk_seq2[i],c("M.N.lcl")]) / 3.92
    female_sd <- (census_sub[chunk_seq[i]:chunk_seq2[i],c("F.N.ucl")] - census_sub[chunk_seq[i]:chunk_seq2[i],c("F.N.lcl")]) / 3.92

    # Use this as our estimate
    out <- data.frame(generation = i,
                      female_mean = mean(census_sub[chunk_seq[i]:chunk_seq2[i],c("F.N.hat")]),
                      female_sd = mean(female_sd),
                      male_mean = mean(census_sub[chunk_seq[i]:chunk_seq2[i],c("M.N.hat")]),
                      male_sd = mean(male_sd))
  })))
  
  # Save
  write.table(chunk_demos,
              paste0("data/",stream,"_simulation_demography.txt"),
              row.names = F,quote = F,sep = "\t")
  
  return(chunk_demos)
})

# Plot all the chunk demos...
demo_plots <- lapply(river_demogs,function(river){
  
  ggplot(river,aes(generation,male_mean+female_mean))+
    geom_line(colour="blue2")+
    geom_ribbon(aes(ymin=(male_mean-2*male_sd)+(female_mean-2*female_sd),
                    ymax=(male_mean+2*male_sd)+(female_mean+2*female_sd)),alpha=0.5)
  
})
cowplot::plot_grid(plotlist=demo_plots,ncol=2)

all_demos <- data.frame(rbindlist(lapply(1:length(river_demogs),function(x){
  river_demogs[[x]]$stream = streams[x]
  return(river_demogs[[x]])
})))
ggplot(all_demos,aes(generation,male_mean+female_mean,colour=stream,fill=stream))+
  geom_line(colour="blue2")+
  geom_ribbon(aes(ymin=(male_mean-2*male_sd)+(female_mean-2*female_sd),
                  ymax=(male_mean+2*male_sd)+(female_mean+2*female_sd)),alpha=0.5)
