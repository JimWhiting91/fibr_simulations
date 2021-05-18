# Read in and visualise mutation frequency trajectories
lib <- c("ggridges","data.table","ggplot2","dplyr")
lapply(lib,library,character.only=T)

res <- data.frame(fread("outputs/LL_test_run.res",fill=T,skip = 15))

# Bin our data into
start_intervals <- cut_interval(res[res$generation==1,"freq"],length=0.1)
names(start_intervals) <- res[res$generation==1,"mutation_id"]

# Add these to our original data
res$start_interval <- NA
for(mut in names(start_intervals)){
  res[res$mutation_id == mut,"start_interval"] <- as.character(start_intervals[mut])
}

# And now visualise...
no_new_res <- na.omit(res)
ggplot(no_new_res,aes(x=generation,y=freq))+
  facet_wrap(~start_interval,ncol=2)+
  geom_line(aes(group = mutation_id),alpha=0.1)

# Add to this, the mean freq for each group through time
mean_freqs <- no_new_res %>% group_by(start_interval,generation) %>% summarise()

# What about allele enrichment
calc_allele_enrich  <- function(x){
  
  # First work out direction
  af_dir <- ifelse(x[x$generation==1,"freq"] < x[x$generation==max(x$generation),"freq"],"up","down")
  final_freq <- ifelse(any(x$generation == max_gen),x[x$generation==max(x$generation),"freq"],0)
  
  # Caculate enrichment
  if(af_dir == "up"){
    af_enrich <- (final_freq - x$freq[which(x$generation == 1)]) / x$freq[which(x$generation == 1)]
  } else {
    af_enrich <- (final_freq - x$freq[which(x$generation == 1)]) / (1-x$freq[which(x$generation == 1)])
  }
  return(data.frame(mut_id=unique(x$mutation_id),
                    af_enrich=af_enrich,
                    start_interval=unique(x$start_interval),
                    start_af=x[x$generation==1,"freq"],
                    end_af=final_freq))
}
max_gen <- max(no_new_res$generation)

# Appy to all mutations
muts <- unique(no_new_res$mutation_id)
enrich_res <- data.frame(rbindlist(lapply(1:length(muts),function(x){
  print(x)
  calc_allele_enrich(no_new_res[no_new_res$mutation_id == muts[x],])
})))

# Plot
hist(enrich_res$af_enrich)
ggplot(enrich_res,)

###################################################
# Effect of selection
# Get all res
res_inputs <- list.files("outputs/",pattern = "LL_test_run_with_selection_rnorm_0.05selected")
res <- data.frame(rbindlist(lapply(res_inputs,function(path){
  
  # Fetch
  tmp <- data.frame(fread(paste0("outputs/",path),fill=T,skip = 15))
  colnames(tmp) <- c("mutation_id","freq","sel_coef","generation")
  tmp$run <- path
  
  # Bin our data into
  start_intervals <- cut_interval(tmp[tmp$generation==1,"freq"],length=0.1)
  start_intervals_df <- data.frame(mutation_id=tmp[tmp$generation==1,"mutation_id"],
                                   start_interval=start_intervals)
  
  # Add these to our original data
  tmp_merge <- merge(tmp,start_intervals_df,by = "mutation_id")
  
  return(tmp_merge)
})))

# Add a mutation_run identifier
res$mutation_id_run <- paste0(res$mutation_id,":",res$run)

# Visualise
ggplot(res,aes(x=generation,y=freq,colour=sel_coef))+
  facet_wrap(~start_interval,ncol=2)+
  geom_line(aes(group = mutation_id_run),alpha=0.1)+
  scale_colour_gradientn(colours = terrain.colors(10))

# Visualise only selected
ggplot(res[res$sel_coef != 0,],aes(x=generation,y=freq,colour=sel_coef))+
  facet_wrap(~start_interval,ncol=2)+
  geom_line(aes(group = mutation_id_run),alpha=0.5)+
  scale_colour_gradient2()

# plot final frequency...
ggplot(res[res$generation==max(res$generation),],aes(x=freq,y=start_interval,fill=sel_coef))+
  #  facet_wrap(~start_interval,ncol=2)+
  # # geom_line(aes(group = mutation_id_run),alpha=0.1)+
  #  geom_density
  #  scale_colour_gradientn(colours = terrain.colors(10))
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Selection Coef")

# plot final frequency...
ggplot(res[res$generation==max(res$generation) & res$sel_coef != 0,],aes(y=freq,x=sel_coef))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~start_interval,ncol=2)

# stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
# scale_fill_viridis_c(name = "Selection Coef")


# With this, we also want to look at the difference between start and end AF
all_muts <- unique(res$mutation_id_run)
max_gen <- max(res$generation)
calc_AF_change <- function(x){
  if(any(x$generation == max_gen)){
    af_change <- x$freq[which(x$generation == max_gen)] - x$freq[which(x$generation == 1)]
  } else {
    af_change <- 0 - x$freq[which(x$generation == 1)]
  }
  return(data.frame(AF_change=af_change,
                    start_interval=unique(x$start_interval),
                    sel_coef=unique(x$sel_coef)))
}

start_end_af <- res %>% group_by(mutation_id_run) %>% group_map(~calc_AF_change(.x))

# Plot these...
plot_af_change <- data.frame(rbindlist(start_end_af))
ggplot(plot_af_change[plot_af_change$sel_coef != 0,],aes(x=sel_coef,y=AF_change))+
  facet_wrap(~start_interval,ncol=2)+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0,linetype="dotted")

# Bin our data into
plot_af_change$sel_coef_interval <- cut_interval(plot_af_change$sel_coef,n=10)


# Average across all of the estimates
plot_af_change_avg <- data.frame(plot_af_change %>% 
                                   group_by(start_interval,sel_coef_interval) %>%
                                   summarise(mean_AF=mean(AF_change),
                                             sd_AF=sd(AF_change),
                                             se_AF=sd(AF_change)/sqrt(length(AF_change))))

# Plot
ggplot(plot_af_change_avg,aes(x=sel_coef_interval,y=mean_AF))+
  facet_wrap(~start_interval,ncol=2,scales = "free_y")+
  geom_point()+
  geom_errorbar(aes(ymin=mean_AF-se_AF,ymax=mean_AF+se_AF))
# geom_errorbar(aes(ymin=mean_AF-2*sd_AF,ymax=mean_AF+2*sd_AF))

##### Allele enrichment and selection ######
#af_enrich_selection <- res %>% group_by(mutation_id_run) %>% group_map(~calc_allele_enrich(.x))
af_enrich_selection <- pbmcapply::pbmclapply(unique(res$mutation_id_run),function(x){
  calc_allele_enrich(res[res$mutation_id_run == x,])
},mc.cores = 6)

# Plot these...
plot_af_enrich <- data.frame(rbindlist(af_enrich_selection))
ggplot(plot_af_change[plot_af_change$sel_coef != 0,],aes(x=sel_coef,y=AF_change))+
  facet_wrap(~start_interval,ncol=2)+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept=0,linetype="dotted")

# Bin our data into
plot_af_change$sel_coef_interval <- cut_interval(plot_af_change$sel_coef,n=10)


