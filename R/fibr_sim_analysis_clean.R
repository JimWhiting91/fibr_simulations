# Derive neutral simulated distributions of selection coefficients...
lib <- c("ggridges","data.table","ggplot2","dplyr","parallel","Rfast","cowplot")
lapply(lib,library,character.only=T)

# results directory
sim_res <- "outputs/21,11,29_FIBR_runs_without_selection_resolveDemo_MeanPopSizes_withGen0_fixGen0/"
fibr_res <- list.files(sim_res)

# Set up colours
fibr_cols <- data.frame(pop=c("IC","IT","IUL","ILL","GHP"),
                        col=c("#756bb1","#bcbddc","#fa9fb5","#dd3497","#088da5"),
                        pop2=c("CA","TY","UL","LL","GHP"))

#### Selection Coef Function ####
calc_sel_coef <- function(x){
  
  # Get parameters
  tau = max(x$generation) - min(x$generation)
  p0 = x[x$generation==min(x$generation),"freq"]
  q0 = 1-p0
  pt = x[x$generation==max(x$generation),"freq"]
  qt = 1-pt
  
  # And Calc
  s = (2/tau)*log((pt*q0)/(qt*p0))
  return(data.frame(mutation_id=unique(x$mutation_id),
                    sel_coef=s))
}

calc_sel_coef_obs <- function(p0,pt,tau){
  
  # Get parameters
  # tau = max(x$generation) - min(x$generation)
  # p0 = x[x$generation==min(x$generation),"freq"]
  q0 = 1-p0
  # pt = x[x$generation==max(x$generation),"freq"]
  qt = 1-pt
  
  # And Calc
  s = (2/tau)*log((pt*q0)/(qt*p0))
  return(s)
}

#### Representative plot of AF trajectories ####
# Just plot a representative sample of 10 runs
C_sample <- grep("CA_",fibr_res,value=T)[2:11]
C_sample_res <- data.frame(rbindlist(lapply(C_sample,function(run){
  
  # Fetch
  tmp <- data.frame(fread(paste0(sim_res,"/",run),fill=T,skip = 16))
  colnames(tmp) <- c("mutation_id","freq","generation")
  tmp$run <- run
  
  # Bin our data into
  start_intervals <- cut_interval(tmp[tmp$generation==0,"freq"],length=0.1)
  start_intervals_df <- data.frame(mutation_id=tmp[tmp$generation==0,"mutation_id"],
                                   start_interval=start_intervals)
  
  # Add these to our original data
  tmp_merge <- merge(tmp,start_intervals_df,by = "mutation_id")
  tmp_merge$mutation_id_full <- paste0(tmp_merge$run,"_",tmp_merge$mutation_id)
  return(tmp_merge)
})))

# Plot em
af_traj_plot <- ggplot(C_sample_res[C_sample_res$generation > 0,],aes(x=generation,y=freq))+
  geom_line(aes(group=mutation_id_full),alpha=0.1)+
  facet_wrap(~start_interval,ncol=2)+
  theme_minimal()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text = element_text(size=14))+
  labs(x="Generations Post Introduction",y="AF")

af_traj_plot

#### Calculate neutral selection coefficients ####
intro_pops <- c("LL","UL","CA","TY")
intro_labs <- c("ILL","IUL","IC","IT")
new_intro <- c("Lower Lallaja","Upper Lallaja","Caigual","Taylor")

# Read in or run if first
if(!(file.exists("outputs/neutral_sim_SC.rds"))){
  pop_neutral_SC <- lapply(intro_pops,function(pop){
    print(paste0("STARTING ",pop))
    # Subset
    pop_runs <- grep(pop,fibr_res,value=T)
    
    # Across all runs, get selection coefficients for all mutations
    run_res <- data.frame(rbindlist(pbmcapply::pbmclapply(pop_runs,function(run){
      
      # Fetch
      tmp <- data.frame(fread(paste0(sim_res,"/",run),fill=T,skip = 16))
      colnames(tmp) <- c("mutation_id","freq","generation")
      tmp$run <- run
      
      # Subset
      min_gen <- 1
      max_gen <- max(tmp$generation)
      tmp <- tmp[tmp$generation %in% c(min_gen,max_gen),]
      
      # Also remove any mutations that are not represented at both extremes
      mut_counts <- table(tmp$mutation_id)
      tmp <- tmp[tmp$mutation_id %in% names(mut_counts[mut_counts==2]),]
      
      # We only want to track minor allele frequencies...
      to_change <- tmp[tmp$generation == 1 & tmp$freq > 0.5,"mutation_id"]
      tmp[tmp$mutation_id %in% to_change,"freq"] <- 1-tmp[tmp$mutation_id %in% to_change,"freq"]
      
      # Now for all mutations, calculate selection coefficients based on a diploid model with no dominance...
      sel_coef_res <- data.frame(rbindlist(lapply(unique(tmp$mutation_id),function(x) {
        #print(x)
        calc_sel_coef(tmp[tmp$mutation_id==x,])
      })))
      
    },mc.cores = 4)))
    
    run_res$pop <- pop
    
    return(run_res)
  })
  saveRDS(pop_neutral_SC,"outputs/neutral_sim_SC.rds")
} else {
  pop_neutral_SC <- readRDS("outputs/neutral_sim_SC.rds")
}

# Merge them all together
all_pop_sims <- data.frame(rbindlist(pop_neutral_SC))
distribution_stats <- data.frame(all_pop_sims[!(all_pop_sims$sel_coef %in% c(Inf,-Inf)),] %>% group_by(pop) %>% summarise(mean=mean(abs(sel_coef),na.rm = T),
                                                                                                                          var=var(abs(sel_coef),na.rm = T)))
distribution_stats

# Visualise by group
to_plot <- all_pop_sims
for(i in 1:length(new_intro)){
  to_plot$pop <- gsub(intro_pops[i],new_intro[i],to_plot$pop)
  distribution_stats$pop <- gsub(intro_pops[i],new_intro[i],distribution_stats$pop)
}
distribution_stats

# Proportion of polymorphisms lost/fixed ----------------------------------
pop_neutral_LP <- data.frame(rbindlist(lapply(intro_pops,function(pop){
  print(paste0("STARTING ",pop))
  # Subset
  pop_runs <- grep(pop,fibr_res,value=T)
  
  # Across all runs, estimate the total number of polymorphisms at start vs end
  run_res <- data.frame(rbindlist(pbmcapply::pbmclapply(pop_runs,function(run){
    
    # Fetch
    tmp <- data.frame(fread(paste0(sim_res,"/",run),fill=T,skip = 16))
    colnames(tmp) <- c("mutation_id","freq","generation")
    tmp$run <- run
    
    # Find the GH minor alleles...
    GH_major_muts <- tmp[tmp$generation==0 & tmp$freq > 0.5,"mutation_id"]
    tmp$maf_freq <- tmp$freq
    tmp[tmp$mutation_id %in% GH_major_muts,"maf_freq"] <- 1 - tmp[tmp$mutation_id %in% GH_major_muts,"maf_freq"]
    
    # Count start and end mutations
    start_muts <- nrow(tmp[tmp$generation==1,])
    start_end_muts <- tmp[tmp$generation %in% c(1,max(tmp$generation)),"mutation_id"]
    start_end_muts_counts <- table(start_end_muts)
    start_end_muts <- length(names(start_end_muts_counts)[start_end_muts_counts == 2])
    
    # Return these
    out <- data.frame(run=run,
                      pop=pop,
                      start=start_muts,
                      end=start_end_muts,
                      fixed_mafs=nrow(tmp[tmp$generation==max(tmp$generation) & tmp$maf_freq == 1,]))
    
  },mc.cores = 4)))
  
  # Get proportion lost/fixed
  run_res$lost_poly <- (run_res$start-run_res$end)/run_res$start
  
  return(run_res)
})))

# Clean output
for(i in 1:4){
  pop_neutral_LP$pop <- gsub(intro_pops[i],intro_labs[i],pop_neutral_LP$pop)
}
pop_neutral_LP$pop_F <- factor(pop_neutral_LP$pop,levels=intro_labs)

# Visualise
prop_lostpoly_fig <- ggplot(pop_neutral_LP,aes(y=lost_poly,x=pop_F,fill=pop_F))+
  geom_violin(draw_quantiles = c(0.025,0.5,0.975),show.legend = F)+
  theme_bw()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14))+
  labs(x="Introduced Population",
       y="Proportion of SNPs\nlost/fixed (per simulation)")+
  scale_fill_manual(breaks = fibr_cols$pop,
                    values = fibr_cols$col)
prop_lostpoly_fig

# Also plot the sum of fixed minor alleles
# Report means and variance
fixed_mafs <- pop_neutral_LP %>% group_by(pop) %>% summarise(mean=mean(lost_poly),
                                                             median=median(lost_poly),
                                                             var=var(lost_poly),
                                                             fixed_mafs=sum(fixed_mafs))
fixed_mafs

# Get max in tails...
max(pop_neutral_LP[pop_neutral_LP$pop=="IC","lost_poly"])
max(pop_neutral_LP[pop_neutral_LP$pop=="IT","lost_poly"])


shared_mut_res <- all_pop_sims[all_pop_sims$sel_coef != Inf,]
shared_muts <- unique(shared_mut_res$mutation_id)
length(shared_muts)

# Count up only mutations where all pops are presented
pop_counts <- sapply(shared_muts,function(x) length(unique(shared_mut_res[shared_mut_res$mutation_id==x,"pop"])))
shared_muts <- shared_muts[which(pop_counts==4)]
shared_mut_res <- shared_mut_res[shared_mut_res$mutation_id %in% shared_muts,]
shared_mut_res$mutation_id_pop <- paste0(shared_mut_res$mutation_id,"_",shared_mut_res$pop)


#### Estimated selection coefs over observed for all pops ####

# Per pop S
GH_AF <- data.frame(fread("data/AF/GHP_AF.txt"))
colnames(GH_AF) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
GH_AF$mutation_id <- paste0(GH_AF$chr,"_",GH_AF$pos)
GH_AF <- GH_AF[grep("chr",GH_AF$chr),]

# Fetch the minor af as lowest row min of both freqs...
GH_maf <- rowMins(as.matrix(GH_AF[,c("REF_FREQ","ALT_FREQ")]),value=T)
GH_maf_index <- rowMins(as.matrix(GH_AF[,c("REF_FREQ","ALT_FREQ")]),value=F)

# These intro pops have new names
intro_pops_af <- c("IC","IT","IUL","ILL")

# How many generations does each run for?
gen_max <- c(8,8,9,9)
names(gen_max) <- intro_pops_af

#if(!(file.exists("outputs/intro_pops_SC.rds"))){
popS <- data.frame(rbindlist(pbmcapply::pbmclapply(intro_pops_af,function(pop){
  
  # Read in AF
  intro_AF <- data.frame(fread(paste0("data/AF/",pop,"_AF.txt")))
  colnames(intro_AF) <- c("chr","pos","N_ALLELES","N_CHR","REF_ALLELE","REF_FREQ","ALT_ALLELE","ALT_FREQ")
  intro_AF$mutation_id <- paste0(intro_AF$chr,"_",intro_AF$pos)
  #intro_AF <- intro_AF[intro_AF$chr == "chr15",]
  intro_AF <- intro_AF[grep("chr",intro_AF$chr),]
  
  # Calculate the selection coefficient
  intro_tmp <- data.frame(GH_start_freq = GH_maf,
                          intro_end_freq =  intro_AF[,c("REF_FREQ","ALT_FREQ")][cbind(seq_along(GH_maf_index), GH_maf_index)])
  intro_tmp$s <- calc_sel_coef_obs(p0 = intro_tmp$GH_start_freq,
                                   pt=intro_tmp$intro_end_freq,
                                   tau=gen_max[pop])
  
  # Build output data.frame
  out <- data.frame(pop=pop,
                    chr=intro_AF$chr,
                    pos=intro_AF$pos,
                    mutation_id=intro_AF$mutation_id,
                    sel_coef=intro_tmp$s)
})))

# Combine with the expecteds
popS_densities <- popS[,c("mutation_id","sel_coef","pop")]
popS_densities$type <- "Observed"
popS_exp <- to_plot
popS_exp$type <- "Simulated"
popS_densities <- rbind(popS_densities,popS_exp)

# Fix the labelling again
for(i in 1:length(intro_pops_af)){
  popS_densities$pop <- gsub(new_intro[i],intro_pops_af[i],popS_densities$pop)
}
popS_densities$pop_F <- factor(popS_densities$pop,levels = c("ILL","IUL","IC","IT"))

# plot
density_cols <- RColorBrewer::brewer.pal(8, "Dark2")[c(1,6)]
popS_densities$type_F <- factor(popS_densities$type,levels=c("Simulated","Observed"))

#pdf("figs/FigSX_simulated_and_observed_sel_coefs_distributions.pdf",width=4,height=8)
sel_coef_density_figs <- ggplot(popS_densities,aes(x=sel_coef,colour=type_F))+
  geom_density(size=1.5)+
  facet_wrap(~pop_F,ncol=1,strip.position = "right",scales = "free_y")+
  #geom_text(data = distribution_stats, aes(label = sd_label), x = 0.3, y = 3, parse = TRUE, hjust = 0)+
  theme_bw()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text = element_text(size=14),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size=14))+
  labs(y="Density",x="Selection Coefficient",colour="")+
  geom_vline(xintercept = 0,linetype="dotted")+
  scale_color_manual(values=density_cols)
sel_coef_density_figs
#dev.off()

# Calculate sum per site
sum_per_site <- data.frame(popS[!(popS$sel_coef %in% c(Inf,-Inf)),] %>% group_by(chr,pos) %>% summarise(sel_coef=sum(sel_coef)))

#### Observed ####
# Run over all observed
obs_muts <- matrix(ncol=4,nrow=nrow(popS_densities[popS_densities$type != "Simulated",])/4)
for(i in 1:length(intro_pops_af)){
  #print(i)
  obs_muts[,i] <- popS_densities[popS_densities$type != "Simulated" & popS_densities$pop == intro_pops_af[i],"sel_coef"]
}
rownames(obs_muts) <-  popS_densities[popS_densities$type != "Simulated","mutation_id"][1:nrow(obs_muts)]

# Apply function over all
# Clean estimates to remove NAs and Infs
obs_muts <- data.frame(obs_muts[is.finite(rowSums(obs_muts)),])
obs_muts$total_sum <- rowSums(obs_muts)

# Plot density
hist(obs_muts$total_sum)

# Extract chr and pos
obs_muts$chr <- sapply(strsplit(rownames(obs_muts),"_"),'[[',1)
obs_muts$pos <- sapply(strsplit(rownames(obs_muts),"_"),'[[',2)

#### Simulated ####
# Can now permute over 1000 of the mutations, repeat each one 10 times...
length(shared_muts)
mut_perms = 1000
mut_iters = 10

simulated_selcoef_sum <- data.frame(rbindlist(pbmcapply::pbmclapply(sample(shared_muts,mut_perms),function(x){
  
  # Subset for only this mutation
  mut_sub <- shared_mut_res[shared_mut_res$mutation_id==x,]
  
  # Get selection coefficients for 10 iterations per mutation
  mut_iter_res <- data.frame(rbindlist(lapply(1:mut_iters,function(y){
    
    # Get a permuted vector to run geom over
    perm_vec <- sapply(unique(mut_sub$pop),function(pop){
      sample(mut_sub[mut_sub$pop == pop,"sel_coef"],1)
    })
    
    out_mat <- matrix(ncol=5,nrow=1)
    out_mat[1,] <- c(perm_vec,sum(perm_vec))
    colnames(out_mat) <- c(names(perm_vec),"rowSum")
    
    data.frame(out_mat)
  })))
  
  mut_iter_res$mutation_id=x
  mut_iter_res
})))

# Plot distribution
hist(simulated_selcoef_sum$rowSum)

# Can now take a cutoff from this distribution in absolute terms and highlight SNPs
sim_cutoff <- quantile(abs(simulated_selcoef_sum$rowSum[is.finite(simulated_selcoef_sum$rowSum)]),0.999)

# Get outlier SNPs
outlier_snps <- obs_muts[abs(obs_muts$total_sum) > sim_cutoff,]
nrow(outlier_snps)/nrow(obs_muts)

# How are these organised among chromosomes
outlier_chr_counts <- table(outlier_snps$chr)
sum(outlier_chr_counts)

# What is the expected based on chr size?
chr_sizes <- data.frame(obs_muts %>% group_by(chr) %>% summarise(max_size=max(as.integer(pos))))
chr_sizes$prop_size <- chr_sizes$max_size/sum(chr_sizes$max_size)

# Make and plot
chr_plot <- data.frame(chr=names(outlier_chr_counts),
                       obs=as.integer(outlier_chr_counts),
                       exp=chr_sizes$prop_size*nrow(outlier_snps))
chr_plot$chr_F <- factor(chr_plot$chr,levels=paste0("chr",1:23))
chr_plot2 <- reshape2::melt(chr_plot[,c("obs","exp","chr_F")])

# Plot as double bar
snp_outlier_by_chrom_fig <- ggplot(chr_plot2,aes(x=chr_F,y=value,fill=variable))+
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+
  scale_fill_brewer(palette="Blues")+
  theme(axis.text.x = element_text(size=14,angle=45),
        axis.text.y = element_text(size=14),
        axis.title = element_text(size=16),
        legend.title = element_blank(),
        legend.text = element_text(size=14))+
  labs(x="",y="SNP Outlier Count")

# And also just plot chr15 given prominence
chr15_muts <- obs_muts[obs_muts$chr == "chr15",]
chr15_muts$outlier <- "No"
chr15_muts[abs(chr15_muts$total_sum) > quantile(abs(simulated_selcoef_sum$rowSum[is.finite(simulated_selcoef_sum$rowSum)]),0.99),"outlier"] <- "99"
chr15_muts[abs(chr15_muts$total_sum) > quantile(abs(simulated_selcoef_sum$rowSum[is.finite(simulated_selcoef_sum$rowSum)]),0.999),"outlier"] <- "99.9"
chr15_muts[abs(chr15_muts$total_sum) > quantile(abs(simulated_selcoef_sum$rowSum[is.finite(simulated_selcoef_sum$rowSum)]),0.9999),"outlier"] <- "99.99"

chr15_snpsum_fig <- ggplot(chr15_muts[is.finite(chr15_muts$total_sum) & chr15_muts$outlier != "No",],aes(as.integer(pos),total_sum,colour=outlier))+
  geom_point()+
  scale_colour_manual(breaks = c("No","99","99.9","99.99"),
                      values = c("gray50","gold2","orange2","red4"))+
  theme_minimal()+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))+
  labs(y="Summed Selection",x="Chr15 Position (Mb)",colour="Simulated\nCutoff (%)")+
  scale_x_continuous(breaks = seq(0,max(as.integer(chr15_muts$pos)),5000000),
                     labels =  seq(0,max(as.integer(chr15_muts$pos)),5000000)/1000000)
chr15_snpsum_fig

# Save these data to RDS
saveRDS(chr15_muts,"outputs/chr15_calculated_selection_coefficients.rds")

# Fetch the maximum lowest pos on chr15 and get individual sel coefs
min_15 <- chr15_muts[is.finite(chr15_muts$total_sum) & chr15_muts$total_sum == min(chr15_muts[is.finite(chr15_muts$total_sum),"total_sum"]),]
max_15 <- chr15_muts[is.finite(chr15_muts$total_sum) & chr15_muts$total_sum == max(chr15_muts[is.finite(chr15_muts$total_sum),"total_sum"]),]



# Add a figure to the supps that includes the simulated generations -------
gen_popsizes <- data.frame(rbindlist(lapply(intro_pops,function(pop){
  tmp <- read.table(paste0("data/",pop,"_simulation_demography.txt"),header=T)
  tmp$pop=fibr_cols[fibr_cols$pop2==pop,"pop"]
  tmp
})))

gen_popsizes$total <- gen_popsizes$female_mean+gen_popsizes$male_mean
gen_popsizes$total_sd <- gen_popsizes$female_sd+gen_popsizes$male_sd

sim_demog_fig <- ggplot(gen_popsizes,aes(generation,y=total,colour=pop,fill=pop))+
  geom_ribbon(aes(ymin=total-2*total_sd,ymax=total+2*total_sd),alpha=0.5,colour=NA)+
  geom_line()+
  theme_bw()+
  scale_fill_manual(breaks = fibr_cols$pop,
                    values = fibr_cols$col)+
  scale_colour_manual(breaks = fibr_cols$pop,
                      values = fibr_cols$col)+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))+
  labs(y="Population Size",x="Sampling Generation",colour="Introduced\nPopulation",fill="Introduced\nPopulation")


# Separate Main and Supp --------------------------------------------------
pdf("figs/FigureX_ObsSNP_simulation_results.pdf",width=14,height=7)
plot_grid(snp_outlier_by_chrom_fig + 
            theme_bw(base_size = 22) +
            theme( 
              plot.title = element_text(hjust = 0.5),
              legend.position="right",
              legend.title = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size=12),
              strip.text = element_text(size=12),
              axis.title = element_text(size=14),
              axis.text.x = element_text(angle=45,hjust=1)),
          
          chr15_snpsum_fig  + 
            theme_bw(base_size = 22) +
            theme( 
              plot.title = element_text(hjust = 0.5),
              legend.position="right",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(size=12),
              strip.text = element_text(size=12),
              axis.title = element_text(size=14)),
          ncol=1,nrow=2,axis = "tblr",align="v",labels=c("A","B"),label_size=24)
dev.off()

pdf("figs/FigureSX_neutral_simulation_results.pdf",width=15,height=10)
plot_grid(sim_demog_fig,
          plot_grid(af_traj_plot,
                    sel_coef_density_figs,
                    prop_lostpoly_fig,
                    ncol=3,nrow=1,axis = "tblr",align="h",labels=c("B","C","D"),label_size=24),
          nrow=2,ncol=1,axis = "tblr",align="v",labels=c("A",""),label_size=24,rel_heights = c(1,1.5))
dev.off()
