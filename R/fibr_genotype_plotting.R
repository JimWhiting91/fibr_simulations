# Genotype Plots for FIBR chr15 region of interest
library(GenotypePlot)
library(cowplot)
library(vcfR)
library(ggtree)
library(ggplot2)

# Fetch files
vcf_path <- "data/chr15_5Mb_region.recode.vcf.gz"
vcf_in <- read.vcfR(vcf_path)
popmap <- read.table("data/fibr_chr15.popmap",header=T)

# Set up all our colours here...
library(RColorBrewer)
fibr_colour_rivers<-data.frame(river = c("Tacarigua","Guanapo","Aripo","Oropouche","Madamas","IC","IT","IUL","ILL"),
                                     colour = c(brewer.pal(n = 8, name = "Dark2")[1:3],
                                                brewer.pal(n = 8, name = "Dark2")[6],
                                                brewer.pal(n = 8, name = "Dark2")[8],
                                                "#756bb1","#bcbddc","#fa9fb5","#dd3497")) 
fibr_colour_rivers$river_F<-factor(fibr_colour_rivers$river,levels= fibr_colour_rivers$river)

# Make base plots
# Raw genotypes
geno_figs <- genotype_plot(vcf_object = vcf_in[is.biallelic(vcf_in),],
                           popmap = popmap,
                           cluster = F,
                           snp_label_size = 50000,invariant_filter = T,missingness = 0.2,polarise_genotypes = "GLP")

# Raw allele frequencies
af_figs <- genotype_plot(vcf_object = vcf_in[is.biallelic(vcf_in),],
                           popmap = popmap,
                           cluster = F,
                           snp_label_size = 50000,
                          invariant_filter = T,
                          missingness = 0.2,
                          polarise_genotypes = "GLP",
                          plot_allele_frequency = T)

# Make a cluster over just the region highlighted
region_to_keep <- which(as.integer(vcf_in@fix[,2]) >= 5066973 & 
                          as.integer(vcf_in@fix[,2]) <= 5135080)
region_figs <- genotype_plot(vcf_object = vcf_in[region_to_keep,],
                           popmap = popmap[popmap$pop %in% c("GHP","GLP","IC","IT","ILL","IUL"),],
                           cluster = T,
                           polarise_genotypes = "GLP",
                           snp_label_size = 5000,invariant_filter = T,missingness = 0.2)

# Modify the dendrogram...
tip_metadata <- data.frame(tip=region_figs$dendro_labels)
tip_metadata$pop <- NA
tip_metadata$river <- NA 
tip_metadata$predation <- NA
for(i in 1:nrow(tip_metadata)){
  tip_metadata$pop[i] <- popmap[popmap$ind == tip_metadata$tip[i],"pop"]
}
river_codes <- c("MAD","G","TAC","O","AP","IC","IT","IUL","ILL")
names(river_codes) <- c("Madamas","Guanapo","Tacarigua","Oropouche","Aripo","IC","IT","IUL","ILL")
for(river in river_codes){
  tip_metadata[grep(river,tip_metadata$pop),"river"] <- names(river_codes)[river_codes == river]
}
tip_metadata$predation <- "Natural LP"
tip_metadata$predation[grep("HP",tip_metadata$pop)] <- "Natural HP"
tip_metadata[tip_metadata$pop %in% c("IC","IT","IUL","ILL"),"predation"] <- "Intro LP"
tip_metadata$tip_N <- 1:nrow(tip_metadata)

# Add tip labels
dendro_with_tree <- region_figs$dendrogram + 
  geom_jitter(aes(y=-2.5,x=tip_metadata$tip_N,colour=tip_metadata$river,shape=tip_metadata$predation),
              size = 3,height = 2.2,width=0,alpha=0.75)+
  scale_colour_manual(breaks=fibr_colour_rivers$river_F,
                      values=fibr_colour_rivers$colour)+
  scale_shape_manual(breaks=c("Natural HP","Natural LP","Intro LP"),
                     values=c(19,17,15))+
  theme(legend.position="left",
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15))+
  labs(colour="River",shape="Predation")+
  guides(colour = guide_legend(override.aes = list(size=8)),
         shape = guide_legend(override.aes = list(size=8)))

# Add dendrogram tips
dendro_with_tips <- region_figs$dendrogram +
  geom_text(aes(x=1:length(region_figs$dendro_labels),
                y=-2.5,
                label=region_figs$dendro_labels))

# Combine
region_figs$dendrogram <- dendro_with_tips
combined_regions <- plot_grid(region_figs$dendrogram,
                              region_figs$genotypes,
                              ncol=2,align="h",axis = "tb")

# Find out where we need to put our annotated box
per_site_missing <- apply(extract.gt(vcf_in), MARGIN = 1, function(x){ sum(is.na(x)) })
per_site_missing <- per_site_missing/(ncol(vcf_in@gt)-1)
keep_missing <- which(per_site_missing < 0.2)
vcf_in_missing <- vcf_in[keep_missing,]

bp <- as.integer(vcf_in_missing@fix[,2])
plotting_pos<-seq(min(bp),max(bp),by=(max(bp)-min(bp))/length(bp))[1:length(bp)]
box_start<-plotting_pos[min(which(as.integer(vcf_in_missing@fix[,2]) >= 5066973))]
box_end<-plotting_pos[max(which(as.integer(vcf_in_missing@fix[,2]) <= 5135080))]
geno_figs$genotypes + annotate("rect", xmin = box_start, xmax = box_end, ymin = -Inf, ymax = Inf,alpha = .5)

# Merge them all together now
pdf("figs/FigureSX_genotype_plots_chr15_fibr_vs_natural.pdf",width=15,height=18)
plot_grid(plot_grid(geno_figs$positions,
                    geno_figs$genotypes + annotate("rect", xmin = box_start, xmax = box_end, ymin = -Inf, ymax = Inf,alpha = .25),
                    af_figs$genotypes + annotate("rect", xmin = box_start, xmax = box_end, ymin = -Inf, ymax = Inf,alpha = .25),
                    align = "v",axis = "tblr",ncol=1,nrow=3,rel_heights = c(0.5,6,5),labels=c("A","","B"),label_size = 36),
          combined_regions,
          ncol=2,align="h",axis = "tblr",nrow=1,rel_widths = c(0.8,1),labels = c("","C"),label_size = 36)
dev.off()
