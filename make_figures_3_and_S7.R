### These scripts produce Figure 3 and Figure S7 in Smadja et al. 2022 BIORXIV-2022-500634 ###


library(gtrellis)
library(circlize)


#### load data ####
## Note: make sure all bed files are in the form chr1 chr2 etc for chromosomeID and not 1, 2, etc


# load local score information with resampling xi=1 - only significant regions:
bed_localsc1Resampl = read.table("Significant_zones0.05_pvalnoFDR_xi1_resampling-chr.bed",header=T,na.string=".")

# load outlier gene information
bed_genesCmax = read.table("outlierGenesCmax.bed",header=T,na.string=".")
bed_genesCmean = read.table("outlierGenesCmean.bed",header=T,na.string=".")

# load significant gene clusters information
bed_clusters = read.table("SigClusters_outliergenes_500kb_50kb_chr.bed",header=T,na.string=".")

# load coldspots of recombination
bed_r_cold = read.table("coldspots.bed",header=F,na.string=".")

# load sterility QTL regions
bed_sterility = read.table("Sterility_regions_mm10.bed",header=T,na.string=".")


require(tidyverse)
bed_localsc1Resampl %>% mutate(start = beg) -> bed_localsc1Resampl
bed_genesCmax %>% glimpse
bed_genesCmean %>% glimpse
bed_clusters %>% glimpse


#### produce Figure 3 ####

pdf("Figure3.pdf",height = 7,width = 12)

gtrellis_layout(n_track=7, 
                equal_width = FALSE,
                track_axis = c(FALSE, TRUE, TRUE, TRUE, FALSE,FALSE,FALSE), 
                #track_ylab = c("", "LC", "C2max", "C2mean", "GC", "R", "Sterility_regions"), 
                track_ylim = c(0, 1, 0, 14000, 0, 8, 0, 8, 0, 1, 0, 1, 0, 1),
                track_height = unit.c(1.5*grobHeight(textGrob("chr1")), unit(2, "null"), unit(2, "null"), unit(2, "null"), unit(2, "mm"), unit(1, "mm"),unit(1, "mm")), 
                nrow=4, 
                byrow = TRUE, 
                category = paste0("chr", c(1:19, "X")),  
                species = "mm10",  
                xlab = "Genomic position", 
                asist_ticks = FALSE,
                compact = TRUE,
                gap = unit(c(5, 2), "mm")
)

# add chromosome tracks
add_track(panel_fun = function(gr) {
  chr = get_cell_meta_data("name")  
  grid.rect(gp = gpar(fill = "#EEEEEE"))
  grid.text(chr)
})

# add local score
add_points_track(bed_localsc1Resampl, bed_localsc1Resampl[[4]], 
                 pch = 16, size = unit(1.5, "mm"), 
                 gp = gpar(col = alpha("#228833",0.99))) 

# add Cmax outlier genes 
add_points_track(bed_genesCmax, bed_genesCmax[[7]],
                 size = unit(1.8, "mm"),
                 pch = 16,
                 gp = gpar(col= ifelse(grepl("Vmn",bed_genesCmax[[5]], fixed=TRUE), 
                                       alpha("#4477AA",0.4), 
                                       ifelse(grepl("Olfr",bed_genesCmax[[5]], fixed=TRUE), alpha("#AA3377",0.4), alpha("#CCBB44",0.4))))) 

# add Cmean outlier genes 
add_points_track(bed_genesCmean, bed_genesCmean[[9]],
                 size = unit(1.8, "mm"), 
                 pch = 16, 
                 gp = gpar(col = ifelse(grepl("Vmn",bed_genesCmean[[5]], fixed=TRUE),
                                        alpha("#4477AA",0.4), 
                                      ifelse(grepl("Olfr",bed_genesCmean[[5]], fixed=TRUE),alpha("#AA3377",0.4), alpha("#999933",0.4)))))

# add significant clusters
add_track(bed_clusters,panel_fun = function(gr) {
  genrang = gr
  grid.rect(genrang[[2]], 
            width=genrang[[3]]-genrang[[2]], height = unit(1, "npc"), 
            default.units = "native", hjust = 0, vjust = 0.5, 
            gp = gpar(fill = "#997700",col="#997700"),draw = TRUE)
})


# add coldspots of recombination
add_track(bed_r_cold,panel_fun = function(gr) {
  genrang = gr
  grid.rect(genrang[[2]], width=genrang[[3]]-genrang[[2]], height = unit(1, "npc"), default.units = "native", hjust = 0, vjust = 0.5, gp = gpar(fill = "#33BBEE",col="#33BBEE", alpha=0.5),draw = TRUE)
})

# add sterility regions
add_track(bed_sterility,panel_fun = function(gr) {
  genrang = gr
  grid.rect(genrang[[2]], width=genrang[[3]]-genrang[[2]], height = unit(1, "npc"), default.units = "native", hjust = 0, vjust = 0.5, gp = gpar(fill = "#666666",col="#666666"),draw = TRUE)
})

dev.off()



#### Produce Figure S7 - zoom on outlier main gene clusters on chromosomes 2, 7, and 9 ####

## chr2

pdf("FigureS7_chr2.pdf",height=6,width =15)
zoom = function(df) {
  gtrellis_layout(data = df, 
                  n_track=4, 
                  nrow=4, 
                  byrow = TRUE, 
                  track_axis = c(FALSE, TRUE, TRUE, TRUE), 
                  track_ylim = c(0, 1, 0, 14000, 0, 8, 0, 8), 
                  track_height = unit.c(1.5*grobHeight(textGrob("chr1")), unit(2, "null"), unit(2, "null"), unit(2, "null")), 
                  xlab = "Genomic position",
                  track_ylab = c("", "local scores", "C2_max values", "C2_mean values"))
  
  # add chromosome tracks
  add_track(panel_fun = function(gr) {
    chr = get_cell_meta_data("name")  
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr)
  })
  
  # add local score
  add_segments_track(bed_localsc1Resampl,bed_localsc1Resampl[[4]], 
                     gp = gpar(col = alpha("#228833",0.99), lwd =4)) 
  
  # add outlier genes Cmax
  add_points_track(bed_genesCmax, bed_genesCmax[[7]],
                   size = unit(3, "mm"),
                   pch = 16,
                   gp = gpar(col= ifelse(grepl("Vmn",bed_genesCmax[[5]], fixed=TRUE), 
                                         alpha("#4477AA",0.99),
                                         ifelse(grepl("Olfr",bed_genesCmax[[5]], fixed=TRUE), alpha("#AA3377",0.99), alpha("#CCBB44",0.99))))) 
  # add outlier genes Cmean
  add_points_track(bed_genesCmean, bed_genesCmean[[9]],
                   size = unit(3, "mm"), 
                   pch = 16, 
                   gp = gpar(col = ifelse(grepl("Vmn",bed_genesCmean[[5]], fixed=TRUE),
                                          alpha("#4477AA",0.99), 
                                          ifelse(grepl("Olfr",bed_genesCmean[[5]], fixed=TRUE),alpha("#AA3377",0.99), alpha("#999933",0.99)))))
}

df = data.frame(chr = c("chr2"),
                start = c(89800000),
                end = c(90400000))

zoom(df)
dev.off()


## chr7

pdf("FigureS7_chr7.pdf",height=6,width =15)
zoom = function(df) {
  gtrellis_layout(data = df, 
                  n_track=4, 
                  nrow=4, 
                  byrow = TRUE, 
                  track_axis = c(FALSE, TRUE, TRUE, TRUE), 
                  track_ylim = c(0, 1, 0, 12000, 0, 6, 0, 4), 
                  track_height = unit.c(1.5*grobHeight(textGrob("chr1")), unit(2, "null"), unit(2, "null"), unit(2, "null")), 
                  xlab = "Genomic position",
                  track_ylab = c("", "local scores", "C2_max values", "C2_mean values"))
  
  # add chromosome tracks
  add_track(panel_fun = function(gr) {
    chr = get_cell_meta_data("name")  
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr)
  })
  
  # add local score
  add_segments_track(bed_localsc1Resampl,bed_localsc1Resampl[[4]], 
                     gp = gpar(col = alpha("#228833",0.99), lwd =4)) 
  
  # add outlier genes Cmax
  add_points_track(bed_genesCmax, bed_genesCmax[[7]],
                   size = unit(3, "mm"),
                   pch = 16,
                   gp = gpar(col= ifelse(grepl("Vmn",bed_genesCmax[[5]], fixed=TRUE), 
                                         alpha("#4477AA",0.99), 
                                         ifelse(grepl("Olfr",bed_genesCmax[[5]], fixed=TRUE), alpha("#AA3377",0.99), alpha("#CCBB44",0.99))))) 
  
  # add outlier genes Cmean
  add_points_track(bed_genesCmean, bed_genesCmean[[9]],
                   size = unit(3, "mm"), 
                   pch = 16, 
                   gp = gpar(col = ifelse(grepl("Vmn",bed_genesCmean[[5]], fixed=TRUE),
                                          alpha("#4477AA",0.99), 
                                          ifelse(grepl("Olfr",bed_genesCmean[[5]], fixed=TRUE),alpha("#AA3377",0.99), alpha("#999933",0.99)))))
}

df = data.frame(chr = c("chr7"),
                start = c(84450000),
                end = c(86750000))

zoom(df)
dev.off()


## chr9

pdf("FigureS7_chr9.pdf",height=6,width =15)
zoom = function(df) {
  gtrellis_layout(data = df, 
                  n_track=4, 
                  nrow=4, 
                  byrow = TRUE, 
                  track_axis = c(FALSE, TRUE, TRUE, TRUE), 
                  track_ylim = c(0, 1, 0, 14000, 0, 8, 0, 8), 
                  track_height = unit.c(1.5*grobHeight(textGrob("chr1")), unit(2, "null"), unit(2, "null"), unit(2, "null")), 
                  xlab = "Genomic position",
                  track_ylab = c("", "local scores", "C2_max values", "C2_mean values"))
  
  # add chromosome tracks
  add_track(panel_fun = function(gr) {
    chr = get_cell_meta_data("name")  
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr)
  })
  
  # add local score
  add_segments_track(bed_localsc1Resampl,bed_localsc1Resampl[[4]], 
                     gp = gpar(col = alpha("#228833",0.99), lwd =4)) 
  
  # add outlier genes Cmax
  add_points_track(bed_genesCmax, bed_genesCmax[[7]],
                   size = unit(3, "mm"),
                   pch = 16,
                   gp = gpar(col= ifelse(grepl("Vmn",bed_genesCmax[[5]], fixed=TRUE), 
                                         alpha("#4477AA",0.99), 
                                         ifelse(grepl("Olfr",bed_genesCmax[[5]], fixed=TRUE), alpha("#AA3377",0.99), alpha("#CCBB44",0.99))))) 
  
  # add outlier genes Cmean
  add_points_track(bed_genesCmean, bed_genesCmean[[9]],
                   size = unit(3, "mm"), 
                   pch = 16, 
                   gp = gpar(col = ifelse(grepl("Vmn",bed_genesCmean[[5]], fixed=TRUE),
                                          alpha("#4477AA",0.99), 
                                          ifelse(grepl("Olfr",bed_genesCmean[[5]], fixed=TRUE),alpha("#AA3377",0.99), alpha("#999933",0.99)))))
}


df = data.frame(chr = c("chr9"),
                start = c(37350000),
                end = c(39150000))
zoom(df)
dev.off()

##### end ####