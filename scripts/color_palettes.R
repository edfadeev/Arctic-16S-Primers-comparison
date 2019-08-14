#load library
library(RColorBrewer)


#####################################
#Color palettes for plots
#####################################
#large colours range
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#4477AA","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tol21rainbow<- c("#771155", 
                 "#AA4488", 
                 "#CC99BB", 
                 "#114477", 
                 "#4477AA", 
                 "#77AADD", 
                 "#117777", 
                 "#44AAAA", 
                 "#77CCCC", 
                 "#117744", 
                 "#44AA77", 
                 "#88CCAA", 
                 "#777711", 
                 "#AAAA44", 
                 "#DDDD77", 
                 "#774411", 
                 "#AA7744", 
                 "#DDAA77", 
                 "#771122", 
                 "#AA4455", 
                 "#DD7788")

phyla.col <- c("Alphaproteobacteria"= "#771155" ,
               "Acidimicrobiia"="#AA4488",
               "Subgroup 21" = "#DDAA77",
               "Actinobacteria"="#CC99BB", 
               "BD2-11 terrestrial group"= "#114477",
               "Oxyphotobacteria" = "#771122",
               "Nitrospinia" = "#AA4455",
               "Deltaproteobacteria"= "#4477AA", 
               "Bacteroidia"= "#77AADD", 
               "Gammaproteobacteria"="#117777", 
               "Epsilonproteobacteria" = "#FFFFFF",
               "Nitrospira" = "#DD7788",
               "Nitrososphaeria" = "#E69F00",
               "Marinimicrobia (SAR406 clade)"= "#34ABAA", 
               "Marinimicrobia (SAR406 clade)_unclassified"= "#34ABAA", 
               "Nitrospinia" = "#77CCCC",  
               "OM190"="#117744", 
               "Anaerolineae"= "#44AA77", 
               "Gemmatimonadetes" = "#F0E442",
               "Phycisphaerae"= "#88CCAA", 
               "Planctomycetacia"= "#777711", 
               "SAR202 clade"= "#AAAA44", 
               "Dehalococcoidia"= "#AAAA44",
               "JG30-KF-CM66" ="#DDDD77", 
               "Thermoanaerobaculia" = "#774411", 
               "Subgroup 22" = "#0072B2", 
               "Verrucomicrobiae" = "#AA7744",
               "Proteobacteria_unclassified" = "#CC79A7",
                "Bacteria_unclassified"= "white",
                "Subgroup 6"= "red",
               "Dadabacteriia" = "blue",
               "Other taxa" = "black")

reg_colours <-c("EGC"= "blue","WSC"= "red")

theme_plot <- theme(axis.title.x = element_text(size = 15), 
                    axis.title.y = element_text(size = 15),
                    axis.text.x = element_text(size = 12),
                    axis.text.y = element_text(size = 12),
                    strip.text = element_text(size = 15),
                    legend.text= element_text(size = 12),
                    plot.title = element_text(size=18, lineheight=.8, face="bold", hjust = 0.5),
                    panel.background = element_rect(fill = 'white', colour = 'white'),
                    panel.border = element_rect(color = "black", fill = NA, size = 1), 
                    strip.background = element_rect(color = "black", size = 1),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank())


taxa_col <- c("Alphaproteobacteria" = "#000000", 
               "Acidimicrobiia"="#E69F00", 
               "Bacteroidia" = "#56B4E9", 
               "Dehalococcoidia"= "#118921", 
               "Deltaproteobacteria"= "#F0E442", 
               "Gammaproteobacteria" = "#0072B2", 
               "Marinimicrobia (SAR406 clade)" = "#D55E00",
              "Marinimicrobia (SAR406 clade)_unclassified"= "#D55E00", 
               "Nitrospinia" = "#4477AA",
               "Planctomycetacia"= "#771155" ,
               "Verrucomicrobiae"="#AA4488",
               "Other"="gray50", 
               "Thiotrichales"= "#114477", 
               "Xanthomonadales"= "#4477AA", 
               "Vibrionales"= "#77AADD")

