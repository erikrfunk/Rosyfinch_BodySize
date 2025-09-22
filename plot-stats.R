library(ggplot2)
library(ggforce)
library(NatParksPalettes)
library(nationalparkcolors)
library(reshape2)
#devtools::install_github("katiejolly/nationalparkcolors")

stats=read.csv("AllObservations_noLowGris_fixedCoords.txt",header = T,sep="\t")
taxon = c('Leucosticte tephrocotis griseonucha','Leucosticte tephrocotis umbrina','Leucosticte tephrocotis tephrocotis',
          'Leucosticte tephrocotis littoralis','Leucosticte tephrocotis wallowa','Leucosticte tephrocotis dawsoni',
          'Leucosticte atrata','Leucosticte australis')
stats = stats[stats$Species %in% taxon,]
stats$Species = factor(stats$Species,levels = taxon) # taxon set in "process-gbif.R"

# Some possible filters
stats_to_plot = stats # No filter
stats_to_plot = stats[stats$Month %in% c(5,6,7,8) | stats$Locality_binary=="Island",] # Summer obs  
stats_to_plot = stats[stats$Month %in% c(1,2,11,12) | stats$Locality_binary=="Island",] # Winter obs

# Some possible palettes
pal=park_palette("GeneralGrant")[8:1]
pal=c("#FFA630","#EBC775","#D7E8BA","#4DA1A9","#2E5077","#483656","#611C35","#6F3147")

# Boxplot or Sina plot
pdf("BodyMass_boxplot.pdf",height = 4,width=7)
ggplot(stats_noStGeorge,aes(x = group, y = Mass_g,color=merged_group)) +  
  geom_sina(outlier.shape = NA)+ 
  #geom_point(position = position_jitterdodge(), alpha=0.2) +
  #scale_x_discrete(limits = rev(levels(stats$Species)))+
  scale_color_manual(values=c("#1A3D82","#B14311")) +
  #scale_color_natparks_d("GeneralGrant")+
  guides (fill = guide_legend(ncol = 1))+
  #xlab("Population")+
  #ylab("Heterozygosity")+
  ggtitle("Body mass by taxon") +
  coord_flip()+
  theme_bw()
dev.off()

# Scatterplot 2 layers
pdf("Bio10_linear_model.pdf",height=4,width=6)
ggplot(stats_to_plot,aes(x=abs(lon),y=Mass_g, color=Species))+
  geom_point() +
  geom_smooth(method="lm",se=T,fullrange=F,alpha=0.25,color="black",linetype="dashed",linewidth=0.75)+
  #geom_point(data=island_melted,aes(x=value,y=Mass_g,color="Mainland.or.Island...0.1."))+
  #facet_wrap(facets=.~variable,scales="free")+
  #scale_color_manual(values = c("#1A3D82","#B14311"))+
  theme_bw()
dev.off()

# hist
plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")+
    theme_bw()
  plt + guides(fill=guide_legend(title=label_column))
}

plot_multi_histogram(stats_noStGeorge,"Mass_g","merged_group")

#Histogram by taxon
data_to_plot = cleaned_obs
pdf("Histogram_all_gbif_bio11_wSummerIslands.pdf.pdf",height=8,width=5)
ggplot(transform(data_to_plot,verbatimScientificName=factor(verbatimScientificName,levels=taxon)),
       aes(x=Mass_g))+
  geom_histogram()+
  facet_wrap(~verbatimScientificName,scales = "free_y",ncol = 1)+
  ggtitle("Histogram of Mean temp, coldest quarter (bio6)\n(May-Aug)")+
  theme_bw()
dev.off()
