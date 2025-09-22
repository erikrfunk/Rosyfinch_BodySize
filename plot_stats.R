library(ggplot2)
library(NatParksPalettes)

master_fn="SampleList_master.csv"
df=read.csv(master_fn,header=T)
head(df)

# A few filter options
df_to_plot = df[which(!(is.na(df$Mass_g))),]
df_to_plot = df[df$Time=="Contemporary" | df$Time=="Historical" & df$Depth>2,]
df_to_plot = df[df$Time=="Historical" & df$Depth>2,]
df_to_plot = df[df$Time=="Contemporary",]
df_to_plot$Species = factor(df_to_plot$Species, levels = c(
  "Leucosticte tephrocotis griseonucha",
  "Leucosticte tephrocotis umbrina",
  "Leucosticte tephrocotis littoralis",
  "Leucosticte tephrocotis tephrocotis",
  "Leucosticte atrata",
  "Leucosticte australis",
  "Leucosticte tephrocotis dawsoni"))
names(df_to_plot)

# Boxplot
ggplot(df_to_plot,aes(x = Species, y = Mass_g,fill=Species)) +  
  geom_boxplot(outlier.shape = NA)+ 
  geom_point(position = position_jitterdodge(), alpha=0.2) +
  scale_fill_natparks_d("Triglav") +
  guides (fill = guide_legend(ncol = 1))+
  #xlab("Population")+
  ylab("Mass (g)")+
  #ggtitle("Change in He after droppinig transitions (and sites)") +
  theme_bw()

# Scatterplot
ggplot(df_to_plot,aes(x=He_sites_noTrans,y=He_NoTrans,color=Pop,shape=Time))+
  geom_point(size=2) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, alpha=0.15)+
  scale_color_natparks_d("Yellowstone")+
  scale_shape_manual(values=c(17,19))+
  ggtitle("He versus He and no transitions")+
  theme_bw()
