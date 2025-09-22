
# Read in the three population Fst files, A is focal
FstAB = read.table("PBS/ASRF_outgroup_aus_only/site_unconstrained/isl_aus_FstSiteUnconstrained_25kbwin.fst", header=F)
FstAC = read.table("PBS/ASRF_outgroup_aus_only/site_unconstrained/isl_arc_FstSiteUnconstrained_25kbwin.fst", header=F)
FstBC = read.table("PBS/ASRF_outgroup_aus_only/site_unconstrained/aus_arc_FstSiteUnconstrained_25kbwin.fst", header=F)
output = "PBS/ASRF_outgroup_aus_only/site_unconstrained/"

#---- For Sites based analyses
# Transform Fst (Cavalli-Sforza 1969): T = -log(1- FST)
FstAB$Tfst = -log(1-FstAB$WEIR_AND_COCKERHAM_FST)
FstAC$Tfst = -log(1-FstAC$WEIR_AND_COCKERHAM_FST)
FstBC$Tfst = -log(1-FstBC$WEIR_AND_COCKERHAM_FST)

Ttmp = merge(FstAB,FstAC,by=c("CHROM","POS"))
TFSTs = merge(Ttmp,FstBC,by=c("CHROM","POS"))
Ts = TFSTs[,c(1:2,4,6,8)]
names(Ts)[3] = "Tab"
names(Ts)[4] = "Tac"
names(Ts)[5] = "Tbc"

# Calculate PBSa given PBSa = (Tab + Tac - Tbc)/2 or, for PBSn1 all three PBS
PBSa = (Ts$Tab + Ts$Tac - Ts$Tbc)/2
PBSb = (Ts$Tab + Ts$Tbc - Ts$Tac)/2
PBSc = (Ts$Tab + Ts$Tac - Ts$Tbc)/2

# Finally, calculate PBEa given PBEa = PBSa - (Tbc * PBSa_med)/Tbc_med
PBSa_med = median(na.omit(PBSa))
Tbc_med = median(na.omit(Ts$Tbc))
PBEa = PBSa - (Ts$Tbc * PBSa_med)/Tbc_med
PBEa_out = cbind(Ts[,1:3],PBEa)
PBEa_out$PBEa[PBEa_out$PBEa < 0] = 0 # Set negative PBS values to 0?
write.table(PBEa_out,output,row.names=F,quote=F)

# Or PBSn1, given PBSn1a = PBSa/(1+PBSa+PBSb+PBSc)
PBSn1a = PBSa/(1+PBSa+PBSb+PBSc)
PBSn1a_out = cbind(Ts[,1:3],PBSn1a)
write.table(PBSn1a_out,paste0(output,"PBSn1a_sites.txt"),row.names=F,quote=F)
quantile(PBSn1a_out$PBSn1a,probs=0.995)

###############################
#---- For window based analyses
# Transform Fst (Cavalli-Sforza 1969): T = -log(1- FST)
FstAB$Tfst = -log(1-FstAB$WEIGHTED_FST)
FstAC$Tfst = -log(1-FstAC$WEIGHTED_FST)
FstBC$Tfst = -log(1-FstBC$WEIGHTED_FST)

Ttmp = merge(FstAB,FstAC,by=c("CHROM","BIN_START"))
TFSTs = merge(Ttmp,FstBC,by=c("CHROM","BIN_START"))
Ts = TFSTs[,c(1:3,7,12,17)]
names(Ts)[4] = "Tab"
names(Ts)[5] = "Tac"
names(Ts)[6] = "Tbc"

# Calculate PBSa given PBSa = (Tab + Tac - Tbc)/2 or, for PBSn1 all three PBS
PBSa = (Ts$Tab + Ts$Tac - Ts$Tbc)/2
PBSb = (Ts$Tab + Ts$Tbc - Ts$Tac)/2
PBSc = (Ts$Tbc + Ts$Tac - Ts$Tab)/2

# Finally, calculate PBEa given PBEa = PBSa - (Tbc * PBSa_med)/Tbc_med
PBSa_med = median(na.omit(PBSa))
Tbc_med = median(na.omit(Ts$Tbc))
PBEa = PBSa - ((Ts$Tbc * PBSa_med)/Tbc_med)

PBEa_out = cbind(Ts[,1:3],PBEa)
PBEa_out$PBEa[PBEa_out$PBEa < 0] = 0 # Set negative PBS values to 0?
write.table(PBEa_out,output,row.names=F,quote=F)
quantile(PBEa_out$PBEa,probs=0.995)

# Or PBSn1, given PBSn1a = PBSa/(1+PBSa+PBSb+PBSc)
PBSn1a = PBSa/(1+PBSa+PBSb+PBSc)
PBSn1a_out = cbind(Ts[,1:3],PBSn1a)
write.table(PBSn1a_out,paste0(output,"_PBSn1a"),row.names=F,quote=F)
quantile(PBSn1a_out$PBSn1a,probs=0.995)

##########################################
# Compare the two sets of PBE calculations and pull out 95 percentile windows overlapping both gris and umb
PBEgris = read.table("PBS/pbs_gris.txt",header=T)
PBEumb = read.table("PBS/pbs_umbrina.txt",header=T)
gris_top = quantile(PBEgris$PBEa,probs=0.95)
umb_top = quantile(PBEumb$PBEa,probs=0.95)

PBEgrisTop = PBEgris[PBEgris$PBEa >= gris_top,]
PBEumbTop = PBEumb[PBEumb$PBEa >= umb_top,]

PBEcommon = merge(PBEgrisTop,PBEumbTop,by=c("CHROM","BIN_START"))
names(PBEcommon)[4] = "PBEgris"
names(PBEcommon)[6] = "PBEumb"
write.table(PBEcommon,"PBS/PBS95percentile_GrisUmbOutliers.txt",quote=F,row.names=F)

# Or get windows overlapping another windowed result
pbs_res = read.table("PBS/ASRF_outgroup_aus_only/site_unconstrained/PBE_ausOnly_SiteUnconstrained_25kbwin.txt",header=T)
gemma_res = read.table("Fst/top_window_overlaps_withGenes.txt",header=T)
names(pbs_res) = c("Chr","WinStart","WinStop","PBEa")
gemma_res$Chr = gsub(";HRSCAF=[0-9]*","",gemma_res$Chr)
pbs_res$Chr = gsub(";HRSCAF=[0-9]*","",pbs_res$Chr)

gemma_res$overlap_res = rep(0,nrow(gemma_res))
for(i in 1:nrow(pbs_res)){
  for(j in 1:nrow(gemma_res)){
    if(pbs_res$Chr[i] == gemma_res$Chr[j] & (max(pbs_res$WinStart[i],gemma_res$WinStart[j]) < min(pbs_res$WinStop[i],gemma_res$WinStop[j]))){
      gemma_res$overlap_res[j] = 1
    }
  }
}
for(i in 1:nrow(gemma_res)){
  if(gemma_res$meanP[i] < 0.304){
  gemma_res$overlap_res[i] = gemma_res$overlap_res[i] + 1
  }
}

# Do a slow loop to pull individual sites that fall within Fst outlier windows
fst_tops = read.table("Fst/top_window_overlaps_withGenes.txt",header=T)
pbs_sites = read.table("PBS/ASRF_outgroup_aus_only/PBE_sites.txt",header=T)
head(pbs_sites)
names(pbs_sites)[4] = "Stat"
pbs_sites$CHROM = as.character(pbs_sites$CHROM)
quantile(pbs_sites$Stat,0.999,na.rm=T)
sig_cutoff = 2.17


sig_sites = c()
window_hits=0
for(i in 1:nrow(fst_tops)){
  win = fst_tops[i,]
  overlap = na.omit(pbs_sites[pbs_sites$CHROM == win$Chromosome & pbs_sites$POS > win$WinStart & pbs_sites$POS < win$WinStop & pbs_sites$Stat >= sig_cutoff,])
  if(nrow(overlap)!=0 & !(is.na(overlap[1,2]))){
    window_hits = window_hits + 1
    overlap$Gene = win$Gene
    sig_sites = rbind(sig_sites,overlap)
  }
}

#write.table(sig_sites,"PBS_sigSites99.9th_FstSigWin_overlap.txt",quote=F,row.names=F)

#---- Maximize either Fst or PBS per window using either the site-constrained or unconstrained approach following Shpak et al. 2025
# Use a windowed Fst file just to set the window coords
fst_wins = read.table("PBS/ASRF_outgroup/ASRF_Island_fst_50kbwins.txt.windowed.weir.fst",header=T)
head(fst_wins)
fst_wins$CHROM = gsub(";HRSCAF=[0-9]*","",fst_wins$CHROM)
coords = fst_wins[,1:3]
# Read in the desired sites-based file, e.g. fst or pbs and set the index column 
sites = read.table("PBS/ASRF_outgroup/sites/PerSite_PBS_Ts.txt",header=T)
sites$CHROM = as.factor(sites$CHROM)
sites_split = split(sites,sites$CHROM)
idx = 4
# Set outdir
outdir = "~/finches/rosyfinches/BodySize_proj/PBS/ASRF_outgroup/sites/"
site_constrained = c()
for(k in 1:length(sites_split)){
  if(names(sites_split)[k] %in% coords$CHROM){
    scaff_df = na.omit(sites_split[[k]])
    scaff_coords = coords[coords$CHROM==names(sites_split)[k],]
    scaff_coords$Stat = NA
    for(i in 1:nrow(scaff_coords)){
      win_max = max(scaff_df[scaff_df$POS>scaff_coords[i,2] & scaff_df$POS<scaff_coords[i,3],idx])
      scaff_coords$Stat[i] = win_max
    }
    write.table(scaff_coords,paste0(outdir,"aus_arc_FstSiteUnconstrained_25kbwin.fst"),row.names=F,quote=F,append=T,col.names = F) # only include the header once
  }
}
pbsstats = read.table("PBS/ASRF_outgroup/sites/PBE_persite_noInf_winMean.txt",header=F)
names(pbsstats) = c("Chromosome","WinStart","WinStop","Stat")
pbsstats$Chromosome = gsub(";HRSCAF=[0-9]*","",pbsstats$Chromosome)
pbsstats = na.omit(pbsstats)
dim(pbsstats)
pbsstats_sigwins = pbsstats[pbsstats$Stat >= quantile(pbsstats$Stat,0.95),]
# merge with top fst windows
fst_tops = read.table("Fst/top_window_overlaps_withGenes.txt",header=T)
overlap = merge(fst_tops,pbsstats_sigwins,by=c("Chromosome","WinStart","WinStop"))

# Instead, first maximize fst within a window for each AB AC and BC comparisons unconstrained by site.
# Then use this unconstrained maximum to calculate a PBE statistic.