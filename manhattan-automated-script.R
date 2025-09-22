library("qqman")

#set list of file names 
fns = list.files("Fst/umbrina/",pattern = "*.windowed.weir.fst",full.names = T)

# Iterate through the files to first get 99.5th percentile, then plot the 
# manhattan with outlier windows highlighted
chr_order<-read.table("~/finches/rosyfinches/scafflist_zfinch_matches",header =F, stringsAsFactors = F)

for(i in 1:length(fns)){
  results<-read.table(fns[i],header = T,stringsAsFactors = F)
  cutoff = quantile(results$WEIGHTED_FST,0.995)
  index_matches<-match(results[,1],chr_order$V2)
  results<-cbind(index_matches,results)
  results<-results[!(is.na(results$index_matches)),]
  f_out = paste0(fns[i],".manhattan.png")
  results$win = seq(1:nrow(results))
  sig_wins = results$win[results$WEIGHTED_FST > cutoff]
  
  png(file=f_out,height=4,width=10,units = "in",res=600)
  manhattan_orange(results, chr = "index_matches", bp="BIN_START",snp = "win",
                   p="WEIGHTED_FST",logp=F,ylab="Mean Fst",main=paste0(fns[i],": 99.5th=",cutoff),
                   highlight = sig_wins)
  dev.off()
}
