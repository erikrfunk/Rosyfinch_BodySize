process_obs = function(obs,taxon,months=c(1:12)){
  obs_cleaned = obs[obs$verbatimScientificName %in% taxon,]
  obs_cleaned = obs_cleaned[obs_cleaned$month %in% months,]
  obs_cleaned$decimalLatitude = as.numeric(obs_cleaned$decimalLatitude)
  obs_cleaned$decimalLongitude = as.numeric(obs_cleaned$decimalLongitude)
  obs_cleaned = obs_cleaned[!is.na(obs_cleaned$decimalLatitude) & !is.na(obs_cleaned$decimalLongitude),]
  obs_cleaned = obs_cleaned[!(obs_cleaned$decimalLatitude==0) & !(obs_cleaned$decimalLongitude==0),]
  return(obs_cleaned)
}

obs = read.csv("Leucosticte_all_gbif.csv",header=T,sep="\t")
taxon = c('Leucosticte tephrocotis griseonucha','Leucosticte tephrocotis umbrina','Leucosticte tephrocotis tephrocotis',
          'Leucosticte tephrocotis littoralis','Leucosticte tephrocotis wallowa','Leucosticte tephrocotis dawsoni','Leucosticte atrata','Leucosticte australis')
months = c(11,12,1,2)
cleaned_obs = process_obs(obs,taxon=taxon,months=months)

#write.table(cleaned_obs,"Leucosticte_gbif_cleaned.csv",row.names = F, quote = F, sep="")
