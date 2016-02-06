data_location = "~/Dropbox/Work/hackathon/nuffield/"
#####
#IoT data
temp = read.csv(paste0(data_location, "IoT/iot_bodytemp_data.csv"))
hr = read.csv(paste0(data_location, "IoT/iot_heartrate_data.csv"))
hr$resting_difference = hr$User_Resting_Heart_Rate - hr$Age_Group_Heart_Rate
sleep = read.csv(paste0(data_location, "IoT/iot_sleepquality_data.csv"))

summary(lm(sleep$Sleep_Quality_Value ~ hr$resting_difference))

#####
#Nuffield
history = read.csv(sep = ',', file = paste0(data_location, "Health_History.rpt"))
metrics = read.csv(sep = ',', file = paste0(data_location, "Metrics.rpt"))
patients = read.csv(sep = ',', file = paste0(data_location, "Patients.rpt"))
habits = read.csv(sep = ',', file = paste0(data_location, "Personal_Habits.rpt"))

#get london postcodes
postcodes = read.csv("~/Dropbox/Work/hackathon/London postcodes.csv")
names(postcodes)
post = as.character(postcodes$Postcode)
post.london = unique(unlist(lapply(strsplit(post," ",fixed=TRUE),"[",1)))
write.csv(post.london, file = "london_postcodes.csv")

names(patients)
sum(as.character(patients$Postal_Area)%in%post.london)

history_columns = c(4:64, 66:85,253,254,259,260)
history_disease = c(4:64, 66:85)
history.cor = history[, history_disease]
for(i in 1:dim(history.cor)[2]){
  history.cor[,i] = as.character(history.cor[,i])
  history.cor[history.cor[i]=='NULL' , i] = NA
  history.cor[history.cor[i]=='No' , i] = F
  history.cor[history.cor[i]=='Yes' , i] = T
  }


correlation = cor(history.cor)
