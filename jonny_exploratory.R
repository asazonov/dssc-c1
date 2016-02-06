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
metrics = read.csv(sep = ',', file = paste0(data_location, "Metrics.rpt")) #dilyana is fixing this
patients = read.csv(sep = ',', file = paste0(data_location, "Patients.rpt"))
habits = read.csv(sep = ',', file = paste0(data_location, "Personal_Habits.rpt"))

sort(unique(patients$Postal_Area))
dim(patients)
