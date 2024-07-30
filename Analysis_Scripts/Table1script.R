#creating table one for analysis
#here we will be creating table one to describ the different cohorts
#that we are using in this study
library(tableone)
library(tidyverse)

#pick out which samples are used for the modelling in the paper in the final set
final_samples_metadata<- read.csv("./pdac_meta.csv")
pdac_metadata_tableone <- final_samples_metadata

#pick which variables we care about for analysis
select_variables <- c("Age.at.Collection", "Race", "Gender.Flag", "Group","Source", 
               "Pre.existing.Diabetes", "Stage", "Grade")
categorical_variables <- c("Race", "Gender.Flag", "Pre.existing.Diabetes","Stage","Group","Grade")
pdac_metadata_tableone <- pdac_metadata_tableone[,select_variables]

replace_str_na <- function(x){
  ifelse(x=="n/a", NA, x)
}
#formatting two of the diagnoses writing types
pdac_metadata_tableone$Diagnosis <- pdac_metadata_tableone$Group %>%str_replace("Benign pancreas", "Benign Pancreas")
pdac_metadata_tableone$Diagnosis <- pdac_metadata_tableone$Group %>%str_replace("Cancer_other", "Cancer Other")
pdac_metadata_tableone$Diagnosis <- pdac_metadata_tableone$Diagnosis %>%str_replace("Islet Cell Tumor", "Cancer Other")

#format all of the variables
pdac_metadata_tableone$Pre.existing.Diabetes  <- replace_str_na(pdac_metadata_tableone$Pre.existing.Diabetes)
pdac_metadata_tableone$Stage  <- replace_str_na(pdac_metadata_tableone$Stage)
pdac_metadata_tableone$Stage <- ifelse(is.na(pdac_metadata_tableone$Stage), "no_info", pdac_metadata_tableone$Stage)
pdac_metadata_tableone$PDAC_Stage <- ifelse(pdac_metadata_tableone$Group !="PDAC", "not_PDAC", pdac_metadata_tableone$Stage)
pdac_metadata_tableone$Other_Stage <- ifelse(pdac_metadata_tableone$Diagnosis !="Cancer Other","not Cancer_Other", pdac_metadata_tableone$Stage)
pdac_metadata_tableone$Grade  <- replace_str_na(pdac_metadata_tableone$Grade)

#collapsing race into only two categories, leaving this out will give you all categories as reported
pdac_metadata_tableone$Race<- ifelse(pdac_metadata_tableone$Race =="white",
                                      "White and Non Hispanic",
                                      "Hispanic or non-White") %>% 
                                      as.factor()
#we will be reporting sex at birth
pdac_metadata_tableone$Gender.Flag <- pdac_metadata_tableone$Gender.Flag %>%str_replace("M\\/F @ birth", "F")


variables <- c("Age.at.Collection", "Race", "Gender.Flag")
categorical_variables <- c("Race", "Gender.Flag")

#now we are going to make the tableone in two parts
#the first part will be the variables that are not the stage variables
#the stage variables will be done separately so that they can be created on a 
#per cohort level and pasted together with the subgroup analyzed variables
pdac_metadata_tableone$Source <- ifelse(pdac_metadata_tableone$Source=="CEDAR_2020", "CEDAR", "BCC")


source_comparison <- CreateTableOne(vars = variables, strata = c("Source", "Diagnosis"),
                                    data = pdac_metadata_tableone,
                                    factorVars = categorical_variables)
source_compare_printed <-print(source_comparison, 
                               quote = FALSE, noSpaces = TRUE,
                               format="f", minMax = T,showAllLevels = T,
                               nonnormal=c("Age.at.Collection") ,
                               printToggle = FALSE)

pdac_metadata_tableone_CEDAR = pdac_metadata_tableone %>% filter(Source=="CEDAR")
pdac_metadata_tableone_BCC = pdac_metadata_tableone %>% filter(Source!="CEDAR")
source_comparison_CEDAR <- CreateTableOne(vars = variables, strata = c("Diagnosis"),
                                    data = pdac_metadata_tableone_CEDAR,
                                    factorVars = categorical_variables)
source_compare_printed_CEDAR <-print(source_comparison_CEDAR, 
                               quote = FALSE, noSpaces = TRUE,
                               format="f", minMax = T,showAllLevels = T,
                               nonnormal=c("Age.at.Collection"),
                               printToggle = FALSE)
#
source_comparison_BCC <- CreateTableOne(vars = variables, strata = c("Diagnosis"),
                                    data = pdac_metadata_tableone_BCC,
                                    factorVars = categorical_variables)
source_compare_printed_BCC <-print(source_comparison_BCC, 
                               quote = FALSE, noSpaces = TRUE,
                               format="f", minMax = T,showAllLevels = T,
                               nonnormal=c("Age.at.Collection"),
                               printToggle = FALSE)
#



just_stage <- CreateTableOne(vars = c("PDAC_Stage", "Other_Stage"), strata = c("Source"),
                                    data = pdac_metadata_tableone,
                                    factorVars = c("Stage"))
just_stage_printed<-print(just_stage, 
                               quote = FALSE, noSpaces = TRUE,minMax=T,
                               showAllLevels = T, nonnormal = c("Stage"),
                               printToggle = FALSE)

#specity where you want the tableone information to be written
#write.table(source_compare_printed, "./pdac_all_group_table1.txt")
#write.table(just_stage_printed,  "./pdac_stage_table1.txt")
color_palette_base <- c('#F6C142','#C2D7EC','#4075B1', '#4075B1', '#DF8244','#68379a', '#4eac5b')
names(color_palette_base) <- c("IPMN", "Acute pancreatitis","Chronic pancreatitis", "Pancreatitis",
                               "Cancer Other", "PDAC", "Benign Pancreas")
color_palette_ordered <- color_palette_base[c("PDAC", "Cancer Other", "IPMN",
                                              "Pancreatitis","Benign Pancreas")]

#make pie chart plot
BCC_patient_numbers <- pdac_metadata_tableone %>% filter(Source=="BCC")
bcc_pie_numbers <- BCC_patient_numbers$Diagnosis %>% table() 
color_palette_bcc <- color_palette_base[names(bcc_pie_numbers)]
pie(as.vector(bcc_pie_numbers), labels =as.vector(bcc_pie_numbers),
    main="BCC", col=color_palette_bcc, cex=1.5) 
legend("right", names(bcc_pie_numbers), cex=0.7, inset=c(-0.2,0),
       fill=color_palette_bcc,bty="n")

CEDAR_patient_numbers <- pdac_metadata_tableone %>% filter(Source=="CEDAR")
cedar_pie_numbers <- CEDAR_patient_numbers$Diagnosis %>% table()
color_palette_cedar <- color_palette_base[names(cedar_pie_numbers)]
pie(as.vector(cedar_pie_numbers),labels = as.vector(cedar_pie_numbers),
    main="CEDAR", col =color_palette_cedar, cex=1.5) 
#legend("right", names(cedar_pie_numbers), cex=0.7, inset=c(-0.2,0), fill=color_palette_cedar, bty="n")



