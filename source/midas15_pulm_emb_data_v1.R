# |----------------------------------------------------------------------------------|
# | Project: Pulmonary Embolism                                                      |
# | Script: Make the data set for the analysis                                       |
# | Authors: Christopher Espana; Davit Sargsyan                                      |   
# | Created: 06/15/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
sink(file = "tmp/log_midas15_pulm_emb_data_v1.txt")
options(scipen = 999)
date()

# Load packages----
require(data.table)
require(knitr)
require(icd)
# require(foreach)
# require(parallel)
# require(doParallel)

# PART I: Load MIDAS----
# NOTE: data was preprocessed. See project 
# 'midas' R script 'export_midas_from_csv_to_rdata_v3.R'
system.time(load("E:/MIDAS/midas15_clean.RData"))
midas15

# Remove unused variables----
midas15[, ZIP := NULL]
midas15[, TOTBIL := NULL]
midas15[, STATUS := NULL]
midas15[, check := NULL]
gc()

# Number of patients
length(unique(midas15$Patient_ID))
# 18,057,028 records of 4,446,438 patients

# Exclude Emergency/Other Outpatient---- 
# 1 = inpatient
# 2 = ER outpatient
# 3 = same day surgery (SDS) outpatient
# 4 = other outpatient (non-ER and non-SDS)
# 5 = non-ER outpatient (3 or 4)
t1 <- table(midas15$YEAR,
            midas15$ADM_TYPE)
kable(format(t1, big.mark = ","))
  # |     |1       |2       |3       |4      |5       |
  # |:----|:-------|:-------|:-------|:------|:-------|
  # |1995 |445,427 |0       |0       |0      |51,019  |
  # |1996 |462,937 |0       |0       |0      |60,475  |
  # |1997 |469,139 |0       |0       |0      |73,844  |
  # |1998 |480,809 |0       |0       |0      |83,075  |
  # |1999 |487,342 |0       |0       |0      |86,042  |
  # |2000 |506,849 |0       |0       |0      |90,692  |
  # |2001 |512,121 |0       |0       |0      |94,236  |
  # |2002 |538,553 |0       |0       |0      |99,614  |
  # |2003 |567,841 |0       |0       |0      |104,101 |
  # |2004 |565,937 |0       |0       |0      |124,838 |
  # |2005 |577,161 |0       |0       |0      |127,869 |
  # |2006 |588,474 |0       |0       |0      |126,024 |
  # |2007 |588,841 |71      |0       |7      |129,072 |
  # |2008 |605,433 |367,436 |78,367  |79,698 |0       |
  # |2009 |605,830 |418,646 |95,288  |78,484 |0       |
  # |2010 |599,233 |456,183 |102,948 |83,270 |0       |
  # |2011 |574,068 |478,511 |105,419 |78,619 |0       |
  # |2012 |568,398 |532,918 |98,152  |84,347 |0       |
  # |2013 |550,449 |540,387 |97,425  |87,278 |0       |
  # |2014 |541,051 |562,339 |102,041 |92,963 |0       |
  # |2015 |530,759 |614,445 |106,091 |98,142 |0       |

# Exclude emergency records (missing before 2008)----
midas15 <- droplevels(subset(midas15,
                             ADM_TYPE != 3))
table(midas15$ADM_TYPE)
midas15[, ADM_TYPE := NULL]
gc()

# Part II: convert diagnoses codes to comorbid conditions----
# Load diagnoses codes----
l1 <- readxl::read_xls("data/Horizon_icd10_to_icd9_cardio_crosswalk_WJK_030918.xls")
l1 <- data.table(l1[!is.na(l1$MGRP),
                    c("MGRP",
                      "MAP_ICD9_CD",
                      "MAP_ICD10_CD")])
l1
length(unique(l1$MGRP))
# 26 comorbidities

# put ICD9 & 10 codes together----
tmp1 <- split(x = l1$MAP_ICD9_CD,
              f = l1$MGRP)
tmp2 <- split(x = l1$MAP_ICD10_CD,
              f = l1$MGRP)
l2 <- list()
for(i in 1:length(tmp1)) {
  l2[[i]] <- c(tmp1[[i]],
               tmp2[[i]])
}

# Remove leading numbers from the names of comorbidities----
names(l2)[1:7] <- substr(x = names(tmp1)[1:7], 
                         start = 2, 
                         stop = nchar(names(tmp1)[1:7]))
names(l2)[8:length(l2)] <- substr(x = names(tmp1)[8:length(tmp1)], 
                                  start = 4, 
                                  stop = nchar(names(tmp1)[8:length(tmp1)]))
data.table(names(l2))
l2 <- as.comorbidity_map(l2)
l2

# Separate diagnostic codes (DX1:DX9)
dx <- data.table(Record_ID = 1:nrow(midas15),
                 Patient_ID = midas15$Patient_ID,
                 midas15[, DX1:DX9])
dx

# cl <- makeCluster(getOption("cl.cores",
#                             detectCores() - 1))
# Option 1
# # On Windows: very slow, noparallelizaion, and huge drain on CPU. 
# # Deleting 'cl' clears he memory. Why doesn't it work?
# registerDoParallel(cl)
# 
# dtt <- foreach(i = 1:9) %dopar% {
#   tmp <- dx[, c(2, i + 2),
#             with = FALSE]
#   tmp[[2]] <- as.icd9(tmp[[2]])
#   comorbid(x = tmp,
#            map = l2,
#            visit_name = "Patient_ID",
#            icd_name = names(dx)[i + 2])
# }
# head(dtt[[1]])
# stopCluster(cl)
# gc()
# Option 2
# # Very slow due to data export to the cluster. How to fix?
# clusterExport(cl = cl,
#               varlist = list("l2",
#                              "dx"))
# dtt <- parLapply(cl = cl,
#                  X = 1:9,
#                  fun = function(i){
#                    require(icd)
#                    comorbid(x = dx,
#                             map = l2,
#                             visit_name = "Patient_ID",
#                             icd_name = names(dx)[i + 2])
#                  })
# Option 3
# Number of patients with each condition----
dtt <- list()
for(i in 1:9){
  dtt[[i]] <- icd9_comorbid(x = dx,
                            map = l2,
                            visit_name = "Patient_ID",
                            icd_name = names(dx)[i + 2])
}
head(dtt[[1]])


comorb <- data.table(Patient_ID = rownames(dtt[[1]]),
                     apply(Reduce("+", dtt),
                           MARGIN = 2,
                           function(a){
                             a > 0
                           }))
head(comorb)
length(unique(comorb$Patient_ID))
kable(format(data.frame(N_Patients = colSums(comorb[, -1])),
             big.mark = ","))
  # |                                      |N_Patients |
  # |:-------------------------------------|:----------|
  # |ATRIAL FIBRILLATION                   |731,976    |
  # |ATRIAL FLUTTER                        |125,400    |
  # |CVA                                   |254,282    |
  # |TIA                                   |184,907    |
  # |HYPERTENSION                          |3,205,789  |
  # |RENAL DISEASE                         |698,103    |
  # |HEART DISEASE                         |853,110    |
  # |RDIOMYOPATHY                          |284,022    |
  # |HYPERLIPIDEMIA                        |1,533,739  |
  # |CORONARY ARTERY DISEASE               |1,338,022  |
  # |PULMONARY EMBOLISM                    |101,592    |
  # |DVT                                   |155,698    |
  # |COPD                                  |730,473    |
  # |DIABETES                              |1,035,601  |
  # |CANCER                                |929,983    |
  # |INTRACRANIAL                          |102,070    |
  # |SPINAL HEMATOMA                       |523        |
  # |INTRAOCULAR VITREOUS HEMORRHAGE       |7,413      |
  # |HEMOPERITONEUM                        |6,278      |
  # |HEMARTHROSIS                          |2,557      |
  # |CARDIAC TAMPONADE                     |3,609      |
  # |INTRAMUSCULAR WITH COMPARTMENT SYNDRO |591        |
  # |GASTROINTESTINAL                      |151,202    |
  # |HEMATURIA                             |56,346     |
  # |NONTRAUMATIC HEMATOMA OF SOFT TISSUE  |3,411      |
  # |TRAUMA                                |1,009      |

# Exclude all cancer patients (~1M)----
dt1 <- subset(midas15,
              !(Patient_ID %in% comorb$Patient_ID[comorb$CANCER]))
comorb <- subset(comorb,
                 !CANCER)
comorb$CANCER <- NULL
length(unique(dt1$Patient_ID)) == nrow(comorb)
# TRUE

rm(midas15)
gc()

# Tables----
out <- list()
for (i in 2:ncol(comorb)) {
  tmp <- addmargins(table(comorb$`PULMONARY EMBOLISM`,
                          comorb[[i]]))
  out[[i - 1]] <- c(PE = tmp[2, 3],
                OTHER = tmp[3, 2],
                BOTH = tmp[2, 2])
  gc()
}
out <- do.call("rbind", out)
rownames(out) <- colnames(comorb)[-1]
out
write.csv(out,
          file = "tmp/pe_counts.csv")

# Save cases and controls, upload to the server
case <- subset(dt1,
               Patient_ID %in% comorb$Patient_ID[comorb$`ATRIAL FIBRILLATION`])
save(case,
     file = file.path("data/case.RData"), 
     compress = FALSE)

ctrl <- subset(dt1,
               Patient_ID %in% comorb$Patient_ID[!comorb$`ATRIAL FIBRILLATION`])
save(ctrl,
     file = file.path("data/ctrl.RData"), 
     compress = FALSE)

sink()
beepr::beep(3)