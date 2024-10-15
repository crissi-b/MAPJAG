

library(readr)

tx_08 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2518-08__20231107__215839/transcripts.csv")
tx_11 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2518-11__20231107__215839/transcripts.csv")
tx_04 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data/wei/Xenium/Xenium_Chris/20231107__215641__110723_lungv46_synovium/output-XETG00150__0013726__2519-04__20231107__215839/transcripts.csv")

files <- list(tx_04, tx_08, tx_11)
donor <- c("D2518_08", "D2518_11", "D2519_04")

for (i in 1:length(donor)){
files[[i]]$donor <- donor[[i]]
files[[i]]$donor_FOV <- paste(files[[i]]$donor, files[[i]]$fov_name, sep="_")
  }
names(files) <- donor

library(data.table)
all_tx <- rbindlist(files)
head(all_tx)

all_tx_f <- all_tx %>% filter(feature_name != "BLANK" & feature_name != "NegControl")

all_tx_f <- all_tx %>% filter(!grepl('BLANK|NegControl', feature_name) & cell_id != "UNASSIGNED")

all_tx_f$cell_id_new <- as.numeric(as.factor(all_tx_f$cell_id))

length(unique(all_tx_f$donor_FOV))

sum(length(unique(tx_04$fov_name)),length(unique(tx_08$fov_name)),length(unique(tx_11$fov_name)))

all_fovs_split <- list()
for (i in 1:length(unique(all_tx_f$donor_FOV))){
  all_fovs_split[[i]] <- filter(all_tx_f, donor_FOV== unique(all_tx_f$donor_FOV)[[i]])
}
names(all_fovs_split) <- unique(all_tx_f$donor_FOV)

out_dirs <- paste("/rds/projects/c/croftap-celldive01/xenium/analysis/tx_counts_FOV/", unique(all_tx_f$donor_FOV), ".csv", sep="")


for (i in 1:length(unique(all_tx_f$donor_FOV))){
  write.csv(all_fovs_split[[i]], out_dirs[[i]])
}


#batch 2
library(readr)

tx_18_24 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data_batch2/Xenium_data_for_Chris_Mahoney/output-XETG00150__0018573__2518-24__20240227__000410/transcripts.csv.gz")
tx_19_15 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/data_batch2/Xenium_data_for_Chris_Mahoney/output-XETG00150__0018573__2519-15__20240227__000410/transcripts.csv.gz")


files2 <- list(tx_18_24, tx_19_15)
donor <- c("D2518_24", "D2519_15")

for (i in 1:length(donor)){
files2[[i]]$donor <- donor[[i]]
files2[[i]]$donor_FOV <- paste(files2[[i]]$donor, files2[[i]]$fov_name, sep="_")
  }
names(files2) <- donor

library(data.table)
all_tx_b2 <- rbindlist(files2)
head(all_tx_b2)

all_tx_b2_f <- all_tx_b2 %>% filter(!grepl('BLANK|NegControl|UnassignedCodeword_', feature_name) & cell_id != "UNASSIGNED")

head(all_tx_b2_f)

all_tx_b2_f$cell_id_new <- as.numeric(as.factor(all_tx_b2_f$cell_id))

length(unique(all_tx_b2_f$donor_FOV))

sum(length(unique(tx_18_24$fov_name)),length(unique(tx_19_15$fov_name)))

all_fovs_split_b2 <- list()
for (i in 1:length(unique(all_tx_b2_f$donor_FOV))){
  all_fovs_split_b2[[i]] <- filter(all_tx_b2_f, donor_FOV== unique(all_tx_b2_f$donor_FOV)[[i]])
}
names(all_fovs_split_b2) <- unique(all_tx_b2_f$donor_FOV)

out_dirs <- paste("/rds/projects/c/croftap-celldive01/xenium/analysis/tx_counts_FOV/", unique(all_tx_b2_f$donor_FOV), ".csv", sep="")


for (i in 1:length(unique(all_tx_b2_f$donor_FOV))){
  write.csv(all_fovs_split_b2[[i]], out_dirs[[i]])
}

#batch 3
library(readr)

tx_18_28 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014095__2518-28__20240331__210931/transcripts.csv.gz")

tx_18_36 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014095__2518-36__20240331__210931/transcripts.csv.gz")
  
tx_18_40 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014097__2518-40__20240331__210931/transcripts.csv.gz")
  
tx_19_18 <- read_csv("/rds/projects/c/croftap-celldive01/xenium/batch3/20240331__210418__BWH_20240331_RA_Kevin_Adam/output-XETG00150__0014097__2519-18__20240331__210931/transcripts.csv.gz")


files3 <- list(tx_18_28, tx_18_36, tx_18_40, tx_19_18)
donor <- c("D2518_28", "D2518_36", "D2518_40", "D2519_18")

for (i in 1:length(donor)){
files3[[i]]$donor <- donor[[i]]
files3[[i]]$donor_FOV <- paste(files3[[i]]$donor, files3[[i]]$fov_name, sep="_")
  }
names(files3) <- donor

library(data.table)
all_tx_b3 <- rbindlist(files3)
head(all_tx_b3)

all_tx_b3_f <- all_tx_b3 %>% filter(!grepl('BLANK|NegControl|UnassignedCodeword_', feature_name) & cell_id != "UNASSIGNED")

all_tx_b3_f$cell_id_new <- as.numeric(as.factor(all_tx_b3_f$cell_id))

length(unique(all_tx_b3_f$donor_FOV))

sum(length(unique(tx_18_28$fov_name)),length(unique(tx_18_36$fov_name)), length(unique(tx_18_40$fov_name)),length(unique(tx_19_18$fov_name)))

all_fovs_split_b3 <- list()
for (i in 1:length(unique(all_tx_b3_f$donor_FOV))){
  all_fovs_split_b3[[i]] <- filter(all_tx_b3_f, donor_FOV== unique(all_tx_b3_f$donor_FOV)[[i]])
}
names(all_fovs_split_b3) <- unique(all_tx_b3_f$donor_FOV)

out_dirs <- paste("/rds/projects/c/croftap-celldive01/xenium/analysis/tx_counts_FOV/", unique(all_tx_b3_f$donor_FOV), ".csv", sep="")


for (i in 1:length(unique(all_tx_b3_f$donor_FOV))){
  write.csv(all_fovs_split_b3[[i]], out_dirs[[i]])
}



