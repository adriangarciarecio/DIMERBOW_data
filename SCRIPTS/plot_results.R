# library
library(tidyverse)
library(dplyr)
library(reshape2)
library("gridExtra")
library(gtools)
library(xlsx)

##INTERACTIONS#########################################################################################
dimers <- list.files(path = '/home/adrian/Documents/bender/PDB/dimers/', pattern = ".pdb")

#COUNTS
l_count_pdb <- c()
l_count_lz <- c()
l_count_ini_1 <- c()
l_count_end_1 <- c()
l_count_ini_2 <- c()
l_count_end_2 <- c()
l_count_ini_3 <- c()
l_count_end_3 <- c()
l_count_ini_4 <- c()
l_count_end_4 <- c()
l_count_ini_5 <- c()
l_count_end_5 <- c()

#COUNTING TOTAL
data_pdb <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/count_inter.txt', skip = 1)
data_lz <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/count_inter.txt', skip = 1)
data_ini_1 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/count_inter_ini.txt', skip = 1)
data_end_1 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/count_inter_end.txt', skip = 1)
data_ini_2 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/count_inter_ini.txt', skip = 1)
data_end_2 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/count_inter_end.txt', skip = 1)
data_ini_3 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/count_inter_ini.txt', skip = 1)
data_end_3 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/count_inter_end.txt', skip = 1)
data_ini_4 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/count_inter_ini.txt', skip = 1)
data_end_4 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/count_inter_end.txt', skip = 1)
data_ini_5 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/count_inter_ini.txt', skip = 1)
data_end_5 <- read.table('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/count_inter_end.txt', skip = 1)

rownames(data_pdb) <- (apply(data_pdb[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_lz) <- (apply(data_lz[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_ini_1) <- (apply(data_ini_1[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_end_1) <- (apply(data_end_1[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_ini_2) <- (apply(data_ini_2[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_end_2) <- (apply(data_end_2[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_ini_3) <- (apply(data_ini_3[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_end_3) <- (apply(data_end_3[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_ini_4) <- (apply(data_ini_4[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_end_4) <- (apply(data_end_4[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_ini_5) <- (apply(data_ini_5[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
rownames(data_end_5) <- (apply(data_end_5[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))

#MERGE DATA 
data_gro_1 <- merge.data.frame(data_ini_1, data_end_1, by="row.names")
rownames(data_gro_1) <- data_gro_1$Row.names
data_gro_1$Row.names <- NULL
data_gro_2 <- merge.data.frame(data_ini_2, data_end_2, by="row.names")
rownames(data_gro_2) <- data_gro_2$Row.names
data_gro_2$Row.names <- NULL
data_gro_3 <- merge.data.frame(data_ini_3, data_end_3, by="row.names")
rownames(data_gro_3) <- data_gro_3$Row.names
data_gro_3$Row.names <- NULL
data_gro_4 <- merge.data.frame(data_ini_4, data_end_4, by="row.names")
rownames(data_gro_4) <- data_gro_4$Row.names
data_gro_4$Row.names <- NULL
data_gro_5 <- merge.data.frame(data_ini_5, data_end_5, by="row.names")
rownames(data_gro_5) <- data_gro_5$Row.names
data_gro_5$Row.names <- NULL

data_pdb_t <- data_pdb[,3:4]
colnames(data_pdb_t) <- c("count_pdb", "%")
data_lz_t <- data_lz[,3:4]
colnames(data_lz_t) <- c("count_lz", "%")

data_gro_1 <- data_gro_1[,c(3,7)] 
data_gro_2 <- data_gro_2[,c(3,7)] 
data_gro_3 <- data_gro_3[,c(3,7)] 
data_gro_4 <- data_gro_4[,c(3,7)] 
data_gro_5 <- data_gro_5[,c(3,7)] 

data_gro <- merge.data.frame(data_gro_1, data_gro_2, by="row.names")
rownames(data_gro) <- data_gro$Row.names
data_gro$Row.names <- NULL
colnames(data_gro) <- c("1.1", "1.2", "2.1", "2.2")

data_gro <- merge.data.frame(data_gro, data_gro_3, by="row.names")
rownames(data_gro) <- data_gro$Row.names
data_gro$Row.names <- NULL
colnames(data_gro) <- c("1.1", "1.2", "2.1", "2.2", "3.1", "3.2")

data_gro <- merge.data.frame(data_gro, data_gro_4, by="row.names")
rownames(data_gro) <- data_gro$Row.names
data_gro$Row.names <- NULL
colnames(data_gro) <- c("1.1", "1.2", "2.1", "2.2", "3.1", "3.2", "4.1", "4.2")

data_gro <- merge.data.frame(data_gro, data_gro_5, by="row.names")
rownames(data_gro) <- data_gro$Row.names
data_gro$Row.names <- NULL
colnames(data_gro) <- c("1.1", "1.2", "2.1", "2.2", "3.1", "3.2", "4.1", "4.2", "5.1", "5.2")

data_gro$mean_ini <- apply(data_gro[,c(1,3,5,7,9)], 1, function(x) mean(x))
data_gro$sd_ini <- apply(data_gro[,c(1,3,5,7,9)], 1, function(x) sd(x))
data_gro$mean_end <- apply(data_gro[,c(2,4,6,8,10)], 1, function(x) mean(x))
data_gro$sd_end <- apply(data_gro[,c(2,4,6,8,10)], 1, function(x) sd(x))

data_gro <- data_gro[,c(11,12,13,14)]

data_all <- merge.data.frame(data_pdb_t, data_lz_t, by="row.names")
rownames(data_all) <- data_all$Row.names
data_all$Row.names <- NULL
data_all <- merge.data.frame(data_all, data_gro, by="row.names")
rownames(data_all) <- data_all$Row.names
data_all$Row.names <- NULL

data_all <- data_all[,c(1,3,5,6,7,8)]

data_all_ord <- data_all[rev(order(data_all$mean_end)),]

#FILTERING
i = 0
for (d in data_all_ord[,dim(data_all_ord)[2]]){
  if (d == 0 & i != 0) {
    data_all_ord_f <- data_all_ord[1:i,c(1,2,dim(data_all_ord)[2]-3,dim(data_all_ord)[2]-2,dim(data_all_ord)[2]-1, dim(data_all_ord)[2])]
    break
  }
  if (d == 0 & i == 0) {
    print ('NO INTERACTIONS!')
  }
  i = i+1
}

data_all_ord_f$interaction <- rownames(data_all_ord_f)
l_group <- c()
for (inter in data_all_ord_f$interaction){
  if (grepl('TM.-', inter)) {
    l_group <- c(l_group, 'TM')
  }
  if (grepl('ECL.-', inter)) {
    l_group <- c(l_group, 'ECL')
  }
  if (grepl('ICL.-', inter)) {
    l_group <- c(l_group, 'ICL')
  }
  if (grepl('CT-', inter)) {
    l_group <- c(l_group, 'ICL')
  }
  if (grepl('NT-', inter)) {
    l_group <- c(l_group, 'ECL')
  }
}

data_all_ord_f$group <- l_group
data_all_ord_f$group <- as.factor(data_all_ord_f$group)

sysm <- read.table(file = "/home/adrian/Documents/bender/RESULTS/symmetry.txt", header = T, sep = '\t')
symmetry <- sysm$SYMMETRY

#COUNTING HH - HT
dimers_sys <- sysm$DIMER[sysm$SYMMETRY == 'YES']
dimers_asys <- sysm$DIMER[sysm$SYMMETRY == 'NO']

matrix_sys <- matrix(0, length(dimers_sys), length(inter_end))
colnames(matrix_sys) <- inter_end
rownames(matrix_sys) <- dimers_sys

l_inter_sys_inter <- vector(mode = "list", length = 1)
l_inter_sys_values <- vector(mode = "list", length = 1)

for (j in 1:length(dimers_sys)){
  dimer <- as.character(dimers_sys[j])
  print(dimer)
  files_pdb <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/', pattern = dimer)
  data_pdb <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/', files_pdb, sep = ""), skip = 1)
  info_pdb <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/', files_pdb, sep = ""))
  info_pdb <- strsplit(info_pdb, "\n")
  count_pdb <- info_pdb[[1]][1]
  l_count_pdb <- c(l_count_pdb, as.integer(count_pdb))
  files_lz <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/', pattern = dimer)
  data_lz <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/', files_lz, sep = ""), skip = 1)
  info_lz <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/', files_lz, sep = ""))
  info_lz <- strsplit(info_lz, "\n")
  count_lz <- info_lz[[1]][1]
  l_count_lz <- c(l_count_lz, as.integer(count_lz))
  files_gro_1 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', pattern = dimer)
  for (file in files_gro_1){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""))
      info_ini <- strsplit(info_ini, "\n")
      count_ini <- info_ini[[1]][1]
      l_count_ini_1 <- c(l_count_ini_1, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""))
      info_end <- strsplit(info_end, "\n")
      count_end <- info_end[[1]][1]
      l_count_end_1 <- c(l_count_end_1, as.integer(count_end))
    }
  }
  #Add data
  inter_pdb <- (apply(data_pdb[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_lz <- (apply(data_lz[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_pdb) <- inter_pdb
  colnames(data_pdb) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_lz) <- inter_lz
  colnames(data_lz) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_pdb_t <- data_pdb[,3:4]
  data_lz_t <- data_lz[,3:4]
  data_ini_t <- data_ini[,3:4] 
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by="row.names")
  data_gro_t <- data_gro[, c(2,4)] 
  colnames(data_gro_t) <- c('counts_ini', 'counts_end')
  rownames(data_gro_t) <- data_gro$Row.names
  
  data_all <- merge.data.frame(data_pdb_t, data_lz_t, by="row.names")
  data_all_t <- data_all[, c(1,2,4,5)] 
  colnames(data_all_t) <- c('interactions', 'counts_pdb', 'counts_ini_1')
  rownames(data_all_t) <- data_all_t$interactions
  
  data_all <- merge.data.frame(data_pdb_t, data_gro_t, by="row.names")
  data_all_t <- data_all[, c(1,2,4,5)] 
  colnames(data_all_t) <- c('interactions', 'counts_pdb', 'counts_ini_1', 'counts_end_1')
  rownames(data_all_t) <- data_all_t$interactions
  
  files_gro_2 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', pattern = dimer)
  for (file in files_gro_2){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""))
      info_ini <- strsplit(info_ini, "\n")
      count_ini <- info_ini[[1]][1]
      l_count_ini_2 <- c(l_count_ini_2, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""))
      info_end <- strsplit(info_end, "\n")
      count_end <- info_end[[1]][1]
      l_count_end_2 <- c(l_count_end_2, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4] 
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by="row.names")
  data_gro_t_2 <- data_gro[, c(2,4)] 
  colnames(data_gro_t_2) <- c('counts_ini_2', 'counts_end_2')
  rownames(data_gro_t_2) <- data_gro$Row.names
  
  files_gro_3 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', pattern = dimer)
  for (file in files_gro_3){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""))
      info_ini <- strsplit(info_ini, "\n")
      count_ini <- info_ini[[1]][1]
      l_count_ini_3 <- c(l_count_ini_3, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""))
      info_end <- strsplit(info_end, "\n")
      count_end <- info_end[[1]][1]
      l_count_end_3 <- c(l_count_end_3, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4] 
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by="row.names")
  data_gro_t_3 <- data_gro[, c(2,4)] 
  colnames(data_gro_t_3) <- c('counts_ini_3', 'counts_end_3')
  rownames(data_gro_t_3) <- data_gro$Row.names
  
  files_gro_4 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', pattern = dimer)
  for (file in files_gro_4){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""))
      info_ini <- strsplit(info_ini, '\n')
      count_ini <- info_ini[[1]][1]
      l_count_ini_4 <- c(l_count_ini_4, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""))
      info_end <- strsplit(info_end, '\n')
      count_end <- info_end[[1]][1]
      l_count_end_4 <- c(l_count_end_4, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4]
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by='row.names')
  data_gro_t_4 <- data_gro[, c(2,4)]
  colnames(data_gro_t_4) <- c('counts_ini_4', 'counts_end_4')
  rownames(data_gro_t_4) <- data_gro$Row.names
  
  files_gro_5 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', pattern = dimer)
  for (file in files_gro_5){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""))
      info_ini <- strsplit(info_ini, '\n')
      count_ini <- info_ini[[1]][1]
      l_count_ini_5 <- c(l_count_ini_5, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""))
      info_end <- strsplit(info_end, '\n')
      count_end <- info_end[[1]][1]
      l_count_end_5 <- c(l_count_end_5, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4]
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by='row.names')
  data_gro_t_5 <- data_gro[, c(2,4)]
  colnames(data_gro_t_5) <- c('counts_ini_5', 'counts_end_5')
  rownames(data_gro_t_5) <- data_gro$Row.names
  
  data_all_sys <- merge.data.frame(data_gro_t, data_gro_t_2, by="row.names")
  rownames(data_all_sys) <- data_all_sys$Row.names
  data_all_sys$Row.names <- NULL
  data_all_sys <- merge.data.frame(data_all_sys, data_gro_t_3, by="row.names")
  rownames(data_all_sys) <- data_all_sys$Row.names
  data_all_sys$Row.names <- NULL
  data_all_sys <- merge.data.frame(data_all_sys, data_gro_t_4, by="row.names")
  rownames(data_all_sys) <- data_all_sys$Row.names
  data_all_sys$Row.names <- NULL
  data_all_sys <- merge.data.frame(data_all_sys, data_gro_t_5, by="row.names")
  rownames(data_all_sys) <- data_all_sys$Row.names
  data_all_sys$Row.names <- NULL
  
  data_all_sys$mean_ini <- apply(data_all_sys[,c(1,3,5,7,9)], 1, function(x) mean(x))
  data_all_sys$mean_end <- apply(data_all_sys[,c(2,4,6,8,10)], 1, function(x) mean(x))
  data_all_sys$sd_ini <- apply(data_all_sys[,c(1,3,5,7,9)], 1, function(x) sd(x))
  data_all_sys$sd_end <- apply(data_all_sys[,c(2,4,6,8,10)], 1, function(x) sd(x))
  data_all_sys_ord <- data_all_sys[rev(order(data_all_sys$mean_end)),]
  i = 0
  for (d in data_all_sys_ord[,dim(data_all_sys_ord)[2]]){
    if (d == 0 & i != 0) {
      data_all_sys_ord_f <- data_all_sys_ord[1:i,c(1,2,dim(data_all_sys_ord)[2]-3,dim(data_all_sys_ord)[2]-2,dim(data_all_sys_ord)[2]-1, dim(data_all_sys_ord)[2])]
      break
    }
    if (d == 0 & i == 0) {
      print ('NO INTERACTIONS!')
    }
    i = i+1
  }
  l_interf <- NULL
  l_valuesf <- NULL
  l_interf <- c(l_interf, dimer)
  interficies <- rownames(data_all_sys_ord_f)
  values <- data_all_sys_ord_f$mean_end #
  for (i in (1:length(values))){
    if (values[i] > 1) { #Value of the number of interactions
      l_interf <- c(l_interf, interficies[i])
      l_valuesf <- c(l_valuesf, values[i])
    }
  }
  if (is_empty(l_valuesf)){
    l_interf <- c(l_interf, NA)
    l_valuesf <- c(l_valuesf, 0)
  }
  # number <- strsplit(dimer, "_")
  l_inter_sys_inter[[as.integer(j)]] <- l_interf 
  l_inter_sys_values[[as.integer(j)]] <- l_valuesf
}

for (i in 1:length(l_inter_sys_inter)){
  l <- l_inter_sys_inter[[i]]
  dime <- l[1]
  rowname <- strsplit(dime, "_")
  rowname <- as.integer(rowname[[1]][3])
  interactions <- l[2:length(l)]
  if (is.na(interactions)){
    next
  }
  print(interactions)
  for (j in 1:length(interactions)){
    print (interactions[j])
    colname <- which(colnames(matrix_sys) == interactions[j])
    print(colname)
    if (is_empty(colname)){
      line <- strsplit(interactions[j], "-")
      if (!is.na(line)){
        int <- paste(line[[1]][2], line[[1]][1], sep = "-")
        colname <- which(colnames(matrix_sys) == int)
        print(colname)
        if (is_empty(colname)){
          next
        }
      }
    }
    matrix_sys[dime, colname] <- matrix_sys[dime, colname] + l_inter_sys_values[[i]][j]
  }
}

data_inter_clas_sys <- as.data.frame(matrix_sys)
# data_inter_clas_sys <- apply(data_inter_clas_sys, 2, function(x) factor(x, levels = c(0,1), labels = c('No', 'Yes')))
data_inter_clas_tm_sys <- data_inter_clas_sys[, which(grepl('TM.-TM.', colnames(data_inter_clas_sys)))]
data_inter_clas_tm_sys <- data_inter_clas_tm_sys[, order(colnames(data_inter_clas_tm_sys))]

matrix_asys <- matrix(0, length(dimers_asys), length(inter_end))
colnames(matrix_asys) <- inter_end
rownames(matrix_asys) <- dimers_asys

l_inter_asys_inter <- vector(mode = "list", length = 1)
l_inter_asys_values <- vector(mode = "list", length = 1)

for (j in 1:length(dimers_asys)){
  dimer <- as.character(dimers_asys[j])
  print(dimer)
  files_pdb <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/', pattern = dimer)
  data_pdb <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/', files_pdb, sep = ""), skip = 1)
  info_pdb <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/', files_pdb, sep = ""))
  info_pdb <- strsplit(info_pdb, "\n")
  count_pdb <- info_pdb[[1]][1]
  l_count_pdb <- c(l_count_pdb, as.integer(count_pdb))
  files_lz <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/', pattern = dimer)
  data_lz <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/', files_lz, sep = ""), skip = 1)
  info_lz <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/', files_lz, sep = ""))
  info_lz <- strsplit(info_lz, "\n")
  count_lz <- info_lz[[1]][1]
  l_count_lz <- c(l_count_lz, as.integer(count_lz))
  files_gro_1 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', pattern = dimer)
  for (file in files_gro_1){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""))
      info_ini <- strsplit(info_ini, "\n")
      count_ini <- info_ini[[1]][1]
      l_count_ini_1 <- c(l_count_ini_1, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""))
      info_end <- strsplit(info_end, "\n")
      count_end <- info_end[[1]][1]
      l_count_end_1 <- c(l_count_end_1, as.integer(count_end))
    }
  }
  #Add data
  inter_pdb <- (apply(data_pdb[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_lz <- (apply(data_lz[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_pdb) <- inter_pdb
  colnames(data_pdb) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_lz) <- inter_lz
  colnames(data_lz) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_pdb_t <- data_pdb[,3:4]
  data_lz_t <- data_lz[,3:4]
  data_ini_t <- data_ini[,3:4] 
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by="row.names")
  data_gro_t <- data_gro[, c(2,4)] 
  colnames(data_gro_t) <- c('counts_ini', 'counts_end')
  rownames(data_gro_t) <- data_gro$Row.names
  
  data_all <- merge.data.frame(data_pdb_t, data_lz_t, by="row.names")
  data_all_t <- data_all[, c(1,2,4,5)] 
  colnames(data_all_t) <- c('interactions', 'counts_pdb', 'counts_ini_1')
  rownames(data_all_t) <- data_all_t$interactions
  
  data_all <- merge.data.frame(data_pdb_t, data_gro_t, by="row.names")
  data_all_t <- data_all[, c(1,2,4,5)] 
  colnames(data_all_t) <- c('interactions', 'counts_pdb', 'counts_ini_1', 'counts_end_1')
  rownames(data_all_t) <- data_all_t$interactions
  
  files_gro_2 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', pattern = dimer)
  for (file in files_gro_2){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""))
      info_ini <- strsplit(info_ini, "\n")
      count_ini <- info_ini[[1]][1]
      l_count_ini_2 <- c(l_count_ini_2, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""))
      info_end <- strsplit(info_end, "\n")
      count_end <- info_end[[1]][1]
      l_count_end_2 <- c(l_count_end_2, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4] 
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by="row.names")
  data_gro_t_2 <- data_gro[, c(2,4)] 
  colnames(data_gro_t_2) <- c('counts_ini_2', 'counts_end_2')
  rownames(data_gro_t_2) <- data_gro$Row.names
  
  files_gro_3 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', pattern = dimer)
  for (file in files_gro_3){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""))
      info_ini <- strsplit(info_ini, "\n")
      count_ini <- info_ini[[1]][1]
      l_count_ini_3 <- c(l_count_ini_3, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""))
      info_end <- strsplit(info_end, "\n")
      count_end <- info_end[[1]][1]
      l_count_end_3 <- c(l_count_end_3, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4] 
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by="row.names")
  data_gro_t_3 <- data_gro[, c(2,4)] 
  colnames(data_gro_t_3) <- c('counts_ini_3', 'counts_end_3')
  rownames(data_gro_t_3) <- data_gro$Row.names
  
  files_gro_4 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', pattern = dimer)
  for (file in files_gro_4){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""))
      info_ini <- strsplit(info_ini, '\n')
      count_ini <- info_ini[[1]][1]
      l_count_ini_4 <- c(l_count_ini_4, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""))
      info_end <- strsplit(info_end, '\n')
      count_end <- info_end[[1]][1]
      l_count_end_4 <- c(l_count_end_4, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4]
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by='row.names')
  data_gro_t_4 <- data_gro[, c(2,4)]
  colnames(data_gro_t_4) <- c('counts_ini_4', 'counts_end_4')
  rownames(data_gro_t_4) <- data_gro$Row.names
  
  files_gro_5 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', pattern = dimer)
  for (file in files_gro_5){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""))
      info_ini <- strsplit(info_ini, '\n')
      count_ini <- info_ini[[1]][1]
      l_count_ini_5 <- c(l_count_ini_5, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""))
      info_end <- strsplit(info_end, '\n')
      count_end <- info_end[[1]][1]
      l_count_end_5 <- c(l_count_end_5, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4]
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by='row.names')
  data_gro_t_5 <- data_gro[, c(2,4)]
  colnames(data_gro_t_5) <- c('counts_ini_5', 'counts_end_5')
  rownames(data_gro_t_5) <- data_gro$Row.names
  
  data_all_asys <- merge.data.frame(data_gro_t, data_gro_t_2, by="row.names")
  rownames(data_all_asys) <- data_all_asys$Row.names
  data_all_asys$Row.names <- NULL
  data_all_asys <- merge.data.frame(data_all_asys, data_gro_t_3, by="row.names")
  rownames(data_all_asys) <- data_all_asys$Row.names
  data_all_asys$Row.names <- NULL
  data_all_asys <- merge.data.frame(data_all_asys, data_gro_t_4, by="row.names")
  rownames(data_all_asys) <- data_all_asys$Row.names
  data_all_asys$Row.names <- NULL
  data_all_asys <- merge.data.frame(data_all_asys, data_gro_t_5, by="row.names")
  rownames(data_all_asys) <- data_all_asys$Row.names
  data_all_asys$Row.names <- NULL
  
  data_all_asys$mean_ini <- apply(data_all_asys[,c(1,3,5,7,9)], 1, function(x) mean(x))
  data_all_asys$mean_end <- apply(data_all_asys[,c(2,4,6,8,10)], 1, function(x) mean(x))
  data_all_asys$sd_ini <- apply(data_all_asys[,c(1,3,5,7,9)], 1, function(x) sd(x))
  data_all_asys$sd_end <- apply(data_all_asys[,c(2,4,6,8,10)], 1, function(x) sd(x))
  data_all_asys_ord <- data_all_asys[rev(order(data_all_asys$mean_end)),]
  i = 0
  for (d in data_all_asys_ord[,dim(data_all_asys_ord)[2]]){
    if (d == 0 & i != 0) {
      data_all_asys_ord_f <- data_all_asys_ord[1:i,c(1,2,dim(data_all_asys_ord)[2]-3,dim(data_all_asys_ord)[2]-2,dim(data_all_asys_ord)[2]-1, dim(data_all_asys_ord)[2])]
      break
    }
    if (d == 0 & i == 0) {
      print ('NO INTERACTIONS!')
    }
    i = i+1
  }
  l_interf <- NULL
  l_valuesf <- NULL
  l_interf <- c(l_interf, dimer)
  interficies <- rownames(data_all_asys_ord_f)
  values <- data_all_asys_ord_f$mean_end #
  for (i in (1:length(values))){
    if (values[i] > 1) { #Value of the number of interactions
      l_interf <- c(l_interf, interficies[i])
      l_valuesf <- c(l_valuesf, values[i])
    }
  }
  if (is_empty(l_valuesf)){
    l_interf <- c(l_interf, NA)
    l_valuesf <- c(l_valuesf, 0)
  }
  # number <- strsplit(dimer, "_")
  l_inter_asys_inter[[as.integer(j)]] <- l_interf 
  l_inter_asys_values[[as.integer(j)]] <- l_valuesf
}

for (i in 1:length(l_inter_asys_inter)){
  l <- l_inter_asys_inter[[i]]
  dime <- l[1]
  rowname <- strsplit(dime, "_")
  rowname <- as.integer(rowname[[1]][3])
  interactions <- l[2:length(l)]
  if (is.na(interactions)){
    next
  }
  for (j in 1:length(interactions)){
    colname <- which(colnames(matrix_asys) == interactions[j])
    if (is_empty(colname)){
      line <- strsplit(interactions[j], "-")
      if (!is.na(line)){
      int <- paste(line[[1]][2], line[[1]][1], sep = "-")
      colname <- which(colnames(matrix_asys) == int)
      if (is_empty(colname)){
        next
      }
      }
    
    }
    matrix_asys[dime, colname] <- matrix_asys[dime, colname] + l_inter_asys_values[[i]][j]
  }
  }

data_inter_clas_asys <- as.data.frame(matrix_asys)
# data_inter_clas_asys <- apply(data_inter_clas_asys, 2, function(x) factor(x, levels = c(0,1), labels = c('No', 'Yes')))
data_inter_clas_tm_asys <- data_inter_clas_asys[, which(grepl('TM.-TM.', colnames(data_inter_clas_asys)))]
data_inter_clas_tm_asys <- data_inter_clas_tm_asys[, order(colnames(data_inter_clas_tm_asys))]

###########COUNTING EACH DIMER
#COUNTS
l_count_pdb <- c()
l_count_lz <- c()
l_count_ini_1 <- c()
l_count_end_1 <- c()
l_count_ini_2 <- c()
l_count_end_2 <- c()
l_count_ini_3 <- c()
l_count_end_3 <- c()
l_count_ini_4 <- c()
l_count_end_4 <- c()
l_count_ini_5 <- c()
l_count_end_5 <- c()

matrix_inter <- matrix(0, length(dimers), length(inter_end))
colnames(matrix_inter) <- inter_end
rownames(matrix_inter) <- dimers

l_inter <- vector(mode = "list", length = 1)

for (dimer in dimers){
  print(dimer)
  print(dimer[:-4])
  files_pdb <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/', pattern = dimer)
  data_pdb <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/', files_pdb, sep = ""), skip = 1)
  info_pdb <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb/', files_pdb, sep = ""))
  info_pdb <- strsplit(info_pdb, "\n")
  count_pdb <- info_pdb[[1]][1]
  l_count_pdb <- c(l_count_pdb, as.integer(count_pdb))
  files_lz <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/', pattern = dimer)
  data_lz <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/', files_lz, sep = ""), skip = 1)
  info_lz <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_pdb_lz/', files_lz, sep = ""))
  info_lz <- strsplit(info_lz, "\n")
  count_lz <- info_lz[[1]][1]
  l_count_lz <- c(l_count_lz, as.integer(count_lz))
  files_gro_1 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', pattern = dimer)
  for (file in files_gro_1){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""))
      info_ini <- strsplit(info_ini, "\n")
      count_ini <- info_ini[[1]][1]
      l_count_ini_1 <- c(l_count_ini_1, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_1/', file, sep = ""))
      info_end <- strsplit(info_end, "\n")
      count_end <- info_end[[1]][1]
      l_count_end_1 <- c(l_count_end_1, as.integer(count_end))
    }
  }
  #Add data
  inter_pdb <- (apply(data_pdb[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_lz <- (apply(data_lz[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_pdb) <- inter_pdb
  colnames(data_pdb) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_lz) <- inter_lz
  colnames(data_lz) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_pdb_t <- data_pdb[,3:4]
  data_lz_t <- data_lz[,3:4]
  data_ini_t <- data_ini[,3:4] 
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by="row.names")
  data_gro_t <- data_gro[, c(2,4)] 
  colnames(data_gro_t) <- c('counts_ini', 'counts_end')
  rownames(data_gro_t) <- data_gro$Row.names
  
  data_all <- merge.data.frame(data_pdb_t, data_lz_t, by="row.names")
  data_all_t <- data_all[, c(1,2,4,5)] 
  colnames(data_all_t) <- c('interactions', 'counts_pdb', 'counts_ini_1')
  rownames(data_all_t) <- data_all_t$interactions
  
  data_all <- merge.data.frame(data_pdb_t, data_gro_t, by="row.names")
  data_all_t <- data_all[, c(1,2,4,5)] 
  colnames(data_all_t) <- c('interactions', 'counts_pdb', 'counts_ini_1', 'counts_end_1')
  rownames(data_all_t) <- data_all_t$interactions
  
  files_gro_2 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', pattern = dimer)
  for (file in files_gro_2){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""))
      info_ini <- strsplit(info_ini, "\n")
      count_ini <- info_ini[[1]][1]
      l_count_ini_2 <- c(l_count_ini_2, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_2/', file, sep = ""))
      info_end <- strsplit(info_end, "\n")
      count_end <- info_end[[1]][1]
      l_count_end_2 <- c(l_count_end_2, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4] 
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by="row.names")
  data_gro_t_2 <- data_gro[, c(2,4)] 
  colnames(data_gro_t_2) <- c('counts_ini_2', 'counts_end_2')
  rownames(data_gro_t_2) <- data_gro$Row.names
  
  files_gro_3 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', pattern = dimer)
  for (file in files_gro_3){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""))
      info_ini <- strsplit(info_ini, "\n")
      count_ini <- info_ini[[1]][1]
      l_count_ini_3 <- c(l_count_ini_3, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_3/', file, sep = ""))
      info_end <- strsplit(info_end, "\n")
      count_end <- info_end[[1]][1]
      l_count_end_3 <- c(l_count_end_3, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4] 
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by="row.names")
  data_gro_t_3 <- data_gro[, c(2,4)] 
  colnames(data_gro_t_3) <- c('counts_ini_3', 'counts_end_3')
  rownames(data_gro_t_3) <- data_gro$Row.names
  
  files_gro_4 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', pattern = dimer)
  for (file in files_gro_4){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""))
      info_ini <- strsplit(info_ini, '\n')
      count_ini <- info_ini[[1]][1]
      l_count_ini_4 <- c(l_count_ini_4, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_4/', file, sep = ""))
      info_end <- strsplit(info_end, '\n')
      count_end <- info_end[[1]][1]
      l_count_end_4 <- c(l_count_end_4, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4]
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by='row.names')
  data_gro_t_4 <- data_gro[, c(2,4)]
  colnames(data_gro_t_4) <- c('counts_ini_4', 'counts_end_4')
  rownames(data_gro_t_4) <- data_gro$Row.names
  
  files_gro_5 <- list.files(path = '/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', pattern = dimer)
  for (file in files_gro_5){
    if (grepl('ini', file)){
      data_ini <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""), skip = 1)
      info_ini <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""))
      info_ini <- strsplit(info_ini, '\n')
      count_ini <- info_ini[[1]][1]
      l_count_ini_5 <- c(l_count_ini_5, as.integer(count_ini))
    }
    if (grepl('_p_', file)){
      data_end <- read.table(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""), skip = 1)
      info_end <- read_file(file = paste('/home/adrian/Documents/bender/RESULTS/COUNTS_gro_5/', file, sep = ""))
      info_end <- strsplit(info_end, '\n')
      count_end <- info_end[[1]][1]
      l_count_end_5 <- c(l_count_end_5, as.integer(count_end))
    }
  }
  #Add data
  inter_ini <- (apply(data_ini[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  inter_end <- (apply(data_end[,c(1,2)], 1, function(x) paste(x[1], x[2], sep = '-')))
  
  rownames(data_ini) <- inter_ini
  colnames(data_ini) <- c('type1', 'type2', 'counts', 'percent')
  rownames(data_end) <- inter_end
  colnames(data_end) <- c('type1', 'type2', 'counts', 'percent')
  
  data_ini_t <- data_ini[,3:4]
  data_end_t <- data_end[,3:4]
  
  data_gro <- merge.data.frame(data_ini_t, data_end_t, by='row.names')
  data_gro_t_5 <- data_gro[, c(2,4)]
  colnames(data_gro_t_5) <- c('counts_ini_5', 'counts_end_5')
  rownames(data_gro_t_5) <- data_gro$Row.names
  
  data_all_t <- merge.data.frame(data_all_t, data_gro_t_2, by="row.names")
  rownames(data_all_t) <- data_all_t$Row.names
  data_all_t$Row.names <- NULL
  data_all_t <- merge.data.frame(data_all_t, data_gro_t_3, by="row.names")
  rownames(data_all_t) <- data_all_t$Row.names
  data_all_t$Row.names <- NULL
  data_all_t <- merge.data.frame(data_all_t, data_gro_t_4, by="row.names")
  rownames(data_all_t) <- data_all_t$Row.names
  data_all_t$Row.names <- NULL
  data_all_t <- merge.data.frame(data_all_t, data_gro_t_5, by="row.names")
  rownames(data_all_t) <- data_all_t$Row.names
  data_all_t$Row.names <- NULL
  
  data_all_t$mean_ini <- apply(data_all_t[,c(3,5,7,9,11)], 1, function(x) mean(x))
  data_all_t$mean_end <- apply(data_all_t[,c(4,6,8,10,12)], 1, function(x) mean(x))
  data_all_t$sd_ini <- apply(data_all_t[,c(3,5,7,9,11)], 1, function(x) sd(x))
  data_all_t$sd_end <- apply(data_all_t[,c(4,6,8,10,12)], 1, function(x) sd(x))
  data_all_t_ord <- data_all_t[rev(order(data_all_t$mean_end)),]
  i = 0
  for (d in data_all_t_ord[,dim(data_all_t_ord)[2]]){
    if (d == 0 & i != 0) {
      data_all_t_ord_f <- data_all_t_ord[1:i,c(1,2,dim(data_all_t_ord)[2]-3,dim(data_all_t_ord)[2]-2,dim(data_all_t_ord)[2]-1, dim(data_all_t_ord)[2])]
      break
    }
    if (d == 0 & i == 0) {
      print ('NO INTERACTIONS!')
    }
    i = i+1
  }
  l_interf <- NULL
  l_interf <- c(l_interf, dimer)
  interficies <- data_all_t_ord_f$interactions
  values <- data_all_t_ord_f$mean_end #
  for (i in (1:length(values))){
    if (values[i] > 1) { #Value of the number of interactions
      l_interf <- c(l_interf, interficies[i])
    }
  }
  number <- strsplit(dimer, "_")
  l_inter[[as.integer(number[[1]][3])]] <- l_interf
  
  # plot_data <- data_all_t_ord_f[,c(1,3,4)]
  # colnames(plot_data) <- c("interactions", "Initial", "End")
  # plot_data <- melt(plot_data, id.vars="interactions")
  # plot_data$interactions <- as.factor(plot_data$interactions)
  # plot <- ggplot(plot_data, aes(interactions, value)) +
  #   xlab('Interactions')  + ylab('Counts') +
  #   geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
  #   scale_x_discrete(breaks = unique(plot_data$interactions)) +
  #   theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.title = element_blank())
  # 
  # plot <- plot +  ggtitle(paste('Interactions of ', dimer, sep = ""))
  # 
  # ggsave(paste('/home/adrian/Documents/bender/RESULTS/INTERFICIES_PLOT/', dimer,sep = ""), device = 'png')
  
}

for (l in l_inter){
  dime <- l[1]
  print (dime)
  rowname <- strsplit(dime, "_")
  rowname <- as.integer(rowname[[1]][3])
  print (rowname)
  interactions <- l[2:length(l)]
  if (is_empty(interactions)){
    next 
  }
  for (i in interactions){
    colname <- which(colnames(matrix_inter) == i )
    if (is_empty(colname)){
      line <- strsplit(i, "-")
      i <- paste(line[[1]][2], line[[1]][1], sep = "-")
      colname <- which(colnames(matrix_inter) == i)
      if (is_empty(colname)){
        next
      }
    }
    print (i)
    matrix_inter[rowname, colname] <- matrix_inter[rowname, colname] + 1 
  }
}

data_inter_clas <- as.data.frame(matrix_inter)
data_inter_clas <- apply(data_inter_clas, 2, function(x) factor(x, levels = c(0,1), labels = c('No', 'Yes')))
data_inter_clas_tm <- data_inter_clas[, which(grepl('TM.-TM.', colnames(data_inter_clas)))]
data_inter_clas_tm <- data_inter_clas_tm[, order(colnames(data_inter_clas_tm))]

max_inter_clas_tm <- matrix_inter[, which(grepl('TM.-TM.', colnames(matrix_inter)))]
max_inter_clas_tm <- max_inter_clas_tm[, order(colnames(max_inter_clas_tm))]

write_csv(data.frame(max_inter_clas_tm), "TM_zone.txt")

interactions <- data.frame(dimers, l_count_pdb, l_count_lz, l_count_ini_1, l_count_end_1, l_count_ini_2, l_count_end_2, l_count_ini_3, 
                           l_count_end_3, l_count_ini_4, l_count_end_4, l_count_ini_5, l_count_end_5)

interactions$MEAN_ini <- apply(interactions[,c(4,6,8,10,12)], 1, function(x) mean(x))
interactions$MEAN_end <- apply(interactions[,c(5,7,9,11,13)], 1, function(x) mean(x))
interactions$SD_ini <- apply(interactions[,c(4,6,8,10,12)], 1, function(x) sd(x))
interactions$SD_end <- apply(interactions[,c(5,7,9,11,13)], 1, function(x) sd(x))

##RMSD#########################################################################################
rmsd_1 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_rmsd_1.txt", header = T)
rmsd_1$mean <- apply(rmsd_1[,3:69], 1, function(x) mean(x))
rmsd_1$sd <- apply(rmsd_1[,3:69], 1, function(x) sd(x))

rmsd_2 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_rmsd_2.txt", header = T)
rmsd_2$mean <- apply(rmsd_2[,3:69], 1, function(x) mean(x))
rmsd_2$sd <- apply(rmsd_2[,3:69], 1, function(x) sd(x))

rmsd_3 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_rmsd_3.txt", header = T)
rmsd_3$mean <- apply(rmsd_3[,3:69], 1, function(x) mean(x))
rmsd_3$sd <- apply(rmsd_3[,3:69], 1, function(x) sd(x))

rmsd_4 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_rmsd_4.txt", header = T)
rmsd_4$mean <- apply(rmsd_4[,3:69], 1, function(x) mean(x))
rmsd_4$sd <- apply(rmsd_4[,3:69], 1, function(x) sd(x))

rmsd_5 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_rmsd_5.txt", header = T)
rmsd_5$mean <- apply(rmsd_5[,3:69], 1, function(x) mean(x))
rmsd_5$sd <- apply(rmsd_5[,3:69], 1, function(x) sd(x))

rmsd_all <- cbind.data.frame(rmsd_1[, c(1,2)], rmsd_1[,c(70,71)], rmsd_2[,c(70,71)], rmsd_3[,c(70,71)], rmsd_4[,c(70,71)], rmsd_5[,c(70,71)])
colnames(rmsd_all) <- c('DIMER', 'FAMILY', 'RMSD_MEAN_1', 'RMSD_SD_1', 'RMSD_MEAN_2' , 'RMSD_SD_2', 'RMSD_MEAN_3' , 'RMSD_SD_3', 'RMSD_MEAN_4' , 'RMSD_SD_4', 'RMSD_MEAN_5' , 'RMSD_SD_5')

rmsd_all$RMSD_MEAN <- apply(rmsd_all[,c(3,5,7,9,11)], 1, function(x) mean(x))
rmsd_all$RMAS_SD <- apply(rmsd_all[,c(4,6,8,10,12)], 1, function(x) sd(x))

rmsd_all$error <- apply(rmsd_all[,c(13,14)], 1, function(x) {
  error <- qt(0.975,df=4-1)*x[[2]]/sqrt(4)
  return(error)
})

# rmsd_all$NormalDistribution <- apply(rmsd_all[,c(3,5,7,9)], 1, function(x) {
#   pv <- shapiro.test(x)$p.value 
#   if (pv <= 0.05) {return('NO')}
#   if (pv > 0.05) {return('YES')}
# })
# nd_rmsd <- rmsd_1[rmsd_1$NormalDistribution == 'YES',]
# nd_rmsd$MeanTest <- apply(nd_rmsd[,3:69], 1, function(x) t.test(x, mu = 1)$p.value)
# nnd_rmsd <- rmsd_1[rmsd_1$NormalDistribution == 'NO',]
# nnd_rmsd$MeanTest <- apply(nnd_rmsd[,3:69], 1, function(x) wilcox.test(x, mu = 1)$p.value)
# rmsd_1$ci_l <- apply(rmsd_1[,3:69], 1, function(x) t.test(x)$conf.int[1])
# rmsd_1$ci_r <- apply(rmsd_1[,3:69], 1, function(x) t.test(x)$conf.int[2])

#DEGREES#########################################################################################
degrees_1 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_angles_1.txt", header = T)
degrees_1$mean <- apply(degrees_1[,4:54], 1, function(x) mean(x))
degrees_1$sd <- apply(degrees_1[,4:54], 1, function(x) sd(x))

degrees_2 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_angles_2.txt", header = T)
degrees_2$mean <- apply(degrees_2[,4:54], 1, function(x) mean(x))
degrees_2$sd <- apply(degrees_2[,4:54], 1, function(x) sd(x))

degrees_3 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_angles_3.txt", header = T)
degrees_3$mean <- apply(degrees_3[,4:54], 1, function(x) mean(x))
degrees_3$sd <- apply(degrees_3[,4:54], 1, function(x) sd(x))

degrees_4 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_angles_4.txt", header = T)
degrees_4$mean <- apply(degrees_4[,4:54], 1, function(x) mean(x))
degrees_4$sd <- apply(degrees_4[,4:54], 1, function(x) sd(x))

degrees_5 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_angles_5.txt", header = T)
degrees_5$mean <- apply(degrees_5[,4:54], 1, function(x) mean(x))
degrees_5$sd <- apply(degrees_5[,4:54], 1, function(x) sd(x))

degrees_all <- cbind.data.frame(degrees_1[, c(1,2)], degrees_1[,c(55,56)], degrees_2[,c(55,56)], degrees_3[,c(55,56)], degrees_4[,c(55,56)], degrees_5[,c(55,56)])
colnames(degrees_all) <- c('DIMER', 'FAMILY', 'DEGREE_MEAN_1', 'DEGREE_SD_1', 'DEGREE_MEAN_2' , 'DEGREE_SD_2', 'DEGREE_MEAN_3' , 'DEGREE_SD_3', 'DEGREE_MEAN_4' , 'DEGREE_SD_4', 'DEGREE_MEAN_5' , 'DEGREE_SD_5')

degrees_all$DEGREE_MEAN <- apply(degrees_all[,c(3,5,7,9,11)], 1, function(x) mean(x))
degrees_all$DEGREE_SD <- apply(degrees_all[,c(4,6,8,10,12)], 1, function(x) sd(x))

degrees_all$error <- apply(degrees_all[,c(13,14)], 1, function(x) {
  error <- qt(0.975,df=4-1)*x[[2]]/sqrt(4)
  return(error)
})

#DIST#########################################################################################
dist_1 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_dist_1.txt", header = T)
dist_1$mean <- apply(dist_1[,4:54], 1, function(x) mean(x))
dist_1$sd <- apply(dist_1[,4:54], 1, function(x) sd(x))

dist_2 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_dist_2.txt", header = T)
dist_2$mean <- apply(dist_2[,4:54], 1, function(x) mean(x))
dist_2$sd <- apply(dist_2[,4:54], 1, function(x) sd(x))

dist_3 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_dist_3.txt", header = T)
dist_3$mean <- apply(dist_3[,4:54], 1, function(x) mean(x))
dist_3$sd <- apply(dist_3[,4:54], 1, function(x) sd(x))

dist_4 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_dist_4.txt", header = T)
dist_4$mean <- apply(dist_4[,4:54], 1, function(x) mean(x))
dist_4$sd <- apply(dist_4[,4:54], 1, function(x) sd(x))

dist_5 <- read.table(file = "/home/adrian/Documents/bender/RESULTS/dimer_dist_5.txt", header = T)
dist_5$mean <- apply(dist_5[,4:54], 1, function(x) mean(x))
dist_5$sd <- apply(dist_5[,4:54], 1, function(x) sd(x))

dist_all <- cbind.data.frame(dist_1[, c(1,2)], dist_1[,c(55,56)], dist_2[,c(55,56)], dist_3[,c(55,56)], dist_4[,c(55,56)], dist_5[,c(55,56)])
colnames(dist_all) <- c('DIMER', 'FAMILY', 'DIST_MEAN_1', 'DIST_SD_1', 'DIST_MEAN_2' , 'DIST_SD_2', 'DIST_MEAN_3' , 'DIST_SD_3', 'DIST_MEAN_4' , 'DIST_SD_4', 'DIST_MEAN_5' , 'DIST_SD_5')

dist_all$DIST_MEAN <- apply(dist_all[,c(3,5,7,9,11)], 1, function(x) mean(x))
dist_all$DIST_SD <- apply(dist_all[,c(4,6,8,10,12)], 1, function(x) sd(x))

dist_all$error <- apply(dist_all[,c(13,14)], 1, function(x) {
  error <- qt(0.975,df=4-1)*x[[2]]/sqrt(4)
  return(error)
})

#FUSION PROT/ ASYS#########################################################################################
fus_prot <- read.table(file = "/home/adrian/Documents/bender/RESULTS/fusion_prot.txt", header = T, sep = '\t')
fus_prot <- fus_prot[order(fus_prot$DIMER),]
fus_prot$FUSION <- as.character(fus_prot$FUSION)
fusion <- c()
for (d in rmsd_all$DIMER){
  d <- strsplit(d, '_')
  d <- d[[1]]
  d <- d[[1]]
  fus <- fus_prot[fus_prot$DIMER == d, ]
  print (str(fus[1,2]))
  fusion <- c(fusion, fus[1,2])
}

info_all <- cbind.data.frame(rmsd_all[, c(1,2)], fusion, symmetry, interactions[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)], rmsd_all[,c(3,4,5,6,7,8,9,10,11,12,13,14,15)], degrees_all[,c(3,4,5,6,7,8,9,10,11,12,13,14,15)], dist_all[,c(3,4,5,6,7,8,9,10,11,12,13,14,15)])
hh <- c('3ODU_CXCR4_6', '4U15_ACM3_24')
info_all$symmetry[info_all$DIMER %in% hh] <- 'YES'
write.csv(info_all, file = '/home/adrian/RESULTS/dimer_general_info.csv', quote = F, row.names = F)

info_all[rev(order(info_all$MEAN_end)),]

resum <- info_all[, c(1,2,4,5,6,17,18,19,20,31,32,33,44,45,46,57,58,59)]

### TOP 20 ###
resum <- resum[rev(order(resum$MEAN_end)),]
resum_20 <- resum[1:20,]

#PDB AND GRO COUNTS##########################################################################################################
pd <- position_dodge(0.1) # move them .05 to the left and right

pdb_data_asys <- resum_asys[,c(1,5,4)]
colnames(pdb_data_asys) <- c("Dimers", "With fusion protein", "Without fusion protein")
pdb_data_asys <- melt(pdb_data_asys, id.vars="Dimers")

pdb_data_sys <- resum_sys[,c(1,5,4)]
colnames(pdb_data_sys) <- c("Dimers", "With fusion protein", "Without fusion protein")
pdb_data_sys <- melt(pdb_data_sys, id.vars="Dimers")
dimers_sys[1:24]
add <- data.frame(Dimers = sapply(dimers_sys[1:24], function(x) paste(x,1000,sep = "")), variable = c(rep("With fusion protein", 24)), value = c(rep(0,24)))
pdb_data_asys <- rbind.data.frame(pdb_data_asys, add)
pdb_asys <- ggplot(pdb_data_asys, aes(Dimers, value)) + 
  ylab('Counts') +
  geom_bar(aes(fill = variable), position = position_dodge2(width = 1, padding = 0), stat="identity", colour="#666666", width = 0.9) +
  scale_color_hue(labels = c("Without lz", "With lz")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14), legend.title = element_blank(), legend.position = "none", 
        legend.direction = "horizontal", legend.text = element_text(size = 14),
        axis.title.x = element_blank(), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16)) +
  scale_y_continuous(limits = c(0,90), expand = c(0,0))

pdb_sys <- ggplot(pdb_data_sys, aes(Dimers, value)) + 
  ylab('Counts') +
  geom_bar(aes(fill = variable), position = position_dodge2(width = 1, padding = 0), stat="identity", colour="#666666", width = 0.9) +
  scale_color_hue(labels = c("Without lz", "With lz")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14), legend.title = element_blank(), legend.position = "none", 
        axis.title.x = element_blank(), axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 16)) +
  scale_y_continuous(limits = c(0,90), expand = c(0,0))  
  


c <- grid.arrange(pdb_sys,                                    # bar plot spaning two columns
                  pdb_asys,                               # box plot and scatter plot
                  ncol = 1, nrow = 2)

ggsave(c, filename = "/home/adrian/Escritorio/RESULTS/counts_pdb.jpg", width = 20, height = 10, dpi = "retina")

#RMSD, DEGREE, DIST PLOT############################################################################################
resum$rmsd_error <- rmsd_all$error
resum$degrees_error <- degrees_all$error
resum$dist_error <- dist_all$error
resum_sys <- resum[resum$symmetry == "YES",]
resum_asys <- resum[resum$symmetry == "NO",]

add <- data.frame(DIMER = sapply(dimers_sys[1:24], function(x) paste(x,1000,sep = "")), RMSD_MEAN = c(rep(0, 24)), rmsd_error = c(rep(0,24)))
time_data_1 <- resum_asys[,c(1,10,19)]
time_data_1_c <- rbind.data.frame(time_data_1, add)
add <- data.frame(DIMER = sapply(dimers_sys[1:24], function(x) paste(x,1000,sep = "")), DEGREE_MEAN = c(rep(0, 24)), degrees_error = c(rep(0,24)))
time_data_2 <- resum_asys[,c(1,13,20)]
time_data_2_c <- rbind.data.frame(time_data_2, add)
add <- data.frame(DIMER = sapply(dimers_sys[1:24], function(x) paste(x,1000,sep = "")), DIST_MEAN = c(rep(0, 24)), dist_error = c(rep(0,24)))
time_data_3 <- resum_asys[,c(1,16,21)]
time_data_3_c <- rbind.data.frame(time_data_3, add)

p1_a <- ggplot(time_data_1_c, aes(DIMER, RMSD_MEAN)) + 
  ylab('RMSD') +
  geom_bar(position = "dodge", stat="identity", fill = '#80aaff', colour="#666666") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank(), axis.text.x.bottom = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank()) + 
  geom_errorbar(aes(ymin=RMSD_MEAN-rmsd_error, ymax=RMSD_MEAN+rmsd_error), width=.1, position=pd, colour="#666666") +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_hline(yintercept = 2.5, linetype = "dashed")+
  scale_y_continuous(limits = c(0,8), expand = c(0,0))


p2_a <- ggplot(time_data_2_c, aes(DIMER, DEGREE_MEAN)) + 
  ylab('Degrees (º)') +
  geom_bar(position = "dodge", stat="identity", fill = '#70db70', colour="#666666") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank(), axis.text.x.bottom = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank()) +
  geom_errorbar(aes(ymin=DEGREE_MEAN-degrees_error, ymax=DEGREE_MEAN+degrees_error), width=.1, position=pd, colour="#666666") +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  scale_y_continuous(limits = c(0,8), expand = c(0,0))


p3_a <- ggplot(time_data_3_c, aes(DIMER, DIST_MEAN)) + 
  coord_cartesian(ylim = c(30,80)) +
  ylab('Distance (A)') +
  geom_bar(position = "dodge", stat="identity", fill = '#ff8080', colour="#666666") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank()) +
  geom_errorbar(aes(ymin=DIST_MEAN-dist_error, ymax=DIST_MEAN+dist_error), width=.1, position=pd, colour="#666666") +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  scale_y_continuous(limits = c(0,80), expand = c(0,0))

gro_data_asys <- resum_asys[,c(1,6,8,7,9)]
colnames(gro_data_asys) <- c("Dimers", "Initial", "SD_ini", "End", "SD_end")
gro_data_asys$error_ini <- apply(gro_data_asys[,c(2,3)], 1, function(x) {
  error <- qt(0.975,df=4-1)*x[[2]]/sqrt(4)
  return(error)
})
gro_data_asys$error_end <- apply(gro_data_asys[,c(4,5)], 1, function(x) {
  error <- qt(0.975,df=4-1)*x[[2]]/sqrt(4)
  return(error)
})

gro_data_asys_f <- melt(gro_data_asys[,c(1,2,4)], id.vars="Dimers")
error <- c(gro_data_asys$error_ini, gro_data_asys$error_end)
gro_data_asys_f$error <- error

add <- data.frame(Dimers = sapply(dimers_sys[1:24], function(x) paste(x,0,sep = "")), variable = c(rep("End", 24)), value = c(rep(0,24)), error = c(rep(0,24)))
gro_data_asys_f_c <- rbind.data.frame(gro_data_asys_f, add)

gro_a <- ggplot(gro_data_asys_f_c, aes(x=as.factor(Dimers), y=value,  fill=variable)) + 
  ylab('Counts') +
  geom_bar(aes(fill = variable), position = "dodge", stat="identity", colour="#666666") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16), legend.title = element_blank(), legend.direction = "horizontal", 
        legend.position = c(0.3,0.85) , axis.title.x = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank()) +
  geom_errorbar(aes(ymin=value-error, ymax=value+error), width=.1, position=position_dodge(.9), colour="#666666") +
  geom_line(position=position_dodge(.9)) +
  geom_point(position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0,80), expand = c(0,0))

time_data_1 <- resum_sys[,c(1,10,19)]
time_data_2 <- resum_sys[,c(1,13,20)]
time_data_3 <- resum_sys[,c(1,16,21)]

p1_s <- ggplot(time_data_1, aes(DIMER, RMSD_MEAN)) + 
  ylab('RMSD') +
  geom_bar(position = "dodge", stat="identity", fill = '#80aaff', colour="#666666") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank(), axis.text.x.bottom = element_blank(),
        axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 16)) + 
  geom_errorbar(aes(ymin=RMSD_MEAN-rmsd_error, ymax=RMSD_MEAN+rmsd_error), width=.1, position=pd, colour="#666666") +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_hline(yintercept = 2.5, linetype = "dashed")+
  scale_y_continuous(limits = c(0,8.5), expand = c(0,0))

p2_s <- ggplot(time_data_2, aes(DIMER, DEGREE_MEAN)) + 
  ylab('Degrees (º)') +
  geom_bar(position = "dodge", stat="identity", fill = '#70db70', colour="#666666") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank(), axis.text.x.bottom = element_blank(),
        axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 16)) +
  geom_errorbar(aes(ymin=DEGREE_MEAN-degrees_error, ymax=DEGREE_MEAN+degrees_error), width=.1, position=pd, colour="#666666") +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_hline(yintercept = 3, linetype = "dashed")+
  scale_y_continuous(limits = c(0,8), expand = c(0,0))


p3_s <- ggplot(time_data_3, aes(DIMER, DIST_MEAN)) + 
  coord_cartesian(ylim = c(30,80)) +
  ylab('Distance (A)') +
  geom_bar(position = "dodge", stat="identity", fill = '#ff8080', colour="#666666") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x = element_blank(), axis.text.x.bottom = element_blank(),
        axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 16)) +
  geom_errorbar(aes(ymin=DIST_MEAN-dist_error, ymax=DIST_MEAN+dist_error), width=.1, position=pd, colour="#666666") +
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  scale_y_continuous(limits = c(0,80), expand = c(0,0))

gro_data_sys <- resum_sys[,c(1,6,8,7,9)]
colnames(gro_data_sys) <- c("Dimers", "Initial", "SD_ini", "End", "SD_end")
gro_data_sys$error_ini <- apply(gro_data_sys[,c(2,3)], 1, function(x) {
  error <- qt(0.975,df=4-1)*x[[2]]/sqrt(4)
  return(error)
})
gro_data_sys$error_end <- apply(gro_data_sys[,c(4,5)], 1, function(x) {
  error <- qt(0.975,df=4-1)*x[[2]]/sqrt(4)
  return(error)
})

gro_data_sys_f <- melt(gro_data_sys[,c(1,2,4)], id.vars="Dimers")
error <- c(gro_data_sys$error_ini, gro_data_sys$error_end)
gro_data_sys_f$error <- error

gro_s <- ggplot(gro_data_sys_f, aes(x=as.factor(Dimers), y=value,  fill=variable)) + 
  ylab('Counts') +
  geom_bar(aes(fill = variable), position = "dodge", stat="identity", colour="#666666") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16), legend.title = element_blank(), legend.direction = "horizontal", 
        legend.position = c(0.85,0.85) , axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 18)) +
  geom_errorbar(aes(ymin=value-error, ymax=value+error), width=.1, position=position_dodge(.9), colour="#666666") +
  geom_line(position=position_dodge(.9)) +
  geom_point(position=position_dodge(.9)) +
  scale_y_continuous(limits = c(0,80), expand = c(0,0))

p_a <- grid.arrange(p1_a,                                    # bar plot spaning two columns
                  p2_a,
                  p3_a,                               # box plot and scatter plot
                  ncol = 1, nrow = 3)

p_s <- grid.arrange(p1_s,                                    # bar plot spaning two columns
                    p2_s,
                    p3_s,                              # box plot and scatter plot
                    ncol = 1, nrow = 3)

p <- grid.arrange(gro_s, gro_a,
                    ncol = 2, nrow = 1)

ggsave(p1_a, filename = "/home/adrian/Escritorio/RESULTS/analysis_asys_rmsd.png", width = 13, height = 5, dpi = "retina")

ggsave(p2_a, filename = "/home/adrian/Escritorio/RESULTS/analysis_asys_degree.png", width = 13, height = 5, dpi = "retina")

ggsave(p3_a, filename = "/home/adrian/Escritorio/RESULTS/analysis_asys_dist.png", width = 13, height = 5, dpi = "retina")

ggsave(p1_s, filename = "/home/adrian/Escritorio/RESULTS/analysis_sys_rmsd.png", width = 13, height = 5, dpi = "retina")

ggsave(p2_s, filename = "/home/adrian/Escritorio/RESULTS/analysis_sys_degree.png", width = 13, height = 5, dpi = "retina")

ggsave(p3_s, filename = "/home/adrian/Escritorio/RESULTS/analysis_sys_dist.png", width = 13, height = 5, dpi = "retina")

ggsave(gro_s, filename = "/home/adrian/Escritorio/RESULTS/analysis_sys_gro.png", width = 13, height = 5, dpi = "retina")

ggsave(gro_a, filename = "/home/adrian/Escritorio/RESULTS/analysis_asys_gro.png", width = 13, height = 5, dpi = "retina")

#INTERFICIES PLOT############################################################################################
values_asys <- apply(t(data_inter_clas_asys), 1, function(x) sum(x))
interactions_asys <- rownames(t(data_inter_clas_asys))
data_asys <- data.frame(values_asys, interactions_asys)
data_asys <- data_asys[rev(order(data_asys$values_asys)),]
l_group <- c()
for (inter in data_asys$interactions_asys){
  if (grepl('TM.-TM.', inter)) {
    l_group <- c(l_group, 'TM')
    next
  }
  if (grepl('TM.-ECL.', inter) | grepl('ECL.-', inter)) {
    l_group <- c(l_group, 'EC')
    next
  }
  if (grepl('ICL.-', inter)| grepl('TM.-ICL.', inter)) {
    l_group <- c(l_group, 'IC')
    next
  }
  if (grepl('TM.-CT', inter) | grepl('CT-', inter)) {
    l_group <- c(l_group, 'IC')
    next
  }
  if (grepl('NT-', inter) | grepl('TM.-NT', inter)) {
    l_group <- c(l_group, 'EC')
    next
  }
}
data_asys$group <- l_group
data_asys$group <- as.factor(data_asys$group)


values_sys <- apply(t(data_inter_clas_sys), 1, function(x) sum(x))
interactions_sys <- rownames(t(data_inter_clas_sys))
data_sys <- data.frame(values_sys, interactions_sys)
data_sys <- data_sys[rev(order(data_sys$values_sys)),]
l_group <- c()
for (inter in data_sys$interactions_sys){
  if (grepl('TM.-TM.', inter)) {
    l_group <- c(l_group, 'TM')
    next
  }
  if (grepl('TM.-ECL.', inter) | grepl('ECL.-', inter)) {
    l_group <- c(l_group, 'EC')
    next
  }
  if (grepl('ICL.-', inter)| grepl('TM.-ICL.', inter)) {
    l_group <- c(l_group, 'IC')
    next
  }
  if (grepl('TM.-CT', inter) | grepl('CT-', inter)) {
    l_group <- c(l_group, 'IC')
    next
  }
  if (grepl('NT-', inter) | grepl('TM.-NT', inter)) {
    l_group <- c(l_group, 'EC')
    next
  }
}
data_sys$group <- l_group
data_sys$group <- as.factor(data_sys$group)

# interactions_sys <- c("TM1-TM1", "TM5-TM5", "ICL2-ICL2", "CT-CT", "TM4-TM5",   "TM4-TM4",   "TM1-TM2",
#                       "ICL3-ICL3", "TM5-TM6", "ECL2-ECL2", "ECL2-TM5",  "ICL2-ICL3", "NT-ECL2",
#                       "ECL3-ECL3", "ECL1-TM1",  "TM1-TM4", "TM2-TM2",  "TM2-ECL1",  "NT-NT", "NT-ECL1",
#                       "ECL1-ECL1", "NT-TM7",  "TM1-TM7", "TM3-TM4", "NT-ECL3", "CT-ICL3", "ECL2-ECL3",
#                       "ICL1-ICL1", "TM5-ICL3",  "CT-TM4",    "TM4-TM7",   "TM4-ICL1", "TM7-CT", "ICL1-CT",
#                       "TM1-CT",    "ECL1-ECL2", "TM3-ICL2",  "TM4-ECL2",  "CT-ICL2",   "TM5-ECL3",
#                       "ICL2-TM4",  "TM7-ECL2",  "ECL1-TM3",  "ECL2-TM6" , "TM4-ECL1",  "TM6-TM6" ,
#                       "TM2-TM4" ,  "TM1-TM5" , "TM5-ICL2",  "ICL1-ICL2", "TM5-TM7",   "TM4-TM6",
#                       "ICL3-TM6",  "TM5-ECL1",  "TM7-TM7",   "ICL1-TM5", "TM3-TM3",   "ECL2-TM2",
#                       "ECL3-TM3",  "TM7-ICL2", "TM7-ECL3",  "ICL3-TM3",  "ICL3-TM4" , "TM2-TM3", "ECL1-ECL3", "ECL1-TM7",
#                       "TM2-TM7",   "TM1-TM3",   "TM1-ICL1", "TM1-ECL2",  "NT-TM5",    "TM6-ECL3", "TM6-TM7", "ICL3-TM7",
#                       "NT-TM3",    "TM4-NT",    "TM3-ECL2",  "CT-TM5",    "TM2-TM5",   "ICL1-ICL3",
#                       "TM1-ICL3",  "TM1-ECL3",  "TM1-TM6",  "NT-TM6",    "NT-TM1" ,   "ECL1-TM6" , "CT-TM3" ,   "CT-TM6",
#                       "ICL2-TM1",  "NT-TM2" ,   "TM3-TM7"  )
# 
# interactions_asys <- c("TM1-TM4" , "TM1-TM5",  "TM1-TM1",  "TM4-TM7",  "TM1-TM2",  "TM1-TM6",  "CT-ICL2",
#                        "ICL1-ICL2", "ECL3-ECL3", "TM2-TM5",  "ECL1-ECL2", "ECL1-ECL3", "TM4-ICL1", "ICL1-ICL3",
#                        "TM1-TM7", "TM5-TM5",  "ICL2-ICL3", "NT-NT", "ICL2-ICL2", "CT-CT",  "CT-TM4", "CT-ICL3",
#                        "TM4-TM5",  "ECL1-TM1",  "TM4-ECL1",  "TM1-ECL3",  "ICL2-TM1",  "NT-ECL3",  "NT-TM1",
#                        "TM4-NT",  "ECL2-TM2", "TM5-ICL2", "TM1-TM3", "NT-TM2",  "ECL2-ECL3", "ECL1-TM7",
#                        "NT-ECL1", "TM1-ECL2", "TM1-ICL3", "TM7-TM7", "TM3-ICL2", "TM4-ECL2", "ICL2-TM4",
#                        "ICL1-TM5",  "TM3-TM3",  "TM1-CT",  "ECL3-TM3",  "TM7-ICL2", "TM7-ECL2",  "TM7-ECL3",
#                        "ICL3-TM3",  "ICL3-TM4", "ECL1-TM3", "TM2-TM3",  "TM2-TM7", "TM1-ICL1", "NT-TM5",
#                        "TM6-ECL3",  "TM6-TM7",  "ICL3-TM7", "TM6-TM6",  "ICL3-TM6", "TM5-TM6", "ECL2-TM6",
#                        "TM2-ECL1",  "TM2-TM2",  "TM7-CT", "ECL1-ECL1", "NT-TM3",  "NT-ECL2", "TM3-ECL2",
#                        "CT-TM5", "TM5-ECL1", "TM3-TM4", "TM2-TM4", "NT-TM6", "NT-TM7", "TM5-ECL3", "ECL2-TM5",
#                        "TM4-TM4", "ICL3-ICL3", "TM5-ICL3",  "ICL1-ICL1", "ICL1-CT", "ECL1-TM6", "CT-TM3",
#                        "CT-TM6", "TM4-TM6", "TM3-TM7", "TM5-TM7",  "ECL2-ECL2")
# # 
interactions_sys <- c("TM1-TM1", "TM5-TM5", "ICL2-ICL2", "CT-CT", "TM4-TM5", "TM5-TM4",   "TM4-TM4",  "TM2-TM1", "TM1-TM2",
                      "ICL3-ICL3", "TM5-TM6", "TM6-TM5", "ECL2-ECL2", "ECL2-TM5", "TM5-ECL2", "ICL3-ICL2", "ICL2-ICL3", "ECL2-NT",
                      "NT-ECL2",
                      "ECL3-ECL3", "ECL1-TM1", "TM1-ECL1",  "TM1-TM4", "TM4-TM1", "TM2-TM2", "ECL1-TM2", "TM2-ECL1",  "NT-NT", "ECL1-NT", "NT-ECL1")
interactions_asys <- c("TM1-TM4" , "TM4-TM1", "TM5-TM1", "TM1-TM5", "TM7-TM4", "TM4-TM7",
                       "TM6-TM1", "TM1-TM6",  "ICL2-CT", "CT-ICL2",
                       "ICL1-ICL2", "ICL2-ICL1", "ECL3-ECL3", "TM5-TM2", "TM2-TM5", "ECL2-ECL1", "ECL1-ECL2", "ECL3-ECL1",
                       "ECL1-ECL3", "TM4-ICL1", "ICL1-TM4", "ICL1-ICL3", "ICL3-ICL1", "TM7-TM1", "ICL3-ICL2",
                       "TM1-TM7", "TM5-TM5",  "ICL2-ICL3", "NT-NT", "ICL2-ICL2", "CT-TM4", "TM4-CT", "CT-ICL3", "TM4-TM5",
                       "ICL3-CT", "TM5-TM4")

data_sys_20_start <- data_sys[data_sys$interactions_sys %in% interactions_sys,]
total <- sum(data_sys_20_start$values_sys)
data_sys_20_start$percent <- sapply(data_sys_20_start$values_sys, function(x) {x/total*100})
data_asys_20_start <- data_asys[data_asys$interactions_asys %in% interactions_asys,]
total <- sum(data_asys_20_start$values_asys)
data_asys_20_start$percent <- sapply(data_asys_20_start$values_asys, function(x) {x/total*100})

data_sys_20_end <- data_sys[data_sys$interactions_sys %in% interactions_sys,]
total <- sum(data_sys_20_end$values_sys)
data_sys_20_end$percent <- sapply(data_sys_20_end$values_sys, function(x) {x/total*100})
data_asys_20_end <- data_asys[data_asys$interactions_asys %in% interactions_asys,]
total <- sum(data_asys_20_end$values_asys)
data_asys_20_end$percent <- sapply(data_asys_20_end$values_asys, function(x) {x/total*100})


labels_sys <- c("TM1-TM1", "ECL2-ECL2", "CT-CT", "ICL3-ICL3", "TM4-TM5", "NT-NT", "TM5-TM5", "ICL2-ICL2", "TM4-TM4",  "ECL2-TM5", 
                      "TM1-TM2", "NT-ECL2", "TM5-TM6", "TM4-TM1", "ECL3-ECL3", "ECL1-TM1", "ICL2-ICL3",  "NT-ECL1", "TM2-ECL1", "TM2-TM2")
labels_asys <- c("ICL1-ICL3", "CT-ICL2", "CT-ICL3", "ICL2-ICL3", "ECL2-ECL1", "TM7-TM4", "ICL2-ICL1", "TM1-TM5", "CT-TM4", "NT-NT", "TM4-TM1", 
                       "TM1-TM6", "ECL3-ECL3", "TM4-TM5", "TM4-ICL1", "ECL1-ECL3", "TM5-TM2", "TM1-TM7", "TM5-TM5", "ICL2-ICL2")

data_sys_20_start$interactions_sys <- factor(data_sys_20_start$interactions_sys, levels=labels_sys)
data_sys_20_start <-data_sys_20_start[order(data_sys_20_start$interactions_sys),]
data_asys_20_start$interactions_asys <- factor(data_asys_20_start$interactions_asys, levels=labels_asys)
data_asys_20_start <-data_asys_20_start[order(data_asys_20_start$interactions_asys),]
data_asys_20_start <- data_asys_20_start[1:20,]
data_sys_20_start <- data_sys_20_start[1:20,]

data_sys_20_end$interactions_sys <- factor(data_sys_20_end$interactions_sys, levels=labels_sys)
data_sys_20_end <-data_sys_20_end[order(data_sys_20_end$interactions_sys),]
data_asys_20_end$interactions_asys <- factor(data_asys_20_end$interactions_asys, levels=labels_asys)
data_asys_20_end <-data_asys_20_end[order(data_asys_20_end$interactions_asys),]
data_asys_20_end <- data_asys_20_end[1:20,]

labels_sys <- c("TM1-TM1", "ECL2-ECL2", "CT-CT", "ICL3-ICL3", "TM4-TM5", "NT-NT", "TM5-TM5", "ICL2-ICL2", "TM4-TM4",  "ECL2-TM5", 
                "TM1-TM2", "ECL2-NT", "TM5-TM6", "TM1-TM4", "ECL3-ECL3", "ECL1-TM1", "ICL2-ICL3",  "ECL1-NT", "ECL1-TM2", "TM2-TM2")

labels_asys <- c("ICL1-ICL3", "ICL2-CT", "ICL3-CT","ICL2-ICL3", "ECL1-ECL2", "TM4-TM7", "ICL1-ICL2", "TM1-TM5", "CT-TM4", "NT-NT", "TM1-TM4", 
                 "TM1-TM6", "ECL3-ECL3", "TM4-TM5", "ICL1-TM4", "ECL1-ECL3", "TM2-TM5", "TM1-TM7", "TM5-TM5", "ICL2-ICL2")
 
data_sys_20_start$interactions_sys <- labels_sys
data_asys_20_start$interactions_asys <- labels_asys 

data_sys_20_end$interactions_sys <- labels_sys
data_asys_20_end$interactions_asys <- labels_asys

# data_all_ord_f$group
# Order data:
data <- data_asys_20_end
# data <- data %>% arrange(group, desc(values_sys))

# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar)
data=rbind(data, to_add)
data=data %>% arrange(group)
data$id=seq(1, nrow(data))

# Get the name and the y position of each label
label_data=data
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data=data %>%
  group_by(group) %>%
  summarize(start=min(id), end=max(id) - empty_bar) %>%
  rowwise() %>%
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

# Make the plot
p = ggplot(data, aes(x=as.factor(id), y=round(percent), fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=round(percent), fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 16, xend = start, yend = 16), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 12, xend = start, yend = 12), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 8, xend = start, yend = 8), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 4, xend = start, yend = 4), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),5), y = c(0,4,8,12,16), label = c("0", "4", "8", "12", "16") , color="black", size=6 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=round(percent), fill=group), stat="identity", alpha=0.5) +
  ylim(-30,30) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=round(percent)+5, label=interactions_asys, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=6, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -1, xend = end, yend = -1), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -5, label=group), hjust=c(1,1,0), colour = "black", alpha=0.8, size=8, fontface="bold", inherit.aes = FALSE)

p

ggsave(p, filename = "/home/adrian/Escritorio/RESULTS/interactions_interfaces_asys_end_20.jpg", width = 10, height = 10, dpi = "retina")


max_inter_clas_tm <- max_inter_clas_tm[, colSums(max_inter_clas_tm != 0) > 0]
data <- melt(max_inter_clas_tm)

data_inter_clas_tm_asys <- data_inter_clas_tm_asys[, colSums(data_inter_clas_tm_asys != 0) > 0]
data_inter_clas_tm_sys <- data_inter_clas_tm_sys[, colSums(data_inter_clas_tm_sys != 0) > 0]

data_inter_clas_tm_asys_c <- apply(data_inter_clas_tm_asys, 1, function(x) replace(x, x>0, 1))
data_inter_clas_tm_sys_c <- apply(data_inter_clas_tm_sys, 1, function(x) replace(x, x>0, 1))

library(pheatmap)
library(RColorBrewer)
# ggplot( data, aes(Var2, Var1) ) +
#   geom_tile()+
#   geom_tile(colour="white", size=0.25, show.legend=F, aes(fill=value) )+
#   scale_fill_gradient2(low="blue", mid="#ff6666", high = "#66ff66") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
# Generte data (modified the mydf slightly)
col1 <- brewer.pal(12, "Set3")

mydf <- data.frame(row.names = rownames(max_inter_clas_tm), dimers =rownames(max_inter_clas_tm), interface = c(rep("SYMMETRIC", 6), rep("ASYMMETRIC", 2), rep("SYMMETRIC", 1), rep("ASYMMETRIC", 1), rep("SYMMETRIC", 1), rep("ASYMMETRIC", 1), 
                                                                          rep("SYMMETRIC", 12), rep("ASYMMETRIC", 1), rep("SYMMETRIC", 2), rep("ASYMMETRIC", 5), rep("SYMMETRIC", 2), rep("ASYMMETRIC", 3),
                                                                          rep("SYMMETRIC", 4), rep("ASYMMETRIC", 1), rep("SYMMETRIC", 8), rep("ASYMMETRIC", 1), rep("SYMMETRIC", 7)))
mydf <- mydf[order(mydf$interface), ]
mydf$dimers <- NULL
data <- t(max_inter_clas_tm)
data <- data[,rownames(mydf)]
s <- pheatmap(data, cluster_cols = F, cluster_rows = F, 
              annotation_col = mydf, gaps_col = c(15), 
              fontsize_row = 8,  color = c("white", "black"), 
              legend_breaks = c(1,0), legend = F)

ggsave(s, filename = "/home/adrian/Escritorio/RESULTS/interficies_distribution_f.jpg", width = 10, height = 10, dpi = "retina")

info_all
plot(resum$RMSD_MEAN, resum$DEGREE_MEAN, col= resum$DIMER)
abline(lm(resum$DEGREE_MEAN~resum$RMSD_MEAN))
text(resum$RMSD_MEAN, resum$DEGREE_MEAN, labels=resum$DIMER, cex= 0.5, pos = 2)

plot(resum$RMSD_MEAN, resum$DIST_MEAN, col= resum$DIMER)
abline(lm(resum$DIST_MEAN~resum$RMSD_MEAN))
text(resum$RMSD_MEAN, resum$DIST_MEAN, labels=resum$DIMER, cex= 0.5, pos = 2)

plot(resum$DIST_MEAN, resum$DEGREE_MEAN, col= resum$DIMER)
abline(lm(resum$DEGREE_MEAN~resum$DIST_MEAN))
text(resum$DIST_MEAN, resum$DEGREE_MEAN, labels=resum$DIMER, cex= 0.5, pos = 2)
