#Path to input directory with images
INFILE_DIR = ""

#Path to output csv for plotting 
OUTFILE = ""


#File # and condition folder
FILE1 = ""
FILE2 = ""
FILE3 = ""
FILE4 = ""
FILE5 = ""
FILE6 = ""

FOLDER1 = ""
FOLDER2 = ""
FOLDER3 = ""
FOLDER4 = ""
FOLDER5 = ""
FOLDER6 = ""




df_1 <- read.csv(INFILE_DIR + "/" + FOLDER1 + "/" + FILE1 +".csv")
df_1 <- read.csv(INFILE_DIR + "/" + FOLDER2 + "/" + FILE2 +".csv")
df_1 <- read.csv(INFILE_DIR + "/" + FOLDER3 + "/" + FILE3 +".csv")
df_1 <- read.csv(INFILE_DIR + "/" + FOLDER4 + "/" + FILE4 +".csv")
df_1 <- read.csv(INFILE_DIR + "/" + FOLDER5 + "/" + FILE5 +".csv")
df_1 <- read.csv(INFILE_DIR + "/" + FOLDER6 + "/" + FILE6 +".csv")


#FILTER OUT PATCHES WITH low # of segmentations or segmentations of > 50000 pixels
remove_patches_1 <- unique(df_1[df_1$Area > 50000,]$Patch)

remove_patches_1_ <- c()
for (patch in 1:max(df_1$Patch)){
  subset_df <- df_1[df_1$Patch == patch,]
  if (nrow(subset_df) < 45) {
    remove_patches_1_ <- c(remove_patches_1_,patch)
  }
}
remove_patches_1 <- c(remove_patches_1,remove_patches_1_)




remove_patches_2 <- unique(df_2[df_2$Area > 50000,]$Patch)

remove_patches_2_ <- c()
for (patch in 1:max(df_2$Patch)){
  subset_df <- df_2[df_2$Patch == patch,]
  if (nrow(subset_df) < 45) {
    remove_patches_2_ <- c(remove_patches_2_,patch)
  }
}
remove_patches_2 <- c(remove_patches_2,remove_patches_2_)



remove_patches_3 <- unique(df_3[df_3$Area > 50000,]$Patch)

remove_patches_3_ <- c()
for (patch in 1:max(df_3$Patch)){
  subset_df <- df_3[df_3$Patch == patch,]
  if (nrow(subset_df) < 45) {
    remove_patches_3_ <- c(remove_patches_3_,patch)
  }
}
remove_patches_3 <- c(remove_patches_3,remove_patches_3_)




remove_patches_4 <- unique(df_4[df_4$Area > 50000,]$Patch)

remove_patches_4_ <- c()
for (patch in 1:max(df_4$Patch)){
  subset_df <- df_4[df_4$Patch == patch,]
  if (nrow(subset_df) < 45) {
    remove_patches_4_ <- c(remove_patches_4_,patch)
  }
}
remove_patches_4 <- c(remove_patches_4,remove_patches_4_)



remove_patches_5 <- unique(df_5[df_5$Area > 50000,]$Patch)

remove_patches_5_ <- c()
for (patch in 1:max(df_5$Patch)){
  subset_df <- df_5[df_5$Patch == patch,]
  if (nrow(subset_df) < 45) {
    remove_patches_5_ <- c(remove_patches_5_,patch)
  }
}
remove_patches_5 <- c(remove_patches_5,remove_patches_5_)



remove_patches_6 <- unique(df_6[df_6$Area > 50000,]$Patch)

remove_patches_6_ <- c()
for (patch in 1:max(df_6$Patch)){
  subset_df <- df_6[df_6$Patch == patch,]
  if (nrow(subset_df) < 45) {
    remove_patches_6_ <- c(remove_patches_6_,patch)
  }
}
remove_patches_6 <- c(remove_patches_6,remove_patches_6_)




df_1 <- df_1[!df_1$Patch %in% remove_patches_1,]
df_2 <- df_2[!df_2$Patch %in% remove_patches_2,]
df_3 <- df_3[!df_3$Patch %in% remove_patches_3,]
df_4 <- df_4[!df_4$Patch %in% remove_patches_4,]
df_5 <- df_5[!df_5$Patch %in% remove_patches_5,]
df_6 <- df_6[!df_6$Patch %in% remove_patches_6,]




df_1$File <- as.factor(FILE1)
df_2$File <- as.factor(FILE2)
df_3$File <- as.factor(FILE3)
df_4$File <- as.factor(FILE4)
df_5$File <- as.factor(FILE5)
df_6$File <- as.factor(FILE6)
df_1$Folder <- as.factor(FOLDER1)
df_2$Folder <- as.factor(FOLDER2)
df_3$Folder <- as.factor(FOLDER3)
df_4$Folder <- as.factor(FOLDER4)
df_5$Folder <- as.factor(FOLDER5)
df_6$Folder <- as.factor(FOLDER6)
df <- rbind(df_1,df_2,df_3,df_4,df_5,df_6)

#FILTER OUT NOISE
df <- df[df$VCAD < 4000,]
df <- df[df$FULL_DAPI > 1000,]


#Scale by area
df$bg_total_MARKER <- df$Area * #df$BG_{MARKER}    *SUBSTITUTE SPECIFIC MARKER BEING TESTED



write.csv(df_all,OUT_FILE)