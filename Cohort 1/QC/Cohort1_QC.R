# Functions  ------------------------------------------------------------

####Load Libraries====
libs_load <- function(x){
  for( i in x ){
    print(paste0("Checking for library: ", i))
    if(require( i , character.only = TRUE ) ){
      print(paste0(i, " already installed. Loading now"))
    }
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      print(paste0(i, " not installed. Trying CRAN for install."))
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
      paste0(i, " installed and loaded successfully")
    }
    if ( ! require(i, character.only=TRUE) ) {
      paste0(i," could not be installed from CRAN. Trying Bionconductor....")
      BiocManager::install(i)
      require( i , character.only = TRUE )
      paste0(i, " installed and loaded successfully")
    }
    if ( ! require(i, character.only=TRUE) ) {
      paste0(i, "could not be installed. Check manually")
    }
    #  Load package after installing
  }
}
#Usage:
#libs_load(c("...", "..."))


####Create PC1 PC2 PC3 PCA Plot Strip (FactoExtra)====
create_pca_plots <- function(pca_result, df, label_col, title) {
  # Convert label column to factor
  label_factor <- as.factor(df[[label_col]])
  
  # Create PCA plots for different axes
  pca_plot1 <- fviz_pca_ind(pca_result, habillage = label_factor,title="", axes = c(1, 2))
  pca_plot2 <- fviz_pca_ind(pca_result, habillage = label_factor,title="", axes = c(1, 3))
  pca_plot3 <- fviz_pca_ind(pca_result, habillage = label_factor,title="", axes = c(2, 3))
  
  # Arrange the plots in a 3x1 grid with a title
  combined_plot <- grid.arrange(pca_plot1, pca_plot2, pca_plot3, nrow = 1, top = textGrob(title, gp = gpar(fontsize = 14, font = 2)))
  
  return(combined_plot)
}

# Example usage:
# Assuming 'NArm_Merged_PCA' is your PCA result and 'NArm_Merged_df' is your data frame
# result <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "Age_bin", "PCA Plot of Age Bins")
# print(result)

####Calculate RLE function====
calcRLE <- function(E, log.data = TRUE, metadata.cols = NULL) {
  if (!inherits(E, what = "data.frame")) {
    stop("E must be a data.frame")
  }
  
  # Separate metadata columns if specified
  if (!is.null(metadata.cols)) {
    metadata <- E[, metadata.cols, drop = FALSE]
    E <- E[, !colnames(E) %in% metadata.cols]
  } else {
    metadata <- NULL
  }
  
  # Log the data if specified
  if (log.data == TRUE) {
    message("Logging the data")
    E <- log2(E + 1)
  }
  
  # Calculate the median value of log2 counts for each gene (column)
  g.medians <- apply(E, 2, FUN = median, na.rm = TRUE)
  
  # For each gene (column), subtract its median
  E.new <- sweep(E, 2, g.medians, FUN = "-")
  
  # Re-add metadata columns if they were separated
  if (!is.null(metadata)) {
    E.new <- cbind(metadata, E.new)
  }
  
  return(E.new)
}

# Example usage:
# Assuming 'gene_expression_data' is your data frame and 'c("SampleID", "Age", "Gender")' are your metadata columns
# result <- calcRLE(gene_expression_data, log.data = TRUE, metadata.cols = c("SampleID", "Age", "Gender"))

####Create RLE Boxplot (GGplot)====
create_rle_boxplot <- function(rle_data, metadata_cols, color_col,title) {
  # Check if the metadata and color columns exist
  if (!all(metadata_cols %in% colnames(rle_data))) {
    stop("One or more specified metadata columns do not exist in the data frame.")
  }
  if (!(color_col %in% metadata_cols)) {
    stop("The specified color column must be one of the metadata columns.")
  }
  
  # Convert RLE data to long format for plotting
  rle_long <- rle_data %>%
    pivot_longer(cols = -all_of(metadata_cols), names_to = "Gene_Symbol", values_to = "NPX")
  
  # Create plotting order variables for each sample, plot in order of global medians
  plot_order <- rle_long %>%
    group_by(SampleID) %>%
    summarise(median_NPX = median(NPX, na.rm = TRUE)) %>%
    arrange(desc(median_NPX)) %>%
    pull(SampleID)
  
  # Calculate median for horizontal line on plot
  overall_median <- median(na.omit(rle_long$NPX))
  
  # Create the boxplot
  plot <- ggplot(rle_long, aes(x = SampleID, y = NPX, fill = !!sym(color_col))) +
    geom_boxplot(outlier.colour = "black", outlier.size = 0.05) +
    scale_x_discrete(limits = plot_order) +
    geom_hline(yintercept = overall_median, linetype = "dashed", colour = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 5, hjust = 1)) +
    ylab("Relative Log Abundance") +
    xlab("") +
    ggtitle(title)
  
  return(plot)
}
# Example usage:
# Assuming 'gene_expression_data' is your data frame and 'c("SampleID", "PlateID")' are your metadata columns
# rle_data <- calcRLE(gene_expression_data, log.data = TRUE, metadata.cols = c("SampleID", "PlateID"))
# plot <- create_rle_boxplot(rle_data, metadata_cols = c("SampleID", "PlateID"), color_col = "PlateID",title="RLE - Inflammation Panel - Plate ID)
# print(plot)

# Load Libraries  ------------------------------------------------------------
libs_load(c("OlinkAnalyze", "tidyverse", "FactoMineR", "factoextra", "gridExtra", "RColorBrewer", "grid", "this.path"))

# Set up file paths & output directories  ------------------------------------------------------------

# Root
setwd(this.path::here())
QC_ROOT= this.path::here()

# Input
INPUT_DIR<-paste0(QC_ROOT, "/Input_Data")

# Output directories
OUTPUT_DIR <- paste0(QC_ROOT, "/Output")
# Check if the Output directory exists (Avoids overwriting)
if (!dir.exists(OUTPUT_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(OUTPUT_DIR, recursive = TRUE)
  message("Directory created: ", OUTPUT_DIR)
} else {
  message("Directory already exists: ", OUTPUT_DIR)
}

#Directory for outputting Detectability Results
DETECT_DIR <- paste0(OUTPUT_DIR, "/Detectability")
if (!dir.exists(DETECT_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(DETECT_DIR, recursive = TRUE)
  message("Directory created: ", DETECT_DIR)
} else {
  message("Directory already exists: ", DETECT_DIR)
}

#Directory for outputting PCAs
PCA_DIR <- paste0(OUTPUT_DIR, "/PCA")
if (!dir.exists(PCA_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(PCA_DIR, recursive = TRUE)
  message("Directory created: ", PCA_DIR)
} else {
  message("Directory already exists: ", PCA_DIR)
}

#Directory for outputting RLE plots
RLE_DIR <- paste0(OUTPUT_DIR, "/RLE")
if (!dir.exists(RLE_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(RLE_DIR, recursive = TRUE)
  message("Directory created: ", RLE_DIR)
} else {
  message("Directory already exists: ", RLE_DIR)
}

#Directory for outputting Post QC data files
POST_QC_DIR <- paste0(OUTPUT_DIR, "/Post_QC_Datasets")
if (!dir.exists(POST_QC_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(POST_QC_DIR, recursive = TRUE)
  message("Directory created: ", POST_QC_DIR)
} else {
  message("Directory already exists: ", POST_QC_DIR)
}

# Read in Raw Data and filter out controls --------------------------------

#Read in wide format NPX Target data using read_NPX function from Olink Analyze. Also read in Olink assay annotation
raw_NPX<-read_NPX(paste0(INPUT_DIR,'/Cohort1_INF&CDM_Raw.xlsx'))

# Read in and Subject Metadata
Cohort1_metadata <- read.csv(paste0(INPUT_DIR,'/Cohort1_metadata.csv'))

#Read in additional assay information and add to data frame
Olink_Assays <- read.csv(paste0(INPUT_DIR,'/Olink_Assays.csv'))
raw_NPX<-raw_NPX %>% left_join(Olink_Assays, by = "Assay")

#Remove unneeded control samples:
filt_NPX <- raw_NPX %>% filter(!(SampleID %in% c("PP1","PP2","PP3","PP4",'IPC 1', 'IPC 2', 'IPC 3', 'Neg 1', 'Neg 2', 'Neg 3'))) #Pooled Plasma, Interplate controls (IPC) and negative controls (Neg) used in Olink QC but not needed here
#Remove control assays
filt_NPX <- filt_NPX %>% filter(!str_detect(Assay, "Ctrl"))

#Check number of samples and assays in filtered data frame
filt_NPX %>% distinct(SampleID) %>% nrow()
filt_NPX%>% distinct(Assay) %>% nrow()
filt_NPX%>% distinct(Gene_Symbol) %>% nrow()

# Check Sample QC warnings (based on Olink internal controls) --------------------------------

#If a sample deviates by more than +/-0.3 NPX from the median of all samples it will be flagged as a warning for that panel
QC_warns<-filt_NPX %>% filter(QC_Warning == 'Warning') #filter rows which have a QC warning
QC_warns %>% group_by(SampleID, Panel,QC_Warning) %>% summarise(count = n(), .groups = 'drop')#Print a summary

#3 samples had QC warnings on the Inflammation assay: HC13, T28, T84. 
#NPX values for T28 are already NA.Need to do this for others:
filt_NPX[filt_NPX$SampleID=="HC13" & filt_NPX$QC_Warning=="Warning" & filt_NPX$Panel=="Olink Inflammation","NPX"]<-NA
filt_NPX[filt_NPX$SampleID=="T84" & filt_NPX$QC_Warning=="Warning" & filt_NPX$Panel=="Olink Inflammation","NPX"]<-NA

# Check detectability % (NPX > LOD) --------------------------------

#Checking protein detection rate (>LLOD) after filtering
filt_NPX <- filt_NPX %>% mutate(Gr_LOD = ifelse(is.na(NPX), "No", ifelse(NPX > LOD, "Yes", "No"))) # Create new col called Gr_LOD: Yes if >LOD, No if NA or <LOD
detection_res<-filt_NPX %>% group_by(Assay, Gr_LOD,Panel,Gene_Symbol) %>% summarise(count = n(), .groups = 'drop') #summarise to protein level (long form)
detection_res <- detection_res %>% spread(key = Gr_LOD, value = count, fill = 0) #swap to wide form for nice summary table
detection_res<-detection_res %>% mutate(Det_Percent = round((Yes / 166) * 100,0)) %>% arrange(Det_Percent) #add detection percent column
detection_res <- detection_res %>% mutate(Detected = ifelse(Det_Percent >= 25, "yes", "no")) 
#return the number of proteins with values <25% detection (threshold used for this study)
non_detected<-detection_res %>% filter(Det_Percent < 25) %>% pull(Assay)
length(non_detected)

####Plot bottom 30 proteins detection % ====
per_detected_plot<-ggplot(detection_res[1:30,],aes(x=Gene_Symbol,y=Det_Percent,fill=Panel))+
  geom_bar(stat="identity")+
  scale_x_discrete(limits = detection_res$Gene_Symbol[1:30]) + 
  scale_y_continuous(limits =c(0,100) ) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=6),
        axis.text.y = element_text(size=8))+
  labs(y = "% Samples Detected")+
  ggtitle("% Assay Detectability - lowest 30 ") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "red")
per_detected_plot
#Save plot and detection results table
ggsave(filename="Bottom 30 Detection Percentage.png", path=DETECT_DIR, plot=per_detected_plot,width=15,height=12,units="cm")
write.csv(detection_res, paste0(DETECT_DIR,"/Cohort 1 Assay Detectability.csv"))

####Delta LOD Analysis====
#Examine pattern of deviation from LOD to check for systematic issues
#Add Delta LOD calculation and detectability status (>=25%) to DF
filt_NPX$Delta_LOD<-filt_NPX$NPX-filt_NPX$LOD
filt_NPX <- filt_NPX %>% left_join(detection_res %>% select(Gene_Symbol, Detected),by ="Gene_Symbol")

#create plotting order variable to plot sample in order of Delta LOD medians
delta_order = filt_NPX %>% group_by(Gene_Symbol) %>% summarise(median_delta = median(na.omit(Delta_LOD))) %>% arrange(desc(median_delta)) %>% pull(Gene_Symbol)

#Create plot and save
Delta_LOD_plot = ggplot(filt_NPX, aes(x = Gene_Symbol, y = Delta_LOD, fill = Detected)) + 
  geom_boxplot(outlier.colour = "black", outlier.size = 0.05, na.rm = TRUE, lwd = 0.5) +  # Set line width of boxes to 0.5
  geom_jitter(width = 0.15, size = 0.15, na.rm = TRUE,color="#6a7178",alpha=0.5) +  # Add jitter to show all points
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") + 
  scale_x_discrete(limits = delta_order) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, size = 5, hjust = 1)) + 
  ylab("NPX deviation from LOD") + 
  xlab("") + 
  ggtitle("Delta LOD - NPX deviation from LOD") +
  labs(fill = ">=25% Detection")
Delta_LOD_plot
#pattern of deviation below LOD appears to have gaussian distribution for marginal proteins that are classed as "detected"
ggsave("Delta LOD Plot.png", plot=Delta_LOD_plot, path=DETECT_DIR, width=40, height=14, units="cm")

#Remove proteins wth <25% detection and save as new DF
detected_NPX<- filt_NPX%>% filter(!Assay %in% non_detected)
detected_NPX %>% distinct(Assay) %>% nrow()

# QC PCAs --------------------------------

####Plate Assignment PCAs####

#split into CDM and Inflamm datasets
INF_df<- detected_NPX %>% select(Gene_Symbol,SampleID, PlateID, NPX) %>% filter(str_detect(PlateID, "INF")) %>% spread(key = Gene_Symbol, value = NPX, fill = NA)
CDM_df<- detected_NPX %>% select(Gene_Symbol,SampleID, PlateID, NPX) %>% filter(str_detect(PlateID, "CDM")) %>% spread(key = Gene_Symbol, value = NPX, fill = NA)

INF_df<-as.data.frame(na.omit(INF_df))#for PCAs, must remove cases with NAs
CDM_df<-as.data.frame(na.omit(CDM_df))#for PCAs, must remove cases with NAs
dim(INF_df)#Check n remaining cases
dim(CDM_df)#check n remaining cases
rownames(INF_df)<-INF_df$SampleID #setrownames for PCA labels
rownames(CDM_df)<-CDM_df$SampleID #setrownames for PCA labels

#Create PCAs and plots
INF_PCA<-PCA(INF_df[,3:ncol(INF_df)],graph=FALSE,scale.unit = TRUE)
CDM_PCA<-PCA(CDM_df[,3:ncol(CDM_df)],graph=FALSE,scale.unit = TRUE)

#create PCA plots and save as strips
INF_plt_PCA <- create_pca_plots(INF_PCA, INF_df, "PlateID", "Inflammation Panel Plate ID PCA")
ggsave("Inflammation Plate ID PCA.png", plot=INF_plt_PCA,path=PCA_DIR , width =80, height=20,units="cm")
CDM_plt_PCA <- create_pca_plots(CDM_PCA, CDM_df, "PlateID", "Cardiometabolic Panel Plate ID PCA")
ggsave("Cardiometabolic Plate ID PCA.png",plot=CDM_plt_PCA, path=PCA_DIR ,width =80, height=20,units="cm")

####Subject Characteristics PCAs####

#Create wide NPX DF for merging with clinical parameters
wide_npx<-detected_NPX %>% select(Gene_Symbol,SampleID, NPX) %>% spread(key = Gene_Symbol, value = NPX, fill = NA)
wide_npx<-Cohort1_metadata %>% left_join(wide_npx,by="SampleID")#Join metadata and wide_npx df
wide_npx<-as.data.frame(wide_npx) #convert tibble to DF
rownames(wide_npx)<-wide_npx$SampleID #set rownames for PCA labels
dim(wide_npx)

wide_npx_na.rm<-na.omit(wide_npx)#DF without NA cases for PCA
dim(wide_npx_na.rm)

#Create PCA object
all_PCA<-PCA(wide_npx_na.rm[,8:ncol(wide_npx_na.rm)],graph=FALSE,scale.unit = TRUE )

#Create PCA plots of axes 1:3 and save

#Age
age_PCA <- create_pca_plots(all_PCA, wide_npx_na.rm, "Age_bin", "INF & CDM Age Bin PCA")
ggsave("Age Bins PCA.png", path=PCA_DIR, plot=age_PCA,width =50, height=20,units="cm")
#Sex
sex_PCA <- create_pca_plots(all_PCA, wide_npx_na.rm, "Sex", "INF & CDM Sex PCA")
ggsave("Sex PCA.png", path=PCA_DIR, plot=sex_PCA,width =50, height=20,units="cm")
#Ethnicity
ethnicity_PCA <- create_pca_plots(all_PCA, wide_npx_na.rm, "Ethnicity", "INF & CDM Ethnicity PCA")
ggsave("Ethnicity PCA.png", path=PCA_DIR, plot=ethnicity_PCA,width =50, height=20,units="cm")
#Group
group_PCA <- create_pca_plots(all_PCA, wide_npx_na.rm, "Group", "INF & CDM Group PCA")
ggsave("Group PCA.png", path=PCA_DIR, plot=group_PCA,width =50, height=20,units="cm")
#Subgroup
subgroup_PCA <- create_pca_plots(all_PCA, wide_npx_na.rm, "Subgroup", "INF & CDM Subgroup PCA")
ggsave("Activity Subgroup PCA.png", path=PCA_DIR, plot=subgroup_PCA,width =50, height=20,units="cm")
#Freezer Time
freezer.time_PCA1<-fviz_pca_ind(all_PCA,col.ind=wide_npx_na.rm$Freezer_Time,axes=c(1,2),title="")+labs(color="Freezer Time")
freezer.time_PCA2<-fviz_pca_ind(all_PCA,col.ind=wide_npx_na.rm$Freezer_Time,axes=c(1,3),title="")+labs(color="Freezer Time")
freezer.time_PCA3<-fviz_pca_ind(all_PCA,col.ind=wide_npx_na.rm$Freezer_Time,axes=c(2,3),title="")+labs(color="Freezer Time")
ggsave("Freezer Time PCA.png", path=PCA_DIR,plot=grid.arrange(freezer.time_PCA1,freezer.time_PCA2,freezer.time_PCA3,nrow=1),width =50, height=20,units="cm")


# QC RLE plots ------------------------------------------------------------

####Plate assignment RLE Plots====
#create RLE DF on INF and CDM dfs separately
INF_RLE<-calcRLE(INF_df, log.data = FALSE, metadata.cols = c("SampleID", "PlateID"))
CDM_RLE<-calcRLE(CDM_df, log.data = FALSE, metadata.cols = c("SampleID", "PlateID"))
#Run the RLE create plot function & Save
INF_plate_RLEp <- create_rle_boxplot(INF_RLE, metadata_cols = c("SampleID", "PlateID"), color_col = "PlateID",title="RLE - Inflammation Panel")
CDM_plate_RLEp <- create_rle_boxplot(CDM_RLE, metadata_cols = c("SampleID", "PlateID"), color_col = "PlateID",title="RLE - Cardiometabolic Panel")
ggsave("Inflammation Panel Plate ID RLE Plot.png", path=RLE_DIR , plot =INF_plate_RLEp, width=40, height=14, units="cm")
ggsave("Cardiometabolic Panel Plate ID RLE Plot.png", path=RLE_DIR, plot =CDM_plate_RLEp, width=40, height=14, units="cm")

####Subject Characteristics RLE Plots====

#Create RLE DF on combined data
#Note NAs were not removed, so some cases will have fewer proteins as part of median
all_RLE<-calcRLE(wide_npx, log.data = FALSE, metadata.cols = colnames(wide_npx[,1:7]))
#Age
age_RLEp <- create_rle_boxplot(all_RLE, metadata_cols = colnames(wide_npx[,1:7]), color_col = "Age_bin",title="RLE - Age Bin")
ggsave("Age Bin RLE Plot.png", path=RLE_DIR, plot =age_RLEp, width=40, height=14, units="cm")
#Sex
sex_RLEp <- create_rle_boxplot(all_RLE, metadata_cols = colnames(wide_npx[,1:7]), color_col = "Sex",title="RLE - Sex")
ggsave("Sex RLE Plot.png", path=RLE_DIR, plot=sex_RLEp, width=40, height=14, units="cm")
#Ethnicity
ethnicity_RLEp <- create_rle_boxplot(all_RLE, metadata_cols = colnames(wide_npx[,1:7]), color_col = "Ethnicity",title="RLE - Ethnicity")
ggsave("Ethnicity RLE Plot.png", path=RLE_DIR, plot =ethnicity_RLEp, width=40, height=14, units="cm")
#Group
group_RLEp <- create_rle_boxplot(all_RLE, metadata_cols = colnames(wide_npx[,1:7]), color_col = "Group",title="RLE - Group")
ggsave("Group RLE Plot.png", path=RLE_DIR, plot =group_RLEp, width=40, height=14, units="cm")
#Subgroup
subgroup_RLEp <- create_rle_boxplot(all_RLE, metadata_cols = colnames(wide_npx[,1:7]), color_col = "Subgroup",title="RLE - Subgroup")
ggsave("Subgroup RLE Plot.png", path=RLE_DIR, plot =subgroup_RLEp, width=40, height=14, units="cm")
#Freezer Time
freezerT_RLEp <- create_rle_boxplot(all_RLE, metadata_cols = colnames(wide_npx[,1:7]), color_col = "Freezer_Time",title="RLE - Freezer Time (yrs)")
ggsave("Freezer Time RLE Plot.png", path=RLE_DIR, plot =freezerT_RLEp, width=40, height=14, units="cm")

#No clear confounding effects observed in PCAs or RLE plots

# Write out a final post-QC long and simpler wide CSV -----------------------------

#Write out raw NPX data
write.csv(wide_npx, paste0(POST_QC_DIR,"/Cohort 1 Post QC Data.csv"))
#Write out Z-score scaled data
Scaled_npx<-wide_npx
Scaled_npx[,8:165]<-scale(wide_npx[,8:165],scale=TRUE,center=TRUE)
write.csv(Scaled_npx, paste0(POST_QC_DIR,"/Cohort 1 Post QC Z-scores.csv"))


