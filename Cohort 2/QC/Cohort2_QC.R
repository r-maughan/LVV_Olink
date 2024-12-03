
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


###Create PC1 PC2 PC3 PCA Plot Strip (FactoExtra)====
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

###Calculate RLE function====
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


###Create RLE Boxplot (GGplot)====
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

###Function to remove the batch effect using a linear model====
adjust_batch_effect <- function(data, batch_col = "Batch", metadata_cols = c("SampleID", "PlateID", "Batch")) {
  # Ensure the data is a data frame
  if (!inherits(data, what = "data.frame")) {
    stop("Data must be a data.frame")
  }
  
  # Ensure the specified batch column exists in the data frame
  if (!(batch_col %in% colnames(data))) {
    stop("The specified batch column does not exist in the data frame.")
  }

  # Initialize a list to store the linear model residuals
  lmfit <- vector('list', ncol(data[, !(colnames(data) %in% metadata_cols)]))
  names(lmfit) <- colnames(data[, !(colnames(data) %in% metadata_cols)])

  # Loop through each protein column and fit a linear model
  for (i in seq_along(lmfit)) {
    lmfit[[i]] <- summary(lm(data[[names(lmfit)[i]]] ~ data[[batch_col]]))$residuals
  }

  # Create an empty data frame to store the adjusted values
  adjusted_data <- data.frame(matrix(nrow = nrow(data), ncol = length(names(lmfit))))
  colnames(adjusted_data) <- names(lmfit)
  rownames(adjusted_data) <- rownames(data)

  # Fill the adjusted data frame with the residuals
  for (i in seq_along(lmfit)) {
    adjusted_data[, i] <- lmfit[[i]]
  }

  # Re-add the metadata columns to the adjusted data frame
  adjusted_data <- cbind(data[, metadata_cols], adjusted_data)

  return(adjusted_data)
}
# Example usage: adjusted_data <- adjust_batch_effect(INF_df, batch_col = "Batch", metadata_cols = c("SampleID", "PlateID", "Batch"))
# Assuming 'INF_df' is your data frame

# Libraries  ------------------------------------------------------------

libs_load(c("OlinkAnalyze", "tidyverse", "FactoMineR", "factoextra", "gridExtra", "RColorBrewer", "ggplot2","grid", "this.path"))

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

#Directory for investigating possible batch Inflammation Batch Effect
INF_PRE_DIR <- paste0(OUTPUT_DIR, "/INF_Pre_BatchAdjust")
if (!dir.exists(INF_PRE_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(INF_PRE_DIR, recursive = TRUE)
  message("Directory created: ", INF_PRE_DIR)
} else {
  message("Directory already exists: ", INF_PRE_DIR)
}

#Directory for confirming nil Cardiometabolic Batch Effect
CDM_PRE_DIR <- paste0(OUTPUT_DIR, "/CDM_Pre_BatchAdjust")
if (!dir.exists(CDM_PRE_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(CDM_PRE_DIR, recursive = TRUE)
  message("Directory created: ", CDM_PRE_DIR)
} else {
  message("Directory already exists: ", CDM_PRE_DIR)
}

#Directory for INF & CDM Combined QC
BOTH_POST_DIR <- paste0(OUTPUT_DIR, "/INF&CDM_post_BatchAdjust")
if (!dir.exists(BOTH_POST_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(BOTH_POST_DIR, recursive = TRUE)
  message("Directory created: ", BOTH_POST_DIR)
} else {
  message("Directory already exists: ", BOTH_POST_DIR)
}

##INF & CDM Combined QC PCAs
POST_PCA_DIR <- paste0(BOTH_POST_DIR, "/PCA")
if (!dir.exists(POST_PCA_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(POST_PCA_DIR, recursive = TRUE)
  message("Directory created: ", POST_PCA_DIR)
} else {
  message("Directory already exists: ", POST_PCA_DIR)
}

##INF & CDM Combined QC RLEs
POST_RLE_DIR <- paste0(BOTH_POST_DIR, "/RLE")
if (!dir.exists(POST_RLE_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(POST_RLE_DIR, recursive = TRUE)
  message("Directory created: ", POST_RLE_DIR)
} else {
  message("Directory already exists: ", POST_RLE_DIR)
}

#Directory for INF & CDM Combined QC
POST_QC_DATA_DIR <- paste0(OUTPUT_DIR, "/Post_QC_Datasets")
if (!dir.exists(POST_QC_DATA_DIR)) {
  # If it doesn't exist, create the directory
  dir.create(POST_QC_DATA_DIR, recursive = TRUE)
  message("Directory created: ", POST_QC_DATA_DIR)
} else {
  message("Directory already exists: ", POST_QC_DATA_DIR)
}

# Inflammation Panel QC and Batch effect --------------------------------

#Because of a known technical issue with INF plates, initially reading in and looking at INF data separately

### Read in INF Data and filter out controls ----
#Read in wide format NPX Target data using read_NPX function from Olink Analyze. Also read in Olink assay annotation
INF_p1<-read_NPX(paste0(INPUT_DIR,"/Cohort2_INF_Plate1_Raw.xlsx"))
INF_p2to4<-read_NPX(paste0(INPUT_DIR,"/Cohort2_INF_Plate2to4_Raw.xlsx"))
Olink_Assays <- read.csv(paste0(INPUT_DIR,"/Olink_Assays.csv"))

#Check sample numbers in each
INF_p1 %>% distinct(SampleID) %>% nrow()
INF_p2to4 %>% distinct(SampleID) %>% nrow()

####Merge P1 & P2-P4 data frames and add gene symbols ====
INF_raw<- bind_rows(INF_p1, INF_p2to4)
INF_raw %>% distinct(SampleID) %>% nrow()
#add gene symbols
INF_raw<-INF_raw %>% left_join(Olink_Assays, by = "Assay")

####Filter out Un-needed control samples & assays ====

#Remove internal PP, IPC and Neg controls (used in Olink QC but not needed here)
INF_filt <- INF_raw %>% filter(!str_detect(SampleID, c("PP|IPC|Neg")))  
#Remove any control assays, these are used by Olink analyze software but are not needed by us
INF_filt <- INF_filt %>% filter(!str_detect(Assay, "Ctrl"))#Remove internal controls

#Check number of samples and assays in filtered data frame
INF_filt %>% distinct(SampleID) %>% nrow()
INF_filt%>% distinct(Assay) %>% nrow()
INF_filt%>% distinct(Gene_Symbol) %>% nrow()


####Examining QC Warns  ====

#If a sample deviates by more than +/-0.3 NPX from the median of all samples it will be flagged as a warning for that panel
INF_warns<-INF_filt %>% filter(QC_Warning == 'Warning') #filter rows which have a QC warning
INF_warns %>% group_by(SampleID, Panel,QC_Warning) %>% summarise(count = n(), .groups = 'drop')#Print a summary

#3 samples had QC warnings on the Inflammation assay: GCA140, GCA142, NGCA88 
INF_filt[INF_filt$SampleID=="GCA140" & INF_filt$QC_Warning=="Warning","NPX"]<-NA
INF_filt[INF_filt$SampleID=="GCA142" & INF_filt$QC_Warning=="Warning","NPX"]<-NA
INF_filt[INF_filt$SampleID=="NGCA88" & INF_filt$QC_Warning=="Warning","NPX"]<-NA

### INF detectability % checks (NPX > LOD) ----

# Create new col called Gr_LOD: Yes if >LOD, No if NA or <LOD
INF_filt <- INF_filt %>% mutate(Gr_LOD = ifelse(is.na(NPX), "No", ifelse(NPX > LOD, "Yes", "No"))) 
INF_detection<-INF_filt %>% group_by(Assay, Gr_LOD,Panel,Gene_Symbol) %>% summarise(count = n(), .groups = 'drop') #summarise to protein level (long form)
INF_detection <- INF_detection %>% spread(key = Gr_LOD, value = count, fill = 0) #swap to wide form for nice summary table
#add detection percent column, sort by detection % and add >=25% detected column
INF_detection<-INF_detection %>% mutate(Det_Percent = round((Yes / 225) * 100,0)) %>% arrange(Det_Percent) %>% mutate(Detected = ifelse(Det_Percent >= 25, "yes", "no")) 

#return the number of proteins with values <25% detection (threshold used for this study)
INF_notDetected<-INF_detection %>% filter(Det_Percent < 25) %>% pull(Assay)
length(INF_notDetected)

#write detectability results to file
write.csv(INF_detection,paste0(INF_PRE_DIR,"/INF Detectability Results.csv"))

####Delta LOD Analysis ====
#visualise pattern of deviation from LOD to check for systematic issues

#Add Delta LOD calculation and detectability status (>=25%) to DF
INF_filt$Delta_LOD<-INF_filt$NPX-INF_filt$LOD
INF_filt <- INF_filt %>% left_join(INF_detection %>% select(Gene_Symbol, Detected),by ="Gene_Symbol")

#create plotting order variable to plot sample in order of Delta LOD medians
INF_D_Order = INF_filt %>% group_by(Gene_Symbol) %>% summarise(median_delta = median(na.omit(Delta_LOD))) %>% arrange(desc(median_delta)) %>% pull(Gene_Symbol)

#Create plot and save
INF_Delta_LOD_p = ggplot(INF_filt, aes(x = Gene_Symbol, y = Delta_LOD, fill = Detected)) + 
  geom_boxplot(outlier.colour = "black", outlier.size = 0.05, na.rm = TRUE, lwd = 0.5) +  # Set line width of boxes to 0.5
  geom_jitter(width = 0.2, size = 0.15, na.rm = TRUE,color="#6a7178",alpha=0.8) +  # Add jitter to show all points
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") + 
  scale_x_discrete(limits = INF_D_Order) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, size = 5, hjust = 1)) + 
  ylab("NPX deviation from LOD") + 
  xlab("") + 
  ggtitle("Delta LOD - NPX deviation from LOD") +
  labs(fill = ">=25% Detection")
INF_Delta_LOD_p
#pattern of deviation below LOD appears to have gaussian distribution for marginal proteins that are classed as "detected"
ggsave("INF Panel Delta LOD Plot.png",path=INF_PRE_DIR,plot=INF_Delta_LOD_p, width=40, height=17, units="cm")

###Output final filtered Inflammation DF as "INF_detected" ====

#Remove proteins wth <25% detection and save as new DF
INF_detected<- INF_filt%>% filter(!Assay %in% INF_notDetected)
INF_detected %>% distinct(Assay) %>% nrow()

### Confirming the batch effect ----

#Select Gene Symbol, SampleID, PlateID and NPX columns and convert to wide format DF
INF_df<- INF_detected %>% select(Gene_Symbol,SampleID, PlateID, NPX) %>% spread(key = Gene_Symbol, value = NPX, fill = NA) %>% as.data.frame()
#Removing samples with NA values for easier processing for PCA etc
INF_df<-INF_df %>% filter(!(SampleID %in% c("GCA140","GCA142","NGCA88")))
rownames(INF_df)<-INF_df$SampleID

####Visualise batch effect by PCA ====
INF_PCA<-PCA(INF_df[,3:81], scale=TRUE,graph=FALSE)
INF_plt_PCA <- create_pca_plots(INF_PCA, INF_df, "PlateID", "PCA: Inflammation Panel Before Batch Correction")
ggsave("Inflammation PlateID PCA before Batch Correction.png", path=INF_PRE_DIR, plot=INF_plt_PCA,width =85, height=25,units="cm")

####Visualise batch effect using RLE plots ====

#Input should be dataframe with rows as samples, cols as proteins, tell it whether to log or not and which are metadata columns
INF_RLE<-calcRLE(INF_df,log.data = FALSE, metadata.cols = c("SampleID","PlateID"))
INF_RLE_p <- create_rle_boxplot(INF_RLE, metadata_cols = c("SampleID", "PlateID"), color_col = "PlateID",title="Inflammation Panel Pre-Correction RLE")
INF_RLE_p
ggsave("Inflammation Plate ID RLE before Batch Correction.png", path=INF_PRE_DIR, plot=INF_RLE_p,width =60, height=22,units="cm")
#Clear batch effect visible whereby Plate 1 NPX values are quite different to P2-4

### Batch effect Adjustment / Removal ----

#Will use linear model based adjustment to correct for batch effect, see function above
#Note that this is not necessary for differential abundance analysis where batch will be added to model as covariate alongisde group

#Create new column in INF_df that designates each sample by batch: Batch1=Plate1, Batch2=Plate2-4
INF_df <- INF_df %>% mutate(Batch = ifelse(str_detect(PlateID, "Plate1"),"Batch1","Batch2")) %>% select(SampleID,PlateID,Batch,everything())
INF_df %>% count(Batch, PlateID) #check levels

INF_adjusted_df <- adjust_batch_effect(INF_df, batch_col = "Batch", metadata_cols = c("SampleID", "PlateID", "Batch"))
rownames(INF_adjusted_df)<-INF_adjusted_df$SampleID
####Check batch removal by PCA  ====
INF_Brm_PCA<-PCA(INF_adjusted_df[,4:82], scale=TRUE,graph=FALSE)
batch<-as.factor(INF_adjusted_df$Batch)
INF_Brm_PCA <- create_pca_plots(INF_Brm_PCA, INF_adjusted_df, "Batch", "PCA: Inflammation Panel After Batch Correction")
ggsave("Inflammation Batch PCA after Batch Correction.png", path=POST_PCA_DIR, plot=INF_Brm_PCA,width =85, height=25,units="cm")

####Create RLE plots and save to view clear batch effect on plate 1 ====

#Input should be dataframe with rows as samples, cols as proteins, tell it whether to log or not and which are metadata columns
INF_Brm_RLE<-calcRLE(INF_adjusted_df,log.data = FALSE, metadata.cols = c("SampleID","PlateID","Batch"))
INF_Brm_RLE_p <- create_rle_boxplot(INF_Brm_RLE, metadata_cols = c("SampleID", "PlateID","Batch"), color_col = "Batch",title="Inflammation Panel Post-Correction RLE")
INF_Brm_RLE_p
ggsave("Inflammation Batch RLE Plot after Batch Correction.png", path=POST_RLE_DIR, plot=INF_Brm_RLE_p,width =60, height=22,units="cm")

#Batch effect no longer detectable
#Therefore, we've confirmed that the addition of batch as covariate to LM used for differential abundance analysis will remove the batch effect
#Original NPX data can be used for DA analysis
#For visualisation and other downstream analyses, will merge data with post-QC CDM dataset and then transform to Z-scores.

# Cardiometabolic Panel QC --------------------------------

### Read in CDM Data and clean ----

#Read in NPX data
CDM_raw<-read_NPX(paste0(INPUT_DIR,"/Cohort2_CDM_Raw.xlsx"))
#Check sample numbers in each
CDM_raw %>% distinct(SampleID) %>% nrow()
#Add gene symbols to dataset
CDM_raw<-CDM_raw %>% left_join(Olink_Assays, by = "Assay")

####Filter out Un-needed control samples & assays ====

#Remove internal PP, IPC and Neg controls (used in Olink QC but not needed here)
CDM_filt <- CDM_raw %>% filter(!str_detect(SampleID, c("PP|IPC|Neg")))  
#Remove any control assays, these are used by Olink analyze software but are not needed by us
CDM_filt <- CDM_filt %>% filter(!str_detect(Assay, "Ctrl"))#Remove internal controls

#Check number of samples and assays in filtered data frame
CDM_filt %>% distinct(SampleID) %>% nrow()
CDM_filt%>% distinct(Assay) %>% nrow()
CDM_filt%>% distinct(Gene_Symbol) %>% nrow()

####Examining QC Warns  ====

#If a sample deviates by more than +/-0.3 NPX from the median of all samples it will be flagged as a warning for that panel
CDM_warns<-CDM_filt %>% filter(QC_Warning == 'Warning') #filter rows which have a QC warning
#no warnings found, proceed


### CDM detectability % checks (NPX > LOD) ----

# Create new col called Gr_LOD: Yes if >LOD, No if NA or <LOD
CDM_filt <- CDM_filt %>% mutate(Gr_LOD = ifelse(is.na(NPX), "No", ifelse(NPX > LOD, "Yes", "No"))) 
CDM_detection<-CDM_filt %>% group_by(Assay, Gr_LOD,Panel,Gene_Symbol) %>% summarise(count = n(), .groups = 'drop') #summarise to protein level (long form)
CDM_detection <- CDM_detection %>% spread(key = Gr_LOD, value = count, fill = 0) #swap to wide form for nice summary table
#add detection percent column, sort by detection % and add >=25% detected column
CDM_detection<-CDM_detection %>% mutate(Det_Percent = round((Yes / 225) * 100,0)) %>% arrange(Det_Percent) %>% mutate(Detected = ifelse(Det_Percent >= 25, "yes", "no")) 

#return the number of proteins with values <25% detection (threshold used for this study)
CDM_notDetected<-CDM_detection %>% filter(Det_Percent < 25) %>% pull(Assay)
length(CDM_notDetected)

#write detectability results to file
write.csv(CDM_detection,paste0(CDM_PRE_DIR,"/CDM Detectability Results.csv"))

####Delta LOD Analysis ====
#visualise pattern of deviation from LOD to check for systematic issues

#Add Delta LOD calculation and detectability status (>=25%) to DF
CDM_filt$Delta_LOD<-CDM_filt$NPX-CDM_filt$LOD
CDM_filt <- CDM_filt %>% left_join(CDM_detection %>% select(Gene_Symbol, Detected),by ="Gene_Symbol")

#create plotting order variable to plot sample in order of Delta LOD medians
CDM_D_order = CDM_filt %>% group_by(Gene_Symbol) %>% summarise(median_delta = median(na.omit(Delta_LOD))) %>% arrange(desc(median_delta)) %>% pull(Gene_Symbol)

#Create plot and save
CDM_Delta_LOD_p = ggplot(CDM_filt, aes(x = Gene_Symbol, y = Delta_LOD, fill = Detected)) + 
  geom_boxplot(outlier.colour = "black", outlier.size = 0.05, na.rm = TRUE, lwd = 0.5) +  # Set line width of boxes to 0.5
  geom_jitter(width = 0.2, size = 0.15, na.rm = TRUE,color="#6a7178",alpha=0.8) +  # Add jitter to show all points
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") + 
  scale_x_discrete(limits = CDM_D_order) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, size = 5, hjust = 1)) + 
  ylab("NPX deviation from LOD") + 
  xlab("") + 
  ggtitle("Delta LOD - NPX deviation from LOD") +
  labs(fill = ">=25% Detection")
CDM_Delta_LOD_p
#pattern of deviation below LOD appears to have gaussian distribution for marginal proteins that are classed as "detected"
ggsave("CDM Panel Delta LOD Plot.png",path=CDM_PRE_DIR,plot=CDM_Delta_LOD_p, width=40, height=17, units="cm")

###Output final filtered Cardiometabolic DF as "CDM_detected" ====

#Remove proteins wth <25% detection and save as new DF
CDM_detected<- CDM_filt%>% filter(!Assay %in% CDM_notDetected)
CDM_detected %>% distinct(Assay) %>% nrow()
#Create wide version too
CDM_df<- CDM_detected %>% select(Gene_Symbol,SampleID, PlateID, NPX) %>% spread(key = Gene_Symbol, value = NPX, fill = NA)

###CDM Panel PCA & RLE vs PlateID====

#PCA
CDM_plt<-as.factor(CDM_df$PlateID)
CDM_PCA<-PCA(CDM_df[,3:91], scale=TRUE,graph=FALSE)
CDM_plt_PCA1<-fviz_pca_ind(CDM_PCA, habillage= CDM_plt,label="none", axes=c(1,2))
CDM_plt_PCA2<-fviz_pca_ind(CDM_PCA, habillage= CDM_plt,label="none", axes=c(1,3))
ggsave("Cardiometabolic Plate ID PCA.png", path=CDM_PRE_DIR, plot=grid.arrange(CDM_plt_PCA1,CDM_plt_PCA2,nrow=1),width =70, height=25,units="cm")

#RLE
CDM_RLE<-calcRLE(CDM_df,log.data = FALSE, metadata.cols = c("SampleID","PlateID"))
CDM_RLE_p <- create_rle_boxplot(CDM_RLE, metadata_cols = c("SampleID", "PlateID"), color_col = "PlateID",title="Cardiometabolic - RLE Plot - PlateID")
CDM_RLE_p
ggsave("Cardiometabolic Plate ID RLE.png", path=CDM_PRE_DIR, plot=CDM_RLE_p,width =60, height=22,units="cm")
#No visible confounding effect of plate ID

# QC on merged dataset --------------------------------

### Read in metadata and merge with panels ----
Cohort2_metadata <- read.csv(paste0(INPUT_DIR,"/Cohort2_metadata.csv"))

#Combine batch adjusted INF and CDM data frames (without plate ID)
Merged_df<- Cohort2_metadata %>% left_join(INF_adjusted_df, by ="SampleID") %>% left_join(CDM_df, by ="SampleID") %>% select(-matches("PlateID|Batch"))
#Check sample and assay numbers
dim(Merged_df)

### PCAs on merged DF ----
#Need NA free version of DF for PCA, give it rownames for PCA labels
NArm_Merged_df <- Merged_df %>% filter(!if_any(13:ncol(Merged_df), is.na))
rownames(NArm_Merged_df)<-NArm_Merged_df$SampleID
dim(NArm_Merged_df)  #check number of rows for number of removed cases

#Create PCA object
NArm_Merged_PCA<-PCA(NArm_Merged_df[,13:180], scale=TRUE,graph=FALSE)
#Age
Merged_Age_PCA <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "Age_bin", "Age Bins PCA plots")
ggsave("Age Bin PCA.png",path=POST_PCA_DIR,plot=Merged_Age_PCA,width =50, height=20,units="cm")
#Sex
Merged_Sex_PCA <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "Sex", "Sex PCA plots")
ggsave("Sex PCA.png",path=POST_PCA_DIR,plot=Merged_Sex_PCA,width =50, height=20,units="cm")
#Ethnicity
Merged_Ethnicity_PCA <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "Ethnicity", "Ethnicity PCA plots")
ggsave("Ethnicity PCA.png",path=POST_PCA_DIR,,plot=Merged_Ethnicity_PCA,width =50, height=20,units="cm")
#Clinical Site
Merged_Site_PCA <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "Site", "Clinical Site PCA plots")
ggsave("Clinical Site PCA.png",path=POST_PCA_DIR,plot=Merged_Site_PCA,width =50, height=20,units="cm")
#Diagnosis
Merged_Diagnosis_PCA <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "Diagnosis", "Diagnosis PCA plots")
ggsave("Diagnosis PCA.png",path=POST_PCA_DIR,plot=Merged_Diagnosis_PCA,width =50, height=20,units="cm")
#Cranial Ischaemia
Merged_Isch_PCA <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "Cranial_Isch_Comp", "Cranial Ischaemic Complications PCA plots")
ggsave("Cranial Ischaemic Complications PCA.png",path=POST_PCA_DIR,plot=Merged_Isch_PCA,width =50, height=20,units="cm")
#PMR
Merged_PMR_PCA <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "PMR", "PMR PCA plots")
ggsave("PMR PCA.png",path=POST_PCA_DIR,plot=Merged_PMR_PCA,width =50, height=20,units="cm")
#TAB
Merged_TAB_PCA <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "TAB", "TAB Result PCA plots")
ggsave("TAB Result PCA.png",path=POST_PCA_DIR,plot=Merged_TAB_PCA,width =50, height=20,units="cm")
#USS
Merged_USS_PCA <- create_pca_plots(NArm_Merged_PCA, NArm_Merged_df, "USS", "USS Result PCA plots")
ggsave("USS Result PCA.png",path=POST_PCA_DIR,plot=Merged_USS_PCA,width =50, height=20,units="cm")
#Freezer Time
Merged_FreezerT_PCA1<-fviz_pca_ind(NArm_Merged_PCA,col.ind=NArm_Merged_df$Freezer_time,axes=c(1,2),title="")+labs(color="Freezer Time")
Merged_FreezerT_PCA2<-fviz_pca_ind(NArm_Merged_PCA,col.ind=NArm_Merged_df$Freezer_time,axes=c(1,3),title="")+labs(color="Freezer Time")
Merged_FreezerT_PCA3<-fviz_pca_ind(NArm_Merged_PCA,col.ind=NArm_Merged_df$Freezer_time,axes=c(2,3),title="")+labs(color="Freezer Time")
ggsave("Freezer Time PCA.png",path=POST_PCA_DIR,plot=grid.arrange(Merged_FreezerT_PCA1,Merged_FreezerT_PCA2, Merged_FreezerT_PCA3, nrow=1,top=textGrob("Freezer Time PCA plots", gp=gpar(fontzise=14,font=2))),width =50, height=20,units="cm")

### RLEs on merged DF ----
#Note that RLE plots have not had NA cases removed - meaning that some cases will have RLE calculated on basis of fewer proteins

#Create RLE object
Merged_RLE<-calcRLE(Merged_df,log.data = FALSE, metadata.cols = colnames(Merged_df[1:14]))
#Age Bins
Merged_RLEp_Age <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "Age_bin",title="All Prots - RLE Plot - Age Bin")
ggsave("Freezer Time RLE.png", path=POST_RLE_DIR,plot=Merged_RLEp_Age, width =55, height=20,units="cm")
#Sex
Merged_RLEp_Sex <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "Sex",title="All Prots - RLE Plot - Sex")
ggsave("Sex RLE.png", path=POST_RLE_DIR,plot=Merged_RLEp_Sex, width =55, height=20,units="cm")
#Ethnicity
Merged_RLEp_Ethnicity <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "Ethnicity",title="All Prots - RLE Plot - Title")
ggsave("Ethnicity RLE.png", path=POST_RLE_DIR,plot=Merged_RLEp_Ethnicity, width =55, height=20,units="cm")
#Clinical Site
Merged_RLEp_Site <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "Site",title="All Prots - RLE Plot - Clinical Site")
ggsave("Clinical Site RLE.png",path=POST_RLE_DIR, plot=Merged_RLEp_Site, width =55, height=20,units="cm")
#Diagnosis
Merged_RLEp_Dx <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "Diagnosis",title="All Prots - RLE Plot - Diagnosis")
ggsave("Diagnosis RLE.png", path=POST_RLE_DIR,plot=Merged_RLEp_Dx, width =55, height=20,units="cm")
#Cranial Ischaemia
Merged_RLEp_Isch <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "Cranial_Isch_Comp",title="All Prots - RLE Plot - Cranial Ischaemic Complication")
ggsave("Cranial Ischaemic Complications RLE.png",path=POST_RLE_DIR,plot=Merged_RLEp_Isch, width =55, height=20,units="cm")
#PMR
Merged_RLEp_PMR <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "PMR",title="All Prots - RLE Plot - PMR Symptoms")
ggsave("PMR Symptoms RLE.png", path=POST_RLE_DIR,plot=Merged_RLEp_PMR, width =55, height=20,units="cm")
#TAB
Merged_RLEp_TAB <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "TAB",title="All Prots - RLE Plot - TAB Result")
ggsave("TAB Result RLE.png",path=POST_RLE_DIR, plot=Merged_RLEp_TAB, width =55, height=20,units="cm")
#USS
Merged_RLEp_USS <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "USS",title="All Prots - RLE Plot - USS Result")
ggsave("USS Result RLE.png", path=POST_RLE_DIR,plot=Merged_RLEp_USS, width =55, height=20,units="cm")
#Freezer Time
Merged_RLEp_FreezerT <- create_rle_boxplot(Merged_RLE, metadata_cols = colnames(Merged_df[1:14]), color_col = "Freezer_time",title="All Prots - RLE Plot - Freezer Time")
ggsave("Freezer Time RLE.png", path=POST_RLE_DIR,plot=Merged_RLEp_FreezerT, width =55, height=20,units="cm")

# Write out a final post-QC long and simpler wide CSV -----------------------------

#First write out Z-score transformed version of batch corrected dataset for visualisations etc
#No sense in reading out adjusted NPX values
Scaled_Merged_df<-Merged_df
Scaled_Merged_df[,13:180]<-scale(Merged_df[,13:180],scale=TRUE,center=TRUE)
write.csv(Scaled_Merged_df,paste0(POST_QC_DATA_DIR,"/Cohort 2 Post QC Batch Adjusted Z-scores.csv"))

#Second write out Unadjusted NPX values for use in differential abundance analysis
#Because of the batch effect, will have to conduct differential abundance analysis on INF and CDM panels separately
#Output as separate datasets
final_INF_df<- Cohort2_metadata %>% left_join(INF_df, by ="SampleID") %>% select(-matches("PlateID"))
write.csv(final_INF_df,paste0(POST_QC_DATA_DIR,"/Cohort 2 Post QC Unadjusted Inflammation Panel Data.csv"))
final_CDM_df<- Cohort2_metadata %>% left_join(CDM_df, by ="SampleID") %>% select(-matches("PlateID|Batch"))
write.csv(final_CDM_df,paste0(POST_QC_DATA_DIR,"/Cohort 2 Post QC Unadjusted Cardiometabolic Panel Data.csv"))

