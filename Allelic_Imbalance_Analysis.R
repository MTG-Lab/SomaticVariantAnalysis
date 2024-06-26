install.packages("ggplot")
library(ggplot2)

if (!require("tidyverse")) {
  install.packages("tidyverse")
}
library(tidyverse)

# FORMAT: CHROM	POS	REF	ALT	GenomicAF	AffectedAF	UnaffectedAF
#Set the input files here
input_files <- c("594R7NAF_out.txt", "621R2NAF_out.txt", "76R6NAF_out.txt", "158R3NAF_out.txt", "198R4NAF_out.txt", "397R4NAF_out.txt", "583R4NAF_out.txt", "650R2NAF_out.txt", "673R5NAF_out.txt", "754R2NAF_out.txt", "806R3NAF_out.txt")
input_files <- c("621R2NAF_out.txt", "76R6NAF_out.txt", "158R3NAF_out.txt")

# Define the output file path
output_file <- "output_both_flow_cells.txt"
output_file <- "output_oldflowcellMP.txt"

result_df <- data.frame(VAF = numeric(), Sample = character(), Patient = character(), stringsAsFactors = FALSE)


for (file in input_files) {
  # Read the input file
  input_df <- read.delim(file, stringsAsFactors = FALSE, header = TRUE)
  
  # Extract the patient ID from the file name
  patient <- gsub("(\\d+).*", "P\\1", file)
  
  # Filter the input dataframe for genomic samples with VAF > 40 or VAF < 60
  genomic_df <- subset(input_df, GenomicAf > 0.4 & GenomicAF < 0.6)
  
  if (nrow(genomic_df) > 0) {
    # Extract the necessary columns for affected samples
    vaf_affected <- genomic_df$AffectedAF * 100
    sample_affected <- paste0(patient, "_affected")
    
    temp_df_affected <- data.frame(VAF = vaf_affected, Sample = sample_affected, Patient = patient, stringsAsFactors = FALSE)
    
    # Append the temporary data frame for affected samples to the result data frame
    result_df <- bind_rows(result_df, temp_df_affected)
    
    # Check if UnaffectedAF column exists
    if ("UnaffectedAF" %in% names(genomic_df)) {
      # Extract the necessary columns for unaffected samples
      vaf_unaffected <- genomic_df$UnaffectedAF * 100
      sample_unaffected <- paste0(patient, "_unaffected")
      
      # Create a temporary data frame for unaffected samples
      temp_df_unaffected <- data.frame(VAF = vaf_unaffected, Sample = sample_unaffected, Patient = patient, stringsAsFactors = FALSE)
      
      # Append the temporary data frame for unaffected samples to the result data frame
      result_df <- bind_rows(result_df, temp_df_unaffected)
    }
    
    # Extract the necessary columns for genomic samples
    vaf_genomic <- genomic_df$GenomicAF * 100
    sample_genomic <- paste0(patient, "_genomic")
    
    temp_df_genomic <- data.frame(VAF = vaf_genomic, Sample = sample_genomic, Patient = patient, stringsAsFactors = FALSE)
    
    # Append the temporary data frame for genomic samples to the result data frame
    result_df <- bind_rows(result_df, temp_df_genomic)
  }
}


# Write the result data frame to the output file
write.table(result_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# plots
#--------------------------MethodsPaper--------------------------------------

library(ggplot2)

# Assuming n is your dataset
n = read.table("output_both_flow_cells.txt", header = TRUE)

# Adjusted theme settings for consistency, including font size
consistent_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                          axis.title = element_text(size = 14),
                          panel.background = element_blank(), 
                          plot.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank())

# Function to create a plot with consistent bar width and theme
create_plot <- function(data, sample_order, filename, width, height) {
  p <- ggplot(data, aes(x = Sample, y = VAF)) +
    geom_boxplot(position = position_dodge(1), width = 0.7) + # Fixed bar width
    scale_x_discrete(limits = sample_order) +
    labs(x = "Sample", y = "Variant Allele Frequency (VAF)") +
    consistent_theme
  
  ggsave(filename, plot = p, device = "pdf", width = width, height = height, dpi = 300)
}

# Define your sample orders for the different plots
sample_order_first <- c("P76_affected", "P76_unaffected", "P76_genomic", "P158_affected", "P158_unaffected", "P158_genomic", "P621_affected", "P621_unaffected", "P621_genomic", "P650_affected", "P650_unaffected", "P650_genomic", "P673_affected", "P673_unaffected", "P673_genomic", "P754_affected", "P754_unaffected", "P754_genomic", "P806_affected", "P806_unaffected", "P806_genomic")
sample_order_second <- c("P198_affected", "P198_unaffected", "P198_genomic", "P397_affected", "P397_unaffected", "P397_genomic", "P583_affected", "P583_unaffected", "P583_genomic", "P650_affected", "P650_unaffected", "P650_genomic", "P673_affected", "P673_unaffected", "P673_genomic", "P754_affected", "P754_unaffected", "P754_genomic", "P806_affected", "P806_unaffected", "P806_genomic")

# Create the plots
create_plot(n, sample_order_first, "plot2.pdf", 7, 7)
create_plot(n, sample_order_second, "plot3.pdf", 7, 7)



#----------------------------Variance for seqround2-------------------------------
library(dplyr)

# Load data from a file
data <- read.table("output_oldflowcellMP_new.txt", header = TRUE, sep = "\t")
data <- read.table("output_both_flow_cells.txt", header = TRUE, sep = "\t")

# Group by Sample and calculate variance and standard deviation
results <- data %>%
  group_by(Sample) %>%
  summarise(
    Variance = var(VAF),
    Std_Deviation = sd(VAF)
  )

# Print results
print(results)

# Optionally, write the results to a file
write.table(results, "output_variance_sd_MP_New.txt", row.names = FALSE, sep = "\t")
write.table(results, "output_both_flow_cells_variance_sd.txt", row.names = FALSE, sep = "\t")


#-------------------Calculate SD from 202 genomic samples---------------------------------
# Load data from a file
data <- read.table("all_AF_OUT.txt", header = TRUE, sep = "\t")

data$GenomicAF <- data$GenomicAF * 100

# Filter for allele frequencies between 40 and 60
filtered_data <- subset(data, GenomicAF >= 0.40 & GenomicAF <= 0.60)

# Calculate standard deviation
std_deviation <- sd(filtered_data$GenomicAF)

# Print standard deviation
print(std_deviation*2)

# Load the data from the results file
results <- read.table("output_variance_sd.txt", header = TRUE, sep = "\t")

# Filter for samples with standard deviation greater than 11.29859
filtered_samples <- subset(results, Std_Deviation > 11.29859)

# Print the samples
print(filtered_samples$Sample)


