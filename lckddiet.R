#Homework 3: BioConductor & Regular Expressions

#The libraries we will use in the study
library(GEOquery)
library(ggplot2)
library(tidyr)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(biomaRt)
library(GenomicFeatures)
library(dplyr)

centered.plot.title = theme(plot.title = element_text(hjust = 0.5))  #for centered plot title

#a)Downloading the GSE115342 Dataset with the Geoquery Package
gse = getGEO("GSE115342", AnnotGPL = TRUE) #downloading dataset with getGEO function.
gse

#for reach numerical data
exp_data = exprs(gse[[1]]) #reaching to expression matrix with exprs function
head(exp_data)

dim(exp_data) #for learn dimension of the expression matrix

#b)Checking the Sample Names and Separating them into Categories
sample_names = phenoData(gse[[1]])$title #controlling sample names
sample_names

#Use regular expression
# Define the categories and patterns
categories = c("chow_cortex", "chow_liver", "lckd_cortex", "lckd_liver")
patterns = c("Ob_cortex_chow_7w_(1|2|3)", 
              "Ob_liver_chow_7w_(1|2|3)", 
              "Ob_cortex_LCKD_chow_7w_(1|2|3)", 
              "Ob_liver_LCKD_7w_(1|2|3)")

# Initialize empty lists to store indices
indices_list = list()

# Loop through each category and pattern
for (i in 1:length(categories)) {
  category = categories[i]
  pattern = patterns[i]
  
  # Find the indices
  number_indices = grep(pattern, sample_names)
  indices = which(grepl(pattern, sample_names))
  
  # Store the indices in the list
  indices_list[[category]] = indices
  
  # Print the number of samples and their indices
  print(paste("Number of samples for", gsub("_", "-", category), ":", length(number_indices)))
  print(paste("Indices of samples for", gsub("_", "-", category), ":", paste(indices, collapse = ", ")))
}

# Assign the indices to the corresponding variables
chow_cortex_ind = indices_list[["chow_cortex"]]
chow_liver_ind = indices_list[["chow_liver"]]
lckd_cortex_ind = indices_list[["lckd_cortex"]]
lckd_liver_ind = indices_list[["lckd_liver"]]

#c)Calculation of p-Values Using t-test

#Create empty vectors to store p-values
p_values_cortex = c()
p_values_liver = c()

#Perform t-test for each gene comparing mRNA levels of regular(chow)-diet and LCKD cases for cortex and liver separately
for (i in 1:nrow(exp_data)) {
  # Perform t-test for cortex
  t_test_cortex = t.test(exp_data[i, chow_cortex_ind], exp_data[i, lckd_cortex_ind])
  # Store p-value
  p_values_cortex = c(p_values_cortex, t_test_cortex$p.value)
  
  # Perform t-test for liver
  t_test_liver = t.test(exp_data[i, chow_liver_ind], exp_data[i, lckd_liver_ind])
  # Store p-value
  p_values_liver = c(p_values_liver, t_test_liver$p.value)
}

# Convert p-values to log10 scale
log_p_values_cortex = -log10(p_values_cortex)
print(paste("logP Values for Cortex:", head(log_p_values_cortex)))
log_p_values_liver = -log10(p_values_liver)
print(paste("logP Values for Liver:", head(log_p_values_liver)))

# Plot scatter
data_for_plot = data.frame(log_p_values_cortex, log_p_values_liver)

plot_1 = ggplot(data = data_for_plot)+ 
       aes(x = log_p_values_cortex, 
           y = log_p_values_liver)+
  geom_point() +
  labs(title = "Relation between Log10 of p-values of Cortex and Liver",
       x = "Log10 of p-values for Cortex", 
       y = "Log10 of p-values for Liver")+
  centered.plot.title
plot_1
  
#d)Achieving Genomic Interval from the Lowest p-value Data (with BioMart and Genomic Features Packages)

# Identify the indices of the five genes with the lowest p-values in the cortex data
lowest_p_values_indices = order(p_values_cortex)[1:5]
head(lowest_p_values_indices)

# Retrieve the data for the genes with the lowest p-values
lowest_genes = fData(gse[[1]])[lowest_p_values_indices, ]

# Extract ENSEMBL IDs and chromosomal locations of these genes
ids = lowest_genes$ENSEMBL_ID

five_chromosome = lowest_genes$CHROMOSOMAL_LOCATION

# Use BioMart to convert ENSEMBL IDs to Entrez IDs
mmusculus = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
entrez_ids = getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'), 
                   filters = 'ensembl_transcript_id', 
                   values = ids, 
                   mart = mmusculus)
# Filter out NA values
entrez_ids = na.omit(entrez_ids)

# Use the TxDb database to determine genomic ranges
txdb = TxDb.Mmusculus.UCSC.mm9.knownGene

# Retrieve genomic ranges based on Entrez IDs
genomic_intervals = genes(txdb, filter = list(gene_id = entrez_ids$entrezgene_id))
# Match genomic intervals to the retrieved Entrez IDs
matched_genomic_intervals = genomic_intervals[genomic_intervals$gene_id %in% entrez_ids$entrezgene_id]
matched_genomic_intervals

#e)Conversion into a Data Frame and its Categorization

exp_data_df = as.data.frame(exp_data)

# Transpose the dataframe
transposed_df = t(exp_data_df)
transposed_df = as.data.frame(transposed_df)

# Use phenoData "title" field to get the sample names
sample_titles = gse[[1]]@phenoData$title

# Initialize an empty vector to store the category names
category_names = c()

# Loop through each sample title and perform pattern matching
for (title in sample_titles) {
  if (grepl("Ob_cortex_chow_7w_(1|2|3)", title)) {
    category_names = c(category_names, "brain_chow")
  } else if (grepl("Ob_liver_chow_7w_(1|2|3)", title)) {
    category_names = c(category_names, "liver_chow")
  } else if (grepl("Ob_cortex_LCKD_chow_7w_(1|2|3)", title)) {
    category_names = c(category_names, "brain_lckd")
  } else if (grepl("Ob_liver_LCKD_7w_(1|2|3)", title)) {
    category_names = c(category_names, "liver_lckd")
  }
}
category_names

# Add category_names as a column to transposed_df
transposed_df = as.data.frame(cbind(category_names, transposed_df))

# Convert the "category_names" column to a factor
transposed_df$category_names = factor(transposed_df$category_names)

# Rename the new column
colnames(transposed_df)[22592] = "Gm2a"

#f) Visualization of the Expression Level of the Gm2a Gene for each Category with Boxplot

which(gse[[1]]@featureData@data$GENE_SYMBOL == "Gm2a") #finding the index of the Gm2a gene
#for validation index information
gse[[1]]@featureData@data$GENE_SYMBOL[22591]

#IMPORTANT: Since the new column is added to the first column, we need to select the 22592nd column for the Gm2a gene from the transposed data.

#for controlling, involves same info
exp_data_df[22591,]

transposed_df[,22592]

#the select function caused problems in the R file, so the select function was used only in the Rmd file
#for select the first column/category information and Gm2a column (for relevant columns)

plot_Gm2a = transposed_df[,c("category_names", "Gm2a")] %>%
  group_by(category_names) %>% #grouping by the category
  ggplot(aes(x = category_names,
             y = Gm2a,
             color = category_names))+ #color scale by category
  geom_boxplot()+ #boxplot have added
  geom_jitter(alpha = 0.9)+
  labs(title = "Boxplot of Gm2a for Category Names",
       x = "Category Names", 
       y = "Gm2a Expression Level")+
  centered.plot.title
plot_Gm2a

#g)Calculation of mRNA Expression Levels according to each Category for each Gene

# Calculate the mean mRNA expression levels for each category
mRNA_df = transposed_df %>%
  group_by(category_names) %>%
  summarise_all(mean, na.rm = TRUE)

# Set row names to the values of the 'category_names' column
rownames = mRNA_df$category_names

# Remove the 'category_names' column from the dataframe
mRNA_df_2 = mRNA_df[,-c(1)] #in .rmd file, have used select function

# Change row names with the 'rownames' vector
mRNA_df_2 = as.data.frame(mRNA_df_2)
rownames(mRNA_df_2) = rownames

# Transpose the dataframe
new_mRNA_df = t(mRNA_df_2)

# Convert to dataframe
new_mRNA_df = as.data.frame(new_mRNA_df)

#h)Correlation of Gene Expression Levels between the Same Diet Type and Different Tissues

#For this task, we'll use the dataframe obtained in part (g).

#Create a scatter plot between brain_lckd and liver_lckd categories.

lckd_brain_liver_plot = new_mRNA_df[,c("brain_lckd", "liver_lckd")] %>%
  ggplot(aes(x = brain_lckd, 
             y = liver_lckd)) +
  geom_point() + # Add points to the plot
  labs(x = "Brain LCKD mRNA Levels", 
       y = "Liver LCKD mRNA Levels", 
       title = "Scatter Plot of mRNA Levels between Brain LCKD and Liver LCKD") + # Set labels and title
  centered.plot.title  # Center the plot title
lckd_brain_liver_plot

#another plot between brain_chow and liver_chow categories.

chow_brain_liver_plot = new_mRNA_df[,c("brain_chow", "liver_chow")] %>%
  ggplot(aes(x = brain_chow, 
             y = liver_chow)) +
  geom_point() + # Add points to the plot
  geom_jitter(alpha = 0.1) + # Add jitter to avoid overlapping points
  labs(x = "Brain chow mRNA Levels", 
       y = "Liver chow mRNA Levels", 
       title = "Scatter Plot of mRNA Levels between Brain chow and Liver chow") + # Set labels and title
  centered.plot.title  # Center the plot title
chow_brain_liver_plot

#i)Converting Data to Tidy Format and Visualizing the Distribution of mRNA Levels for 4 Different Categories with Boxplot (>1000: mRNA expression level)

tidy_expr_data = new_mRNA_df %>%
  pivot_longer(cols = c("brain_chow", "brain_lckd", "liver_chow", "liver_lckd"), 
               names_to = "Category", values_to = "Expression_Level") %>%
  filter(Expression_Level > 1000)

tidy_expr_plot = ggplot(tidy_expr_data, 
       aes(x = Category, 
           y = Expression_Level, 
           color = Category)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.1) +
  theme_minimal() +
  labs(x = "Sample Category",
       y = "Expression Level", 
       title = "Boxplot of mRNA levels across categories")+
  centered.plot.title
tidy_expr_plot
