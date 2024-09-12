# TrinotateR-PfamFigs: A Workflow for Visualizing Pfam Domain Frequencies and Relative Expression from Trinotate Annotation

## Authors
# -	Dany Domínguez Pérez (danydguezperez@gmail.com)
# -	Maria Vittoria Modica (mariavittoria.modica@szn.it)

'General Description
This script provides a reproducible workflow for the visualization of Pfam domain frequencies and relative expression using Trinotate annotation outputs. While the TrinotateR package offers tools to summarize annotation reports, it lacks built-in support for plotting Pfam domain data. This script extends that functionality, allowing users to generate publication-quality visualizations from RNAseq data.

Extended Features
This workflow extends the functionality of TrinotateR by incorporating several key enhancements:

Pfam Domain Parsing: Parses Pfam domain annotations from Trinotate outputs, ensuring proper handling of multiple hits and fields.
TPM Expression Integration: Merges TPM (Transcripts Per Million) values from Salmon quantification into the parsed Pfam data, matching transcript identifiers. This enables the addition of expression data to functional annotations.
Top 20 Pfam Domains: Sorts and summarizes the top 20 Pfam domains based on total occurrences or TPM values, providing users with an immediate snapshot of the most prevalent or highly expressed domains.
Customizable Visualizations:
Overlay Bar Plot: Displays different annotation types (genes, transcripts, proteins, total) in a single bar for each domain, with transparency for easy comparison.
Stacked Bar Plot: Stacks annotation counts for clearer insight into the relative contributions of genes, transcripts, and proteins to each Pfam domain.
Pie Charts: Illustrates the relative expression (based on TPM) of the top 20 Pfam domains, offering a concise visual summary of highly expressed domains.
Reproducibility and Flexibility: The script uses orginal data from the single-end and paired-end RNAseq datasets from the false black coral S. savaglia, providing the step-by-step worflow generated ready-to-publish figures.'

#To begin, install the necessary dependencies:
# Install devtools and other packages below if needed
install.packages("devtools")

# Install TrinotateR from GitHub
devtools::install_github("cstubben/trinotateR")

# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("forcats")
# install.packages("ggplot2")

# Load required libraries
library(trinotateR)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)

# Step-by-Step Workflow
'1. Import Trinotate Output
First, load the Trinotate annotation report into R using the read_trinotate function from the TrinotateR package. 
This function reads the data into a fast and efficient data.table format for easy manipulation.'

'Example 1 - Trinotate_Ss_SE_report.tsv: Trinotate output from from the Single-end transcriptome outputs from the forward unpaired reads of Savalia savaglia RNAseq, containing annotations including Pfam protein domains.
Example Data is freely available at Mendeley Data (http://dx.doi.org/10.17632/3rtbr7c9s8.1)'

# Load Trinotate annotation report
x <- read_trinotate("Trinotate_Ss_SE_report.tsv", header = TRUE)

# Check the imported data
head(x)
colnames(x)
nrow(x)

#summary_trinotate returns the number of unique and total annotations in the table.
summary_trinotate(x)

'Most of the annotations contain mutliple hits in a backtick-delimited list and each hit contains multiple fields in a caret-delimited list.
For example, the second Pfam annotation below contains two hits and each hit contains a pfam id, symbol, name, alignment and e-value. The 'split_pfam' function from the 'TrinotateR' package splits multiple hits and fields, so the second Pfam annotation is now printed in rows 2 and 3 below.'
na.omit(x$Pfam)[1:2]

'2. Parse Pfam Domains
The split_pfam function splits multiple Pfam domain annotations into separate rows for better readability and manipulation.'

# Parse Pfam domains from the annotation report
x1 <- split_pfam(x)

# Check the parsed data
head(x1, 3)

'3. Summarize Pfam Domains
Summarize the Pfam domains to get the total number of unique genes, transcripts, proteins, and annotations associated with each Pfam domain.'

# Summarize Pfam domains
x2 <- summary_pfam(x1)

# View the summary
head(x2)
attr(x2, "count")

'4. Sort and Export Top 20 Pfam Domains
Sort the Pfam domains by total counts in descending order, and extract the top 20 most frequent domains.'

# Sort by total counts and extract top 20
x2_sorted <- x2[order(-x2$total), ]
x2_top20 <- head(x2_sorted, 20)

# Export top 20 Pfam domains
write.table(x2_top20, file = "x2_top20.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(x2_top20, file = "x2_top20.csv", row.names = FALSE)

'5. Reshape Data for Plotting
Reshape the dataframe from wide to long format for easy plotting using pivot_longer from the tidyr package.'

# Reshape data for plotting
x2_top20_long <- x2_top20 %>%
  pivot_longer(cols = c(genes, transcripts, proteins, total), 
               names_to = "annotations", 
               values_to = "counts")

# View reshaped data
print(x2_top20_long)

'6. Visualize Pfam Domain Frequencies
6.1 Overlay Plot with Transparency
This plot displays the top 20 Pfam domains with overlaid bars for each annotation type (genes, transcripts, proteins, and total). The transparency (alpha = 0.8) allows you to visualize how different annotation types contribute to the total counts.'

# Reorder and plot with transparency
x2_top20_long_overlay <- x2_top20_long %>%
  mutate(name = fct_reorder(name, counts, .fun = max, .desc = FALSE))

x2_top20_long_overlay %>%
  mutate(annotations = fct_relevel(annotations, "total", "transcripts", "proteins", "genes")) %>%
  arrange(name, annotations) %>%
  ggplot(aes(x = name, y = counts, fill = annotations)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.8) +  # Overlay bars with transparency
  coord_flip() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) +
  labs(title = "Top 20 Pfam Domain Counts by Annotation Type",
       x = "Pfam Domains", 
       y = "Counts")

'6.2 Stacked Bar Plot
The stacked bar plot displays the top 20 Pfam domains, stacking the counts of each annotation type (genes, transcripts, proteins) within each domain.'

# Calculate total counts and reorder for stacked plot
x2_top20_long_stacked <- x2_top20_long %>%
  group_by(name) %>%
  mutate(total_count = sum(counts[annotations == "total"])) %>%
  ungroup() %>%
  mutate(name = fct_reorder(name, total_count, .desc = FALSE))

# Plot stacked bar chart
x2_top20_long_stacked %>%
  mutate(annotations = fct_relevel(annotations, "transcripts", "genes", "proteins")) %>%
  ggplot(aes(x = name, y = counts, fill = annotations)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) +
  labs(title = "Top 20 Pfam Domain Counts by Annotation Type",
       x = "Pfam Domains", 
       y = "Counts")

'Example 2 - Trinotate_Ss_PE_report.tsv: - Trinotate output from from the paired-end de novo assembly of Savalia savaglia, containing annotations including Pfam protein domains.
Trinotate output from from the Single-end transcriptome outputs from the forward unpaired reads of Savalia savaglia RNAseq, containing annotations including Pfam protein domains.
Example Data is freely available at Mendeley Data repository (http://dx.doi.org/10.17632/ycrmntp78w.3)'

# Load the Trinotate report for paired-end RNAseq data
y <- read_trinotate("Trinotate_Ss_PE_report.tsv", header=TRUE)

# Check the first few rows and column names
head(y)
colnames(y)
# Optionally view the data (uncomment for large dataframes)
# View(y)

# Check the total number of rows
nrow(y)

# Summarize Trinotate annotations
summary_trinotate(y)

# Example of checking Pfam annotations in the dataset
na.omit(y$Pfam)[1:2]

# Parse Pfam domains using split_pfam from TrinotateR
y1 <- split_pfam(y)

# Check the parsed Pfam domain data
head(y1, 3)

# Summarize Pfam domain data (unique genes, transcripts, proteins, etc.)
y2 <- summary_pfam(y1)

# View the summarized Pfam domain data
head(y2)
View(y2)

# Check the total counts of Pfam annotations
attr(y2, "count")

# Step 1: Sort by the 'total' column in descending order
y2_sorted <- y2[order(-y2$total), ]

# Export table containing all Pfam domains sorted by total counts
write.table(y2_sorted, file = "y2_sorted.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Step 2: Extract the top 20 Pfam domains by Total Counts
y2_top20 <- head(y2_sorted, 20)

# View the top 20 Pfam domains
print(y2_top20)

# Export top 20 Pfam domains as .txt and .csv
write.table(y2_top20, file = "y2_top20.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(y2_top20, file = "y2_top20.csv", row.names = FALSE)

#Reshape the Data for Plotting
# Reshape the dataframe from wide to long format for plotting
y2_top20_long <- y2_top20 %>%
  pivot_longer(cols = c(genes, transcripts, proteins, total), 
               names_to = "annotations", 
               values_to = "counts")

# View the reshaped data
print(y2_top20_long)

# Export reshaped data
write.table(y2_top20_long, file = "y2_top20_long.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(y2_top20_long, file = "y2_top20_long.csv", row.names = FALSE)

#Visualize Pfam Domain Frequencies
#Overlay Plot with Transparency

# Reorder based on maximum counts for overlay plot
y2_top20_long_overlay <- y2_top20_long %>%
  mutate(name = fct_reorder(name, counts, .fun = max, .desc = FALSE))

# Create overlay plot with transparency
y2_top20_long_overlay %>%
  mutate(annotations = fct_relevel(annotations, "total", "transcripts", "proteins", "genes")) %>%
  arrange(name, annotations) %>%
  ggplot(aes(x = name, y = counts, fill = annotations)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.8) +  # Overlay bars with transparency
  coord_flip() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) +
  labs(title = "Top 20 Pfam Domain Counts by Annotation Type",
       x = "Pfam Domains", 
       y = "Counts")

#Stacked Bar Plot

# Step 1: Calculate total counts for each Pfam domain
y2_top20_long_stacked <- y2_top20_long %>%
  group_by(name) %>%
  mutate(total_count = sum(counts[annotations == "total"])) %>%
  ungroup() %>%
  mutate(name = fct_reorder(name, total_count, .desc = FALSE))

# Step 2: Create stacked bar plot
y2_top20_long_stacked %>%
  mutate(annotations = fct_relevel(annotations, "transcripts", "genes", "proteins")) %>%
  ggplot(aes(x = name, y = counts, fill = annotations)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) +
  labs(title = "Top 20 Pfam Domain Counts by Annotation Type",
       x = "Pfam Domains", 
       y = "Counts")

#Pfam Domain Expression Analysis: Visualizing TPM for Pfam Annotations

'This section of the workflow focuses on integrating TPM (Transcripts Per Million) expression data with Pfam domain annotations. Using the results from Salmon quantification, we merge TPM values with parsed Pfam domains to visualize relative Pfam domain expression for both single-end and paired-end RNAseq data.

The approach will correctly add TPM values from the Salmon output to the corresponding transcripts in the parsed Trinotate data (i.e., x1 for Ss_SE and y1 for Ss_PE), matching by the transcript identifier.

Workflow Overview:
Merge Salmon TPM values with Trinotate Pfam annotations.
Sum TPM values for each unique Pfam domain.
Sort and export the data for further analysis.
Visualize top 20 Pfam domains by total TPM using pie charts.'

#Example 1: Pfam Expression for Single-End RNAseq Data

'Step 1: Merge TPM Values with Pfam Annotations
First, read the Salmon quantification file and merge TPM values into the x1 dataframe, which contains Pfam annotations.'

# Read the quantification file (Salmon output was renamed from quant.sf to quant_Ss_SE for the single-end assembly)
quant_data <- read.delim("quant_Ss_SE.sf", header = TRUE)

# Merge TPM values into x1 based on transcript name
x1_quant <- x1 %>%
  left_join(quant_data %>% select(Name, TPM), by = c("transcript" = "Name"))

# Inspect the merged data
head(x1_quant)

'Step 2: Summing TPM by Pfam Domain
Next, group the data by the Pfam domain (name) and sum the TPM values. Keep the first occurrence of non-numeric columns (such as gene, symbol, etc.).'

# Sum TPM for each Pfam domain while keeping other columns
x1_pfam_expression <- x1_quant %>%
  group_by(name) %>%
  summarize(across(.cols = c(gene, transcript, protein, pfam, symbol, align, evalue), .fns = first),  
            total_TPM = sum(TPM, na.rm = TRUE)) %>%
  ungroup()  # Return to a normal dataframe

# Inspect the summarized data
head(x1_pfam_expression)

'Step 3: Sort and Export the Data
Now, sort the data by total TPM in descending order and export it for further use.'

# Sort by total_TPM in descending order
x1_pfam_expression_sorted <- x1_pfam_expression %>%
  arrange(desc(total_TPM))

# Inspect the sorted dataframe
head(x1_pfam_expression_sorted)

# Export as .txt (tab-delimited)
write.table(x1_pfam_expression_sorted, file = "x1_pfam_expression_sorted.txt", sep = "\t", row.names = FALSE, quote = FALSE)

'Step 4: Extract Top 20 Pfam Domains
We will now extract the top 20 Pfam domains based on total TPM for visualization.'

# Select the top 20 Pfam domains
x1_pfam_selected_20 <- x1_pfam_expression_sorted %>%
  select(name, symbol, total_TPM) %>%
  arrange(desc(total_TPM)) %>%
  slice(1:20)

# Inspect the top 20 Pfam domains
head(x1_pfam_selected_20)

'Step 5: Create Pie Chart for Top 20 Pfam Domains
Finally, visualize the top 20 Pfam domains by total TPM using a pie chart.'

library(ggplot2)

# Create pie chart of Pfam domains by total TPM
ggplot(x1_pfam_selected_20, aes(x = "", y = total_TPM, fill = factor(name, levels = x1_pfam_selected_20$name))) +
  geom_bar(width = 1, stat = "identity", color = "black") +  # Add black border
  coord_polar(theta = "y") +
  theme_minimal() +
  labs(title = "Pie Chart of Pfam Domains by Total TPM",
       x = NULL, y = NULL, fill = "Pfam Domain") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")  # Show the legend with Pfam domain names

# Create pie chart of all Pfam domains by total TPM on Ss_SE
# Extracting total_TPM, symbol, and name into a New Object
# Extract the desired columns into a new object

x1_pfam_selected <- x1_pfam_expression_sorted %>%
  select(name, symbol, total_TPM)

# Inspect the new object
head(x1_pfam_selected)

'Step 1: Prepare the Data
First, ensure that the data is sorted and then select the top 20 values:'

# Sort by total_TPM in descending order (if not already done)
x1_pfam_selected <- x1_pfam_selected %>%
  arrange(desc(total_TPM))

# Extract the top 20 entries for labeling
top20_pfam <- x1_pfam_selected %>%
  slice(1:20)

'Step 2: Create the Pie Chart
Now, create the pie chart using ggplot2, label the top 20 slices, and display the name labels on the side.'

library(ggplot2)

nrow(x1_pfam_selected)
#3564
# Create the pie chart without symbols on the slices, with names in the legend
ggplot(x1_pfam_selected, aes(x = "", y = total_TPM, fill = factor(name, levels = top20_pfam$name))) +
  geom_bar(width = 1, stat = "identity", color = "black") +  # Add color = "black" for borders
  coord_polar(theta = "y") +
  theme_minimal() +
  labs(title = "Pie Chart of Pfam Domains by Total TPM",
       x = NULL, y = NULL, fill = "Pfam Domain") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")  # Show the legend with names


'Example 2: Pfam Expression for Paired-End RNAseq Data
Step 1: Merge TPM Values with Pfam Annotations
Repeat the process for the paired-end data, merging TPM values from Salmon with the parsed Pfam domains in y1.'

# Read the quantification file (Salmon output)
quant_data <- read.delim("quant_Ss_PE.sf", header = TRUE)

# Merge TPM values into y1 based on transcript name
y1_quant <- y1 %>%
  left_join(quant_data %>% select(Name, TPM), by = c("transcript" = "Name"))

# Inspect the merged data
head(y1_quant)

'Step 2: Summing TPM by Pfam Domain
Next, group the data by Pfam domain (name) and sum the TPM values for each domain.'

# Sum TPM for each Pfam domain while keeping other columns
y1_pfam_expression <- y1_quant %>%
  group_by(name) %>%
  summarize(across(.cols = c(gene, transcript, protein, pfam, symbol, align, evalue), .fns = first),  
            total_TPM = sum(TPM, na.rm = TRUE)) %>%
  ungroup()  # Return to a normal dataframe

# Inspect the summarized data
head(y1_pfam_expression)

'Step 3: Sort and Export the Data
Sort the data by total TPM and export it for further use.'

# Sort by total_TPM in descending order
y1_pfam_expression_sorted <- y1_pfam_expression %>%
  arrange(desc(total_TPM))

# Inspect the sorted dataframe
head(y1_pfam_expression_sorted)

# Export as .txt (tab-delimited)
write.table(y1_pfam_expression_sorted, file = "y1_pfam_expression_sorted.txt", sep = "\t", row.names = FALSE, quote = FALSE)

'Step 4: Extract Top 20 Pfam Domains
Now, extract the top 20 Pfam domains based on total TPM for visualization.'

# Select the top 20 Pfam domains
y1_pfam_selected_20 <- y1_pfam_expression_sorted %>%
  select(name, symbol, total_TPM) %>%
  arrange(desc(total_TPM)) %>%
  slice(1:20)

# Inspect the top 20 Pfam domains
head(y1_pfam_selected_20)

'Step 5: Create Pie Chart for Top 20 Pfam Domains
Finally, visualize the top 20 Pfam domains by total TPM using a pie chart.'

# Create pie chart of Pfam domains by total TPM
ggplot(y1_pfam_selected_20, aes(x = "", y = total_TPM, fill = factor(name, levels = y1_pfam_selected_20$name))) +
  geom_bar(width = 1, stat = "identity", color = "black") +  # Add black border
  coord_polar(theta = "y") +
  theme_minimal() +
  labs(title = "Pie Chart of Pfam Domains by Total TPM",
       x = NULL, y = NULL, fill = "Pfam Domain") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")  # Show the legend with Pfam domain names

# Create pie chart including the expression of all Pfam domains by total TPM on Ss_PE
# Extracting total_TPM, symbol, and name into a New Object
# Extract the desired columns into a new object
y1_pfam_selected <- y1_pfam_expression_sorted %>%
  select(name, symbol, total_TPM)

# Inspect the new object
head(y1_pfam_selected)

'Step 1: Prepare the Data
First, ensure that the data is sorted and then select the top 20 values:'

# Sort by total_TPM in descending order (if not already done)
y1_pfam_selected <- y1_pfam_selected %>%
  arrange(desc(total_TPM))

# Extract the top 20 entries for labeling
top20_pfam <- y1_pfam_selected %>%
  slice(1:20)

'Step 2: Create the Pie Chart
Now, create the pie chart using ggplot2, label the top 20 slices, and display the name labels on the side.'

library(ggplot2)

nrow(y1_pfam_selected)
#5797
# Create the pie chart without symbols on the slices, with names in the legend
ggplot(y1_pfam_selected, aes(x = "", y = total_TPM, fill = factor(name, levels = top20_pfam$name))) +
  geom_bar(width = 1, stat = "identity", color = "black") +  # Add color = "black" for borders
  coord_polar(theta = "y") +
  theme_minimal() +
  labs(title = "Pie Chart of Pfam Domains by Total TPM",
       x = NULL, y = NULL, fill = "Pfam Domain") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")  # Show the legend with names

#Conclusion
'Remarks
This extended workflow offers an efficient approach to parse, summarize, and visualize both the frequency and expression of Pfam domains. 
By integrating Salmon TPM values, users gain the ability to visualize not only the prevalence of Pfam domains but also their relative expression in different datasets.

With customizable visualizations, the workflow provides a straightforward way to generate publication-ready figures tailored to specific needs. 
Whether working with single-end or paired-end RNAseq data, this script ensures that users can present complex transcriptomic annotations clearly and effectively.'
