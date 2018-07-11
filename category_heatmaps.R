library(ComplexHeatmap)
library(ggplot2)
library(ensemblRestWrapper)

### If using RStudio set this to your working directory
# setwd("/home/user/category_heatmaps")

# Reads in a vector of TSV files. It is expected that they have headers containing at least
# ProteinID and logFC fields. Any other fields are ignored. The ProteinID should be in a 
# format supported by UniProt.
get_data <- function(files=c("avium_up.csv", "avium_down.csv"))
{
	data_files <- lapply(files, read.csv, sep="\t", stringsAsFactors=FALSE)
	combined_df <- do.call(rbind, data_files)	
	
	# Drop the any fields we don't actually need
	combined_df <- data.frame(combined_df$ProteinID, combined_df$logFC)
	
	names(combined_df) <- c("name", "logFC")

	return(combined_df)
}

# Decides what data is "important" enough to bother plotting and returns it.
# Currenlty this is defined as abs(logFC) > 1.0
get_significant_data <- function(data)
{
	lower_cut <- -1.0
	upper_cut <- 1.0
	data <- data[data$logFC < lower_cut | data$logFC > upper_cut,]		
	return(data)
}

# Plots simple heatmap with a category label if given in category parameter.
ggplot_heatmap <- function(data, category="")
{
	p <- ggplot(data, aes(1,name)) + geom_tile(aes(fill=logFC),color="white") + 
									scale_fill_gradient(low="steelblue", high="red") +  
									theme(axis.title.x=element_blank(),
									axis.text.x=element_blank(),
									axis.ticks.x=element_blank()) +
									labs(y=category) 
	return(p)
}

# Plots a heatmap from data with categories specified in the tabulated_categories_df
plot_categories <- function(data, tabulated_categories_df, GO_list) {	
	plot_list = list()
	for (i in 1:length(tabulated_categories_df$GO)) {
		print(i)
		category = tabulated_categories_df$GO[[i]]
		wanted_genes <- c()
		for (j in 1:length(GO_list)) {
			if (category %in% GO_list[[j]]) {
				wanted_genes <- c(wanted_genes, names(GO_list)[[j]])
			}	
		}
		plot_data <- data[data$name %in% wanted_genes,]
		print(plot_data)
		plot_list[[i]] <- ggplot_heatmap(plot_data, category=category)
	}
	return(plot_list)
}

# Determines what categories are "important" and returns the selection.
# Currently this is any category with at least 3 members.
select_important_categories <- function(tabulated_df) {
	tabulated_df <- tabulated_df[tabulated_df$Count >= 3,]
	return (tabulated_df)
}

# Determines the "main" GO categories to use for the final heatmap
# Takes a list of GO terms as input
get_main_categories <- function(GO_list) {
	# get a list of all unique GO terms
	unique_GOs <- unique(unlist(GO_list))
	# tabulate the numbers of genes each belongs to
	tabulated_df <- data.frame(unique_GOs, rep(0, length(unique_GOs)))	
	names(tabulated_df) <- c("GO", "Count")
	
	for (i in 1:length(GO_list)) {
		gene_GOs <- GO_list[[i]]
		if (is.null(gene_GOs))
			next
		for (j in 1:length(gene_GOs)) {
			
			tabulated_df[tabulated_df$GO == gene_GOs[[j]],]$Count <- 1 + tabulated_df[tabulated_df$GO == gene_GOs[[j]],]$Count
				
		}
	}
	tabulated_df <- tabulated_df[order(tabulated_df$Count, decreasing=TRUE),]
	tabulated_df <- select_important_categories(tabulated_df)
	
	return(tabulated_df)
					
}

get_genes_from_proteins <- function(proteins)
{
	gene_list <- lapply(as.character(proteins), get_gene_id)
	return (as.vector(gene_list))
}

remove_unknown_genes <- function(gene_list, proteins)
{
	known_entries = list()
	unknown_genes <- sapply(gene_list, identical, character(0))
	unknown_gene_names <- proteins[unknown_genes]
	known_entries$proteins <- proteins[!unknown_genes]
	
	if (length(unknown_gene_names) > 0) {
		print("Warning: Could not convert the following protein IDs to genes:") 
		print(unknown_gene_names)
	}
	
	known_entries$gene_list <- gene_list[! unknown_genes]
	return (known_entries)
	
}

# Assigns category to proteins and returns a list of proteins with their respective categories.
get_gene_GOs <- function(proteins)
{
	gene_list <- get_genes_from_proteins(proteins)
	
	# If could not convert to gene, let's let the user know, remove from vectors, and continue
	known_entries <- remove_unknown_genes(gene_list, proteins)

	GO_list <- lapply(as.character(known_entries$gene_list), get_GO_info)	
	names(GO_list) <- known_entries$proteins
	return(GO_list)	
}

save_plots <- function(plots, path="./plots/")
{
	dir.create(path, showWarnings = FALSE)
	for (i in 1:length(plots)) {
		ggsave(paste(path, i, ".pdf", sep=""), plots[[i]], device="pdf")
	}	
}

test_example <- function()
{
	data <- get_data()
	data <- get_significant_data(data)
	ggplot_heatmap(data)
	
	GO_list <- get_gene_GOs(data$name)
	tabulated_df <- get_main_categories(GO_list)
	plots <- plot_categories(data,tabulated_df, GO_list)
	save_plots(plots)
}

#test_example()