#!/usr/bin/env Rscript

library(ComplexHeatmap)
library(ggplot2)
library(ensemblRestWrapper)
require(gridExtra)

### If using RStudio set this to your working directory
setwd("/home/matthew/lab_root/mycobacterium/category_heatmap")

# Reads in a vector of TSV files. It is expected that they have headers containing at least
# ProteinID and logFC fields. Any other fields are ignored. The ProteinID should be in a 
# format supported by UniProt.
get_data <- function(files=c("avium_up.tsv", "avium_down.tsv"))
{
	data_files <- lapply(files, read.csv, sep="\t", stringsAsFactors=FALSE)
	combined_df <- do.call(rbind, data_files)	
	
	# Drop the any fields we don't actually need
	combined_df <- data.frame(combined_df$ProteinID, combined_df$logFC)
	
	names(combined_df) <- c("name", "logFC")

	return(combined_df)
}

# Decides what data is "important" enough to bother plotting and returns it.
# Currently this is defined as abs(logFC) > 1.0
get_significant_data <- function(data)
{
	lower_cut <- -1.0
	upper_cut <- 1.0
	data <- data[data$logFC < lower_cut | data$logFC > upper_cut,]		
	return(data)
}

# Plots simple heatmap with a category label if given in category parameter.
ggplot_heatmap <- function(data, category="", do_use_legend=TRUE, plot_range=c(-1,1))
{
	p <- ggplot(data, aes(1,name)) + geom_tile(aes(fill=logFC),color='white') + 
									scale_fill_gradient(low='steelblue', high='red', limits=plot_range) +  
									theme(axis.title.x=element_blank(),
									axis.text.x=element_blank(),
									axis.ticks.x=element_blank()) +
									labs(y=category) 
	if (!do_use_legend) {
		p <- p + theme(legend.position="")	
	}


	return(p)
}

# Plots a heatmap from data with categories specified in the tabulated_categories_df
plot_categories <- function(data, tabulated_categories_df, GO_list, do_use_category_label=FALSE, do_use_legend=FALSE, plot_range=c(-1,1)) {	
	plot_list = list()
	for (i in 1:length(tabulated_categories_df$GO)) {
		category = tabulated_categories_df$GO[[i]]
	
		if (do_use_category_label) {
			category_label = category
		} else {
			category_label = ""
		}
		
		wanted_genes <- c()
		for (j in 1:length(GO_list)) {
			if (category %in% GO_list[[j]]) {
				wanted_genes <- c(wanted_genes, names(GO_list)[[j]])
			}	
		}
		plot_data <- data[data$name %in% wanted_genes,]
		plot_list[[i]] <- ggplot_heatmap(plot_data, category=category_label, do_use_legend=do_use_legend, plot_range=plot_range)
	}
	names(plot_list) <- tabulated_categories_df$GO
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

# Checks for protein ID that could not be converted to gene IDs (e.g. if data not available)
# It returns a named list consisting of successfully converted genes(known_entries$gene_list) and 
# the corresponding proteins(known_entries$proteins)
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

# Saves a list of single plots
save_plots <- function(plots, path="./plots/")
{
	dir.create(path, showWarnings = FALSE)
	
	for (i in 1:length(plots)) {
		ggsave(paste(path, i, ".pdf", sep=""), plots[[i]], device="pdf")
	}	
}

# Saves plots horizontally side-by-side
save_grid_plots <- function(plots, name, path="plots")
{
	dir.create(path, showWarnings = FALSE)
	p <- arrangeGrob(grobs=plots, ncol=length(plots))
	cwd = getwd()
	filename <- paste(name, ".png", sep="")
	filename <- gsub(":", "_", filename)
	fp <- file.path(cwd,path,filename)
	ggsave(fp, p, device="png")
}


get_plots <- function(data, plot_range, is_left=TRUE, is_right=FALSE)
{
	data <- get_significant_data(data)
	GO_list <- get_gene_GOs(data$name)
	tabulated_df <- get_main_categories(GO_list)
	
	if (is_left) {
		plots <- plot_categories(data, tabulated_df, GO_list, do_use_category_label=TRUE, do_use_legend=FALSE, plot_range=plot_range)
	} else if (is_right) {
		plots <- plot_categories(data, tabulated_df, GO_list, do_use_category_label=FALSE, do_use_legend=TRUE, plot_range=plot_range)
	} else {
		plots <- plot_categories(data, tabulated_df, GO_list, do_use_category_label=FALSE, do_use_legend=FALSE, plot_range=plot_range)
	}
	
	return (plots)
}

# Plots and saves heatmaps for common categories between two conditions
# The first parameter(input_files) should be a list of  lists containing
# TSV files for each condition.
# You can specify where they should be saved by save_location. 
plot_comparison_heatmap <- function(input_files, save_location="./plots/")
{
	# First determine min/max for all data sets for use in defining color range for plots
	data_min <- Inf
	data_max <- -Inf
	data_files <- list()
	for (i in 1:length(input_files)) {
		data_file <- get_data(files=input_files[[i]])
		data_files[[i]] <- data_file
		
		if (data_min > min(data_file$logFC)) {
			data_min <- min(data_file$logFC)
		}

		if (data_max < max(data_file$logFC)) {
			data_max <- max(data_file$logFC)
		}
		
	}

	plots <- list()	
	for (i in 1:length(data_files)) {
		if (i == 1) {
			new_plot <- get_plots(data_files[[i]], is_left=TRUE, plot_range=c(data_min,data_max))
		} else if (i == length(data_files)) {
			new_plot <- get_plots(data_files[[i]], is_left=FALSE, is_right=TRUE, plot_range=c(data_min,data_max))
		} else {
			new_plot <- get_plots(data_files[[i]], is_left=FALSE, plot_range=c(data_min,data_max))
		}
		plots[[i]] <- new_plot
	}

	dir.create(save_location, showWarnings = FALSE)

	common_categories <- c(names(plots[[1]]))
	for (i in 1:(length(plots)-1)) {
		next_plots_intersect <- intersect(names(plots[[i]]), names(plots[[i+1]]))
		common_categories <- intersect(common_categories, next_plots_intersect)
	}	
	
	for (i in 1:length(common_categories)) {
		# get plots corresponding to this cateogry from all plots
		category_plots <- list()
		for (j in 1:length(plots)) {
			category_plots[[j]] <- plots[[j]][[common_categories[[i]]]]
		}
		save_grid_plots(category_plots, name=common_categories[[i]])
	}
}

# This example shows how to create side-by-side heatmaps.
# It creates one plot per common GO term between the two treatments
test_example_comparison_heatmaps <- function()
{
	# These should be the TSV files corresponding to each condition you want to compare
	# The files should have the protein IDs in one column and the logFC in another
	experimental_treatment1 <- list("avium1.tsv")
	experimental_treatment2 <- list("avium2.tsv")
	experimental_treatment3 <- list("avium2.tsv")
	input_files <- list(experimental_treatment1, experimental_treatment2, experimental_treatment3)

	# This will create and save side-by-side comparisons
	# They will be saved in whatever location you specify with save_location
	plot_comparison_heatmap(input_files, save_location="./plots/")	
}


test_example_comparison_heatmaps()
