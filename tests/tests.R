#!/usr/bin/env Rscript

source("../category_heatmaps.R")


context("Testing data loading")

test_data_location <- "./test_data.csv"
test_data <- get_data(files=list(test_data_location))

test_that('data dimensions correct', {
		expect_equal(ncol(test_data), 2)
})

test_that('data types correct', {
	expect_is(test_data, 'data.frame')	
	expect_is(test_data$name, 'factor')
	expect_is(test_data$logFC, 'numeric')
})


context('Testing GO retrieval')
proteins <- head(test_data$name, n=5)
genes <- get_genes_from_proteins(proteins)

test_that('get back correct genes from protein IDs', {
	expect_equal(genes[[1]], 'MAV_4335')	
})

test_that('removed unknown genes', {
	known_entries <- remove_unknown_genes(genes, proteins)
	print('Length: ')
	expect_equal(length(known_entries$gene_list), 4)
	expect_equal('A0A0H3A355' %in% proteins, TRUE)
	expect_equal('A0A0H3A355' %in% known_entries$proteins, FALSE)
}) 

context('Testing data selection and subsetting')
GO_list <- get_gene_GOs(proteins)

test_that('get back GOs', {
	expect_equal(length(GO_list) > 0, TRUE)
})

test_that('GO list has names', {
	expect_equal(names(GO_list)[[1]], as.character(proteins[[1]]))
})

test_that('get back correct GOs', {
	expect_equal('GO:0003677' %in% GO_list[[1]], TRUE)
})

