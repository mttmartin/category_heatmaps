#!/usr/bin/env Rscript

library(testthat)
source("../category_heatmaps.R")


context("Testing data loading")

test_data_location <- "./tests/test_data.csv"
test_data <- get_data(files=list(test_data_location))

test_that('data dimensions correct', {
		expect_equal(ncol(test_data), 2)
})

test_that('data types correct', {
	expect_is(test_data, 'data.frame')	
	expect_is(test_data$name, 'factor')
	expect_is(test_data$logFC, 'numeric')
})




context("Testing GO retrieval")
proteins <- head(test_data$name, n=5)
genes <- get_genes_from_proteins(proteins)

test_that('get back correct genes from protein IDs', {
	expect_equal(genes[[1]], "MAV_4335")	
})

test_that('removed unknown genes', {
	removed_genes <- remove_unknown_genes(genes, proteins)
	expect_equal(length(removed_genes), 4)
	expect_equal("A0A0H3A355" %in% removed_genes, FALSE)
	expect_equal("A0A0H3A355" %in% proteins, TRUE)
}) 
context("Testing data selection and subsetting")

