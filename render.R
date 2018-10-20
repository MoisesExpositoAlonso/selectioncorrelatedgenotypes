library(rmarkdown)
render("gwasandsimulations.Rmd",
       output_format=c("html_document","pdf_document","doc_document"))