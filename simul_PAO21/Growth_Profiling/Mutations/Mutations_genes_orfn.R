library(readr)

# Import the data
Genes_list <- read_csv("Genes.txt", col_names = FALSE)
iPau21_react <- read.delim("iPau21_react.tsv")
Mutgenes_Caz5d <- read_csv("Mutgenes_Caz5d.txt",col_names = FALSE)

# KO genes of Caz5d founded in the model
Mutgenes_model <- Mutgenes_Caz5d[Mutgenes_Caz5d$X1 %in% Genes_list$X1,]
write.csv(Mutgenes_model, "mutgenes_Caz5d_model.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)

# Assign TRUE as default to all the genes
for (gene in Genes_list$X1) {
  assign(gene,TRUE)
}

# Assign FALSE to the genes KO
for (gene in Mutgenes_Caz5d$X1) {
  assign(gene,FALSE)
}

# Transforms the rules in the table as conditions in r
iPau21_react$rule <- gsub("or", "|", iPau21_react$rule)
iPau21_react$rule <- gsub("and", "&", iPau21_react$rule)

# Read the string as a condition, if is FALSE the bounds change to 0
for (reaction in 1:nrow(iPau21_react)) {
  info_reaction <- iPau21_react[reaction,]
  rule <- info_reaction$rule
  if (rule == "") {
    cat(info_reaction$abbreviation, "\t", info_reaction$name, "\n",
        file = "Reactions_without_rules.txt", append = TRUE)
  } else {
    if (eval(parse(text = rule)) == TRUE) {
    } else {
      iPau21_react[reaction,]$lowbnd <- 0
      iPau21_react[reaction,]$uppbnd <- 0
      cat(info_reaction$abbreviation, "\t", info_reaction$name, "\t", rule, "\n",
          file = "reactions_inactivated.txt", append = TRUE)
    }
  }
}

# Save the table in the correct format
iPau21_react$rule <- gsub("\\|", "or", iPau21_react$rule)
iPau21_react$rule <- gsub("&", "and", iPau21_react$rule)
write.table(iPau21_react, file = "iPau21Caz5d_react.tsv", sep = "\t", 
            row.names = FALSE ,quote = FALSE, na = "")
