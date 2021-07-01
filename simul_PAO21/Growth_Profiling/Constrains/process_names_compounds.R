library(readxl)
library(dplyr)

# Import the data
cp_well <- read_excel("Compuestos_pozos.xlsx")
info_exreact <- read_delim("info_exreact.txt", 
                           "\t", escape_double = FALSE, col_names = FALSE, 
                           trim_ws = TRUE)
# Process the table
cp_well_model <- cp_well[cp_well$source != "Negative Control",]
cp_well_model <- cp_well_model[cp_well_model$model == "YES",]
cp_well_model <- cp_well_model %>% select(-model)

# Check that the compounds of the plates marked as "YES" in the column "model", 
# are in the model. If the result of the sum is 0, means that we have the 
# same compounds in the plates sources that in the model
sum(!unique(cp_well_model$id) %in% unique(info_exreact$compound))

colnames(info_exreact_well)[2] <- "id"
process_table <- inner_join(cp_well_model,info_exreact_well)

# Save the table
write_csv(process_table, "exreacts_sourcesPM.csv", quote_escape = FALSE)

################################################################################
# Create a file the minimuns compounds and the compound of the well for each well

# Import the minimun compounds
compuestos_min <- read_excel("compuestos_min.xlsx")

# Create a folder to save the files created
outdir_fig <- "constrains_files"
if (dir.exists(outdir_fig) == FALSE) {
  dir.create(outdir_fig)
}

react_glucose <- "EX_cpd00027(e)"

for (env in 1:nrow(process_table)) {
  plate <- process_table[env,]$plate
  well <- process_table[env,]$well
  source <- process_table[env,]$abbreviation
  reaction <- process_table[env,]$reaction
  
  file_name <- paste0("PA021_constrains_", plate, "_", well, ".tab")
  file_dir <- paste0("constrains_files/", file_name)
  file.copy("PAO21_constrains.tab", file_dir)
  
  cat("#", source, "\n", file = file_dir, append = TRUE)
  cat(reaction,"\t-10\t10\n", file = file_dir, append = TRUE)
  
  if (plate == "PM3" | plate == "PM4") {
    cat("# Glucose", "\n", file = file_dir, append = TRUE)
    cat(react_glucose,"\t-10\t10", file = file_dir, append = TRUE)
  }
}

