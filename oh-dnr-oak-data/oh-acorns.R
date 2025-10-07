library(pdftools)
library(tidyverse)

pdfoh <- "oh-dnr-oak-data/Ohio-acorns.pdf"

# Convert pdf to text, one string for each page
oak_text <- pdf_text(pdfoh)

# Loop through pages:
for (i in 1:length(oak_text)) {
  
  p <- oak_text[i]
  # Separate text for each line
  rows <- str_split(p, "\n") %>% unlist()
  
  # Specify which rows correspond to table:
  # cbind(1:length(rows), nchar(rows))
  header <- which(nchar(rows) == 0)[1] + 1
  tablerows <- (header + 1):(header + 38)
  pretitle <- grep("TABLE", rows)
  titlerows <- (pretitle + 1):(header - 2)
  
  # Extract table title
  p_title <- rows[titlerows] %>%
    str_trim() %>%
    paste(collapse = " ")
  
  # Grab data for table
  p_table <- rows[tablerows] %>% str_trim()
  
  # Convert charater strings to matrix
  mat <- matrix(" ", nrow = length(p_table), ncol = max(nchar(p_table)))
  for (j in 1:length(p_table)) {
    mat[j, 1:nchar(p_table[j])] <- p_table[j] %>% str_split_1(pattern = "")
  }
  
  # Find column start/ends
  colvals <- apply(mat, 2, function(x) all(x == " "))
  colvals <- !colvals
  
  # Compute start/endpoints of runs
  rle_cols <- rle(colvals)
  ends <- cumsum(rle_cols$lengths)
  starts <- c(1, lag(ends)[-1] + 1)
  cols_info <- data.frame(column = rle_cols$values,
                          start = starts,
                          end = ends) %>%
    filter(column == TRUE)

  # Create empty matrix
  tabb <- matrix(NA, nrow = length(p_table), ncol = nrow(cols_info))
  # Fill matrix with table values:
  for (rr in 1:length(p_table)) {
    for (cc in 1:ncol(tabb)) {
      tabb[rr, cc] <- str_sub(p_table[rr],
                              cols_info$start[cc],
                              cols_info$end[cc])
    }
  }
  tabb <- data.frame(tabb) %>%
    mutate(across(where(is.character), str_trim)) %>%
    mutate(X1 = str_remove(X1, "¹")) %>%
    mutate(across(X2:X12, as.numeric))
  colnames(tabb) <- c("wildlife_area", paste0("y", 2014:2024))
  
  assign(paste0("table", i), tabb)
  assign(paste0("title", i), p_title)
}

# Write to file
write.csv(table1,
          "oh-dnr-oak-data/percent-white-oak-trees-with-acorns.csv",
          row.names = FALSE)
write.csv(table2,
          "oh-dnr-oak-data/percent-red-oak-trees-with-acorns.csv",
          row.names = FALSE)
write.csv(table3,
          "oh-dnr-oak-data/average-percent-of-white-crowns-with-acorns.csv",
          row.names = FALSE)
write.csv(table4,
          "oh-dnr-oak-data/average-percent-of-red-crowns-with-acorns.csv",
          row.names = FALSE)