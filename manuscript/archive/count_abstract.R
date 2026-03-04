txt <- readLines(file.path(here::here(), "manuscript", "CorePAM_manuscript_submission.md"), warn = FALSE)
abstract_start <- grep("^# Abstract", txt)
trial_reg <- grep("Trial registration", txt)
abstract_lines <- txt[(abstract_start+1):(trial_reg[1]-1)]
abstract_text <- paste(abstract_lines, collapse = " ")
abstract_clean <- gsub("[[:punct:]]", " ", abstract_text)
abstract_clean <- gsub("[[:space:]]+", " ", abstract_clean)
abstract_clean <- trimws(abstract_clean)
words <- strsplit(abstract_clean, "[[:space:]]+")[[1]]
words <- words[nchar(words) > 0]
cat("Abstract word count:", length(words), "\n")
cat("Limit: 350 words\n")
cat("Status:", if(length(words) <= 350) "OK" else "OVER LIMIT!", "\n")
