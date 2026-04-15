# Check if the correct number of arguments is provided
args <- commandArgs(trailingOnly = TRUE)
if ("-h" %in% args || "--help" %in% args) {
  cat(
    "Reverse Translation Script\n",
    "Usage: Rscript SHMART.R [MUTATION_CALLER_OUTPUT_VREGION] [MUTATION_CALLER_OUTPUT_JREGION] [POOLED_REFERENCE_NUC_FASTA] [QUERY_AA_FASTA] [POOLED_REFERENCE_AA_FASTA] [FIVEMER_SUB_CSV] [OUTPUT_TRUE_TRANSLATIONS] [OUTPUT_FALSE_TRANSLATIONS] [OUTPUT_FASTA]\n\n",
    "Parameters:\n",
    "  MUTATION_CALLER_OUTPUT_VREGION: V region mutation caller output file\n",
    "  MUTATION_CALLER_OUTPUT_JREGION: J region mutation caller output file\n",
    "  POOLED_REFERENCE_NUC_FASTA: FASTA of pooled reference nucleotide sequences\n",
    "  QUERY_AA_FASTA: FASTA of query amino acid sequences\n",
    "  POOLED_REFERENCE_AA_FASTA: FASTA of pooled reference amino acid sequences\n",
    "  FIVEMER_SUB_CSV: CSV of five-mer substitution probabilities\n",
    "  OUTPUT_TRUE_TRANSLATIONS: Output file for valid reverse translations\n",
    "  OUTPUT_FALSE_TRANSLATIONS: Output file for invalid reverse translations\n",
    "  OUTPUT_FASTA: Output FASTA file\n"
  )
  quit(status = 0)
}
if (length(args) < 9) {
  stop("Missing Input Files. Use -h or --help for usage information.")
}

# Function to load required packages
pkgLoad <- function(packages = "requirements") {
  if (length(packages) == 1L && packages == "requirements") {
    packages <- c("Biostrings", "plyr", "dplyr")
  }
  
  message("Packages to be checked: ", paste(packages, collapse = ", "))
  
  packagecheck <- match(packages, utils::installed.packages()[, "Package"])
  
  message("Installed packages check: ", paste(packagecheck, collapse = ", "))
  
  packagestoinstall <- packages[is.na(packagecheck)]
  
  if (length(packagestoinstall) > 0L) {
    message("Packages to install: ", paste(packagestoinstall, collapse = ", "))
  } else {
    message("All requested packages are already installed")
  }
  
  if (length(packagestoinstall) > 0L) {
    # Check if default library path is writable
    default_lib_path <- .libPaths()[1]
    if (!file.access(default_lib_path, 2) == 0) {
      # Prompt user for an alternative library path
      message("Default library path is not writable. Please specify an alternative library path.")
      new_lib_path <- readline(prompt = "Enter writable library path: ")
      if (new_lib_path != "") {
        .libPaths(new_lib_path)
      }
    }
    
    utils::install.packages(packagestoinstall, repos = "http://cran.us.r-project.org")
  }
  
  for (package in packages) {
    message("Loading package: ", package)
    suppressPackageStartupMessages(
      library(package, character.only = TRUE, quietly = TRUE)
    )
  }
  message("All packages loaded successfully")
}

pkgLoad()

print(args)
options(show.error.locations = TRUE)

### args:
## light
# args <- c("VRegion.MutationCaller_Output.txt",                                                                                 
# "JRegion.MutationCaller_Output.txt",                                                                                 
# "/datacommons/dhvi/scripts/ReverseTranslation/REVERSE_TRANSLATE_FIX/databases/Pooled.Homo_sapiens.light.VJ.fasta",   
# "CH235_UCA4_U1788NPAG0_light.fasta",                                                                                 
# "/datacommons/dhvi/scripts/ReverseTranslation/REVERSE_TRANSLATE_FIX/databases/Pooled.Homo_sapiens.light.VJ.aa.fasta",
# "/datacommons/dhvi/scripts/ReverseTranslation/REVERSE_TRANSLATE_FIX/databases/FivemerMutationData.csv",              
# "CH235_UCA4_U1788NPAG0_light_ReverseTranslated.True.txt",                                                            
# "CH235_UCA4_U1788NPAG0_light_ReverseTranslated.FAILED.txt",                                                          
# "CH235_UCA4_U1788NPAG0_light_ReverseTranslated.True.fasta")

# heavy
args <- c("heavy/VRegion.MutationCaller_Output.txt",                                                                                 
           "heavy/JRegion.MutationCaller_Output.txt",                                                                                 
           "/datacommons/dhvi/scripts/ReverseTranslation/REVERSE_TRANSLATE_FIX/databases/Pooled.Homo_sapiens.heavy.VJ.fasta",   
           "heavy/CH235_UCA4_U1788NPAG0_heavy.fasta",                                                                                 
           "/datacommons/dhvi/scripts/ReverseTranslation/REVERSE_TRANSLATE_FIX/databases/Pooled.Homo_sapiens.heavy.VJ.aa.fasta",
           "/datacommons/dhvi/scripts/ReverseTranslation/REVERSE_TRANSLATE_FIX/databases/FivemerMutationData.csv",              
           "heavy/CH235_UCA4_U1788NPAG0_heavy_ReverseTranslated.True.txt",                                                            
           "heavy/CH235_UCA4_U1788NPAG0_heavy_ReverseTranslated.FAILED.txt",                                                          
           "heavy/CH235_UCA4_U1788NPAG0_heavy_ReverseTranslated.True.fasta")

# Function to read FASTA files
ReadFasta <- function(file) {
  fasta <- readLines(file)
  ind <- grep(">", fasta)
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 1)[-1], length(fasta)))
  seqs <- rep(NA, length(ind))
  for (i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  DF <- data.frame(V1 = gsub(">", "", fasta[ind]), V2 = seqs)
  return(DF)
}

# Function to write FASTA file
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

#Codon table
codon_table <- list(
  "A" = c("GCT","GCC","GCA","GCG"),
  "C" = c("TGT","TGC"),
  "D" = c("GAT","GAC"),
  "E" = c("GAA","GAG"),
  "F" = c("TTT","TTC"),
  "G" = c("GGT","GGC","GGA","GGG"),
  "H" = c("CAT","CAC"),
  "I" = c("ATT","ATC","ATA"),
  "K" = c("AAA","AAG"),
  "L" = c("TTA","TTG","CTT","CTC","CTA","CTG"),
  "M" = c("ATG"),
  "N" = c("AAT","AAC"),
  "P" = c("CCT","CCC","CCA","CCG"),
  "Q" = c("CAA","CAG"),
  "R" = c("CGT","CGC","CGA","CGG","AGA","AGG"),
  "S" = c("TCT","TCC","TCA","TCG","AGT","AGC"),
  "T" = c("ACT","ACC","ACA","ACG"),
  "V" = c("GTT","GTC","GTA","GTG"),
  "W" = c("TGG"),
  "Y" = c("TAT","TAC"),
  "X" = c("NNN")
)

# Function to apply mutation logic
apply_mutation_logic <-function(row) {
  # 1) Check for missing values
  if (any(is.na(row))) {
    return(NA)
  }
  
  #Extract row fields
  qseq            <- as.character(row["qseq"])            # Aligned query amino acids
  sseq            <- as.character(row["sseq"])            # Reference/germline amino acids
  sseq_nucleotide <- as.character(row["sseq_nucleotide"]) # Full reference/germline nucleotides
  
  if (is.na(qseq) || is.na(sseq) || is.na(sseq_nucleotide)) {
    return(NA)
  }
  
  # Function to Minimal Change
  calculate_minimal_change <- function(current_codon, target_codons) {
    min_changes <- Inf
    closest_codon <- NULL
    
    for (target_codon in target_codons) {
      changes <- sum(strsplit(current_codon, "")[[1]] != strsplit(target_codon, "")[[1]])
      if (changes < min_changes) {
        min_changes    <- changes
        closest_codon  <- target_codon
      }
    }
    closest_codon
  }
  
  # Minimal change - Tie-Break 
  calculate_minimal_change_tiebreak <- function(current_codon, target_codons) {
    min_changes   <- Inf
    closest_codon <- NULL
    all_paths     <- NULL
    tied_codons   <- character()
    
    for (target_codon in target_codons) {
      changes <- sum(strsplit(current_codon, "")[[1]] != strsplit(target_codon, "")[[1]])
      if (changes < min_changes) {
        min_changes   <- changes
        closest_codon <- target_codon
        tied_codons   <- c(target_codon)
        all_paths     <- generate_all_paths(current_codon, target_codons)
      } else if (changes == min_changes) {
        tied_codons <- c(tied_codons, target_codon)
        paths       <- generate_all_paths(current_codon, target_codons)
        # Tiebreak: pick the codon with fewer steps in 'paths'
        if (length(paths[[target_codon]]) < length(all_paths[[closest_codon]])) {
          closest_codon <- target_codon
          all_paths     <- paths
        }
      }
    }
    list("codon" = closest_codon, "paths" = all_paths, "tied_codons" = tied_codons)
  }
  
  #Path Generation for Ties
  generate_all_paths <- function(current_codon, target_codons) {
    all_paths <- list()
    for (tcodon in target_codons) {
      paths <- list()
      generate_paths_recursive(current_codon, tcodon, paths, "")
      all_paths[[tcodon]] <- paths
    }
    all_paths
  }
  
  generate_paths_recursive <- function(current_codon, target_codon, paths, current_path) {
    if (current_codon == target_codon) {
      paths[[length(paths) + 1]] <- trimws(current_path, whitespace = "->")
      return()
    }
    for (pos in seq_len(nchar(current_codon))) {
      if (substring(current_codon, pos, pos) != substring(target_codon, pos, pos)) {
        mutated_codon <- current_codon
        substring(mutated_codon, pos, pos) <- substring(target_codon, pos, pos)
        generate_paths_recursive(
          mutated_codon, 
          target_codon, 
          paths, 
          paste(current_path, mutated_codon, sep = " -> ")
        )
      }
    }
  }
  
  seq_length <- min(nchar(qseq), nchar(sseq))
  
  # Tied Codons (Positions)
  all_tied_codons <- vector("list", length = seq_length)
  for (i in seq_len(seq_length)) {
    all_tied_codons[[i]] <- character()
  }
  
  # Print positions with ties
  positions_with_ties <- which(sapply(all_tied_codons, length) > 1)
  if (length(positions_with_ties) > 0) {
    cat("All tied codons with minimal changes:\n")
    for (pos_i in positions_with_ties) {
      cat("Position", pos_i, ":", all_tied_codons[[pos_i]], "\n")
    }
  }
  
  # Generate Mutational Path 
  generate_mutational_paths <- function(current_codon, target_codon) {
    differing_positions <- which(strsplit(current_codon, "")[[1]] != strsplit(target_codon, "")[[1]])
    if (length(differing_positions) == 0) {
      return(list("No mutational paths found."))
    }
    paths <- list()
    for (pos in differing_positions) {
      mutated_codon <- current_codon
      substring(mutated_codon, pos, pos) <- substring(target_codon, pos, pos)
      paths[[length(paths) + 1]] <- paste(current_codon, " -> ", mutated_codon, " -> ", target_codon, sep="")
    }
    paths
  }
  
  generate_mutational_paths_pos <- function(current_codon, target_codon, mutated_position) {
    differing_positions <- which(strsplit(current_codon, "")[[1]] != strsplit(target_codon, "")[[1]])
    if (length(differing_positions) == 0) {
      return(list("No mutational paths found."))
    }
    paths <- list()
    for (pos in differing_positions) {
      mutated_codon <- current_codon
      substring(mutated_codon, pos, pos) <- substring(target_codon, pos, pos)
      paths[[length(paths) + 1]] <- list(
        "Mutation Path"       = paste(current_codon, "->", mutated_codon, "->", target_codon),
        "Amino Acid Position" = mutated_position,
        "Nucleotide Positions"= c(mutated_position*3 - 2,
                                  mutated_position*3 - 1,
                                  mutated_position*3)
      )
    }
    paths
  }
  
  process_mutational_path <- function(sseq_nucleotide, mutational_path, nucleotide_positions, probabilities) {
    s1 <- sseq_nucleotide
    current_nucleotide_positions <- nucleotide_positions
    total_probability <- 1
    
    steps <- unlist(strsplit(mutational_path, " -> "))
    for (ii in seq_len(length(steps) - 1)) {
      current_codon <- steps[ii]
      next_codon    <- steps[ii + 1]
      
      diffpos <- which(strsplit(current_codon, "")[[1]] != strsplit(next_codon, "")[[1]])[1]
      # build the 5-mer positions
      nucleotide_positions <- c(
        current_nucleotide_positions[diffpos]-2,
        current_nucleotide_positions[diffpos]-1,
        current_nucleotide_positions[diffpos],
        current_nucleotide_positions[diffpos]+1,
        current_nucleotide_positions[diffpos]+2
      )
      base_string <- substring(s1, nucleotide_positions[1], nucleotide_positions[5])
      
      # Probability
      current_mutability <- probabilities[probabilities$Fivemer == base_string, "Mutability"]
      sum_mutability <- sum(sapply(1:(nchar(s1)-4), function(x) {
        mut_5mer <- substring(s1, x, x+4)
        rowp     <- probabilities[probabilities$Fivemer == mut_5mer, "Mutability"]
        if (length(rowp) > 0) rowp else 0
      }))
      mutated_base <- substring(next_codon, diffpos, diffpos)
      mutability_by_nucleotide <- probabilities[probabilities$Fivemer == base_string, mutated_base]
      step_probability <- (current_mutability / sum_mutability) * mutability_by_nucleotide
      total_probability <- total_probability * step_probability
      
      cat("Probability for step", ii, ":", step_probability, "\n")
      
      # Update s1
      s1 <- paste0(
        substring(s1, 1, nucleotide_positions[3]-1),
        mutated_base,
        substring(s1, nucleotide_positions[3]+1)
      )
    }
    cat("Total probability:", total_probability, "\n")
  }
  
  # Read probabilities FivemerMutation data
  probabilities <- read.csv(args[6], stringsAsFactors = FALSE)
  
  # Track mutational paths
  mutational_paths_list     <- list()
  nucleotide_positions_list <- list()
  
  # For each position with ties
  for (pos_i in positions_with_ties) {
    cat("Tied Codons for Position", pos_i, "with minimal changes:\n")
    for (tied_codon in all_tied_codons[[pos_i]]) {
      mutated_position <- pos_i
      paths <- generate_mutational_paths_pos(
        substring(sseq_nucleotide, (pos_i - 1)*3 + 1, pos_i*3),
        tied_codon,
        mutated_position
      )
      cat("  - Mutational Paths for", tied_codon, ":\n")
      for (p in paths) {
        cat("    - Mutation Path:", p[["Mutation Path"]], "\n")
        cat("      Amino Acid Position:", p[["Amino Acid Position"]], "\n")
        cat("      Nucleotide Positions:", p[["Nucleotide Positions"]], "\n")
        
        mutational_paths_list     <- c(mutational_paths_list, list(mutational_path = p[["Mutation Path"]]))
        nucleotide_positions_list <- c(nucleotide_positions_list, list(nucleotide_positions = p[["Nucleotide Positions"]]))
      }
    }
    cat("\n")
  }
  
  for (xx in seq_along(mutational_paths_list)) {
    cat("Processing Mutational Path:", mutational_paths_list[[xx]], "\n")
    process_mutational_path(
      sseq_nucleotide,
      mutational_paths_list[[xx]],
      nucleotide_positions_list[[xx]],
      probabilities
    )
  }
  
  # Final output sequence
  output_sequence <- ""
  
  for (i in seq_len(seq_length)) {
    q_aa <- substring(qseq, i, i)
    s_aa <- substring(sseq, i, i)
    
    if (q_aa != s_aa) {
      # mismatch => either tie scenario or minimal-change
      if (length(all_tied_codons[[i]]) > 0) {
        # Tied scenario
        current_codons <- codon_table[[q_aa]]
        if (is.null(current_codons)) current_codons <- "NNN"
        
        tie_res <- calculate_minimal_change_tiebreak(
          substring(sseq_nucleotide, (i - 1)*3 + 1, i*3),
          current_codons
        )
        output_sequence <- paste0(output_sequence, tie_res$codon)
      } else {
        # Normal minimal-change
        current_codons <- codon_table[[q_aa]]
        if (is.null(current_codons)) current_codons <- "NNN"
        
        current_codon <- substring(sseq_nucleotide, (i - 1)*3 + 1, i*3)
        best_codon    <- calculate_minimal_change(current_codon, current_codons)
        output_sequence <- paste0(output_sequence, best_codon)
      }
    } else {
      # match => copy reference codon
      keep_codon <- substring(sseq_nucleotide, (i - 1)*3 + 1, i*3)
      output_sequence <- paste0(output_sequence, keep_codon)
    }
  }
  
  # Return final sequence
  return(output_sequence)
}

# Function for trimming sequences
trim_sequence <- function(start, end, seq) {
  if (start <= nchar(seq) & end <= nchar(seq)) {
    return(substring(seq, start, end))
  } else {
    return(NA)
  }
}

process_sequences <- function(df) {
  # Reverse-translate each row using apply_mutation_logic
  df$qseq_nucleotide <- apply(df, 1, apply_mutation_logic)
  
  df$trans_op <- translate(
    DNAStringSet(df$qseq_nucleotide),
    no.init.codon   = TRUE, 
    if.fuzzy.codon  = "X"
  )
  df$trans_op <- as.character(df$trans_op)
  
  df$qstart_aa <- as.numeric(df$qstart)
  df$qend_aa   <- as.numeric(df$qend)
  
  # Remove leading/trailing '-,X,N' from qseq
  df$qseq <- gsub("^-+|-+$", "", df$qseq)
  df$trans_op <- gsub("^X+|X+$", "", df$trans_op)
  df$trans_op_trimmed <- df$trans_op
  df$query_sequence    <- df$qseq
  
  # Crosscheck translations with query
  df$check <- (df$qseq == df$trans_op)
  
  df$qstart_nucleotide <- (df$qstart - 1) * 3 + 1
  df$qend_nucleotide   <- df$qend * 3
  
  # finalreverse-translated sequence
  df$reverse_translated_sequence <- gsub("^N+|N+$", "",df$qseq_nucleotide)
  
  # if multiple GML hits - choose the best alignment per qseqid (highest pident)
  df <- df %>%
    distinct(qseqid, .keep_all = TRUE) %>%
    group_by(qseqid) %>%
    slice_max(pident)
  
  return(df)
}

# Read Align Mutation Caller Output
v <- read.csv(args[1], sep="\t", header=TRUE)
j <- read.csv(args[2], sep="\t", header=TRUE)

# Read reference nucleotide sequences
seqs <- ReadFasta(args[3])
seqs$sseqid <- sub("\\|.*$", "", seqs$V1)
seqs <- seqs[, c(3, 2)]
colnames(seqs) <- c("sseqid", "sseq_nucleotide")
V_seqs <- plyr::join(v, seqs, type="left")
J_seqs <- plyr::join(j, seqs, type="left")

# Read query aa sequences
InputSeqs <- ReadFasta(args[4])
InputSeqs$V1 <- gsub(InputSeqs$V1, pattern=": [H/L]", replacement=":")

colnames(InputSeqs) <- c("qseqid", "Input_sequence")

# Read reference AA sequences
GermlineSeqs <- ReadFasta(args[5])
GermlineSeqs$sseqid <- sub("\\|.*$", "", GermlineSeqs$V1)
GermlineSeqs <- GermlineSeqs[, c(3, 2)]
colnames(GermlineSeqs) <- c("sseqid", "Germline_sequence")

Vseqs_Input <- plyr::join(V_seqs, InputSeqs, type="left")
VFin <- plyr::join(Vseqs_Input, GermlineSeqs, type="left")
Jseqs_Input <- plyr::join(J_seqs, InputSeqs, type="left")
JFin <- plyr::join(Jseqs_Input, GermlineSeqs, type="left")

# Process sequences for V and J
VFin <- process_sequences(VFin)
JFin <- process_sequences(JFin)

colnames(VFin)[2:ncol(VFin)] <- paste("VRegion", colnames(VFin)[2:ncol(VFin)], sep = "_")
colnames(JFin)[2:ncol(JFin)] <- paste("JRegion", colnames(JFin)[2:ncol(JFin)], sep = "_")

mergedVJ <- plyr::join(VFin, JFin, type = "left")

pick_random_codon <- function(aa, codon_table) {
  codons <- codon_table[[aa]]
  if (is.null(codons)) {
    return("NNN")
  }
  sample(codons, 1)
}

create_reverse_translated_sequence_full <- function(row) {
  
  v_nuc      <- as.character(row[["VRegion_reverse_translated_sequence"]])
  v_aa_start <- as.numeric(row[["VRegion_qstart"]])
  v_aa_end   <- as.numeric(row[["VRegion_qend"]])
  
  j_nuc      <- as.character(row[["JRegion_reverse_translated_sequence"]])
  j_aa_start <- as.numeric(row[["JRegion_qstart"]])
  j_aa_end   <- as.numeric(row[["JRegion_qend"]])
  
  full_aa    <- as.character(row[["VRegion_Input_sequence"]])
  if (is.na(v_aa_end) || is.na(j_aa_start) || is.na(full_aa)) {
    return(NA)
  }
  
  # If V region starts after AA position 1 => fill [1..(v_aa_start-1)] with random codons
  random_prepend_nuc <- ""
  if (v_aa_start > 1) {
    missing_aa_seq <- substring(full_aa, 1, v_aa_start - 1)
    for (i in seq_len(nchar(missing_aa_seq))) {
      aa_char <- substring(missing_aa_seq, i, i)
      random_prepend_nuc <- paste0(random_prepend_nuc, pick_random_codon(aa_char, codon_table))
    }
  }
  
  # Check Overlap or Gap in VJ
  overlap_len <- 0
  if (j_aa_start <= v_aa_end) {
    overlap_len <- (v_aa_end - j_aa_start + 1)
    if (overlap_len < 0) {
      overlap_len <- 0
    }
  }
  
  # Identify nDn-part in the chain (the gap between V end and J start)
  d_aa_start <- v_aa_end + 1
  d_aa_end   <- j_aa_start - 1
  
  # Build final sequence
  
  final_nuc <- paste0(random_prepend_nuc, v_nuc)
  
  # random codons for the D region
  if (d_aa_start <= d_aa_end && j_aa_start > v_aa_end) {
    d_aa_seq <- substring(full_aa, d_aa_start, d_aa_end)
    d_nucleotides <- ""
    for (i in seq_len(nchar(d_aa_seq))) {
      aa_char <- substring(d_aa_seq, i, i)
      d_nucleotides <- paste0(d_nucleotides, pick_random_codon(aa_char, codon_table))
    }
    final_nuc <- paste0(final_nuc, d_nucleotides)
  }
  
  if (!is.na(j_nuc) && j_nuc != "") {
    if (overlap_len > 0) {
      overlap_nuc_len <- overlap_len * 3
      if (overlap_nuc_len > nchar(j_nuc)) {
        overlap_nuc_len <- nchar(j_nuc)
      }
      leftover_j_nuc <- substring(j_nuc, overlap_nuc_len + 1, nchar(j_nuc))
      final_nuc <- paste0(final_nuc, leftover_j_nuc)
    } else {
      final_nuc <- paste0(final_nuc, j_nuc)
    }
  }
  
  # Check if j_aa_end < Length IP Sequence
  
  leftover_aa_start <- j_aa_end + 1
  leftover_aa_end   <- nchar(full_aa)
  if (leftover_aa_start <= leftover_aa_end) {
    leftover_aa_seq <- substring(full_aa, leftover_aa_start, leftover_aa_end)
    leftover_nuc    <- ""
    for (i in seq_len(nchar(leftover_aa_seq))) {
      aa_char <- substring(leftover_aa_seq, i, i)
      if (aa_char == "S") {
        # you said: "S is TCA for J"
        leftover_nuc <- paste0(leftover_nuc, "TCA")
      } else {
        leftover_nuc <- paste0(leftover_nuc, pick_random_codon(aa_char, codon_table))
      }
    }
    final_nuc <- paste0(final_nuc, leftover_nuc)
  }
  
  return(final_nuc)
}


mergedVJ$Full_Reverse_Translated_Sequence <- apply(mergedVJ, 1, create_reverse_translated_sequence_full)
mergedVJ_bkup <- mergedVJ

NAs <- mergedVJ %>%
    filter(
        is.na(Full_Reverse_Translated_Sequence) |
            JRegion_qstart < 0 |
            abs(JRegion_qend - nchar(JRegion_Input_sequence)) > 5 |
            VRegion_pident < 50 |
            JRegion_pident < 50
    )


NAs$Discard <- ifelse(
    (is.na(NAs$VRegion_sseqid) | NAs$VRegion_sseqid == "") &
        (is.na(NAs$JRegion_sseqid) | NAs$JRegion_sseqid == ""),
    "MISSING VJ CALL",
    ifelse(
        (is.na(NAs$VRegion_sseqid) | NAs$VRegion_sseqid == ""),
        "MISSING V CALL",
        ifelse(
            (is.na(NAs$JRegion_sseqid) | NAs$JRegion_sseqid == ""),
            "MISSING J CALL",
            ifelse(
                NAs$VRegion_pident < 50 | NAs$JRegion_pident < 50,
                ifelse(
                    NAs$VRegion_pident < 50 & NAs$JRegion_pident < 50,
                    "VJ region %ID to GML <50",
                    ifelse(
                        NAs$VRegion_pident < 50,
                        "V region %ID to GML <50",
                        "J region %ID to GML <50"
                    )
                ),
                ifelse(
                    NAs$JRegion_qstart < 0 |
                        abs(NAs$JRegion_qend - nchar(NAs$JRegion_Input_sequence)) > 5,
                    "NON-SENSICAL J REGION",
                    NA
                )
            )
        )
    )
)

NAs <- NAs[,c(1,51)]
write.table(NAs, file = "Dropout_Seqeunces.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

mergedVJ <- subset(
    mergedVJ_bkup,
    !(
        is.na(Full_Reverse_Translated_Sequence) |
            JRegion_qstart < 0 |
            abs(JRegion_qend - nchar(JRegion_Input_sequence)) > 5 |
            VRegion_pident < 50 |
            JRegion_pident < 50
    )
)

# Check FullSeq 
mergedVJ$trans_Reverse_Translated_Sequence_full <- as.character(
  translate(DNAStringSet(mergedVJ$Full_Reverse_Translated_Sequence),
            no.init.codon = TRUE, if.fuzzy.codon = "X")
)

mergedVJ$VRegion_Input_sequence <- sub("^X", "", mergedVJ$VRegion_Input_sequence)
mergedVJ$CHECK_FULL <- (mergedVJ$trans_Reverse_Translated_Sequence_full == mergedVJ$VRegion_Input_sequence)
true_rows  <- subset(mergedVJ, CHECK_FULL == TRUE)
false_rows <- subset(mergedVJ, CHECK_FULL == FALSE)


FinalVJ <- true_rows[,c(1,2,26,14,50,25,13,15,12,49,37,39,36,6:9,30:33)]

columnnamesFin <- c(
  "InputSeqID", "V.GML.ID", "J.GML.ID", "InputSequence", 
  "Reverse_Translated_Sequence_full", "VRegion_RevTrans_Sequence", 
  "VRegion_GML_Seq", "VRegion_GML_AA_Seq", "VRegion_Mutations", 
  "JRegion_RevTrans_Sequence", "JRegion_GML_Seq", "JRegion_GML_AA_Seq","JRegion_Mutations", 
  "V_input_start_pos", "V_input_end_pos", "V_GML_start_pos", "V_GML_end_pos", 
  "J_input_start_pos", "J_input_end_pos", "J_GML_start_pos", "J_GML_end_pos"
)

colnames(FinalVJ) <- columnnamesFin
write.table(FinalVJ, file = args[7], col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

FinalVJ_failed <- false_rows[,c(1,2,26,14,50,25,13,15,12,49,37,39,36,6:9,30:33)]
colnames(FinalVJ_failed) <-columnnamesFin
write.table(FinalVJ_failed, file = args[8], col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

out_fas <- FinalVJ[,c(1,5)]
colnames(out_fas) <- c("name","seq")
writeFasta(out_fas, args[9])

