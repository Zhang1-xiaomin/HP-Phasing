
library("dplyr")
library("stringr")
source("00.workflow_utils.R")

write_fasta <- function(sequences, file_name) {
  lines <- unlist(lapply(sequences, function(seq) c(seq$header, seq$sequence)))
  writeLines(lines, con = file_name)
}

cfg <- parse_config()
samp <- cfg$sample_id
ref.samp <- cfg$reference_sample_id
person <- cfg$patient_id
Hap.OP <- file.path(cfg$block_output_dir, samp)
ensure_dir(Hap.OP)

readHap.P <- cfg$output_dir
relaibLocus.P <- cfg$isnv_table_dir
readHap.F <- file.path(readHap.P, paste0("readHap.", samp, "_2_", ref.samp), "2.read-haplotype.Stat.txt")
readHaps <- read.table(readHap.F, sep = "\t", quote = "", header = FALSE)
colnames(readHaps) <- c("count", "S.E", "haplotype")
readHaps$id <- rownames(readHaps)
readHaps <- readHaps[!is.na(readHaps$haplotype), ]
readHaps$S <- str_split_fixed(readHaps$S.E, '-', 3)[, 2]
readHaps$E <- str_split_fixed(readHaps$S.E, '-', 3)[, 3]
readHaps$count <- as.numeric(readHaps$count)
readHaps$S <- as.numeric(readHaps$S)
readHaps$E <- as.numeric(readHaps$E)

SNV.table <- read.table(file.path(relaibLocus.P, paste0(person, ".reliable.locus.CCS.table.txt")),
                        sep = "\t", quote = "", header = TRUE)
block.table <- read.table(file.path(relaibLocus.P, paste0(samp, ".blocks.txt")),
                          sep = "\t", quote = "", header = TRUE)

blocks.phasing.df <- list()
for (n.block in seq(1, nrow(block.table))) {
  print(paste0("block : ", n.block))
  B <- block.table[n.block, "block"]
  B.S <- block.table[n.block, "Block.Start"]
  B.E <- block.table[n.block, "Block.End"]
  SNV.table.Cites <- SNV.table$X.Posi[SNV.table$X.Posi >= B.S & SNV.table$X.Posi <= B.E]
  SNV.table.Cites.df <- SNV.table.Cites %>% as.data.frame()
  colnames(SNV.table.Cites.df) <- "locus"
  Fit.df <- readHaps %>% filter(S >= B.S, E <= B.E, count >= 3)
  Fit.df <- Fit.df[order(Fit.df$count, decreasing = TRUE), ]

  n <- 0
  for (hap.n in seq(1, nrow(Fit.df))) {
    n <- n + 1
    hap <- Fit.df[hap.n, "haplotype"]
    Ct <- Fit.df[hap.n, "count"]
    hap.Nm <- paste0("Hap_", n, ".", Ct)
    bases <- str_split(hap, "_")[[1]]
    pos.c <- c(); base.c <- c()
    for (Pbase in bases) {
      base <- gsub('[0-9]', '', Pbase[1])
      Pos <- gsub('[A-Z]', '', Pbase[1])
      pos.c <- c(pos.c, Pos)
      base.c <- c(base.c, base)
    }
    hap.df <- data.frame("locus" = pos.c, "base" = base.c)
    colnames(hap.df)[2] <- paste0(hap.Nm)
    SNV.table.Cites.df <- merge(SNV.table.Cites.df, hap.df, all.x = TRUE)
  }
  rownames(SNV.table.Cites.df) <- SNV.table.Cites.df$locus
  SNV.table.Cites.df$locus <- NULL

  dfFormg <- SNV.table.Cites.df
  dfFormg.t <- t(dfFormg) %>% as.data.frame()
  dfFormg.t[is.na(dfFormg.t)] <- "Z"
  dfFormg.t$combined_column <- apply(dfFormg.t, 1, function(row) paste(row, collapse = ""))
  dfFormg.t <- dfFormg.t[order(dfFormg.t$combined_column), ]

  Out.list <- list(); O.index <- 0
  for (blockHaps in rownames(dfFormg.t)) {
    O.index <- O.index + 1
    Out.list[[O.index]] <- list(header = paste0(">", blockHaps), sequence = paste(dfFormg.t[blockHaps, "combined_column"], collapse = ""))
  }
  write_fasta(Out.list, file.path(Hap.OP, paste0(samp, ".", B, ".", B.S, "-", B.E, ".ReadHaps.fasta")))

  dfFormg.t.2 <- dfFormg.t
  dfFormg.t.2$combined_column <- NULL
  dfFormg.t.2 <- t(dfFormg.t.2) %>% as.data.frame()
  dfFormg.t.2$locus <- rownames(dfFormg.t.2)
  REFs <- dfFormg.t.2[, c("locus", colnames(dfFormg.t.2)[1])]
  n.otherHaps <- length(colnames(dfFormg.t.2)) - 1

  if (n.otherHaps > 1) {
    for (comCol in colnames(dfFormg.t.2)[2:n.otherHaps]) {
      match <- c()
      for (refCol in colnames(REFs)[2:length(colnames(REFs))]) {
        DF <- merge(REFs[, c("locus", refCol)], dfFormg.t.2[, c("locus", comCol)])
        C1 <- colnames(DF)[2]; C2 <- colnames(DF)[3]
        DF <- DF %>% mutate(comp = case_when(.data[[C1]] == "Z" | .data[[C2]] == "Z" ~ "Z", .data[[C1]] == .data[[C2]] ~ "Same", TRUE ~ "Not-same"))
        Diff.Ct <- nrow(DF[DF$comp == "Not-same", ])
        if (nrow(DF[DF$comp == "N", ]) != nrow(DF) & Diff.Ct == 0) match <- c(match, refCol)
      }
      if (length(match) == 0) REFs <- merge(REFs, dfFormg.t.2[, c("locus", comCol)], all.x = TRUE)
      if (length(match) == 1) {
        DF <- merge(REFs[, c("locus", match[1])], dfFormg.t.2[, c("locus", comCol)])
        C1 <- colnames(DF)[2]; C2 <- colnames(DF)[3]
        DF <- DF %>% mutate(comp = case_when(.data[[C1]] == "Z" | .data[[C2]] == "Z" ~ "N", .data[[C1]] == .data[[C2]] ~ "Same", TRUE ~ "Not-same"))
        DF <- DF[order(DF$locus), ]
        C1.readN <- strsplit(C1, "[.]")[[1]][2] %>% as.numeric()
        C2.readN <- strsplit(C2, "[.]")[[1]][2] %>% as.numeric()
        maxC <- if (C1.readN >= C2.readN) C1 else C2
        DF <- DF %>% mutate(mergeBase = case_when(.data[[C1]] == "Z" & .data[[C2]] == "Z" ~ "Z", !.data[[C1]] == "Z" & .data[[C2]] == "Z" ~ .data[[C1]], .data[[C1]] == "Z" & !.data[[C2]] == "Z" ~ .data[[C2]], comp == "Same" ~ .data[[C1]], comp == "Not-same" ~ .data[[maxC]], TRUE ~ ""))
        mg.DF <- DF[, c("locus", "mergeBase")]
        colnames(mg.DF)[2] <- paste0(strsplit(C1, "[.]")[[1]][1], ".", C1.readN + C2.readN)
        REFs <- merge(REFs, mg.DF, all.x = TRUE)
        REFs <- REFs %>% dplyr::select(-c(C1))
      }
      if (length(match) > 1) {
        sum_reads <- 0
        for (match.Ref in match) sum_reads <- sum_reads + (strsplit(match.Ref, "[.]")[[1]][2] %>% as.numeric())
        for (match.Ref in match) {
          DF <- merge(REFs[, c("locus", match.Ref)], dfFormg.t.2[, c("locus", comCol)])
          C1 <- colnames(DF)[2]; C2 <- colnames(DF)[3]
          C1.readN <- strsplit(C1, "[.]")[[1]][2] %>% as.numeric()
          C2.readN <- strsplit(C2, "[.]")[[1]][2] %>% as.numeric()
          DF <- DF %>% mutate(mergeBase = case_when(.data[[C1]] == "Z" & .data[[C2]] == "Z" ~ "Z", !.data[[C1]] == "Z" & .data[[C2]] == "Z" ~ .data[[C1]], .data[[C1]] == "Z" & !.data[[C2]] == "Z" ~ .data[[C2]], .data[[C1]] == .data[[C2]] ~ .data[[C1]], TRUE ~ ""))
          mg.DF <- DF[, c("locus", "mergeBase")]
          colnames(mg.DF)[2] <- paste0(strsplit(C1, "[.]")[[1]][1], ".", C1.readN + ceiling(C2.readN / sum_reads))
          REFs <- merge(REFs, mg.DF, all.x = TRUE)
          REFs <- REFs %>% dplyr::select(-c(C1))
        }
      }
    }
    REFs <- REFs[, !apply(REFs, 2, function(col) any(grepl("Z", col)))]
    all.Haps <- colnames(REFs)[2:length(colnames(REFs))]
    Len.c <- c(); for (Hap in all.Haps) Len.c <- c(Len.c, strsplit(Hap, "[.]")[[1]][length(strsplit(Hap, "[.]")[[1]])] %>% as.numeric())
    Len.sumNum <- sum(Len.c); Len.cPerct <- Len.c / Len.sumNum
    OK.index <- setdiff(seq(1, length(Len.cPerct)), which(Len.cPerct < 0.02))
    Ok.Haps <- all.Haps[OK.index]; OK.lensSum <- Len.c[OK.index] %>% sum()

    new.Hap.ID <- c(); III <- 0
    for (Hap in Ok.Haps) {
      III <- III + 1
      Num <- strsplit(Hap, "[.]")[[1]][length(strsplit(Hap, "[.]")[[1]])] %>% as.numeric()
      freqt <- round(Num / OK.lensSum, 2)
      new.Hap.ID <- c(new.Hap.ID, paste0("Phased.Hap_", III, "_Freq:", freqt, "_", Num))
    }

    REFs.O <- REFs[, c("locus", Ok.Haps)]
    colnames(REFs.O) <- c("locus", new.Hap.ID)
    rownames(REFs.O) <- REFs.O$locus; REFs.O$locus <- NULL
    REFs.O <- REFs.O[as.character(SNV.table.Cites), ]
    REFs.t <- t(REFs.O) %>% as.data.frame(); REFs.t$combined_column <- apply(REFs.t, 1, function(row) paste(row, collapse = ""))
    if (length(Ok.Haps) == 1) rownames(REFs.t) <- new.Hap.ID
    Out.Haps <- list(); O.index <- 0
    for (blockHaps in rownames(REFs.t)) {
      O.index <- O.index + 1
      Out.Haps[[O.index]] <- list(header = paste0(">", samp, ".", blockHaps), sequence = paste(REFs.t[blockHaps, "combined_column"], collapse = ""))
    }
    write_fasta(Out.Haps, file.path(Hap.OP, paste0(samp, ".", B, ".", B.S, "-", B.E, ".fasta")))
  }
}
