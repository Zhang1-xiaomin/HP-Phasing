
rm(list = ls())
library("Biostrings")
library("dplyr")
library("ggplot2")
library("stringr")
source("00.workflow_utils.R")

read_fasta <- function(file_path) {
  lines <- readLines(file_path)
  sequences <- list(); header <- NULL
  for (line in lines) {
    if (startsWith(line, ">")) { header <- sub("^>", "", line); sequences[[header]] <- "" } else { sequences[[header]] <- paste0(sequences[[header]], line) }
  }
  sequences
}

cfg <- parse_config()
meta_info <- load_workflow_metadata(cfg$metadata_file)
Pers.refs <- meta_info$Pers.refs
samp_list <- meta_info$samp_list
O.path <- cfg$plot_output_dir
ensure_dir(O.path)

for (person in names(Pers.refs)) {
  print(person)
  block.color <- data.frame(); block.color.2 <- data.frame()
  ref.samp <- Pers.refs[[person]]
  SNV.table <- read.table(file.path(cfg$isnv_table_dir, paste0(person, ".reliable.locus.CCS.table.txt")), sep = "\t", quote = "", header = TRUE)
  sort.Samp <- c(ref.samp, setdiff(samp_list[[person]], ref.samp))

  Hap.OP <- file.path(cfg$block_output_dir, ref.samp)
  block.table <- read.table(file.path(cfg$isnv_table_dir, paste0(ref.samp, ".blocks.txt")), sep = "\t", quote = "", header = TRUE)
  Ref.majorHap <- data.frame()
  for (n.block in seq(1, nrow(block.table))) {
    B <- block.table[n.block, "block"]; B.S <- block.table[n.block, "Block.Start"]; B.E <- block.table[n.block, "Block.End"]
    SNVs.Freq.smp <- SNV.table[, c("X.Posi", gsub("-", ".", paste0(ref.samp, "_2_", ref.samp)))]
    block.SNV.sites <- SNVs.Freq.smp %>% filter(X.Posi >= B.S, X.Posi <= B.E)
    block.sites <- block.SNV.sites$X.Posi
    PhasedHaps.F1 <- file.path(Hap.OP, paste0(ref.samp, ".", B, ".", B.S, "-", B.E, ".fasta"))
    PhasedHaps.F2 <- file.path(Hap.OP, paste0(ref.samp, ".", B, ".", B.S, "-", B.E, ".ReadHaps.fasta"))
    if (file.exists(PhasedHaps.F1)) {
      fasta_data <- read_fasta(PhasedHaps.F1)
      majorHap <- names(fasta_data)[1]
      majorHap.Freq <- strsplit(strsplit(majorHap, ":")[[1]][2], "_")[[1]][1] %>% as.numeric()
      if (length(names(fasta_data)) >= 2) {
        for (OtheHap in seq(2, length(names(fasta_data)))) {
          otherHap <- names(fasta_data)[OtheHap]
          otherHap.Freq <- strsplit(strsplit(otherHap, ":")[[1]][2], "_")[[1]][1] %>% as.numeric()
          if (otherHap.Freq > majorHap.Freq) { majorHap.Freq <- otherHap.Freq; majorHap <- otherHap }
        }
      }
      MajorBase <- fasta_data[[majorHap]]
    } else {
      fasta_data <- read_fasta(PhasedHaps.F2)
      majorHap <- names(fasta_data)[1]
      MajorBase <- fasta_data[[majorHap]]
    }
    Majordf <- data.frame(posi = block.sites, base = unlist(strsplit(MajorBase, "")))
    colnames(Majordf)[2] <- "Major.Hap"
    Ref.majorHap <- if (nrow(Ref.majorHap) == 0) Majordf else rbind(Ref.majorHap, Majordf)
  }

  Samp.y.Start <- 1
  plt.Haps <- data.frame()
  for (samp in sort.Samp) {
    Hap.OP <- file.path(cfg$block_output_dir, samp)
    block.table <- read.table(file.path(cfg$isnv_table_dir, paste0(samp, ".blocks.txt")), sep = "\t", quote = "", header = TRUE)
    Phased.site.Base <- list(); y.val.min <- Samp.y.Start; xuni.breks <- c(); block.HapOlst <- list()
    for (n.block in seq(1, nrow(block.table))) {
      B <- block.table[n.block, "block"]; B.S <- block.table[n.block, "Block.Start"]; B.E <- block.table[n.block, "Block.End"]
      SNVs.Freq.smp <- SNV.table[, c("X.Posi", gsub("-", ".", paste0(samp, "_2_", ref.samp)))]
      SNVs.Freq.smp$xuniLocus <- rownames(SNVs.Freq.smp)
      block.SNV.sites <- SNVs.Freq.smp %>% filter(X.Posi >= B.S, X.Posi <= B.E)
      block.sites <- block.SNV.sites$X.Posi
      for (block.site in block.sites) Phased.site.Base[[as.character(block.site)]] <- list()
      PhasedHaps.F1 <- file.path(Hap.OP, paste0(samp, ".", B, ".", B.S, "-", B.E, ".fasta"))
      PhasedHaps.F2 <- file.path(Hap.OP, paste0(samp, ".", B, ".", B.S, "-", B.E, ".ReadHaps.fasta"))
      PhasedHaps.F <- if (file.exists(PhasedHaps.F1)) PhasedHaps.F1 else PhasedHaps.F2
      fasta_data <- read_fasta(PhasedHaps.F)
      majorHap <- names(fasta_data)[1]
      majorHap.Freq <- strsplit(strsplit(majorHap, ":")[[1]][2], "_")[[1]][1] %>% as.numeric()
      Hap.Nmlst <- c(majorHap); Hap.Freq <- list(); Hap.Freq[[majorHap]] <- majorHap.Freq
      if (length(names(fasta_data)) >= 2 && file.exists(PhasedHaps.F1)) {
        for (OtheHap in seq(2, length(names(fasta_data)))) {
          otherHap <- names(fasta_data)[OtheHap]
          otherHap.Freq <- strsplit(strsplit(otherHap, ":")[[1]][2], "_")[[1]][1] %>% as.numeric()
          Hap.Nmlst <- c(Hap.Nmlst, otherHap); Hap.Freq[[otherHap]] <- otherHap.Freq
          if (otherHap.Freq > majorHap.Freq) { majorHap.Freq <- otherHap.Freq; majorHap <- otherHap }
        }
      }
      Hap.df <- data.frame(); Freqs <- list(); Freq.c <- c()
      for (Hhp in names(Hap.Freq)) {
        fq <- Hap.Freq[[Hhp]]; Freq.c <- unique(c(Freq.c, fq))
        if (!as.character(fq) %in% names(Freqs)) Freqs[[as.character(fq)]] <- c(Hhp) else Freqs[[as.character(fq)]] <- c(Freqs[[as.character(fq)]], Hhp)
      }
      Freq.c <- sort(Freq.c, decreasing = TRUE)
      for (Hap in names(fasta_data)) {
        baseSeq <- unlist(strsplit(fasta_data[[Hap]], ""))
        df <- data.frame(citeIdx = block.sites, base = baseSeq)
        colnames(df) <- c("citeIdx", Hap)
        Hap.df <- if (nrow(Hap.df) == 0) df else merge(Hap.df, df)
      }
      Hap.df.s <- Hap.df
      Hap.df.s$A_count <- rowSums(data.frame(Hap.df.s[, Hap.Nmlst]) == "A")
      Hap.df.s$G_count <- rowSums(data.frame(Hap.df.s[, Hap.Nmlst]) == "G")
      Hap.df.s$C_count <- rowSums(data.frame(Hap.df.s[, Hap.Nmlst]) == "C")
      Hap.df.s$T_count <- rowSums(data.frame(Hap.df.s[, Hap.Nmlst]) == "T")
      stat.cols <- c("A_count", "G_count", "C_count", "T_count")
      Hap.df.s$zero.Count <- rowSums(Hap.df.s[, stat.cols] == 0)
      Hap.df.s <- Hap.df.s %>% mutate(Mut_or_not = case_when(zero.Count == 3 ~ "No", TRUE ~ "Mut")) %>% dplyr::select(-c("A_count", "G_count", "C_count", "T_count"))
      Hap.df.s <- merge(Hap.df.s, Ref.majorHap, by.x = "citeIdx", by.y = "posi", all.x = TRUE)
      for (block.site in block.sites) {
        MutNot <- Hap.df.s[Hap.df.s$citeIdx == as.numeric(block.site), "Mut_or_not"]
        BASE <- Hap.df.s[Hap.df.s$citeIdx == as.numeric(block.site), "Major.Hap"]
        if (MutNot == "No") Phased.site.Base[[as.character(block.site)]][[BASE]] <- 1
      }
      Hap.df.s <- merge(Hap.df.s, SNVs.Freq.smp, by.x = "citeIdx", by.y = "X.Posi", all.x = TRUE)
      xuni.brek <- Hap.df.s[1, "xuniLocus"] %>% as.numeric(); xuni.breks <- c(xuni.breks, xuni.brek)
      Phased.block.plt <- data.frame(); y.val <- Samp.y.Start
      for (Freq.idx in seq(1, length(Freq.c))) {
        Freq <- Freq.c[Freq.idx]; Freq.haps <- Freqs[[as.character(Freq)]]
        for (Freq.hap in Freq.haps) {
          y.val <- y.val - 1; if (y.val < y.val.min) y.val.min <- y.val
          baseSeq <- unlist(strsplit(fasta_data[[Freq.hap]], ""))
          df <- data.frame(Start = seq(1, length(baseSeq)), End = seq(1, length(baseSeq)), Base = baseSeq)
          df$line.whith.Freq <- Freq; df$sample <- samp; df$block <- B; df$Start.Posi <- block.sites; df$End.Posi <- block.sites; df$y <- y.val
          Phased.block.plt <- rbind(Phased.block.plt, df)
        }
      }
      block.HapOlst[[B]] <- Phased.block.plt
      Samp.y.Start <- y.val.min - 2
    }
    if (length(block.HapOlst) > 0) plt.Haps <- rbind(plt.Haps, bind_rows(block.HapOlst))
  }
  if (nrow(plt.Haps) > 0) openxlsx::write.xlsx(plt.Haps, file.path(O.path, paste0(person, ".phased.Hap.forPlot.xlsx")), overwrite = TRUE)
}
