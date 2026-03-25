
library("dplyr")
library("stringr")
source("00.workflow_utils.R")

cfg <- parse_config()
meta_info <- load_workflow_metadata(cfg$metadata_file)
Pers.refs <- meta_info$Pers.refs
samp_list <- meta_info$samp_list
readHap.P <- cfg$output_dir
relaibLocus.P <- cfg$isnv_table_dir
ensure_dir(cfg$isnv_table_dir)

persons_to_run <- if (!is.null(cfg$patient_id) && cfg$patient_id != "example_patient") cfg$patient_id else names(Pers.refs)

for (person in persons_to_run) {
  print(person)
  ref.samp <- Pers.refs[[person]]
  for (samp in samp_list[[person]]) {
    print(samp)
    readHap.F <- file.path(readHap.P, paste0("readHap.", samp, "_2_", ref.samp), "2.read-haplotype.Stat.txt")
    readHaps <- read.table(readHap.F, sep = "\t", quote = "", header = FALSE)
    colnames(readHaps) <- c("count", "S.E", "haplotype")
    readHaps$id <- rownames(readHaps)
    readHaps <- readHaps[!is.na(readHaps$haplotype), ]

    SNV.table <- read.table(file.path(relaibLocus.P, paste0(person, ".reliable.locus.CCS.table.txt")),
                            sep = "\t", quote = "", header = TRUE)

    readHaps$id <- rownames(readHaps)
    readHaps$S <- str_split_fixed(readHaps$S.E, '-', 3)[, 2]
    readHaps$E <- str_split_fixed(readHaps$S.E, '-', 3)[, 3]
    readHaps$count <- as.numeric(readHaps$count)
    readHaps$S <- as.numeric(readHaps$S)
    readHaps$E <- as.numeric(readHaps$E)

    locus.list <- list()
    locus.depth <- list()
    for (Locus in SNV.table$X.Posi) {
      Locus <- as.numeric(Locus)
      right.df <- readHaps %>% filter(S <= Locus & E >= Locus)
      locus.list[[paste0("P.", Locus)]] <- c(right.df$id)
      locus.depth[[paste0("P.", Locus)]] <- sum(right.df$count)
    }

    blockLst <- list()
    for (n in seq(1, (nrow(SNV.table) - 1))) {
      locus <- SNV.table[n, "X.Posi"]
      Next.locus <- SNV.table[n + 1, "X.Posi"]
      Flg.locus <- paste0("P.", locus)
      Flg.Nextlocus <- paste0("P.", Next.locus)
      readHaps.locus <- locus.list[[Flg.locus]]
      readHaps.Nextlocus <- locus.list[[Flg.Nextlocus]]
      shared.Read <- intersect(readHaps.locus, readHaps.Nextlocus)

      shared.df <- readHaps %>% filter(id %in% shared.Read)
      shared.dep <- sum(shared.df$count)
      depth.locus <- locus.depth[[Flg.locus]]
      depth.Nextlocus <- locus.depth[[Flg.Nextlocus]]
      deapMean <- mean(c(depth.locus, depth.Nextlocus))

      if (shared.dep >= deapMean / 2) {
        if (length(names(blockLst)) == 0) {
          BlockFlag <- locus
        }
      } else {
        BlockFlag <- Next.locus
      }
      blockLst[[paste0("P.", BlockFlag)]] <- Next.locus
    }

    Block.Start <- c()
    Block.End <- c()
    for (block in names(blockLst)) {
      Block.Start <- c(Block.Start, block)
      Block.End <- c(Block.End, blockLst[[block]])
    }
    blocks <- data.frame("Sample" = samp,
                         "Block.Start.1" = Block.Start,
                         "Block.End" = Block.End)
    blocks$Block.Start <- str_split_fixed(blocks$Block.Start.1, '[.]', 2)[, 2]
    blocks$block.No <- rownames(blocks)
    blocks <- blocks %>% mutate(block = paste0("block.", block.No)) %>% mutate(Block.length = as.numeric(Block.End) - as.numeric(Block.Start) + 1)
    blocks$Block.Start.1 <- NULL
    blocks$block.No <- NULL

    O.df <- blocks[, c("Sample", "block", "Block.Start", "Block.End", "Block.length")]
    write.table(O.df, file.path(relaibLocus.P, paste0(samp, ".blocks.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
