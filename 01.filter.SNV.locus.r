
library("dplyr")
source("00.workflow_utils.R")

cfg <- parse_config()
meta_info <- load_workflow_metadata(cfg$metadata_file)
Persons <- names(meta_info$Pers.refs)
ensure_dir(cfg$isnv_table_dir)

for (person in Persons) {
  print(person)

  CCSF <- file.path(cfg$ccs_table_dir, person, "all.iSNV_with_SNP.pyResults.txt")
  NGSF <- file.path(cfg$ngs_table_dir, person, "all.iSNV_with_SNP.pyResults.txt")
  CCS.df <- read.table(CCSF, sep = "\t", comment.char = "", header = TRUE)
  NGS.df <- read.table(NGSF, sep = "\t", comment.char = "", header = TRUE)
  samps.CCS <- colnames(CCS.df)[3:ncol(CCS.df)]

  kaopu.locus <- list()
  share.locus <- c()
  for (CCS.smp in samps.CCS) {
    CCS.subdf <- CCS.df[, c("X.Posi", CCS.smp)]
    NGS.subdf <- NGS.df[, c("X.Posi", CCS.smp)]
    colnames(CCS.subdf) <- c("X.Posi", "CCS")
    colnames(NGS.subdf) <- c("X.Posi", "NGS")
    mg.df <- merge(CCS.subdf, NGS.subdf, all = TRUE) %>% na.omit()
    mg.df[mg.df == "NO"] <- "0"
    mg.df$CCS <- as.numeric(mg.df$CCS)
    mg.df$NGS <- as.numeric(mg.df$NGS)
    mg.df <- mg.df %>% mutate(diff = CCS - NGS) %>% filter(abs(diff) < 0.2)
    if (length(share.locus) == 0) {
      share.locus <- mg.df$X.Posi
    } else {
      share.locus <- intersect(share.locus, mg.df$X.Posi)
    }
    kaopu.locus[[CCS.smp]] <- mg.df
  }

  filter.table <- CCS.df %>% filter(X.Posi %in% share.locus)
  filter.table$new.locus <- rownames(filter.table)

  write.table(filter.table,
              file.path(cfg$isnv_table_dir, paste0(person, ".reliable.locus.CCS.table.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
