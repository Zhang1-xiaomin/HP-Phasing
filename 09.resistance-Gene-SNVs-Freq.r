
rm(list = ls())
library("Biostrings")
library("dplyr")
library("ggplot2")
source("00.workflow_utils.R")

cfg <- parse_config()
meta_info <- load_workflow_metadata(cfg$metadata_file)
Pers.refs <- meta_info$Pers.refs
samp_list <- meta_info$samp_list
ensure_dir(cfg$resistance_output_dir)

resist.Genes <- list()
resist.Genes[["Cla"]] <- c("23S ribosomal RNA")
resist.Genes[["Lef"]] <- c("DNA gyrase subunit A")

O.list <- list()
for (person in names(Pers.refs)) {
  print(person)
  ref.samp <- Pers.refs[[person]]
  sort.Samp <- c(ref.samp, setdiff(samp_list[[person]], ref.samp))
  SNV.table <- read.table(file.path(cfg$ngs_table_dir, person, "all.iSNV_with_SNP.pyResults.txt"), sep = "\t", quote = "", header = TRUE, comment.char = "")
  RefSamp.annoF <- file.path(cfg$annotation_dir, paste0(ref.samp, ".gff"), paste0(ref.samp, ".gff"))
  RefSamp.anno <- read.delim(RefSamp.annoF, header = FALSE, comment.char = "#", sep = "\t", col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
  all.Ps.SNV <- data.frame()
  for (Res.durg in names(resist.Genes)) {
    print(Res.durg)
    Res.Genes <- resist.Genes[[Res.durg]]
    for (Res.Gene in Res.Genes) {
      print(Res.Gene)
      Res.Gene.l <- RefSamp.anno %>% filter(grepl(Res.Gene, attributes))
      for (cpy in seq(1, nrow(Res.Gene.l))) {
        Res.S <- Res.Gene.l$start[cpy]; Res.E <- Res.Gene.l$end[cpy]
        SNV.table.Res <- SNV.table %>% filter(X.Posi >= Res.S, X.Posi <= Res.E)
        if (nrow(SNV.table.Res) >= 1) {
          SNV.table.Res$Gene <- Res.Gene
          SNV.table.Res$Copy <- cpy
          all.Ps.SNV <- rbind(all.Ps.SNV, SNV.table.Res)
        }
      }
    }
  }
  O.list[[person]] <- all.Ps.SNV
}
openxlsx::write.xlsx(O.list, file.path(cfg$resistance_output_dir, "All.resis.Gene.SNV.iSNVtable.xlsx"))
