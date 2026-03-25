
library("Biostrings")
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
ensure_dir(cfg$block_output_dir)

P <- c(); S <- c(); Bc <- c(); N <- c()
for (person in names(Pers.refs)) {
  print(person)
  for (samp in samp_list[[person]]) {
    Hap.OP <- file.path(cfg$block_output_dir, samp)
    block.table <- read.table(file.path(cfg$isnv_table_dir, paste0(samp, ".blocks.txt")), sep = "\t", quote = "", header = TRUE)
    for (n.block in seq(1, nrow(block.table))) {
      print(paste0("block : ", n.block))
      B <- block.table[n.block, "block"]; B.S <- block.table[n.block, "Block.Start"]; B.E <- block.table[n.block, "Block.End"]
      PhasedHaps.F1 <- file.path(Hap.OP, paste0(samp, ".", B, ".", B.S, "-", B.E, ".fasta"))
      PhasedHaps.F2 <- file.path(Hap.OP, paste0(samp, ".", B, ".", B.S, "-", B.E, ".ReadHaps.fasta"))
      PhasedHaps.F <- if (file.exists(PhasedHaps.F1)) PhasedHaps.F1 else PhasedHaps.F2
      fasta_data <- read_fasta(PhasedHaps.F)
      P <- c(P, person); S <- c(S, samp); Bc <- c(Bc, B); N <- c(N, length(names(fasta_data)))
    }
  }
}

Stat.df <- data.frame(Person = P, Sample = S, Block = Bc, Hap.n = N)
Stat.df.plt <- dplyr::group_by(Stat.df, Sample, Hap.n) |> dplyr::summarise(count = dplyr::n(), .groups = "drop")

p.bar <- ggplot2::ggplot(Stat.df.plt, ggplot2::aes(x = Hap.n, y = count)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge") +
  ggplot2::labs(title = "Histogram of haplotype number", x = "Hap number", y = "Frequency") +
  ggplot2::theme_minimal() + ggplot2::facet_wrap(~ Sample)

ggplot2::ggsave(plot = p.bar, filename = file.path(cfg$plot_output_dir, "phased.Haps.Num.jpg"), width = 10, height = 10)
