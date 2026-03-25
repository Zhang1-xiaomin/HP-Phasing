
rm(list = ls())
library("Biostrings")
library("dplyr")
library("ggplot2")
source("00.workflow_utils.R")

cfg <- parse_config()
meta_info <- load_workflow_metadata(cfg$metadata_file)
Pers.refs <- meta_info$Pers.refs
samp_list <- meta_info$samp_list
Happlot.path <- cfg$plot_output_dir
ensure_dir(Happlot.path)

chunck.n <- 6
all.longChunck.df <- data.frame()
for (person in names(samp_list)) {
  print(person)
  phased.Hap.f <- file.path(Happlot.path, paste0(person, ".phased.Hap.forPlot.xlsx"))
  phased.Hapspltdf <- openxlsx::read.xlsx(phased.Hap.f)
  phased.Hapspltdf[, c("Start", "End", "Start.Posi", "End.Posi", "line.whith.Freq")] <- lapply(phased.Hapspltdf[, c("Start", "End", "Start.Posi", "End.Posi", "line.whith.Freq")], as.numeric)
  phased.Hapspltdf <- phased.Hapspltdf %>% mutate(Length = End.Posi - Start.Posi) %>% mutate(Locus.n = End - Start + 1)
  Chuck.long <- phased.Hapspltdf %>% filter(Locus.n >= chunck.n, Comp_to_RefMajorHap == "NS")
  all.longChunck.df <- rbind(all.longChunck.df, Chuck.long)
}

p.hist.Len <- ggplot(all.longChunck.df, aes(x = Length)) +
  geom_histogram(binwidth = 1000, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "", x = "Length of chunk in haplotype with recombination", y = "No. of chunk") +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, size = 7), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"), legend.position = "none")
ggsave(plot = p.hist.Len, filename = file.path(Happlot.path, paste0("phased.Hap-recomb.", chunck.n, ".allPersons.jpg")), width = 4, height = 2.3)
ggsave(plot = p.hist.Len, filename = file.path(Happlot.path, paste0("phased.Hap-recomb.", chunck.n, ".allPersons.pdf")), width = 4, height = 2.3)

all.snv.c <- c(); chruk.snv.c <- c(); samp.c <- c(); pers.c <- c()
for (person in names(samp_list)) {
  print(person)
  ref.samp <- Pers.refs[[person]]
  sort.Samp <- c(ref.samp, setdiff(samp_list[[person]], ref.samp))
  block.table <- read.table(file.path(cfg$isnv_table_dir, paste0(ref.samp, ".blocks.txt")), sep = "\t", quote = "", header = TRUE)
  SNV.table <- read.table(file.path(cfg$isnv_table_dir, paste0(person, ".reliable.locus.CCS.table.txt")), sep = "\t", quote = "", header = TRUE)
  for (samp in sort.Samp) {
    print(samp)
    SNVs.Freq.smp <- SNV.table[, c("X.Posi", gsub("-", ".", paste0(samp, "_2_", ref.samp)))]
    colnames(SNVs.Freq.smp) <- c("X.Posi", "Freq")
    SNV.locus <- SNVs.Freq.smp %>% filter(Freq != "NO")
    SNV.locus[, colnames(SNV.locus)] <- lapply(SNV.locus[, colnames(SNV.locus)], as.numeric)
    Tot.SNV <- SNV.locus %>% filter(Freq >= 0.05, Freq <= 0.95)
    Chuck.l <- all.longChunck.df %>% filter(sample == samp)
    if (nrow(Chuck.l) > 0) {
      all.Chuck.SNVsites <- c()
      for (Chuck.n in seq(1, nrow(Chuck.l))) {
        Chuck.S <- Chuck.l[Chuck.n, "Start.Posi"]; Chuck.E <- Chuck.l[Chuck.n, "End.Posi"]
        Chuck.SNV.sites <- Tot.SNV %>% filter(X.Posi >= Chuck.S, X.Posi <= Chuck.E)
        if (length(Chuck.SNV.sites$X.Posi) >= chunck.n) all.Chuck.SNVsites <- c(all.Chuck.SNVsites, Chuck.SNV.sites$X.Posi)
      }
      n.all.SNV <- length(Tot.SNV$X.Posi); n.Chuck.SNV <- length(unique(all.Chuck.SNVsites))
    } else {
      n.all.SNV <- length(Tot.SNV$X.Posi); n.Chuck.SNV <- 0
    }
    all.snv.c <- c(all.snv.c, n.all.SNV); chruk.snv.c <- c(chruk.snv.c, n.Chuck.SNV); samp.c <- c(samp.c, samp); pers.c <- c(pers.c, person)
  }
}

chunk.snv.df <- data.frame(person = pers.c, sample = samp.c, all.snv.No = all.snv.c, chruk.snv.No = chruk.snv.c)
chunk.snv.df <- chunk.snv.df %>% mutate(perent.chruk.snv = round(chruk.snv.No * 100 / all.snv.No, 2)) %>% mutate(X.kilo = round(all.snv.No / 1000, 4))

p.point <- ggplot(chunk.snv.df, aes(x = X.kilo, y = perent.chruk.snv)) +
  geom_point(alpha = 0.8, size = 1, shape = 1) +
  geom_smooth(method = "loess", color = "red", size = 0.5) +
  scale_color_manual(values = c("blue", "black")) +
  geom_vline(xintercept = c(0.5), col = "lightgrey", lwd = 0.5) +
  scale_x_continuous(breaks = seq(0, 9)) +
  labs(x = "SNV no. per sample(kilo)", y = "SNVs within chrunks(%)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.title = element_blank(), panel.grid.major = element_blank()) + ylim(0, 100)

ggsave(plot = p.point, filename = file.path(Happlot.path, paste0("phased.Hap-chrunk.SNV.", chunck.n, ".jpg")), width = 4, height = 2.3)
ggsave(plot = p.point, filename = file.path(Happlot.path, paste0("phased.Hap-chrunk.SNV.", chunck.n, ".pdf")), width = 4, height = 2.3)
