
suppressPackageStartupMessages({
  library(yaml)
})

parse_config <- function(default_config = "config/config.yaml") {
  args <- commandArgs(trailingOnly = TRUE)
  config_file <- if (length(args) >= 1 && nzchar(args[1])) args[1] else default_config
  if (!file.exists(config_file)) {
    stop("Config file not found: ", config_file,
         "\nCopy config/config.example.yaml to config/config.yaml and edit it before running.")
  }
  cfg <- yaml::read_yaml(config_file)
  cfg$config_file <- config_file
  cfg
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

load_workflow_metadata <- function(metadata_file) {
  if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)
  meta <- read.delim(metadata_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  req <- c("patient_id", "sample_id", "reference_sample_id")
  miss <- setdiff(req, colnames(meta))
  if (length(miss) > 0) stop("Missing required metadata columns: ", paste(miss, collapse = ", "))
  meta <- meta[!is.na(meta$patient_id) & !is.na(meta$sample_id) & !is.na(meta$reference_sample_id), ]
  Pers.refs <- tapply(meta$reference_sample_id, meta$patient_id, function(x) unique(x)[1])
  samp_list <- split(meta$sample_id, meta$patient_id)
  list(metadata = meta, Pers.refs = Pers.refs, samp_list = samp_list)
}
