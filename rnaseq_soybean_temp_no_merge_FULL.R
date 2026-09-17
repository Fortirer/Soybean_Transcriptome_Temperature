## =====================================================================
##  Soybean_Temp  -  RNA-seq completo, 20 amostras, SEM merge/sort
##  Rsubread align (por lane)  ->  featureCounts (4 lanes -> 1 coluna)
##
##  Resolve o "No space left on device": nunca existem BAMs de mais de
##  uma amostra no disco ao mesmo tempo. Pico de uso ~= 4 BAMs de lane.
##
##  Retomavel: se parar no meio, rode de novo e ele continua de onde
##  estava (pula amostras que ja tem o .rds de contagem).
## =====================================================================

suppressPackageStartupMessages({
  library(Rsubread)
  library(Rsamtools)
  library(data.table)
})

## --------------------------- CONFIG ---------------------------------
raw_dir    <- "D:/AnalisesPreteritas/analysis/190624_NB501279_0080_AHKTTNBGX7/190624_NB501279_0080_AHKTTNBGX7_fastq"
ref_dir    <- "D:/AnalisesPreteritas/analysis/190624_NB501279_0080_AHKTTNBGX7/reference"
out_dir    <- "D:/AnalisesPreteritas/analysis/rna_rsubread_results_no_merge_temperature"
ref_fasta  <- file.path(ref_dir, "Gmax_880_v6.0.fa")
annot_gff3 <- file.path(ref_dir, "Gmax_880_Wm82.a6.v1.gene_exons.gff3")

threads <- 4

## Apaga os BAMs de lane logo depois de contar a amostra.
## TRUE e obrigatorio enquanto o D: estiver cheio. Os FASTQ ficam
## intactos, entao qualquer amostra pode ser realinhada depois.
delete_bams_after_count <- TRUE

## Limpa 02_merged_bam (chunks temporarios + _merged.BAM + _sorted.bam
## das amostras 150-153). Nada disso e usado nesta rota.
clean_merged_dir <- TRUE

## strandSpecific: 0 unstranded, 1 forward, 2 reverse.
## "auto" testa os tres na primeira amostra e escolhe o melhor.
strand_setting <- "auto"

lane_bam_dir <- file.path(out_dir, "01_lane_bam")
merge_dir    <- file.path(out_dir, "02_merged_bam")
count_dir    <- file.path(out_dir, "03_counts")
log_dir      <- file.path(out_dir, "04_logs")
index_base   <- file.path(out_dir, "05_index", "soybean_index")
per_sample   <- file.path(count_dir, "per_sample")
fc_tmp       <- file.path(out_dir, "tmp_fc")

for (d in c(lane_bam_dir, count_dir, log_dir, per_sample, fc_tmp))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

logfile <- file.path(log_dir, format(Sys.time(), "run_%Y%m%d_%H%M%S.log"))
say <- function(...) {
  msg <- paste0(format(Sys.time(), "[%H:%M:%S] "), ..., "\n")
  cat(msg); cat(msg, file = logfile, append = TRUE)
}


## ------------------------ DISCO / LIMPEZA ---------------------------
free_gb <- function(drive = "D") {
  x <- try(system(sprintf(
    'powershell -NoProfile -Command "(Get-PSDrive %s).Free/1GB"', drive),
    intern = TRUE), silent = TRUE)
  if (inherits(x, "try-error") || !length(x)) return(NA_real_)
  round(suppressWarnings(as.numeric(x[1])), 1)
}
dir_gb <- function(p) {
  if (!dir.exists(p)) return(0)
  f <- list.files(p, full.names = TRUE, recursive = TRUE)
  round(sum(file.info(f)$size, na.rm = TRUE) / 1024^3, 1)
}
drop <- function(files, what) {
  files <- files[file.exists(files)]
  if (!length(files)) return(invisible(NULL))
  gb <- round(sum(file.info(files)$size, na.rm = TRUE) / 1024^3, 1)
  file.remove(files)
  say("Removidos ", length(files), " ", what, " (", gb, " GB liberados)")
}

say("=== INICIO === Livre em D: ", free_gb(), " GB")
say("01_lane_bam: ", dir_gb(lane_bam_dir), " GB | 02_merged_bam: ",
    dir_gb(merge_dir), " GB")

if (clean_merged_dir && dir.exists(merge_dir)) {
  ## chunks temporarios do samtools sort (causam o erro "File exists")
  drop(list.files(merge_dir, pattern = "_sorted\\.[0-9]{3,4}(-[0-9]+)?\\.bam$",
                  full.names = TRUE), "chunks temporarios do sort")
  drop(list.files(merge_dir, pattern = "_merged\\.BAM$", full.names = TRUE),
       "BAMs merged")
  drop(list.files(merge_dir, pattern = "_sorted\\.bam(\\.bai)?$",
                  full.names = TRUE), "BAMs sorted/indices")
  say("Livre em D: apos limpeza: ", free_gb(), " GB")
}

if (!file.exists(paste0(index_base, ".00.b.array")) &&
    !length(list.files(dirname(index_base), pattern = "^soybean_index")))
  stop("Indice do genoma nao encontrado em ", dirname(index_base))


## ---------------- ANOTACAO: GFF3 Phytozome -> SAF -------------------
saf_rds <- file.path(count_dir, "Gmax_880_v6_exon_saf.rds")
if (file.exists(saf_rds)) {
  saf <- readRDS(saf_rds)
  say("SAF em cache: ", nrow(saf), " exons / ",
      length(unique(saf$GeneID)), " genes")
} else {
  say("Construindo SAF a partir do GFF3 (exon -> mRNA -> gene via Parent)...")
  ln  <- readLines(annot_gff3, warn = FALSE)
  ln  <- ln[!startsWith(ln, "#") & nzchar(ln)]
  gff <- fread(text = paste(ln, collapse = "\n"), sep = "\t", header = FALSE,
               quote = "", col.names = c("seqid","source","type","start","end",
                                         "score","strand","phase","attr"))
  rm(ln); gc()
  ga <- function(x, k) sub(paste0(".*", k, "=([^;]+).*"), "\\1", x)

  mr <- gff[type == "mRNA"]
  tx2gene <- data.table(tx = ga(mr$attr, "ID"), gene = ga(mr$attr, "Parent"))
  ex <- gff[type == "exon"]
  ex[, tx := ga(attr, "Parent")]
  ex <- merge(ex, tx2gene, by = "tx", all.x = TRUE)

  saf <- as.data.frame(unique(ex[!is.na(gene),
           .(GeneID = gene, Chr = seqid, Start = start,
             End = end, Strand = strand)]))
  saveRDS(saf, saf_rds)
  say("SAF pronto: ", nrow(saf), " exons / ",
      length(unique(saf$GeneID)), " genes")
  rm(gff, ex, mr, tx2gene); gc()
}


## ------------------------- FASTQ -> amostras -------------------------
pat <- "^(.+)_L00([1-4])_(R[12])_001\\.fastq\\.gz$"
ff  <- list.files(raw_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)
fq  <- data.frame(filepath = ff, filename = basename(ff),
                  stringsAsFactors = FALSE)
fq  <- fq[grepl(pat, fq$filename), , drop = FALSE]
fq$sample <- sub(pat, "\\1", fq$filename)
fq$lane   <- paste0("L00", sub(pat, "\\2", fq$filename))
fq$read   <- sub(pat, "\\3", fq$filename)
fq <- fq[fq$sample != "Undetermined_S0", , drop = FALSE]

## ordena 150_S1, 151_S2, ... 454_S20 pelo numero do S
ord <- order(as.integer(sub(".*_S", "", unique(fq$sample))))
samples <- unique(fq$sample)[ord]
say("Amostras: ", length(samples), " -> ", paste(samples, collapse = ", "))


## ----------------------- funcoes de trabalho -------------------------
align_sample <- function(s) {
  sub_fq <- fq[fq$sample == s, , drop = FALSE]
  bams <- character(0)
  for (ln in sort(unique(sub_fq$lane))) {
    d  <- sub_fq[sub_fq$lane == ln, , drop = FALSE]
    r1 <- d$filepath[d$read == "R1"]; r2 <- d$filepath[d$read == "R2"]
    if (length(r1) != 1 || length(r2) != 1) {
      warning("Lane incompleta: ", s, " ", ln); next
    }
    bam <- file.path(lane_bam_dir, paste0(s, "_", ln, ".BAM"))
    if (!file.exists(bam)) {
      say("  alinhando ", s, " ", ln, " (livre: ", free_gb(), " GB)")
      if (!is.na(free_gb()) && free_gb() < 15)
        stop("Menos de 15 GB livres em D:. Libere espaco antes de continuar.")
      Rsubread::align(index = index_base, readfile1 = r1, readfile2 = r2,
                      input_format = "gzFASTQ", output_file = bam,
                      output_format = "BAM", nthreads = threads,
                      phredOffset = 33, unique = FALSE, indels = 5, TH1 = 2)
    } else {
      say("  BAM ja existe: ", basename(bam))
    }
    bams <- c(bams, bam)
  }
  bams
}

count_bams <- function(bams, strand) {
  featureCounts(files = bams, annot.ext = saf, isGTFAnnotationFile = FALSE,
                useMetaFeatures = TRUE, allowMultiOverlap = FALSE,
                isPairedEnd = TRUE, countReadPairs = TRUE,
                requireBothEndsMapped = FALSE, checkFragLength = FALSE,
                strandSpecific = strand, autosort = TRUE, tmpDir = fc_tmp,
                nthreads = threads, verbose = FALSE)
}

assigned_pct <- function(fc) {
  st <- fc$stat
  round(100 * sum(st[st$Status == "Assigned", -1]) / sum(st[, -1]), 1)
}


## ------------- deteccao automatica de strandedness -------------------
strand_file <- file.path(count_dir, "strandedness.rds")
if (identical(strand_setting, "auto")) {
  if (file.exists(strand_file)) {
    strand_setting <- readRDS(strand_file)$best
    say("Strandedness em cache: ", strand_setting)
  } else {
    say("Detectando strandedness na amostra ", samples[1], " (1 lane)...")
    probe <- align_sample(samples[1])[1]
    res <- sapply(0:2, function(k) assigned_pct(count_bams(probe, k)))
    names(res) <- c("unstranded(0)", "forward(1)", "reverse(2)")
    say("  % Assigned -> ", paste(names(res), res, sep = ": ", collapse = " | "))
    strand_setting <- as.integer(which.max(res) - 1)
    saveRDS(list(best = strand_setting, pct = res), strand_file)
    say("  escolhido strandSpecific = ", strand_setting)
  }
}


## --------------------- LOOP PRINCIPAL --------------------------------
for (s in samples) {
  rds <- file.path(per_sample, paste0(s, "_counts.rds"))
  if (file.exists(rds)) { say("[skip] ", s, " ja contado"); next }

  say("---- ", s, " ----  livre em D: ", free_gb(), " GB")
  bams <- align_sample(s)
  if (!length(bams)) { warning("sem BAMs para ", s); next }

  say("  featureCounts (", length(bams), " lanes, strand=", strand_setting, ")")
  fc <- count_bams(bams, strand_setting)
  say("  Assigned: ", assigned_pct(fc), "%")

  saveRDS(list(counts = rowSums(fc$counts),
               genes  = fc$annotation$GeneID,
               length = fc$annotation$Length,
               stat   = fc$stat,
               nlanes = length(bams)), rds)
  rm(fc); gc()

  if (delete_bams_after_count) drop(bams, paste("BAMs de lane de", s))
  unlink(list.files(fc_tmp, full.names = TRUE), recursive = TRUE)
}


## --------------------- MATRIZ FINAL ----------------------------------
rds_files <- file.path(per_sample, paste0(samples, "_counts.rds"))
rds_files <- rds_files[file.exists(rds_files)]
lst <- lapply(rds_files, readRDS)
snames <- sub("_counts\\.rds$", "", basename(rds_files))

counts <- do.call(cbind, lapply(lst, `[[`, "counts"))
dimnames(counts) <- list(lst[[1]]$genes, snames)

stats <- do.call(cbind, lapply(lst, function(z) rowSums(z$stat[, -1, drop = FALSE])))
dimnames(stats) <- list(lst[[1]]$stat$Status, snames)
stats <- stats[rowSums(stats) > 0, , drop = FALSE]

## metadados: 1xx=Amb, 2xx=Elev, 3xx=Temp, 4xx=ElevTemp (CONFIRME)
grp <- substr(snames, 1, 1)
meta <- data.frame(
  sample    = snames,
  code      = sub("_S[0-9]+$", "", snames),
  treatment = factor(c("1"="Amb","2"="Elev","3"="Temp","4"="ElevTemp")[grp],
                     levels = c("Amb","Elev","Temp","ElevTemp")),
  rep       = as.integer(substr(sub("_S[0-9]+$", "", snames), 3, 3)) + 1L,
  stringsAsFactors = FALSE
)

write.csv(counts, file.path(count_dir, "counts_gene_by_sample.csv"), quote = FALSE)
write.csv(data.frame(GeneID = lst[[1]]$genes, Length = lst[[1]]$length),
          file.path(count_dir, "gene_lengths.csv"), row.names = FALSE, quote = FALSE)
write.csv(stats, file.path(count_dir, "featureCounts_summary.csv"), quote = FALSE)
write.csv(meta,  file.path(count_dir, "sample_metadata.csv"),
          row.names = FALSE, quote = FALSE)
saveRDS(list(counts = counts, length = lst[[1]]$length, stat = stats,
             meta = meta, strand = strand_setting),
        file.path(count_dir, "featureCounts_Soybean_Temp.rds"))

say("=== FIM ===  genes: ", nrow(counts), " | amostras: ", ncol(counts))
say("Livre em D: ", free_gb(), " GB")
print(meta)
cat("\n% por categoria (featureCounts):\n")
print(round(prop.table(stats, 2) * 100, 1))
cat("\nArquivos em: ", count_dir, "\n")

## Proximo passo (DE): DESeq2 ou edgeR usando counts + meta$treatment.
