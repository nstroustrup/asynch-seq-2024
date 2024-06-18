loadGenomeCov <- function(annots) {
  # load genome coverage
  dir <- paste0("../data/genomecov/", annots$Directory[1], "/")
  files <- paste0(annots$Basename, ".txt.gz")
  genomecov <- pblapply(files, function(x) {
    d <- fread(paste0(dir, x))
    colnames(d) <- c("Chromosome", "Position", "Count")
    d <- d[Chromosome == "I" &
             Position >= 15064838 - 20000 &
             Position <= 15068346 + 20000]
    d[, Basename := gsub("[.]txt[.]gz", "", x)]
    d <- merge(d, annots, all.x = TRUE, by = c("Basename"))
    d
  })
  genomecov <- rbindlist(genomecov)
  setkey(genomecov, Chromosome, Position)
  genomecov
}

plotCov <- function(x) {
  cols <- c("#0b5394ff", "grey5", "#990000ff")
  gp_cov <- ggplot() +
    geom_line(aes(Position, Coverage), x) +
    theme_classic() +
    xlab("Position in Chromosome I")  + 
    ylab("Normalized Read Counts") +
    scale_x_continuous(labels = scales::comma, breaks = seq(15063000, 15068000, by = 2000)) +
    scale_y_continuous(labels = scales::comma) +
    annotate("text", 
             x = c(mean(rrn1), mean(rrn2), mean(rrn3)), 
             y = max(x$Coverage)*1.2, 
             color = cols,
             label = c("rrn-1.1", "rrn-2.1", "rrn-3.1")) +
    geom_segment(aes(x = rrn1[1], y = 1, 
                     xend = rrn1[1], yend = max(x$Coverage)/0.9), 
                 color = cols[1]) +
    geom_segment(aes(x = rrn1[2], y = 1, 
                     xend = rrn1[2], yend = max(x$Coverage)/0.9), 
                 color = cols[1]) +
    geom_segment(aes(x = rrn2[1], y = 1, 
                     xend = rrn2[1], yend = max(x$Coverage)/0.9), 
                 color = cols[2]) +
    geom_segment(aes(x = rrn2[2], y = 1, 
                     xend = rrn2[2], yend = max(x$Coverage)/0.9), 
                 color = cols[2]) +
    geom_segment(aes(x = rrn3[1], y = 1, 
                     xend = rrn3[1], yend = max(x$Coverage)/0.9), 
                 color = cols[3]) +
    geom_segment(aes(x = rrn3[2], y = 1, 
                     xend = rrn3[2], yend = max(x$Coverage)/0.9), 
                 color = cols[3]) +
    theme(text = element_text(size = 15), legend.position = "top")
}
