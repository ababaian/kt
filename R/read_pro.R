# read_pro
#' Reads a .pro file created by 'diamond'
#' @param pro.path relative system path to .fev file
#' @return A diamond-pro data.frame object
#' @keywords palmid diamond pro
#' @examples
#' # palmDB Alignment file (.pro)
#' pro.path <- system.file( "extdata", "waxsys.pro", package = "palmid")
#' pro.df <- read.pro(pro.path)
#'
#' @import dplyr ggplot2
#' @export

# Diamond Command Used
# time diamond blastx -q $filename_noz -d /$outdate.dmnd -p $THREADS \
# -c 1 --masking 0\
# --tmpdir tmp_$accession \
# --target-indexed \
# --sensitive -s 1 \
# --evalue 1e-8 \
# -k 1\
# -f 6 qseqid qstart qend qlen qstrand \
# sseqid  sstart send slen \
# pident evalue cigar \
# qseq_translated full_qseq \
# > $diamond_output_file 

read_pro <- function(pro.path) {
  # read fev as tsv
  pro.df <- utils::read.csv2(pro.path, header = F, sep = "\t",
                             stringsAsFactors=FALSE)
  
  pro.cols <- c("qseqid", "qstart", "qend", "qlen", "qstrand",
                "sseqid", "sstart", "send", "slen",
                "pident", "evalue", "cigar",
                "qseq_translated", "full_qseq")
  colnames(pro.df) <- pro.cols
  
  # set df.types
  pro.df$qseqid <- as.factor(pro.df$qseqid)
  pro.df$qstart <- as.numeric(pro.df$qstart)
  pro.df$qend   <- as.numeric(pro.df$qend)
  pro.df$qlen   <- as.numeric(pro.df$qlen)
  pro.df$qstrand<- as.factor(pro.df$qstrand)
  pro.df$sseqid <- as.character(pro.df$sseqid)
  pro.df$sstart <- as.numeric(pro.df$sstart)
  pro.df$send   <- as.numeric(pro.df$send)
  pro.df$slen   <- as.numeric(pro.df$slen)
  pro.df$pident <- as.numeric(pro.df$pident)
  pro.df$evalue <- as.numeric(pro.df$evalue)
  pro.df$cigar  <- as.character(pro.df$cigar)
  pro.df$qseq_translated  <- as.character(pro.df$qseq_translated)
  pro.df$full_qseq  <- as.character(pro.df$full_qseq)
  
  # # Initialize empty taxonomy columns (use get.palmTax)
  # pro.df$tspe <- as.character(NA)
  # pro.df$tfam <- as.character(NA)
  # pro.df$tphy <- as.character(NA)
  
  return(pro.df)
}
