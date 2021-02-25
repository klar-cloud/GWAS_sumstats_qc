### DESCRIPTION

# This script finds abnormal differences in allelic freqencies
# between two traits and removes these alleles

# For traits that were already processed with this script the version of 
# data after correction will be taken

# Directories structure for this script:
# -<main data directory>
#   - sumstats_qc (f) (open/save data)
#     - plots (f) (save plots)

### PACKAGES

suppressPackageStartupMessages(library(rlist))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(optparse))


######################## PARAMETERS #########################
option_list = list(  
  make_option("--sumstats", type='character', default=NULL, help='A path to summary statistics after first part of preprocessing; example: "path/to/file/trait_prep1.csv"'),
  make_option("--reference", type='character',default=NULL, help='A path to (unzipped) reference file; file needs to have A1_REF, A2_REF, FREQ_REF and COORDINATE column (or CHR&BP columns). Additionally if you plan to use it for adding rs numbers it must have RSID column); only space/tab-separated values; example: "path/to/file/reference.csv"'),
  make_option("--out", type='character', default=NULL, dest='save_dir', help='A name of the trait you are processing. Will be used to generate save name.'),
  make_option("--density-plot", action="store_true", default=FALSE, dest='density_plot', help='Do you want to create a density plot to visualize freqency outliers?'),
  make_option("--n-nonoutliers", type='integer', default=1000000, dest='n_nonoutliers',  help='How many frequency non-outliers SNPs do you want in visualisation (default=1000000)'),
  make_option("--sd-outliers", type='integer', default=10, dest='sd_outliers', help='How far SNPs have to be in distribution to be cathegorized as ouliers? (in SD; default=10SD)'),
  make_option("--check-rsid", action="store_true", dest='check_rsid', help='Do you want to examine rs numbers column, remove duplicates or create one is mission?'),
  make_option("--method-rsid", type='character', dest='method_rsid', default='bioconductor', help='What method do you want to use to add rs numbers? Default is package "bioconductor", but you can also use your "reference" if it has RSID and COORDINATE(or CHR&BP) columns'),
  make_option("--infere-beta", action="store_true", dest='infere_beta', help='Do you want to use frequency, n and z-score to infere missing beta/sd?'),
  make_option("--add-zscore", action="store_true", dest='add_zscore', help='Do you want to calculate z-score from beta and sd if missing?'),
  make_option("--test", action="store_true", default=FALSE, help='Do you want to load only first 10,000 SNPs from sumstats and max 100,000 SNPs from reference for testing?'))
  
opt <- parse_args(OptionParser(option_list=option_list))

load_sumstats       <- opt$sumstats
load_reference      <- opt$reference
save_dir            <- opt$save_dir

save_density_plot   <- opt$density_plot
params_n_outl       <- opt$n_nonoutliers
params_sd_outl      <- opt$sd_outliers

params_check_rsid   <- opt$check_rsid
params_method_rsid  <- opt$method_rsid
params_infere_beta  <- opt$infere_beta
params_add_zscore   <- opt$add_zscore
params_test         <- opt$test

# CREATE LOG FILE
save_log = paste0(save_dir, '.qc2.log')
sink(save_log, append=FALSE, split=TRUE)

# PRINT PARSERS
cat('Argument       :', 'Value\n')
cat('.......................\n')
cat('sumstats       :', opt$sumstats, '\n')
cat('reference      :', opt$reference, '\n')
cat('save_dir       :', opt$save_dir, '\n')
cat('density_plot   :', opt$density_plot, '\n')
cat('n_nonoutliers  :', opt$n_nonoutliers, '\n')
cat('sd_outliers    :', opt$sd_outliers, '\n')
cat('check_rsid     :', opt$check_rsid, '\n')
cat('method_rsid    :', opt$method_rsid, '\n')
cat('infere_beta    :', opt$infere_beta, '\n')
cat('add_zscore     :', opt$add_zscore, '\n')
cat('test           :', opt$test, '\n\n')

# CHECK IF LOAD FILES EXIST
if (is.null(load_sumstats)){
  stop('Please specify pathway to summary statistics\n') } else
    if (!(file.exists(load_sumstats))){
      stop('The specified load path for summary statistics does not exist\n')
    }
if (is.null(load_reference)){
  stop('Please specify pathway to reference\n') } else
    if (!(file.exists(load_reference))){
      stop('The specified load path for reference does not exist\n')
    }


if(params_method_rsid=='bioconductor'){
  suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP144.GRCh37))
}

#############################################################

### FUNCTIONS
# function to fix alleles:
# flip: allele on opposite strand A <-> T (same effect allele)
# swap: A1 <-> A2 (different effect allele, so effect estimate is in the opposite direction, so is the allele frequency)
# palindromic SNPs, A/T and C/G cannot be solved this way

fix_alleles = function(a1_ref, a2_ref, a1_match, a2_match){
  
  m_alleles = as.matrix(cbind(a1_ref, a2_ref, a1_match, a2_match))
  m_num = matrix(NA, nrow=nrow(m_alleles), ncol=ncol(m_alleles))
  m_num[m_alleles == "A"] = -1
  m_num[m_alleles == "C"] = -2
  m_num[m_alleles == "T"] = 1
  m_num[m_alleles == "G"] = 2
  m_flip = m_num * -1
  palin = ifelse(m_num[,3] == m_flip[,4], 1, 0)
  asis  = ifelse(m_num[,1] == m_num[,3] &
                   m_num[,2] == m_num[,4], 1, 0)
  swap  = ifelse((m_num[,1] == m_num[,4] &
                    m_num[,2] == m_num[,3]) |
                   (m_flip[,1] == m_num[,4] &
                      m_flip[,2] == m_num[,3]), 1, 0)
  flip  = ifelse((m_flip[,1] == m_num[,3] &
                    m_flip[,2] == m_num[,4]) |
                   (m_flip[,1] == m_num[,4] &
                      m_flip[,2] == m_num[,3]), 1, 0)
  action = ifelse(palin == 1 | is.na(m_num[,3]) | is.na(m_num[,4]), "excl",
                  ifelse(asis == 1, "asis",
                         ifelse(swap == 1, "swap",
                                ifelse(swap == 1 & flip == 1, "flipswap",
                                       ifelse(flip == 1, "flip", "error")))))
  
  return(action)
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#FUNCTIONS
# Create a new column with rsid numnber based on chromosome and basepair column
add_rsid <- function(trait){
  trait[,'RSID'] <- NA
  snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
  
  for(chr in 1:22){
    chr_str = paste0(chr)
    chr_snps <- snpsBySeqname(snps, chr_str)
    trait <- trait %>% mutate(RSID = ifelse(CHR==chr, mcols(chr_snps)$RefSNP_id[match((BP), start(chr_snps))], RSID))
    cat('Added RS numbers for chromosom ', chr, '/', '22\n')
    if (chr %% 5 == 0){
      rm(chr_snps)
      invisible(gc())
    }
  }
  
  return(trait)  
}

distinct_withNA <- function(data, colname, quocolname){
  na <- data[is.na(data[[colname]]),]
  not_na <- data[!is.na(data[[colname]]),] %>%
    distinct(!!quocolname, .keep_all = TRUE)
  
  distinct_data <- rbindlist(list(not_na, na))
  return(distinct_data)
}

### SCRIPT
cat("Quality check of summary statistics dataset (part 2)", "\n\n")

# Load summary 
if (params_test){
  load_nrows_sumstats = 10000
  load_nrows_ref      = 100000

} else {
  load_nrows_sumstats = Inf
  load_nrows_ref      = Inf
}

cat('Loading sumstats\n')
trait <- fread(load_sumstats, fill=TRUE, nrow=load_nrows_sumstats)
num_row <- nrow(trait)
cat('Loaded', num_row, 'SNPs from sumstats\n')

# Load reference panel
cat('Loading reference panel\n')
ref = fread(load_reference, nrow=load_nrows_ref)
num_row <- nrow(ref)
cat('Loaded', num_row, 'SNPs from reference panel\n')

if (!("COORDINATE" %in% colnames(ref))) { 
  if (!("CHR" %in% colnames(ref)) & "BP" %in% colnames(ref)) {
    stop('The refrence file lacks necessary columns (COORDINATE or CHR&BP)')
    sink()
  } else {
    ref <- mutate(ref, COORDINATE = paste0(as.character(ref$CHR), ':', as.character(ref$BP))) 
    }
}

if (!("FREQ_REF" %in% colnames(ref))){
  stop('The refrence file lacks necessary column (FREQ)')
  sink()
}
if (!("A1_REF" %in% colnames(ref)) | !("A2_REF" %in% colnames(ref))){
  stop('The refrence file lacks necessary columns (A1, A2)')
  sink()
}
if (params_method_rsid=='reference' & !("RSID" %in% colnames(ref))){
  stop('The refrence file lacks necessary column (RSID)')
  sink()
}

# Merge reference panel with sumstats
cat('Merging sumstats with reference panel \n')
num_row <- nrow(trait)
trait_ref <- inner_join(ref, trait, by='COORDINATE')

rm(trait)
if (!(params_method_rsid=='reference')){
  rm(ref)
}
invisible(gc()) 

if (nrow(trait_ref)-num_row == 0){
  cat('All SNPs from dataset are present in the reference panel\n')
} else {
  cat("Removed", (nrow(trait_ref)-num_row)*-1, "SNPs during merging of datset with reference panel\n")
}

# Check how alleles match reference alleles
cat('\nFixing mismatching alleles\n')
trait_ref <- trait_ref %>% mutate(ACTION = fix_alleles(A1_REF, A2_REF, A1, A2)) 

if ('ZSCORE' %in% colnames(trait_ref)){ 
  trait_ref <- trait_ref %>% 
    mutate(ZSCORE = ifelse(ACTION %in% c("swap", "flipswap"), -1 * ZSCORE,
                           ifelse(ACTION %in% c("excl", "error"), NA, ZSCORE)))}


if ("BETA" %in% colnames(trait_ref)){ 
  trait_ref <- trait_ref %>% 
    mutate(BETA = ifelse(ACTION %in% c("swap", "flipswap"), -1 * BETA,
                         ifelse(ACTION %in% c("excl", "error"), NA, BETA)))}

trait_ref <- trait_ref %>%
  mutate(FREQ = ifelse(ACTION %in% c("swap", "flipswap"), 1 - FREQ,
                       ifelse(ACTION %in% c("excl", "error"), NA, FREQ)))

trait_ref <- select(trait_ref, -c("ACTION", "A1", "A2"))

num_row <- nrow(trait_ref)
trait_ref <- trait_ref[complete.cases(trait_ref), ]
removed_row = (nrow(trait_ref) - num_row)*-1

cat("Removed", removed_row, "palindromic SNPs")

invisible(gc())

# FIND OUTLIERS
cat('\nSearching for frequency outliers\n')

# filter SNPs that are shared between traits, have frequency higher than 0.01, and have abnormal frequency difference
trait_ref <- trait_ref %>%
  #filter(! is.na(FREQ) & ! is.na(FREQ_REF)) %>%
  #filter(FREQ > 0.01 | FREQ_REF > 0.01) %>%
  mutate(DIFF1 = abs(FREQ - FREQ_REF) / sqrt(FREQ * (1-FREQ))) %>%
  mutate(DIFF2 = abs(FREQ - FREQ_REF) / sqrt(FREQ_REF * (1-FREQ_REF))) %>%
  mutate(OUT = ifelse(DIFF1 > (mean(DIFF1) + 10 * sd(DIFF1)), 1,
                      ifelse(DIFF2 > (mean(DIFF2) + 10 * sd(DIFF2)), 1,  0)))


# Check whether differences between traits are not extremaly small
# If yes, skip the rest of the script
if (max(trait_ref$DIFF1) < 0.1 &
    max(trait_ref$DIFF2) < 0.1){
  invisible(gc())
  cat('ATTENTION: Differences between frequencies datset and reference are very small.\nSkipped removing freq outliers.\n\n' )
  
  trait <- select(trait_ref, -c("FREQ_REF", "DIFF1", "DIFF2", "OUT"))
  rm(trait_ref)
  invisible(gc())
  
  names(trait)[names(trait) == 'A1_REF'] <- 'A1'
  names(trait)[names(trait) == 'A2_REF'] <- 'A2'
  
} else {
  
  trait_ref <- select(trait_ref, -c("DIFF1", "DIFF2"))
  
  invisible(gc())
  
  outliers <- trait_ref %>% 
    select(c("COORDINATE", "FREQ", "FREQ_REF", "OUT")) %>%
    filter(OUT==1) %>%
    mutate(DENSITY = 0)

  # PLOTTING
  if (!is.null(save_density_plot)){
    # Check whether if number of alleles for plotting is not bigger than
    # the total on non-ouliers; if yes, take the total_number of outliers
    total_nonoutliers <- nrow(trait_ref %>% 
                                filter(OUT==0)) 
    
    if (total_nonoutliers<params_n_outl){
      params_n_outl<-total_nonoutliers
    } 
    
    # Find non-outliers
    non_outliers <- trait_ref %>%
      select(c("COORDINATE", "FREQ", "FREQ_REF", "OUT")) %>% 
      filter(OUT==0) %>% 
      sample_n(params_n_outl) %>%
      mutate(DENSITY = get_density(FREQ, FREQ_REF, n=100)) 
    
    cat('Plotting frequency density\n')
    save_name_plot = paste0(save_dir, '.outliers_plot.png')
    
    png(save_name_plot, width=630, h=480)
    print({
      p <- ggplot() +
        geom_point(d=outliers, aes(FREQ, FREQ_REF), color="red") +
        geom_point(d=non_outliers, aes(FREQ, FREQ_REF, color = DENSITY)) +
        scale_color_viridis() +
        theme_classic()
      p + xlab(paste('trait frequency')) + ylab(paste('reference frequency'))
    })
    dev.off() 
    
    cat('Succesfully ploted frequency density\n')
    
    rm(non_outliers)
    invisible(gc())
    
  }else{
    invisible(gc())
  }

  
  # REMOVE OUTLIERS from sumstats
  trait <- select(trait_ref, -c("OUT", "FREQ_REF"))
  
  rm(trait_ref)
  invisible(gc())
  
  names(trait)[names(trait) == 'A1_REF'] <- 'A1'
  names(trait)[names(trait) == 'A2_REF'] <- 'A2'
  
  trait <- trait %>% filter(!(COORDINATE %in% as.vector(outliers$COORDINATE)))
  cat('Removed', nrow(outliers), 'freqency outliers\n\n')
  
  rm(outliers)
}

### DESCRIPTION

# This script can: 
# - load summary statistics
# - check rs identificators column, creates one if necessary and removes duplicates
# - calculates z-score from beta and se
# - infere beta and se from zscore, freq and n
# - remove SNPs with low MAF
# - remove rows with missing values
# - remove unecessary columns
# - rename column coordinate to id

# Directories structure for this script:
# -<main data directory>
#   - sumstats_raw (f) (open data)
#   - sumstats_postqc (f) (save data)

### PACKAGES

params_sort=TRUE

#Check rs id column
if(params_check_rsid){
  if ("RSID" %in% colnames(trait)) {
    cat('Detected rs numbers column \n')
    
    nrow_trait = nrow(trait)
    trait <- distinct_withNA(trait, colname="RSID", quocolname=quo(RSID))
    removed_row = (nrow(trait)-nrow_trait)*-1
    
    cat("Removed", removed_row, "SNPs with repeated rs numbers\n")
    
    # Add rs if rsid is just number (because of the bug)
    if (is.numeric(trait$RSID)){
      trait <- trait %>% mutate(RSID = ifelse(is.na(trait$RSID), NA, paste0('rs', as.character(trait$RSID))))
      cat('Corrected RSID column\n')
    }
  }
  # ADD missing rs id colum

  if (!('RSID' %in% colnames(trait))) {
    cat("No rs numbers column, creating one.\n")

        # using external bioconductor package
    if (params_method_rsid=='bioconductor'){
      cat("Adding rs numbers using bioconductor...\n")
      trait <- add_rsid(trait)
      invisible(gc())
      cat("Created rs number column\n")
    }
    
    # using reference panel (it must have at least two columns: RSID, COORDINATE)
    if (params_method_rsid=='reference'){
      cat("Searching rs numbers in reference file\n")
    
      ref = ref[, c("RSID", "COORDINATE")]
        
      num_row <- nrow(trait)
      trait <- inner_join(ref, trait, by='COORDINATE')
      rm(ref)
      invisible(gc())

      cat("\nRemoved", (nrow(trait)-num_row)*-1, "SNPs during merging of" , trait_name, "with reference panel\n")
      cat("Created rs number column\n")
    }
    
    nrow_trait = nrow(trait)
    trait<- distinct_withNA(trait, colname="RSID", quocolname = quo(RSID))
    removed_row = (nrow(trait)-nrow_trait)*-1
    cat("Removed", removed_row, "SNPs with repeated rs numbers\n")
  }
  # Change any 'weird' rs numbers into NA
  trait <- trait %>% mutate(RSID = ifelse(grepl("^(rs)\\d+", trait$RSID), trait$RSID, NA))
  
  cat("Detected", sum(is.na(trait$RSID)), 'SNPs with missing rs number\n')
}

# Calculate z-score
if (params_add_zscore &
    'BETA' %in% colnames(trait) &
    'SE' %in% colnames(trait) &
    !('ZSCORE' %in% colnames(trait))){
  trait <- trait %>% mutate(ZSCORE = BETA/SE)
  cat('Calculated ZSCORE\n')
}

# Infere beta/se
if (params_infere_beta &
    !('BETA' %in% colnames(trait)) &
    'ZSCORE' %in% colnames(trait) & 
    'N' %in% colnames(trait) &
    'FREQ' %in% colnames(trait) ){
  
  trait <- trait %>% mutate(SE = 1 / sqrt(2 * FREQ * (1 - FREQ) * (N + ZSCORE^2))) %>% 
    mutate(BETA = ZSCORE * SE)
  cat('Inferred BETA and SE from ZSCORE\n')
}

# Sort columns
sort_order <- c("RSID", "COORDINATE", "CHR", "BP", "A1", "A2", "N", "FREQ", "PVAL", "ZSCORE", "BETA", "SE" )
miss_colnames <- setdiff(sort_order, colnames(trait))

if (!(length(miss_colnames)==0)){
  for (col in 1:length(miss_colnames)){
    new_sort_order <- sort_order[sort_order!=miss_colnames[col]] } 
} else { new_sort_order <- sort_order }
trait <- subset(trait, select=c(sort_order))

# Saving
save_trait <- paste0(save_dir, '.qc2.csv')
fwrite(trait, save_trait, col.names=T, row.names=F, quote=F, sep="\t", na="NA")

cat('\nSaved', nrow(trait), 'quality checked SNPs\n')

cat('Sucessfully finished quality check\n\n')

head(trait, n=1)

rm(trait)
invisible(gc())

sink()