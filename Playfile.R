#Helper function
snp_in_gene <- function(pos, gene_starts, gene_ends, tol) {
  return(any(pos >= gene_starts - tol & pos <= gene_ends + tol))
}
#Helper function
snp_in_gene_vectorized <- Vectorize(snp_in_gene, vectorize.args = "pos")

# Function for determining if a SNP is in a gene
# Or generally if all rows in a dataframe with columns Chromosome and position is within
# a window in a dataframe defined by columns Chromosome, start and end +- tol
# Where tol is the length outside the window that is included in the window
snp_in_gene_allchromosomes <- function(snps, genes, tol) {
  log_vector <- c()
  snps %>%
    arrange(Chromosome, Position) -> snps
  for (chrom in unique(snps$Chromosome)) {
    gene_starts <- genes$start[genes$Chromosome == chrom]
    gene_ends <- genes$end[genes$Chromosome == chrom]
    log_vector <- c(log_vector, snp_in_gene_vectorized(snps$Position[snps$Chromosome == chrom],
                                                       gene_starts, gene_ends, tol))
  }
  snps %>%
    mutate(in_gene = log_vector) %>%
    return()
}

#Helper function
window_contains_SNP <- function(pos, window_starts, window_ends, tol) {
  return(any(pos >= window_starts - tol & pos <= window_ends + tol))
}
# Helper function
window_contains_SNP_vectorized <- Vectorize(window_contains_SNP, vectorize.args = c("window_starts", "window_ends"))

# Function for determining if a window contains a SNP
# Or generally if all rows in a dataframe windows with columns Chromosome, start, end
# contains position in a dataframe defined by columns Chromosome, Position +- tol
# Where tol is the length outside the window that is included in the window
window_contains_SNP_allchromosomes <- function(snps, windows, tol) {
  log_vector <- c()
  
  windows %>%
    arrange(Chromosome, start) -> windows
  for (chrom in unique(windows$Chromosome)) {
    positions <- snps$Position[snps$Chromosome == chrom]
    window_starts <- windows$start[windows$Chromosome == chrom]
    window_ends <- windows$end[windows$Chromosome == chrom]
    log_vector <- c(log_vector, window_contains_SNP_vectorized(snps$Position[snps$Chromosome == chrom],
                                                              window_starts, window_ends, tol))
  }
  windows %>%
    mutate(contains_snp = log_vector) %>%
    return()
}

#Helper function
find_windows <- function(pos, chrom, chroms, starts, ends, tol) {
  return(which(pos >= starts - tol & pos <= ends + tol & chroms == chrom))
}
#Helper function
find_windows_vectorized <- Vectorize(find_windows, vectorize.args = c("chroms", "starts", "ends"))
#Helper function
find_windows_vectorized_snps <- Vectorize(find_windows, vectorize.args = c("pos", "chrom"))


# Function for mutate a dataframe with SNPs to include the values 
# Found in column column in a dataframe windows with columns Chromosome, start, end +- tol
get_window_vals_snpwise <- function(snps, windows, column, tol) {
  window_indices <- find_windows_vectorized_snps(snps$Position, snps$Chromosome,
                                                 windows$Chromosome, windows$start, windows$end, tol)
  my_vals <- windows[column][window_indices]
  snps %>%
    mutate(window_vals = my_vals) %>%
    return()
}

# Function for making a manhattan plot for multiple traits
# Takes a dataframe with columns Chromosome, Position and P.value
# Will include roomie position as a dashed line if include_roomie is TRUE
myManhattan_multitrait <- function(x, include_roomie = TRUE) {
  
  
  if (!all(c("Chromosome", "P.value", "Position") %in% colnames(x))) {
    stop("Input must be data frame with colnames Chromosome, P.value and Position")
  }
  
  x$Chromosome <- as.factor(x$Chromosome)
  x$P.value <- -log10(x$P.value)
  bonferroni_cutoff <- -log10(0.05/nrow(x))
  
  
  
  
  x %>% 
    # Compute chromosome size
    group_by(Chromosome) %>% 
    summarise(chr_len=max(Position)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len + 5*10^7 * as.numeric(Chromosome)) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(x, by = "Chromosome") %>%
    
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate(BPcum=Position+tot) %>%
    mutate( significant=ifelse(P.value < bonferroni_cutoff, T, F)) -> x
  
  if (include_roomie) {
    x %>%
      filter(Chromosome == 2) %>%
      mutate(pos = Position - 78606912) %>%
      filter(pos == min(pos)) %>%
      pull(BPcum) -> roomie_pos
  }
  
  axisdf = x %>%
    group_by(Chromosome) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  breaksdf = x %>%
    group_by(Chromosome) %>%
    summarize(center=( max(BPcum) ) )
  
  label <- unique(breaksdf$Chromosome)
  
  p <- ggplot(x, aes(x=BPcum, y=P.value)) +
    
    # Show all points
    geom_point( aes(color=Trait), alpha=0.8, size=3, show.legend = F) +
    
    # custom X axis:
    scale_x_continuous(label = label, breaks= axisdf$center, expand = c(0, 0)) +
    viridis::scale_color_viridis(discrete = T) +
    geom_vline(xintercept = roomie_pos, linetype = "dashed", color = "red", linewidth = 1) +
    # Add thresholds
    # Labels
    ylab("-log10(P.value)") +
    xlab("Chromosome") +
    # Custom the theme:
    theme_bw() +
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 12),
      title = element_text(size = 10),
      text = element_text(face = "bold")
    )
  return(p)
}

# Function for making a manhattanplot for a single trait
# Takes a dataframe with columns Chromosome, Position and P.value
# Will include roomie position as a dashed line if include_roomie is TRUE
myManhattan <- function(x, include_roomie = TRUE) {
  
  
  if (!all(c("Chromosome", "P.value", "Position") %in% colnames(x))) {
    stop("Input must be data frame with colnames Chromosome, P.value and Position")
  }
  
  x$Chromosome <- as.factor(x$Chromosome)
  x$P.value <- -log10(x$P.value)
  bonferroni_cutoff <- -log10(0.05/nrow(x))
  
  
  
  
  x %>% 
    # Compute chromosome size
    group_by(Chromosome) %>% 
    summarise(chr_len=max(Position)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len + 10^7 * as.numeric(Chromosome)) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(x, by = "Chromosome") %>%
    
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate(BPcum=Position+tot) %>%
    mutate( significant=ifelse(P.value < bonferroni_cutoff, T, F)) -> x
  
  if (include_roomie) {
    x %>%
      filter(Chromosome == 2) %>%
      mutate(pos = Position - 78606912) %>%
      filter(pos == min(pos)) %>%
      pull(BPcum) -> roomie_pos
  }
  
  axisdf = x %>%
    group_by(Chromosome) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  breaksdf = x %>%
    group_by(Chromosome) %>%
    summarize(center=( max(BPcum) ) )
  
  label <- unique(breaksdf$Chromosome)
  
  p <- ggplot(x, aes(x=BPcum, y=P.value)) +
    
    # Show all points
    geom_point( aes(color=Chromosome), alpha=0.8, size=3, show.legend = FALSE) +
    
    # custom X axis:
    scale_x_continuous(label = label, breaks= axisdf$center, expand = c(0, 0)) +
    viridis::scale_color_viridis(discrete = T) +
    geom_vline(xintercept = roomie_pos, linetype = "dashed", color = "red", linewidth = 1) +
    # Add thresholds
    # Labels
    ylab("-log10(P.value)") +
    xlab("Chromosome") +
    # Custom the theme:
    theme_bw() +
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 12),
      title = element_text(size = 10),
      text = element_text(face = "bold")
    )
  return(p)
}

# Function for making a manhattan-like plot for a single trait (no -log10 transformation)
# Takes a dataframe with columns Chromosome, Position and P.value
# Will include roomie position as a dashed line if include_roomie is TRUE
myManhattan_like <- function(x, include_roomie = TRUE) {
  
  
  if (!all(c("Chromosome", "P.value", "Position") %in% colnames(x))) {
    stop("Input must be data frame with colnames Chromosome, P.value and Position")
  }
  
  x$Chromosome <- as.factor(x$Chromosome)
  
  
  
  
  x %>% 
    # Compute chromosome size
    group_by(Chromosome) %>% 
    summarise(chr_len=max(Position)) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len + 10^7 * as.numeric(Chromosome)) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(x, by = "Chromosome") %>%
    
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate(BPcum=Position+tot) -> x
  
  if (include_roomie) {
    x %>%
      filter(Chromosome == 2) %>%
      mutate(pos = Position - 78606912) %>%
      filter(pos == min(pos)) %>%
      pull(BPcum) -> roomie_pos
  }
  
  axisdf = x %>%
    group_by(Chromosome) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  breaksdf = x %>%
    group_by(Chromosome) %>%
    summarize(center=( max(BPcum) ) )
  
  label <- unique(breaksdf$Chromosome)
  
  p <- ggplot(x, aes(x=BPcum, y=P.value)) +
    
    # Show all points
    geom_point( aes(color=Chromosome), alpha=0.8, size=3, show.legend = FALSE) +
    
    # custom X axis:
    scale_x_continuous(label = label, breaks= axisdf$center, expand = c(0, 0)) +
    viridis::scale_color_viridis(discrete = T) +
    geom_vline(xintercept = roomie_pos, linetype = "dashed", color = "red", linewidth = 1) +
    # Add thresholds
    # Labels
    ylab("Score") +
    xlab("Chromosome") +
    # Custom the theme:
    theme_bw() +
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 12),
      title = element_text(size = 10),
      text = element_text(face = "bold")
    )
  return(p)
}

################################################################################
# Load packages
pacman::p_load(tidyverse, data.table, viridis)
################################################################################
# GWAS filtering steps

genes <- read_delim("20210713_Lj_Gifu_v1.3_predictedGenes.gff3", delim = "\t",
                    col_names = F, skip = 10) %>%
  filter(X3 == "gene") %>%
  select(X1, X4, X5) %>%
  rename(Chromosome = X1, start = X4, end = X5) %>%
  mutate(Chromosome = as.numeric(str_remove(Chromosome, "LjG1.1_chr"))) %>%
  na.omit()

# Filters:
Methods <- c("FarmCPU", "Gemma", "ISIS EM-BLASSO") #Which methods to include
p_value <- 10^-6 #P.value filter
filter_genes <- TRUE #Whether to only include SNPs within x bp of gene
gene_range = 1000 #Range to include genes (in both directions)
snps_per_trait <- 200 #Max number of significant unique SNPs per trait
round_position <- -4 #How much to round position of snps of later filters.
                     # 0, -1, -2, -3, -4 corresponds to rounding to nearest 1, 10, 100, 1000, 10000 bp
reps_per_trait <- 2 #How many methods/trials/snps per rounded position per trait to keep
Trials_plus_traits_across_traits <- 3 #How many trials plus traits per rounded position to keep (minimum is 2)

# Remove snps below p value of 10^-6 and BLINK significants
sig <- read_csv("20240508_Lotus_rarefied_GWAS_Results_sig.csv") %>%
  filter(Method %in% Methods) %>%
  filter(P.value < p_value)

# Remove SNPs not within 1000 bp of a gene
sig %>%
  snp_in_gene_allchromosomes(genes, gene_range) %>%
  filter(!filter_genes & in_gene) %>%
  select(!in_gene) -> sig_gene

# Split SNPs into bacterial gwas and non_bacterial_gwas
sig_gene %>%
  mutate(Trial = str_extract(Trait, "(?<=_)(ave)|[1-3]$"),
         Trait = str_remove(Trait, "_((ave)|[1-3])$")) -> sig_gene
sig_gene %>%
  filter(is.na(Trial) & Trait != "OW_2017") -> other_gwas
sig_gene %>%
  filter(!is.na(Trial)) -> sig_bact

# Remove traits with more than 200 significant SNPs
sig_bact %>%
  count(Trait, SNP) %>%
  count(Trait) %>%
  filter(n < snps_per_trait) %>%
  select(!n) %>%
  left_join(sig_bact) -> sig_bact_traits

# Rounds positions to nearest 10000
sig_bact_traits %>%
  mutate(new_position = round(Position, round_positions)) -> try
# Remove SNPs with no significant SNPs at the same rounded position for the same trait
# (i.e. remove SNPs that are not significant for another replicate or method or closeby snp)
try %>%
  count(Chromosome, new_position, Trait) %>%
  filter(n >= reps_per_trait) %>%
  select(!n) %>%
  left_join(try) -> sig_bact_traits_snps

# Remove SNPs that are not significant at the same rounded position for all traits
# for a different trial or method
sig_bact_traits_snps %>%
  group_by(Chromosome, new_position) %>%
  summarize(n = length(unique(Trait)) + length(unique(Trial))) %>%
  filter(n >= reps_across_traits) %>%
  left_join(sig_bact_traits_snps) %>%
  ungroup() -> sig_bact_traits_snps_consistent

# Rounds position to nearest 10000 and groups traits
other_gwas %>%
  mutate(new_pos = round(Position, -4),
         old_trait = Trait, 
         Trait = case_when(Trait %like% "temp" | Trait %like% "OW" ~ "Temperature",
                           Trait %like% "Seed" ~ "Seed",
                           Trait %like% "ft" | Trait %like% "FT" | Trait %like% "FP" ~ "Flowering",
                           Trait %like% "contr" | Trait %like% "salt" ~ "Salt",
                           .default = Trait)) -> other_gwas

# Filter SNPs that are not significant twice for the same trait in the same rounded position
# Can be 2 methods, 2 replicates or nearby snp
other_gwas %>%
  count(Chromosome, new_pos, Trait) %>%
  filter(n > 1) %>%
  left_join(other_gwas) %>%
  mutate(Trait = old_trait) %>%
  select(!c(in_gene, old_trait, Trial)) -> other_gwas_filtered

################################################################################
# Make differential density plots

# Expression for values to compare
# Needs to evaluate to a boolean for logistic regression
My_expression = expr(My expression here)

# Connect SNPs to windows
SNPs %>%
  window_contains_SNP_allchromosomes(windows, 0) -> windows_withsnps

windows_withsnps %>%
  mutate(expression = eval(My_expression)) -> windows_withdifferential

# Make density plots for K2-K8

window_withdifferential %>%
  select(c(V2:V8, expression)) %>%
  pivot_longer(cols = V2:V8, names_to = "K", values_to = "Variance") %>%
  mutate(K = str_replace(K, "V", "K = ")) %>%
  ggplot(aes(x = Variance, fill = expression, color = expression)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~K) +
  theme_minimal(base_size = 15) +
  scale_fill_viridis(discrete = T) +
  scale_color_viridis(discrete = T) +
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75)) +
  scale_y_continuous(breaks = c(0, 2, 4)) +
  theme(panel.grid = element_blank())
  
# Make normal density plot

# Which column to compare?
my_column = expr(My column goes here)
  
windows_withdifferential %>%
  ggplot(aes(x = my_column)) +
  geom_density(alpha = 1, color = expression, fill = expression) +
  theme_minimal(base_size = 30) +
  scale_fill_viridis(discrete = T)
  scale_color_viridis(discrete = T) +
  theme(panel.grid = element_blank())

################################################################################
# Logistic regression

windows_withdifferential %>%
  mutate(Score = ifelse(expression, 1, 0)) -> check

your_model = as.formula("Score ~ Your model here")
your_null_model = as.formula("Score ~ Your null model here")
    
glm(your_model, data = check, family = "binomial") -> model
glm(your_null_model, data = check, family = "binomial") -> null
  
anova(model, null, test = "Chisq")