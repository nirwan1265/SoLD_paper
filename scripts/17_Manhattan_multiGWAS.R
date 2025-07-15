library(vroom)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(CMplot)

# 1) point to your folder
dir_path <- "/Users/nirwantandukar/Documents/Research/results/SAP/GWAS_results/raw_gwas"

# 2) grab all .txt files, name them by the basename sans “.txt”
files <- list.files(dir_path, pattern="\\.txt$", full.names=TRUE)
names(files) <- basename(files) %>% str_remove("\\.txt$")

# 3) read the first file as “master” for rs/chr/ps
master <- vroom(files[[1]]) %>% 
  select(rs, chr, ps)

# 4) for each file, read just rs + p_wald, rename p_wald → file‑name
pval_list <- imap(files, ~ vroom(.x) %>%
                    dplyr::select(rs, p_wald) %>%
                    rename(!!.y := p_wald)
)

# 5) find the SNPs common to *all* files
common_snps <- reduce(pval_list, inner_join, by = "rs") %>% pull(rs)

# 6) build your final table:
final_df <- master %>%
  filter(rs %in% common_snps) %>%      # keep only shared SNPs
  inner_join(reduce(pval_list, inner_join, by = "rs"), by = "rs") %>%
  select(rs, chr, ps, everything())

# 7) sanity check
dim(final_df)      # should be ~27k rows × (3 + number_of_files) cols
head(final_df)

final_df <- na.omit(final_df)
colnames(final_df) <- c("SNP","Chromosome","Position",
                         "Sum_TG","TG_54_6","TG_54_7","TG_56_4","TG_56_6")


final_df <- final_df %>%
  mutate(
    SNP        = factor(SNP),
    Chromosome = factor(Chromosome),
    Position   = as.integer(Position)
  )

# verify
str(final_df)


# 1) Define your threshold and the trait column names
thr    <- 1e-7
traits <- c("Sum_TG", "TG_54_6", "TG_54_7","TG_56_4","TG_56_6")


# 2) Build a list where each element is the vector of SNPs passing that trait’s cutoff
highlight_list <- lapply(traits, function(tr) {
  final_df[final_df[[tr]] < thr, "SNP"]
})

# 3) Name each element after its trait
names(highlight_list) <- traits


# 1) Pull each element’s SNP vector as character
snp_vecs <- lapply(highlight_list, function(df) as.character(df$SNP))

# 2) Intersect them all to get the common SNPs
common_SNPs <- Reduce(intersect, snp_vecs)

# 3) (Optional) If you want to keep only those in each sub‑list:
highlight_list_common <- lapply(highlight_list, function(df) {
  df %>% filter(as.character(SNP) %in% common_SNPs)
})



# SNPs <-  final_df[
#   final_df$Sum_TG < 1e-7 |
#     final_df$TG_54_6 < 1e-7 |
#     final_df$TG_54_7 < 1e-7 |
#     final_df$TG_56_4 < 1e-7 |
#     final_df$TG_56_6 < 1e-7, 1]
# 

CMplot(final_df,type="p",plot.type="m",LOG10=TRUE,highlight=highlight_list_common,highlight.type="l",
         threshold=1e-7,threshold.col="black",threshold.lty=1,col=c("grey80","black"),
         signal.cex=0.5, signal.col="red", highlight.col="green",highlight.cex=0.7,
         file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

