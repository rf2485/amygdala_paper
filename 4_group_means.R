source("1_data_preparation.R")
library(arsenal)
library(ggpubr)

#define remove_outliers function
remove_outliers <- function(x, na.rm = TRUE) 
{
  ## Find 25% and 75% Quantiles using inbuild function
  quant <- quantile(x, probs=c(.25, .75), na.rm = na.rm)
  
  ## Find Interquantile range and multiply it by 1.5 
  ## to derive factor for range calculation
  H <- 1.5 * IQR(x, na.rm = na.rm)
  
  y <- x
  
  ## fill the outlier elements with NA
  y[x < (quant[1] - H)] <- NA
  y[x > (quant[2] + H)] <- NA
  
  y
}

#import stats tables
mtr_aparc_aseg <- read_tsv("freesurfer/mtr_aparc+aseg.tsv") %>%
  rename(participant_id=`Measure:mean`)
mtr_wm <- read_tsv("freesurfer/mtr_wm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_wm=Seg0001)
mtr_gm <- read_tsv("freesurfer/mtr_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_gm=Seg0001)
mtr_ctx_wm <- read_tsv("freesurfer/mtr_ctx_wm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_ctx_wm=Seg0001)
mtr_ctx_gm <- read_tsv("freesurfer/mtr_ctx_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_ctx_gm=Seg0001)
mtr_lh_ctx_gm <- read_tsv("freesurfer/mtr_lh_ctx_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_lh_ctx_gm=Seg0001)
mtr_rh_ctx_gm <- read_tsv("freesurfer/mtr_rh_ctx_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_rh_ctx_gm=Seg0001)
mtr_subcort_gm <- read_tsv("freesurfer/mtr_subcort_gm.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_subcort_gm=Seg0001)
mtr_ad_sig <- read_tsv("freesurfer/mtr_AD_sig.tsv") %>%
  rename(participant_id=`Measure:mean`, mtr_ad_sig=Seg0001)

failed_qc <- c('sub-CC510255', #abnormality in left temporal pole
               'sub-CC620821', #segmentation errors
               'sub-CC621011', #segmentation errors 
               'sub-CC711027', #severe motion artifacts, segmentation errors
               'sub-CC710551', #segmentation errors, weird contrast
               'sub-CC721434' #segmentation errors
)

#combine demographics and mtr tables
mtr_tr30 <- left_join(mti_over_55_tr30, mtr_aparc_aseg) %>%
  left_join(., mtr_wm) %>%
  left_join(., mtr_gm) %>%
  left_join(., mtr_ctx_wm) %>%
  left_join(., mtr_ctx_gm) %>%
  left_join(., mtr_lh_ctx_gm) %>%
  left_join(., mtr_rh_ctx_gm) %>%
  left_join(., mtr_subcort_gm) %>%
  left_join(., mtr_ad_sig) %>%
  filter(!participant_id %in% failed_qc) %>% #remove rows that failed QC
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)
mtr_tr50 <- left_join(mti_over_55_tr50, mtr_aparc_aseg) %>%
  left_join(., mtr_wm) %>%
  left_join(., mtr_gm) %>%
  left_join(., mtr_ctx_wm) %>%
  left_join(., mtr_ctx_gm) %>%
  left_join(., mtr_lh_ctx_gm) %>%
  left_join(., mtr_rh_ctx_gm) %>%
  left_join(., mtr_subcort_gm) %>%
  left_join(., mtr_ad_sig) %>%
  filter(!participant_id %in% failed_qc) %>% #remove rows that failed QC
  mutate(across(where(is.double), remove_outliers)) #remove outliers (change to NA)

#save to file
write_csv(mtr_tr30, 'final_mtr_tr30.csv')
write_csv(mtr_tr50, 'final_mtr_tr50.csv')

mtr_tr30$SCD <- factor(mtr_tr30$SCD, #change SCD column into factor type
                                levels = c(1, 0),
                                labels = c('SCD', 'Control'))
mtr_tr50$SCD <- factor(mtr_tr50$SCD, #change SCD column into factor type
                       levels = c(1, 0),
                       labels = c('SCD', 'Control'))

#generate t test tables
results <- "asis"
mtr_table_tr30 <- tableby(SCD ~ mtr_ad_sig + `Left-Hippocampus` + `Right-Hippocampus` + 
                            mtr_gm +  mtr_wm + mtr_ctx_gm + mtr_ctx_wm +
                            mtr_lh_ctx_gm + mtr_rh_ctx_gm + mtr_subcort_gm +
                            `Left-Cerebral-White-Matter` + `Right-Cerebral-White-Matter`,
                      data = mtr_tr30, total = FALSE)
summary(mtr_table_tr30, text = TRUE) #view
mtr_table_tr50 <- tableby(SCD ~ mtr_ad_sig + `Left-Hippocampus` + `Right-Hippocampus` + 
                            mtr_gm +  mtr_wm + mtr_ctx_gm + mtr_ctx_wm +
                            mtr_lh_ctx_gm + mtr_rh_ctx_gm + mtr_subcort_gm +
                            `Left-Cerebral-White-Matter` + `Right-Cerebral-White-Matter`,
                          data = mtr_tr50, total = FALSE)
summary(mtr_table_tr50, text = TRUE) #view

ggboxplot(mtr_tr30, x = "SCD", y = "mtr_ad_sig",
          error.plot = "errorbar",add = "jitter", add.params = list(color = "SCD"), 
          legend = "none", combine = T) +
  stat_compare_means(label.x.npc = "center", hide.ns = T,
                     aes(label = paste0("p = ", after_stat(p.format)))) +
  labs(y = "MTR AD Sig" 
  ) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank()
  )
