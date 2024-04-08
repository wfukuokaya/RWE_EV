pacman::p_load(missForest, lubridate, gdata, survminer, gtsummary, WeightIt, cobalt, smd, survival, ggsurvfit, patchwork, contsurvplot, ggsci)
library(tidyverse)
set.cobalt.options(binary = "std")

# define ggplot2 theme...
theme_jikei <- function(base_size = 10, 
                        dark_text = "#1A242F") {
  mid_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 7)[2]
  
  light_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 7)[4]
  
  theme_grey(base_size = base_size) +
    theme(text = element_text(color = mid_text, lineheight = 1.1),
          plot.title = element_text(color = dark_text, size = rel(1.2)),
          plot.subtitle = element_text(size = rel(1.1)),
          axis.text.y = element_text(color = mid_text, size = rel(1)),
          axis.title.y = element_text(size = rel(1)),
          axis.text.x = element_text(color = mid_text, size = rel(1)),
          axis.title.x = element_text(size = rel(1)),
          legend.position = "top",
          legend.title = element_blank(),
          panel.grid = element_line(color = "#F3F4F5"),
          plot.caption = element_text(size = rel(1))) 
}

d0 # raw data
dae # data for adverse events
dae_grade # data for adverse events with grade
djoined # dae + dae_grade

d # d0 + djoined = full analysis set

# stacked bar chart for objective response rate
d %>% count(bor) # calculate number of patients according to the responses
prop.res <- d %>% group_by(bor) %>% summarise(n = n()) %>% mutate(freq = n / sum(n)*100)
prop.res$bor <- casefold(prop.res$bor, upper = TRUE)

prop.res$res <- factor(prop.res$bor, levels = c("PD", "SD", "PR", "CR")) # change order

bar <- ggplot(aes(x = 1, y = freq, fill = res, label = round(freq, 1)), data = prop.res) + 
  geom_bar(stat = "identity") + 
  geom_text(position = position_stack(vjust = 0.5), color = "white", size = 5) + 
  scale_fill_npg() + 
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x=element_blank(),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    legend.title = element_blank()
  ) + 
  labs(y = "Proportion of patients with the best overall response")
bar # show stacked bar chart

# baseline characteristics
tbl1 <- d %>%
  select(age, sex, bmi, ps, smoking, mi, chd, cvd, pvd, dementia, copd, ld, dm, malg, m_diag, local, num, cpi, prim, hist_cat, symptom, chemo_cisplatin, metcat, metvol_score, liver, hb, bellmunt_cat) %>%
  mutate_at(.vars = c("sex", "smoking", "m_diag", "local", "cpi", "prim", "hist_cat", "symptom", "chemo_cisplatin", "metcat", "liver"),
            .funs = str_to_sentence) %>% 
  mutate(
    sex = sex %>% factor(levels = c("Male", "Female")),
    smoking = smoking %>% factor(levels = c("Current", "Former", "Never")),
    cpi = cpi %>% factor(levels = c("Pembrolizumab", "Avelumab")),
    hist_cat = hist_cat %>% factor(levels = c("Urothelial carcinoma", "Uc with variant histology", "Pure non-uc")),
    metcat = metcat %>% factor(levels = c("Ln only", "Visceral metastasis", "Ln and visceral metastasis"))
  ) %>% 
  tbl_summary(
    label = list(
      age ~ "Age, year",
      sex ~ "Sex",
      bmi ~ "Body mass index",
      ps ~ "ECOG performance status",
      smoking ~ "Smoking status",
      mi ~ "Previous history of myocardial infarction",
      chd ~ "Previous history of conjestive heart disease",
      cvd ~ "Previous history of cerebrovascular disease",
      pvd ~ "Previous history of peripheral vascular disease",
      dementia ~ "Previous history of dementia",
      copd ~ "Previous history of COPD",
      ld ~ "Previous history of liver disease",
      dm ~ "Previous history of diabetes mellitus",
      malg ~ "Presence of other malignancy",
      m_diag ~ "de novo metastatic disease",
      local ~ "Previous history of local treatment",
      num ~ "Number of previous systemic therapy",
      cpi ~ "Previous checkpoint inhibitor",
      prim ~ "Primary location of tumor",
      hist_cat ~ "Tumor histology",
      symptom ~ "Symptomatic disease",
      metcat ~ "Location of metastasis",
      chemo_cisplatin ~ "Previous cisplatin-based chemotherapy",
      metvol_score ~ "Number of viscera with metastasis",
      liver ~ "Liver metastasis",
      hb ~ "Baseline hemoglobin concentration",
      bellmunt_cat ~ "Bellmunt risk score"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}-{p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") 

# patients died within 1 month after treatment initiation
tbl1_early <- d %>% 
  filter(death == 1) %>% 
  filter(fu < 1) %>% 
  select(age, sex, bmi, ps, smoking, mi, chd, cvd, pvd, dementia, copd, ld, dm, malg, m_diag, local, num, cpi, prim, hist_cat, symptom, chemo_cisplatin, metcat, metvol_score, liver, hb, bellmunt_cat) %>%
  mutate_at(.vars = c("sex", "smoking", "m_diag", "local", "cpi", "prim", "hist_cat", "symptom", "chemo_cisplatin", "metcat", "liver"),
            .funs = str_to_sentence) %>% 
  mutate(
    sex = sex %>% factor(levels = c("Male", "Female")),
    smoking = smoking %>% factor(levels = c("Current", "Former", "Never")),
    cpi = cpi %>% factor(levels = c("Pembrolizumab", "Avelumab")),
    hist_cat = hist_cat %>% factor(levels = c("Urothelial carcinoma", "Uc with variant histology", "Pure non-uc")),
    metcat = metcat %>% factor(levels = c("Ln only", "Visceral metastasis", "Ln and visceral metastasis"))
  ) %>% 
  tbl_summary(
    label = list(
      age ~ "Age, year",
      sex ~ "Sex",
      bmi ~ "Body mass index",
      ps ~ "ECOG performance status",
      smoking ~ "Smoking status",
      mi ~ "Previous history of myocardial infarction",
      chd ~ "Previous history of conjestive heart disease",
      cvd ~ "Previous history of cerebrovascular disease",
      pvd ~ "Previous history of peripheral vascular disease",
      dementia ~ "Previous history of dementia",
      copd ~ "Previous history of COPD",
      ld ~ "Previous history of liver disease",
      dm ~ "Previous history of diabetes mellitus",
      malg ~ "Presence of other malignancy",
      m_diag ~ "de novo metastatic disease",
      local ~ "Previous history of local treatment",
      num ~ "Number of previous systemic therapy",
      cpi ~ "Previous checkpoint inhibitor",
      prim ~ "Primary location of tumor",
      hist_cat ~ "Tumor histology",
      symptom ~ "Symptomatic disease",
      metcat ~ "Location of metastasis",
      chemo_cisplatin ~ "Previous cisplatin-based chemotherapy",
      metvol_score ~ "Number of viscera with metastasis",
      liver ~ "Liver metastasis",
      hb ~ "Baseline hemoglobin concentration",
      bellmunt_cat ~ "Bellmunt risk score"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}-{p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") 

# orr
theme_gtsummary_journal(journal = "jama")
oruni <- tbl_uvregression(
  d[c("age", "sex", "bmi", "ps_cat", "smoking_cat", 
      "local", "num", "cpi", "prim", "hist_cat", 
      "symptom", "metcat", "chemo_cisplatin", "metvol_score", "liver",
      "hb", "bellmunt", "or_bin")],
  method = glm,
  y = or_bin, 
  method.args = list(family = binomial),
  exponentiate = TRUE
)

# treatment-related adverse events
dae_all <- djoined %>% 
  select(trae) %>% 
  mutate_all(str_to_sentence) %>% 
  tbl_summary(
    statistic = list(all_continuous() ~ "{median} ({p25}-{p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)), 
    sort = list(everything() ~ "frequency")
  )

dae_g3 <- djoined %>% 
  filter(grade == "grade 3" | grade == "grade 4" | grade == "grade 5") %>% 
  select(trae) %>% 
  mutate_all(str_to_sentence) %>% 
  tbl_summary(
    statistic = list(all_continuous() ~ "{median} ({p25}-{p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)), 
    sort = list(everything() ~ "frequency")
  )

tbl_ae <- tbl_merge(tbls = list(dae_all, dae_g3), 
                    tab_spanner = c("**All grade**", "**Grade 3 or higher**"))

# dose reductions
d_dose <- d %>% filter(!is.na(dose_rdctn))

# kaplan-meier estimations
pfs_fit <- survfit(Surv(pfs, prog) ~ 1, data = d)
pfs_plot_crude <- ggsurvplot(pfs_fit, 
                             risk.table = TRUE, 
                             title = "Progression-free survival",
                             xlab = "Months since treatment initiation", 
                             ylab = "PFS probability (95% CI)",
                             legend = "none",
                             censor = TRUE, censor.shape = "O", censor.size = 2.2,
                             palette = "lancet", size = 0.5,  break.time.by = 2,
                             ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                             risk.table.col = "strata",
                             tables.height = 0.12, risk.table.fontsize = 3.0,
                             tables.theme = survminer::theme_cleantable(), conf.int = TRUE, conf.int.fill = "#00468BFF") 

pfs_plot_crude$table <- pfs_plot_crude$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

os_fit <- survfit(Surv(fu, death) ~ 1, data = d)
os_plot_crude <- ggsurvplot(os_fit, 
                            risk.table = TRUE, 
                            title = "Overall survival",
                            xlab = "Months since treatment initiation", 
                            ylab = "OS probability (95% CI)",
                            legend = "none",
                            censor = TRUE, censor.shape = "O", censor.size = 2.2,
                            palette = "lancet", size = 0.5,  break.time.by = 2,
                            ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                            risk.table.col = "strata",
                            tables.height = 0.12, risk.table.fontsize = 3.0,
                            tables.theme = survminer::theme_cleantable(), conf.int = TRUE, conf.int.fill = "#00468BFF") 

os_plot_crude$table <- os_plot_crude$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# merge...
merge <- arrange_ggsurvplots(list(pfs_plot_crude, os_plot_crude),
                             nrow = 1, ncol = 2, print = FALSE)

# time on treatment
tot_fit <- survfit(Surv(tot, discont_bin) ~ 1, data = d)
tot_plot_crude <- ggsurvplot(tot_fit, 
                             risk.table = TRUE, 
                             title = "Time on treatment",
                             xlab = "Months since treatment initiation", 
                             ylab = "TOT probability",
                             legend = "none", 
                             censor = TRUE, censor.shape = "O", censor.size = 2.2,
                             palette = "lancet", size = 0.5,  break.time.by = 2,
                             ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                             risk.table.col = "strata",
                             tables.height = 0.12, risk.table.fontsize = 3.0,
                             tables.theme = survminer::theme_cleantable(), conf.int = TRUE, conf.int.fill = "#00468BFF") 

tot_plot_crude$table <- tot_plot_crude$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# survival information
library(pec)
fu_time <- quantile(prodlim::prodlim(Hist(fu, death) ~ 1, data = d, reverse = TRUE))

theme_gtsummary_journal(journal = "jama")
surv_data <- tbl_survfit(list(pfs_fit, os_fit, tot_fit), prob = 0.5) # median survival probability

# kaplan-meier estimations
pfs_fit <- survfit(Surv(pfs, prog) ~ or, data = d)
pfs_plot_crude <- ggsurvplot(pfs_fit, 
                             risk.table = TRUE, 
                             title = "Progression-free survival",
                             xlab = "Months since treatment initiation", 
                             ylab = "PFS probability",
                             legend = "top", legend.labs = c("SD-PD", "PR-CR"),
                             legend.title = "", 
                             censor = TRUE, censor.shape = "O", censor.size = 2.2,
                             palette = "lancet", size = 0.5,  break.time.by = 2,
                             ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                             risk.table.col = "strata",
                             tables.height = 0.12, risk.table.fontsize = 3.0,
                             tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

pfs_plot_crude$table <- pfs_plot_crude$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

os_fit <- survfit(Surv(fu, death) ~ or, data = d)
os_plot_crude <- ggsurvplot(os_fit, 
                            risk.table = TRUE, 
                            title = "Overall survival",
                            xlab = "Months since treatment initiation", 
                            ylab = "OS probability",
                            legend = "top", legend.labs = c("SD-PD", "PR-CR"),
                            legend.title = "", 
                            censor = TRUE, censor.shape = "O", censor.size = 2.2,
                            palette = "lancet", size = 0.5,  break.time.by = 2,
                            ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                            risk.table.col = "strata",
                            tables.height = 0.12, risk.table.fontsize = 3.0,
                            tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

os_plot_crude$table <- os_plot_crude$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# merge...
merge <- arrange_ggsurvplots(list(pfs_plot_crude, os_plot_crude), 
                             nrow = 1, ncol = 2, print = FALSE)

# log-rank test
pval_pfs <- survdiff(Surv(pfs, prog) ~ or, data = d) %>% print(digit = 4)
pval_os <- survdiff(Surv(fu, death) ~ or, data = d) %>% print(digit = 4)