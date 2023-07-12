library(tidyverse)
library(survival)
library(survminer)
library(broom)

# Figure 6, survival analysis 
## I: 90-day survival analysis between low, mediam vs high diversity in all patients (initial samples of each patient were used)
## J: 90-day survival analysis between Bifidobacterium expansion vs no expansion in patients exposed to lactulose

# Table S7: adjusted hazard ratio for Fig 6, I, from Cox proportional hazard model
# Table S8: adjusted hazard ratio for Fig 6, J, from Cox proportional hazard model

# load in data ------------------------------------------------------------

first_meta <- read_csv("./data/LD847_HD22.quant.meta.csv") %>% 
  filter(initial_sample == "Yes", cohort != "Healthy\nDonor") %>% 
  mutate(vital_status = if_else(survival_90_days == "alive",
                                0, 1))

first_lactulose <- read_csv("./data/90day_survival_demographics_230518.csv") 

# Fig 6, I --------------------------------------------------------------------

fit <- survfit(Surv(time_to_event_death_90, vital_status) ~
                 diversity,
               data = first_meta)

gsur <- ggsurvplot(fit, 
                   conf.int = F, 
                   legend.labs=c("High", "Medium", "Low"),
                   pval = TRUE,
                   risk.table = TRUE,
                   pval.coord = c(5, 0.25),
                   tables.height = 0.2,
                   tables.theme = theme_cleantable(),
                   xlim = c(0,90), 
                   break.x.by = 10,
                   size=0.7,                    
                   xlab="time (days since enrollment)",
                   ylab=" probability of transplant-free survival",
                   ggtheme = theme_bw(), 
                   palette = c("#0072B5FF",  "#E18727FF","#BC3C29FF"),
                   ylim=c(0,1),
                   surv.scale="percent",
                   tables.col="strata",
                   risk.table.col = "strata",
                   risk.table.y.text = FALSE,
                   tables.y.text = FALSE, 
                   legend.title="") 


pdf("results/6.I.LD262.90day_survival.diversity.pdf", width = 6.25, height = 5.5)
print(gsur, newpage = FALSE)
dev.off()

# Fig 6, J --------------------------------------------------------------------

fit <- survfit(Surv(time_to_event_death_90, vital_status) ~
                 lactulose_bifido_status,
               data = first_lactulose)

gsur <- ggsurvplot(fit, 
                   conf.int = F, 
                   pval = TRUE,
                   risk.table = TRUE,
                   pval.coord = c(5, 0.25),
                   tables.height = 0.2,
                   tables.theme = theme_cleantable(),
                   xlim = c(0,90), 
                   break.x.by = 10,
                   size=0.7,                    
                   xlab="time (days since enrollment)",
                   ylab=" probability of transplant-free survival",
                   ggtheme = theme_bw(), 
                   palette = c("#BC3C29FF","#0072B5FF"),
                   ylim=c(0,1),
                   surv.scale="percent",
                   tables.col="strata",
                   risk.table.col = "strata",
                   risk.table.y.text = FALSE,
                   tables.y.text = FALSE, 
                   legend.title="") 


pdf("results/6.J.LD262.90day_survival.bifido.pdf", width = 6.25, height = 5.5)
print(gsur, newpage = FALSE)
dev.off()

# Table S7 ----------------------------------------------------------------

res.cox <- coxph(Surv(time_to_event_death_90, vital_status) ~
                   diversity + 
                   meld_na_final +
                   age +
                   sex, data = first_meta)

summary(res.cox)

tidy(res.cox) %>% 
  mutate(`exp(coef)` = exp(estimate)) %>% 
  write_csv("results/7S.LD262.90day_coxph.diversity.csv")

# Table S8 ----------------------------------------------------------------

res.cox <- coxph(Surv(time_to_event_death_90, vital_status) ~
                   lactulose_bifido_status + 
                   meld_na_final +
                   age +
                   sex, data = first_lactulose)

summary(res.cox)

tidy(res.cox) %>% 
  mutate(`exp(coef)` = exp(estimate)) %>% 
  write_csv("results/8S.LD262.90day_coxph.bifido.csv")
