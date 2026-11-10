# Install once
install.packages(c("tidyverse", "lme4", "lmerTest", "emmeans", "performance",
                   "broom.mixed", "afex", "DHARMa", "psych", "TOSTER"))

# Load every session
library(tidyverse)
library(lme4)
library(emmeans)
library(performance)
library(broom.mixed)
library(DHARMa)

set.seed(2024)  # reproducibility for any bootstrapping

# Target structure:
# participant   : factor, unique ID (e.g., "FR_01", "SR_23")
# language      : factor, "French" / "Serbian"
# condition     : factor, "pitch_accent" / "segmental"
# seq_length    : integer, 3 / 4 / 5
# block_order   : integer, position of block in experiment (1, 2, 3, ...)
# trial_num     : integer, trial position within experiment
# target_word   : factor, "proneja" / "dremina" / "fotava_fokava"
# speaker       : factor, "male" / "female"
# accent_type   : factor, "LF" / "LR" / "SF" / "SR" / NA for segmental
# target_seq    : character, e.g., "ALA"
# response_seq  : character, e.g., "ALL"

SRT_data <- read_excel("Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/EXPERIMENT2/RESULTS/AllConcat.xlsx")

SRT_data <- rename(
  SRT_data,
  participant = Subject,
  speaker = TalkerGender,
  language = LanguageGroup,
  seq_length = SequenceNumber,
  target_word = Stimulus,
  condition = Contrast,
  response_seq = Response,
  target_seq = Correct_response,
)

names(SRT_data)
n_distinct(dat$participant)

dat <- SRT_data
names(dat)
dat$condition <- as.factor(dat$condition)

# Your condition column contains: "both", "consonant", "length", "pitch"
# -> "consonant" is the segmental control; the other three are pitch-accent sub-types
dat <- dat %>%
  mutate(
    contrast_type = condition,  # preserve fine-grained: pitch / length / both / consonant
    condition = if_else(condition == "consonant", "segmental", "pitch_accent"),
    condition = factor(condition, levels = c("segmental", "pitch_accent")),
    language  = factor(language, levels = c("French", "Serbian"))
  )

# ============ ============
dat <- dat %>%
  mutate(
    response_seq = toupper(trimws(as.character(response_seq))),
    target_seq   = toupper(trimws(as.character(target_seq)))
  ) %>%
  filter(!is.na(response_seq), nchar(target_seq) > 0)

# Check: every string should now consist only of A/L (Serbian) or W/P (French)
dat %>% distinct(language, target_seq) %>% count(language)
unique(unlist(strsplit(dat$response_seq, "")))


# Sanity checks
dat %>% count(condition, contrast_type)
dat %>% distinct(participant, language) %>% count(language)
# NOTE: you previously found 52 French (Methods says 51) — resolve this:
dat %>% count(participant) %>% arrange(desc(n)) %>% head(10)   # duplicates?

# Length-mismatched responses: score as incorrect (keep them in the data)
swap <- function(x) chartr("ALWP", "LAPW", x)

dat <- dat %>%
  mutate(
    len_mismatch    = nchar(response_seq) != nchar(target_seq),
    correct_literal = as.integer(response_seq == target_seq),
    reversed_target = swap(target_seq),
    is_reversal     = as.integer(response_seq == reversed_target),
    correct_primary = pmax(correct_literal, is_reversal),
    correct_strict  = correct_literal
  )

mean(dat$len_mismatch)
dat %>% group_by(language, condition) %>% summarise(acc = mean(correct_primary), .groups = "drop")



# Descriptive table (reproduce Table 2)
desc <- dat %>%
  group_by(language, condition, seq_length) %>%
  summarise(acc = mean(correct_primary),
            n_trials = n(), .groups = "drop")
print(desc)

desc %>%
  pivot_wider(names_from = seq_length, values_from = acc, id_cols = c(language, condition)) %>%
  write_csv("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/table2_descriptives.csv")

# Whole-sequence chance: (1/2)^n  — keep seq_length NUMERIC here
chance <- tibble(seq_length = c(3, 4, 5),
                 chance     = 2 * 0.5^c(3, 4, 5))   # .25, .125, .0625

# Above-chance t-tests, pitch-accent condition (filter FIRST)
above_chance <- dat %>%
  filter(condition == "pitch_accent") %>%
  group_by(language, seq_length, participant) %>%
  summarise(acc = mean(correct_primary), .groups = "drop") %>%
  left_join(chance, by = "seq_length") %>%
  group_by(language, seq_length) %>%
  summarise(t = list(t.test(acc, mu = unique(chance), alternative = "greater")),
            .groups = "drop") %>%
  mutate(stats = map(t, broom::tidy)) %>%
  unnest(stats) %>%
  mutate(p_bonf = p.adjust(p.value, method = "bonferroni")) %>%
  dplyr::select(-t)

above_chance
write_csv(above_chance, "/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/table_above_chance.csv")

dat <- dat %>%
  mutate(seq_length_f = factor(seq_length, levels = c(3, 4, 5)))

contrasts(dat$language)     <- c(-0.5, 0.5)          # French = -0.5, Serbian = +0.5
contrasts(dat$condition)    <- c(-0.5, 0.5)          # segmental = -0.5, pitch = +0.5
contrasts(dat$seq_length_f) <- MASS::contr.sdif(3)   # 4-3, 5-4 (no library(MASS)!)

# Record these codings in the paper.

# 1) Check every factor has 2+ levels in the CURRENT dat
dat <- dat %>% droplevels()
sapply(dat[c("language", "condition", "seq_length_f", "target_word")], nlevels)
dat %>% count(language, condition, seq_length_f)   # no empty cells, no single-level factor
# If language or condition shows 1 level, the problem is upstream:
# check the Language/Contrast labels (spelling/case) and the exclusion step.

# 2) Re-apply contrasts (they must be set on the current dat)
contrasts(dat$language)     <- c(-0.5, 0.5)
contrasts(dat$condition)    <- c(-0.5, 0.5)
contrasts(dat$seq_length_f) <- MASS::contr.sdif(3)

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# 3) Refit all models explicitly on the same data (no update())
m_full <- glmer(
  correct_primary ~ language * condition * seq_length_f +
    (1 + condition | participant) + (1 | target_word),
  data = dat, family = binomial, control = ctrl
)

m_no3way <- glmer(
  correct_primary ~ language * condition + language * seq_length_f +
    condition * seq_length_f +
    (1 + condition | participant) + (1 | target_word),
  data = dat, family = binomial, control = ctrl
)

m_no2way <- glmer(
  correct_primary ~ language * seq_length_f + condition * seq_length_f +
    (1 + condition | participant) + (1 | target_word),
  data = dat, family = binomial, control = ctrl
)

isSingular(m_full)
summary(m_full)
anova(m_full, m_no3way)      # 3-way interaction
anova(m_no3way, m_no2way)    # CRITICAL Language x Condition interaction

# Coefficient table with odds ratios
tidy(m_full, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  write_csv("table_omnibus_model.csv")

emm_lc <- emmeans(m_full, ~ condition | language, type = "response")
emm_lc
pairs(emm_lc)   # pitch vs segmental WITHIN each group

emm_cl <- emmeans(m_full, ~ language | condition, type = "response")
pairs(emm_cl)   # French vs Serbian WITHIN each condition

# Interaction contrast (difference of differences)
contrast(emmeans(m_full, ~ language * condition), interaction = "pairwise")

as_tibble(emm_lc) %>% write_csv("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/table_emmeans_lang_by_cond.csv")


dat_pa <- dat %>% filter(condition == "pitch_accent") %>% droplevels()

contrasts(dat_pa$language)     <- c(-0.5, 0.5)
contrasts(dat_pa$seq_length_f) <- MASS::contr.sdif(3)

m_pa <- glmer(
  correct_primary ~ language * seq_length_f + speaker +
    (1 + seq_length_f | participant) + (1 | target_word),
  data = dat_pa, family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
# If singular: (1 | participant) + (1 | target_word), then (1 | participant)

summary(m_pa)
emmeans(m_pa, ~ language | seq_length_f, type = "response")

dat_seg <- dat %>% filter(condition == "segmental") %>% droplevels()
contrasts(dat_seg$language)     <- c(-0.5, 0.5)
contrasts(dat_seg$seq_length_f) <- MASS::contr.sdif(3)

m_seg <- glmer(
  correct_primary ~ language * seq_length_f +
    (1 | participant),
  data = dat_seg, family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

summary(m_seg)
emmeans(m_seg, ~ language | seq_length_f, type = "response")


dat_pa <- dat_pa %>%
  mutate(contrast_type = factor(contrast_type, levels = c("pitch", "length", "both")))

m_acc <- glmer(
  correct_primary ~ language * contrast_type + speaker +
    (1 | participant) + (1 | target_word),
  data = dat_pa, family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

summary(m_acc)

# Interpret via emmeans, NOT raw coefficients
emm_acc <- emmeans(m_acc, ~ contrast_type | language, type = "response")
emm_acc
pairs(emm_acc)
contrast(emmeans(m_acc, ~ language | contrast_type), method = "pairwise")

as_tibble(emm_acc) %>% write_csv("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/table_emmeans_contrast_types.csv")

# (a) Speaker as moderator
m_spk <- update(m_pa, . ~ . + language:speaker)
emmeans(m_spk, ~ language | speaker, type = "response")

# (b) Exclude female rising-accent trials (rising = contrast_type involving pitch rise;
#     adjust the filter to however rising trials are identifiable in your data)
dat_norise_f <- dat_pa %>%
  filter(!(speaker == "female" & contrast_type %in% c("pitch", "both")))
m_pa_sens1 <- update(m_pa, data = dat_norise_f)

# (c) Male-speaker trials only
m_pa_sens2 <- update(m_pa,
                     . ~ . - speaker,
                     data = dat_pa %>%
                       filter(speaker == "male") %>%
                       droplevels())

# Compare the language effect across all three fits
map_dfr(list(full        = m_pa,
             no_fem_rise = m_pa_sens1,
             male_only   = m_pa_sens2),
        ~tidy(.x, effects = "fixed", conf.int = TRUE), .id = "model") %>%
  filter(str_detect(term, "language")) %>%
  write_csv("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/table_speaker_sensitivity.csv")

m_full_strict <- update(m_full, correct_strict ~ .)

m_strict_no3way <- update(m_full_strict, . ~ . - language:condition:seq_length_f)
m_strict_no2way <- update(m_strict_no3way, . ~ . - language:condition)
anova(m_strict_no3way, m_strict_no2way)   # does Language x Condition survive?

# Side-by-side comparison of the critical interaction
bind_rows(
  tidy(m_full,        effects = "fixed") %>% mutate(scoring = "reversals_correct"),
  tidy(m_full_strict, effects = "fixed") %>% mutate(scoring = "strict")
) %>%
  filter(str_detect(term, "language.*condition")) %>%
  write_csv("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/table_scoring_sensitivity.csv")

library(TOSTER)

pa_subj <- dat_pa %>%
  group_by(participant, language) %>%
  summarise(acc = mean(correct_primary), .groups = "drop")

tsum <- pa_subj %>% group_by(language) %>%
  summarise(m = mean(acc), s = sd(acc), n = n())

tsum_TOST(m1 = tsum$m[tsum$language == "French"],
          sd1 = tsum$s[tsum$language == "French"],
          n1  = tsum$n[tsum$language == "French"],
          m2  = tsum$m[tsum$language == "Serbian"],
          sd2 = tsum$s[tsum$language == "Serbian"],
          n2  = tsum$n[tsum$language == "Serbian"],
          low_eqbound = -0.4, high_eqbound = 0.4, eqbound_type = "SMD")

rel <- dat %>%
  group_by(participant, condition) %>%
  mutate(half = if_else(row_number() %% 2 == 0, "even", "odd")) %>%
  group_by(participant, language, half, condition) %>%
  summarise(acc = mean(correct_primary), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = acc) %>%
  mutate(diff = segmental - pitch_accent) %>%
  dplyr::select(participant, half, diff) %>%
  pivot_wider(names_from = half, values_from = diff)

r_raw <- cor(rel$odd, rel$even, use = "complete.obs")
r_sb  <- (2 * r_raw) / (1 + r_raw)   # Spearman–Brown
r_raw; r_sb

sim <- simulateResiduals(m_full, n = 1000)
plot(sim)
testDispersion(sim)
check_collinearity(m_full)

sim <- simulateResiduals(m_full, n = 1000)
plot(sim)
testDispersion(sim)
check_collinearity(m_full)

# Figure 1: accuracy by group x condition x length (bars + 95% CI + individual points)
subj_means <- dat %>%
  group_by(participant, language, condition, seq_length) %>%
  summarise(acc = mean(correct_primary), .groups = "drop")

ggplot(subj_means, aes(x = factor(seq_length), y = acc, fill = language)) +
  stat_summary(fun = mean, geom = "col",
               position = position_dodge(0.9), width = 0.8) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               position = position_dodge(0.9), width = 0.2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             alpha = 0.25, size = 0.8) +
  geom_segment(data = chance,
               aes(x = as.numeric(factor(seq_length)) - 0.45,
                   xend = as.numeric(factor(seq_length)) + 0.45,
                   y = chance, yend = chance),
               inherit.aes = FALSE, linetype = "dashed") +
  facet_wrap(~condition, labeller = as_labeller(
    c(pitch_accent = "Lexical pitch accent", segmental = "Segmental control"))) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Sequence length",
       y = "Proportion correct (whole sequence; reversals scored correct)",
       fill = "Group") +
  theme_bw(base_size = 12)

ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/figure1_accuracy.pdf", width = 8, height = 4.5)

# Figure 2: pitch-accent sub-types by group
subj_acc_types <- dat_pa %>%
  group_by(participant, language, contrast_type) %>%
  summarise(acc = mean(correct_primary), .groups = "drop")

ggplot(subj_acc_types, aes(x = contrast_type, y = acc, fill = language)) +
  stat_summary(fun = mean, geom = "col",
               position = position_dodge(0.9), width = 0.8) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               position = position_dodge(0.9), width = 0.2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15,
                                             dodge.width = 0.9),
             alpha = 0.25, size = 0.8) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Pitch-accent contrast type", y = "Proportion correct", fill = "Group") +
  theme_bw(base_size = 12)

ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/figure2_accent_types.pdf", width = 7, height = 4.5)

sessionInfo()
save(dat, m_full, m_full_strict, m_pa, m_seg, m_acc,
     file = "/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/models_reversals_correct.RData")

