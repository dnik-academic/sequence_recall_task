# Packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)
library(multcomp)
library(DHARMa)
library(performance)
library(effectsize)
library(patchwork)

# ---------------------------
# Data assumptions and setup
# ---------------------------
# Data frame: df with columns
# correct (0/1), correct_strict (0/1),
# Language (French/Serbian),
# Condition (segmental/pitch-accent),
# SequenceLength (3/4/5),
# SpeakerCategory (e.g., female/male),
# ContrastType (pitch/length/both),
# Participant (ID), TargetWord (ID)
dat <- read_excel("Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/EXPERIMENT2/RESULTS/AllConcat.xlsx")
out_dir <- "/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/"

dat <- rename(
  dat,
  participant = Subject,
  speaker = TalkerGender,
  language = LanguageGroup,
  seq_length = SequenceNumber,
  target_word = Stimulus,
  condition = Contrast,
  response_seq = Response,
  target_seq = Correct_response,
)
dat <- dat %>%
  mutate(
    contrast_type = condition,  # preserve fine-grained: pitch / length / both / consonant
    condition = if_else(condition == "consonant", "segmental", "pitch_accent"),
    condition = factor(condition, levels = c("segmental", "pitch_accent")),
    language  = factor(language, levels = c("French", "Serbian"))
  )

df <- dat %>%
  mutate(
    Language = factor(language, levels = c("French", "Serbian")),
    Condition = factor(Condition, levels = c("segmental", "pitch-accent")),
    SequenceLength = factor(seq_length, levels = c(3, 4, 5), ordered = TRUE),
    SpeakerCategory = factor(speaker),
    ContrastType = factor(contrast_type, levels = c("pitch", "length", "both")),
    Participant = factor(participant),
    TargetWord = factor(target_word)
  )

# Effect coding
contrasts(df$Language) <- contr.sum(2) / 2
contrasts(df$Condition) <- contr.sum(2) / 2
# Successive differences for sequence length
contrasts(df$SequenceLength) <- MASS::contr.sdif(3)

# ---------------------------
# Omnibus model (2.1)
# ---------------------------
m_full <- glmer(
  Accuracy ~ language * condition * seq_length +
    (1 + condition | participant) + (1 | target_word),
  data = df, family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
)

# Nested models for LR tests (already reported in text, kept here for emmeans/figures)
m_no3way <- glmer(
  Accuracy ~ language * condition + language * seq_length + condition * seq_length +
    (1 + condition | participant) + (1 | target_word),
  data = df, family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
)

# ---------------------------
# Figure 1: Language × Condition (marginal over length)
# ---------------------------
emm_lc <- emmeans(m_full, ~ language * condition, type = "response")
fig1 <- as.data.frame(emm_lc) %>%
  ggplot(aes(x = condition, y = prob,
             color = language, group = language)) +
  geom_line(position = position_dodge(width = 0.25),
            linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.25),
             size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.25),
                width = 0.12,
                linewidth = 0.7) +
  scale_y_continuous(
    name = "Predicted accuracy (probability)",
    limits = c(0.5, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(x = "Condition", color = "Language") +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

fig1
ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/fig1_language_condition.png", fig1, width = 6.5, height = 4.2, dpi = 300)

# ---------------------------
# Figure 2: Language × Condition × Sequence Length
# ---------------------------

out_dir <- "/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/"

# ---------------------------
# Data (lowercase names used by all models; levels match the figures)
# ---------------------------
df <- dat %>%
  mutate(
    language    = factor(language, levels = c("French", "Serbian")),
    condition   = factor(condition, levels = c("segmental", "pitch_accent")),
    seq_length  = factor(seq_length, levels = c(3, 4, 5)),
    speaker     = factor(speaker),
    participant = factor(participant),
    target_word = factor(target_word)
  )

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

df_pitch <- df %>% filter(condition == "pitch_accent") %>% droplevels()
df_segm  <- df %>% filter(condition == "segmental")    %>% droplevels()

# ---------------------------
# Models
# ---------------------------
m_full <- glmer(
  Accuracy ~ language * condition * seq_length +
    (1 + condition | participant) + (1 | target_word),
  data = df, family = binomial(link = "logit"), control = ctrl
)

m_pitch <- glmer(
  Accuracy ~ language * seq_length + speaker +
    (1 + seq_length | participant) + (1 | target_word),
  data = df_pitch, family = binomial(link = "logit"), control = ctrl
)

m_segm <- glmer(
  Accuracy ~ language * seq_length + (1 | participant),
  data = df_segm, family = binomial(link = "logit"), control = ctrl
)

# ---------------------------
# Shared theme
# ---------------------------
theme_fig <- theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# ---------------------------
# Fig 2: Language x Condition, faceted by sequence length
# ---------------------------
emm_lcs <- emmeans(m_full, ~ language * condition | seq_length, type = "response")

fig2 <- as.data.frame(emm_lcs) %>%
  ggplot(aes(x = condition, y = prob, color = language, group = language)) +
  geom_line(position = position_dodge(width = 0.25), linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.25), size = 2.8) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position = position_dodge(width = 0.25),
                width = 0.12, linewidth = 0.7) +
  scale_y_continuous(name = "Predicted accuracy (probability)",
                     breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "Condition", color = "Language") +
  facet_wrap(~ seq_length, nrow = 1) +
  theme_fig +
  theme(strip.text = element_text(size = 11, face = "bold"))

fig2
ggsave(paste0(out_dir, "fig2_three_way_interaction.png"), fig2,
       width = 8.5, height = 3.6, dpi = 300)
# ---------------------------
# Figure 3: Forest plot of key pairwise contrasts (odds ratios)
# ---------------------------

pairs_lc      <- contrast(emm_lc, method = "revpairwise", by = "language",  adjust = "none") # within-language
pairs_between <- contrast(emm_lc, method = "revpairwise", by = "condition", adjust = "none") # between-language within condition

tidy_or <- function(contr, contrast_set) {
  as.data.frame(summary(contr, infer = TRUE, type = "link")) %>%
    mutate(
      OR       = exp(estimate),
      OR_lower = exp(asymp.LCL),
      OR_upper = exp(asymp.UCL),
      contrast_set = contrast_set
    )
}

df_within  <- tidy_or(pairs_lc,      "Within-language (Pitch-accent vs Segmental)")
df_between <- tidy_or(pairs_between, "Between-language (French vs Serbian)")

print(df_within)
print(df_between)

forest_dat <- bind_rows(df_within, df_between) %>%
  mutate(
    stratum = coalesce(
      if ("language"  %in% names(.)) as.character(language)  else NA_character_,
      if ("condition" %in% names(.)) as.character(condition) else NA_character_
    ),
    label = paste0(contrast_set, " - ", stratum),
    label = factor(label, levels = rev(unique(label)))
  )

fig3 <- forest_dat %>%
  ggplot(aes(y = label, x = OR)) +
  geom_point(size = 2.8) +
  geom_errorbar(
    aes(xmin = OR_lower, xmax = OR_upper),
    orientation = "y",
    width = 0.2,
    linewidth = 0.7
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10("Odds ratio (log scale)") +
  ylab(NULL) +
  theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

fig3
ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/fig3_forest_OR.png", fig3, width = 7, height = 4.8, dpi = 300)

# ---------------------------
# Fig 4: Pitch-accent only, sequence length x language
# ---------------------------
emm_pitch <- emmeans(m_pitch, ~ language | seq_length, type = "response")

fig4 <- as.data.frame(emm_pitch) %>%
  ggplot(aes(x = seq_length, y = prob, color = language, group = language)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.12, linewidth = 0.7) +
  scale_y_continuous(name = "Predicted accuracy (probability)",
                     breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "Sequence length", color = "Language") +
  theme_fig

fig4
ggsave(paste0(out_dir, "fig4_pitch_seq_by_lang.png"), fig4,
       width = 6.5, height = 4.2, dpi = 300)

# ---------------------------
# Fig 5: Segmental only, sequence length x language
# ---------------------------
emm_segm <- emmeans(m_segm, ~ language | seq_length, type = "response")

fig5 <- as.data.frame(emm_segm) %>%
  ggplot(aes(x = seq_length, y = prob, color = language, group = language)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.12, linewidth = 0.7) +
  scale_y_continuous(name = "Predicted accuracy (probability)",
                     breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "Sequence length", color = "Language") +
  theme_fig

fig5
ggsave(paste0(out_dir, "fig5_segmental_seq_by_lang.png"), fig5,
       width = 6.5, height = 4.2, dpi = 300)
# ---------------------------
# Contrast Type in pitch-accent (3.4)
# ---------------------------
m_pitch_ct <- glmer(
  Accuracy ~ language * contrast_type + speaker +
    (1 | participant) + (1 | target_word),
  data = df_pitch, family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
)

emm_ct <- emmeans(m_pitch_ct, ~ contrast_type | language, type = "response")
fig6 <- as.data.frame(emm_ct) %>%
  ggplot(aes(x = contrast_type, y = prob, fill = language)) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.65
  ) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.7),
    width = 0.15,
    linewidth = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "French" = "#4477AA",
      "Serbian" = "#EE9944"
    )
  ) +
  scale_y_continuous(
    name = "Predicted accuracy (probability)",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(x = "Contrast type", fill = "Language") +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

fig6
ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/fig6_contrast_type_by_language.png", fig6, width = 7, height = 4.4, dpi = 300)

# ---------------------------
# Speaker Category sensitivity (3.5.1)
# ---------------------------
emm_speak <- emmeans(m_pitch, ~ language | speaker, type = "response")
fig7 <- as.data.frame(emm_speak) %>%
  ggplot(aes(x = speaker, y = prob,
             color = language, group = language)) +
  geom_line(
    position = position_dodge(width = 0.25),
    linewidth = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.25),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.25),
    width = 0.12,
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = c(
      "French" = "#4477AA",
      "Serbian" = "#EE9944"
    )
  ) +
  scale_y_continuous(
    name = "Predicted accuracy (probability)",
    limits = c(0.5, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(
    x = "Speaker category",
    color = "Language"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

fig7
ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/fig7_speaker_category.png", fig7, width = 6.5, height = 4.2, dpi = 300)

# ---------------------------
# Chance performance visualization (S1)
# ---------------------------
# Participant-level accuracy by Language × SequenceLength, all trials
df_part_acc <- df %>%
  group_by(participant, language, seq_length) %>%
  summarise(acc = mean(Accuracy), .groups = "drop")

figS1 <- df_part_acc %>%
  ggplot(aes(seq_length, acc, color = language)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
               position = position_dodge(width = 0.4), width = 0.12) +
  stat_summary(fun = mean, geom = "point",
               position = position_dodge(width = 0.4), size = 2.8) +
  scale_y_continuous("Accuracy", limits = c(0.5, 1)) +
  xlab("Sequence length") +
  theme_bw()
figS1
ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/figS1_chance_performance.png", figS1, width = 6.5, height = 4.2, dpi = 300)

# ---------------------------
# Equivalence test visualization (3.5.3)
# ---------------------------
# Participant-level pitch-accent accuracy
df_pitch_part <- df_pitch %>%
  group_by(participant, language) %>%
  summarise(acc = mean(Accuracy, na.rm = TRUE), .groups = "drop")

acc_stats <- df_pitch_part %>%
  group_by(language) %>%
  summarise(
    mu = mean(acc, na.rm = TRUE),
    sd = sd(acc, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = language,
    values_from = c(mu, sd, n)
  )

mu_diff <- acc_stats$mu_French - acc_stats$mu_Serbian
sd_pooled <- sqrt(((acc_stats$n_French - 1)*acc_stats$sd_French^2 + (acc_stats$n_Serbian - 1)*acc_stats$sd_Serbian^2) /
                    (acc_stats$n_French + acc_stats$n_Serbian - 2))
se_diff <- sd_pooled * sqrt(1/acc_stats$n_French + 1/acc_stats$n_Serbian)
# 90% CI for TOST visualization
t_90 <- qt(0.95, df = acc_stats$n_French + acc_stats$n_Serbian - 2)
ci90_low <- mu_diff - t_90 * se_diff
ci90_high <- mu_diff + t_90 * se_diff
eq_bound <- 0.40 * sd_pooled

fig8 <- tibble(
  diff = mu_diff,
  ci_low = ci90_low,
  ci_high = ci90_high,
  lower_bound = -eq_bound,
  upper_bound = eq_bound
) %>%
  ggplot(aes(x = 1, y = diff)) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = c(-eq_bound, eq_bound),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.7
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.05,
    linewidth = 0.9,
    color = "#4477AA"
  ) +
  geom_point(
    size = 3.5,
    color = "#4477AA"
  ) +
  scale_x_continuous(
    breaks = 1,
    labels = "French − Serbian",
    limits = c(0.7, 1.3)
  ) +
  ylab("Difference in mean accuracy (90% CI)") +
  xlab("Language comparison") +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

fig8
ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/fig8_equivalence_tost.png", fig8, width = 5.2, height = 4.2, dpi = 300)

# ---------------------------
# Split-half reliability (3.5.4)
# ---------------------------
# Create odd/even halves within Participant × Condition
df_halves <- df %>%
  arrange(participant, condition) %>%
  group_by(participant, condition) %>%
  mutate(
    idx = row_number(),
    half = if_else(idx %% 2 == 1, "odd", "even")
  ) %>%
  ungroup()

half_scores <- df_halves %>%
  group_by(participant, condition, half) %>%
  summarise(
    acc = mean(Accuracy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = acc
  ) %>%
  mutate(
    diff = `segmental` - `pitch_accent`
  ) %>%
  dplyr::select(participant, half, diff) %>%
  pivot_wider(
    names_from = half,
    values_from = diff
  )

# Spearman–Brown corrected reliability
r_half <- cor(half_scores$odd, half_scores$even, use = "complete.obs")
r_sb <- (2 * r_half) / (1 + r_half)

fig9 <- half_scores %>%
  ggplot(aes(odd, even)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  annotate("text", x = min(half_scores$odd, na.rm = TRUE),
           y = max(half_scores$even, na.rm = TRUE),
           hjust = 0, vjust = 1,
           label = paste0("r = ", sprintf("%.3f", r_half),
                          "\nSpearman–Brown = ", sprintf("%.3f", r_sb))) +
  xlab("Odd-half (segmental − pitch-accent)") +
  ylab("Even-half (segmental − pitch-accent)") +
  theme_bw()
fig9
ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/fig9_split_half_reliability.png", fig9, width = 5.2, height = 4.2, dpi = 300)

# ---------------------------
# Model diagnostics (S2)
# ---------------------------
simres <- simulateResiduals(m_full, n = 1000)
p1 <- plotSimulatedResiduals(simres) # base plotting
# Save component plots separately using DHARMa helpers
png("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/figS2_dharma_residuals.png", width = 1500, height = 1000, res = 170)
par(mfrow = c(2, 2))
plot(simres)                      # residuals vs. predicted, QQ, etc.
testUniformity(simres, plot = TRUE)
testDispersion(simres, plot = TRUE)
testZeroInflation(simres, plot = TRUE)
dev.off()

# ---------------------------
# Stricter scoring: omnibus model for sensitivity (3.5.2)
# ---------------------------
m_full_strict <- glmer(
  correct_strict ~ language * condition * seq_length +
    (1 + condition | participant) + (1 | target_word),
  data = df, family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))
)

# Optional: a quick figure analogous to Fig 1 under strict scoring
emm_lc_strict <- emmeans(m_full_strict, ~ language * condition, type = "response")
fig1_strict <- as.data.frame(emm_lc_strict) %>%
  ggplot(aes(x = condition, y = prob,
             color = language, group = language)) +
  geom_line(
    position = position_dodge(width = 0.25),
    linewidth = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.25),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.25),
    width = 0.12,
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = c(
      "French" = "#4477AA",
      "Serbian" = "#EE9944"
    )
  ) +
  scale_y_continuous(
    name = "Predicted accuracy (strict scoring)",
    limits = c(0.5, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(x = "Condition", color = "Language") +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

fig1_strict
ggsave("/Users/dusannikolic/Documents/PAPERS/2026/SEQUENCE_RECALL_TASK/RESULTS_NEW/fig1s_language_condition_strict.png", fig1_strict, width = 6.5, height = 4.2, dpi = 300)

# ---------------------------
# Save all figures to a single directory if desired
# ---------------------------
# dir.create("figures", showWarnings = FALSE)
# file.copy(list.files(pattern = "png$"), "figures", overwrite = TRUE)

####################################################
#########################################################