library(tidyverse)
library(ggbeeswarm)
library(cowplot)
library(ggtext)
library(ZamanianLabThemes)
library(here)

data <- read_csv(here('fig3', 'data', 'dim_bpa_stacks_treated.csv'), col_types = 'nccccccnnnccc')

ci <- data %>% 
  mutate(
    ci_total = (top_worms - bottom_worms) / total_worms,
    ci = (top_worms - bottom_worms) / (top_worms + bottom_worms),
    species = case_when(
      species == 'Bpa' ~ '*B. pahangi*',
      species == 'Dim' ~ '*D. immitis*'
    )) %>%  
  pivot_longer(c('ci', 'ci_total'), names_to = 'measure', values_to = 'value')

stat_layer <- ci %>% 
  filter(
    top_worms + bottom_worms > 0,
    # remove the first two reps that didn't have controls
    experiment_date > 20220621,
    # measure == 'ci',
  ) %>%
  mutate(treatment = case_when(
           treatment == 'AZS' ~ 'Albendazole sulfoxide',
           treatment == 'IVM' ~ 'Ivermectin',
           TRUE ~ treatment
         ),
         conc = case_when(
           conc == '1p' ~ '1%',
           conc == '10uM' ~ '10 µM',
           conc == '100uM' ~ '100 µM'
         ),
         treatment = str_c(conc, '<br>', treatment)) %>%
  mutate(treatment = fct_relevel(treatment, "1%<br>DMSO")) %>% 
  group_by(measure, species) %>% 
  group_nest() %>% 
  mutate(otm = map(data, ~ PMCMRplus::kwManyOneDunnTest(.x$value ~ as.factor(.x$treatment), p.adjust.method = "bonferroni"))) %>% 
  mutate(summary = map(otm, ~ summary(.x))) %>%
  mutate(pval = map(summary, ~ as_tibble(pluck(.x, 3), rownames = "treatment"))) %>%
  mutate(estimate = map(summary, ~ as_tibble(pluck(.x, 4), rownames = "treatment"))) %>%
  mutate(tidy = map2(pval, estimate, ~ left_join(.x, .y, by = "treatment"))) %>%
  unnest(tidy) %>%
  select(-data, -otm, -summary, -pval, -estimate, pval = `1%<br>DMSO.x`, estimate = `1%<br>DMSO.y`) %>%
  mutate(sig = case_when(
    pval <= 0.0001 ~ "****",
    pval <= 0.001 ~ "***",
    pval <= 0.01 ~ "**",
    pval <= 0.05 ~ "*",
    pval > 0.05 ~ ""
  )) %>% 
  mutate(value = 1.1,
         measure = case_when(
           measure == 'ci' ~ '(Cue - Control) / (Cue + Control)',
           measure == 'ci_total' ~ '(Cue - Control) / Total'
         ))

(treated <- ci %>% 
    filter(
      top_worms + bottom_worms > 0,
      # remove the first two reps that didn't have controls
      experiment_date > 20220621,
      # measure == 'ci',
    ) %>%
    mutate(
      measure = case_when(
        measure == 'ci' ~ '(Cue - Control) / (Cue + Control)',
        measure == 'ci_total' ~ '(Cue - Control) / Total'
      ),
      treatment = case_when(
        treatment == 'AZS' ~ 'Albendazole sulfoxide',
        treatment == 'IVM' ~ 'Ivermectin',
        TRUE ~ treatment
      ),
      conc = case_when(
        conc == '1p' ~ '1%',
        conc == '10uM' ~ '10 µM',
        conc == '100uM' ~ '100 µM'
      ),
      treatment = str_c(conc, '<br>', treatment)
    ) %>% 
    # mutate(treatment = factor(treatment, levels = c('DMSO', 'Albendazole sulfoxide', 'Ivermectin'))) %>% 
    ggplot(aes(x = treatment, y = value)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_quasirandom(
      shape = 21, fill = alpha('black', 0.65), color = alpha('black', 0)
    ) +
    stat_summary(
      fun.data = mean_cl_normal, 
      # size = 0.6,
      color = 'indianred'
    ) +
    geom_text(data = stat_layer,
              aes(label = sig)) +
    annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf)+
    labs(x = 'Treatment', y = 'Chemotaxis index') +
    scale_y_continuous(limits = c(-1.1, 1.1)) +
    facet_grid(cols = vars(measure), rows = vars(species), scales = 'free_x') +
    theme_minimal() +
    theme(
      strip.text.x = element_markdown(size = 8),
      strip.text.y = element_markdown(size = 8),
      axis.text.x = element_markdown(size = 8, angle = 0, hjust = 0.5),
      axis.text.y = element_markdown(size = 8),
      axis.title.x = element_markdown(size = 9, face = 'plain'),
      axis.title.y = element_markdown(size = 9, face = 'plain'),
      axis.line = element_line(size = 0.25, color = "black"),
      axis.ticks = element_line(size = 0.25, color = "black"),
      legend.position = 'empty'
    ) +
    NULL
)

save_plot(here('Fig3/plots/Fig3.pdf'), treated, base_width = 7.08661, base_height = 4)
 