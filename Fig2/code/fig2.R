library(tidyverse)
library(ggbeeswarm)
library(cowplot)
library(ggtext)
library(ZamanianLabThemes)
library(here)


data <- read_csv(here('fig2', 'data', 'Bpa_stacks_untreated.csv'), col_types = 'nccccccnnnccc')

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
  group_by(measure) %>% 
  group_nest() %>% 
  mutate(t = map(data, ~t.test(.x$value ~ .x$cue, alternative = 'greater')),
         tidy = map(t, broom::tidy)) %>% 
  unnest(tidy) %>% 
  mutate(sig = case_when(
    p.value <= 0.0001 ~ "****",
    p.value <= 0.001 ~ "***",
    p.value <= 0.01 ~ "**",
    p.value <= 0.05 ~ "*",
    p.value > 0.05 ~ ""
  )) %>% 
  mutate(value = 1.1,
         measure = case_when(
           measure == 'ci' ~ '(Cue - Control) / (Cue + Control)',
           measure == 'ci_total' ~ '(Cue - Control) / Total'
         ),
         cue = 'FBS')
  
(cue_control <- ci %>% 
    # only keep data where strains were done side-by-side
    filter(
      top_worms + bottom_worms > 0,
      cue %in% c('FBS', 'Water'),
      # measure == 'ci'
      ) %>%
    mutate(
      measure = case_when(
        measure == 'ci' ~ '(Cue - Control) / (Cue + Control)',
        measure == 'ci_total' ~ '(Cue - Control) / Total'
      )
    ) %>%
    ggplot(aes(x = cue, y = value)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_quasirandom(shape = 21, fill = alpha('black', 0.65), color = alpha('black', 0)) +
    geom_text(data = stat_layer,
              aes(label = sig)) +
    stat_summary(fun.data = mean_cl_normal, color = 'indianred') +
    labs(x = 'Cue', y = 'Chemotaxis index', color = 'Concentration', shape = 'Date') +
    scale_y_continuous(limits = c(-1.1, 1.1)) +
    scale_x_discrete(limits = c('Water', 'FBS')) +
    facet_grid(cols = vars(measure)) +
    theme_minimal() +
    theme(
      strip.text = element_markdown(size = 8),
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

save_plot(here('Fig2/plots/Fig2.pdf'), cue_control, base_width = 4.5, base_height = 2)
