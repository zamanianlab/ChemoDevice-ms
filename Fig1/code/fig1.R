library(tidyverse)
library(ggbeeswarm)
library(cowplot)
library(ggtext)
library(ZamanianLabThemes)
library(magick)
library(here)

data <- read_csv(here('fig1', 'data', 'cel24stacks.csv'), col_types = 'nccccccnnnccc') %>% 
  drop_na(total_worms) %>% 
  mutate(wells = '24') %>% 
  bind_rows(read_csv(here('Fig1', 'data', 'cel96stacks.csv'), col_types = 'nccccccnnnccc') %>% 
              drop_na(total_worms) %>% 
              filter(# for some reason there are wells where the initial count is lower than the post-assay count
                top_worms + bottom_worms <= total_worms,
                cue %in% c('NaCl', 'Quinine')) %>% 
              mutate(wells = '96'))
 

# fig1b -------------------------------------------------------------------

gravity_data <- data %>% 
  filter(wells == '24',
         experiment_date < 20211201)

gravity_ci <- gravity_data %>% 
  mutate(
    ci_total = case_when(
      cue_side == 'top' ~ (top_worms - bottom_worms) / (total_worms),
      cue_side == 'bottom' ~ (bottom_worms - top_worms) / (total_worms)),
    ci = case_when(
      cue_side == 'top' ~ (top_worms - bottom_worms) / (top_worms + bottom_worms),
      cue_side == 'bottom' ~ (bottom_worms - top_worms) / (top_worms + bottom_worms)),
    genotype = 'N2',
    cue = case_when(
      cue == 'quinine' ~ 'Quinine',
      TRUE ~ cue
    )) %>%
  pivot_longer(c('ci', 'ci_total'), names_to = 'measure', values_to = 'value')

gravity_stat_layer <- gravity_ci %>% 
  group_by(cue_conc, measure) %>% 
  group_nest() %>% 
  mutate(t = map(data, ~t.test(.x$value ~ .x$cue_side, alternative = 'two.sided'))) %>% 
  mutate(tidy = map(t, broom::tidy)) %>% 
  unnest(tidy) %>%
  mutate(sig = case_when(
    p.value <= 0.0001 ~ "****",
    p.value <= 0.001 ~ "***",
    p.value <= 0.01 ~ "**",
    p.value <= 0.05 ~ "*",
    p.value > 0.05 ~ ""
  )) %>% 
  mutate(value = 1.1,
         cue_side = 'bottom',
         measure = case_when(
           measure == 'ci' ~ '(Cue - Control) / (Cue + Control)',
           measure == 'ci_total' ~ '(Cue - Control) / Total'
         ))

gravity_reps <- gravity_ci %>% 
  distinct(strain, genotype, treatment, cue_conc, cue, experiment_date) %>% 
  group_by(strain, genotype, treatment, cue_conc, cue) %>% 
  arrange(desc(experiment_date)) %>% 
  mutate(rep = row_number()) %>% 
  ungroup()

(fig1b <- gravity_ci %>% 
    mutate(measure = case_when(
      measure == 'ci' ~ '(Cue - Control) / (Cue + Control)',
      measure == 'ci_total' ~ '(Cue - Control) / Total'
    )) %>% 
    filter(measure == '(Cue - Control) / (Cue + Control)',
           top_worms + bottom_worms > 2) %>% 
    left_join(gravity_reps) %>% 
    ggplot(aes(x = cue_side, y = value)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_quasirandom(size = 1, shape = 21, fill = alpha('black', 0.5), color = alpha('black', 0)) +
    geom_text(data = gravity_stat_layer %>% filter(measure == '(Cue - Control) / (Cue + Control)'),
              aes(label = sig)) +
    stat_summary(fun.data = mean_cl_normal, color = 'indianred', size = 0.3) +
    labs(x = 'Side', y = 'Chemotaxis index', color = 'Replicate') +
    scale_y_continuous(limits = c(-1, 1.3)) +
    scale_x_discrete(limits = c('top', 'bottom'), labels = c('Top', 'Bottom')) +
    MetBrewer::scale_color_met_d('Juarez') +
    facet_grid(rows = vars(cue_conc), cols = vars(measure)) +
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


# fig1c -------------------------------------------------------------------

split_data <- read_csv(here('Fig1', 'data', 'cel24stacks_sliced.csv'),
                       col_types = 'ncccnc')

bubbles <- data %>% 
  filter(wells == '24') %>% 
  distinct(experiment_date, strain, cell, total_worms, .keep_all = TRUE) %>% 
  select(experiment_date, strain, cell, top_worms, bottom_worms, cue, cue_conc, cue_side, -total_worms) %>%
  filter(experiment_date %in% split_data$experiment_date) %>% 
  rename(top_bubble = top_worms, bottom_bubble = bottom_worms) %>% 
  pivot_longer(cols = contains('bubble'), names_to = 'slice', values_to = 'worms')

all_split <- bind_rows(split_data, bubbles) %>% 
  mutate(slice = factor(slice, 
                        levels = c('bottom_bubble', 'bottom', 'middle', 'top', 'top_bubble')))


all_sum <- all_split %>% 
  group_by(experiment_date, strain, cell, cue, cue_conc, cue_side) %>% 
  summarize(total = sum(worms))

all_join <- all_split %>% 
  left_join(all_sum) %>% 
  # filter(total > 4) %>% 
  mutate(percent_worms = worms / total) %>% 
  mutate(
    cue = case_when(
      cue == 'quinine' ~ 'Quinine',
      TRUE ~ cue),
    slice = case_when(
      cue == 'NaCl' & slice == 'top_bubble' ~ 'Attractant',
      cue == 'Quinine' & slice == 'top_bubble' ~ 'Repellant',
      slice == 'top' ~ 'Top',
      slice == 'middle' ~ 'Middle',
      slice == 'bottom' ~ 'Bottom',
      slice == 'bottom_bubble' ~ 'Control'),
    genotype = case_when(
      strain == 'N2' ~ 'N2',
      strain == 'PR678' ~ '*tax-4(p678)*')) %>% 
  mutate(slice = factor(slice, levels = rev(c('Attractant', 'Repellant', 'Top', 'Middle', 'Bottom', 'Control'))))

split_line <- all_join %>% 
  mutate(dummy = case_when(
           slice %in% c('Attractant', 'Repellant') ~ 1,
           slice == 'Top' ~ 2,
           slice == 'Middle' ~ 3,
           slice == 'Bottom' ~ 4,
           slice == 'Control' ~ 5
         )) %>% 
  group_by(cue, strain, genotype) %>% 
  group_nest() %>% 
  mutate(lm = map(data, ~lm(.x$percent_worms ~ .x$dummy))) %>% 
  mutate(tidy = map(lm, broom::tidy)) %>% 
  unnest(tidy) 

slope_p <- tibble(split_line %>% 
                    filter(term == '.x$dummy') %>% 
                    select(cue, strain, genotype, p.value) %>% 
                    rename(p.slope = p.value))

split_line <- split_line %>% 
  select(cue:genotype, term, estimate) %>% 
  pivot_wider(values_from = estimate, names_from = term) %>% 
  rename(intercept = 4, slope = 5) %>% 
  left_join(slope_p)

(fig1c <- all_join %>%
    mutate(
      genotype = factor(genotype, levels = c('*tax-4(p678)*', 'N2')),
      dummy = case_when(
        slice %in% c('Attractant', 'Repellant') ~ 1,
        slice == 'Top' ~ 2,
        slice == 'Middle' ~ 3,
        slice == 'Bottom' ~ 4,
        slice == 'Control' ~ 5
      )) %>% 
    ggplot(aes(x = dummy, y = percent_worms)) +
    geom_abline(data = split_line,
                aes(slope = slope, intercept = intercept),
                size = 1) +
    geom_quasirandom(aes(fill = interaction(dummy, cue)),
                     size = 1, alpha = 0.75, shape = 21, color = 'grey70') +
    geom_label(data = split_line,
               aes(label = str_glue('m = {round(slope, digits = 2)}\np = {formatC(p.slope, format = "e", digits = 2)}')),
               x = 3.7, y = 0.82, size = 2.5, hjust = 0, fill = alpha('white', 0.25)) +
    facet_grid(cols = vars(genotype), rows = vars(cue)) +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4, 5),
      labels = c('Cue', 'Top', 'Middle', 'Bottom', 'Control')) +
    scale_fill_manual(values = rev(c('#F7F0F0', '#BACCCE', '#7EA8AD', '#41848B', '#046069',
                                     '#F7F0F0', '#D4B5B5', '#B07A7A', '#8D3F3F', '#690404'))) +
    labs(x = 'Slice', y = 'Percent in slice', color = 'Strain', fill = 'Strain') +
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
    NULL)


# fig1d -------------------------------------------------------------------
 
top_data <- data %>% 
  filter(experiment_date >= 20211201,
         # wells == '24',
         # for some reason there are wells where the initial count is lower than the post-assay count
         top_worms + bottom_worms <= total_worms)

ci <- top_data %>% 
  filter(top_worms + bottom_worms != 0) %>%
  mutate(
    ci_total = (top_worms - bottom_worms) / (total_worms),
    ci = (top_worms - bottom_worms) / (top_worms + bottom_worms),
    genotype = case_when(
      strain == 'N2' ~ 'N2',
      strain == 'PR678' ~ '*tax-4<br>(p678)*',
      strain == 'PR672' ~ '*che-1<br>(p672)*',
      strain == 'SP1234' ~ '*dyf-2<br>(m160)*',
      strain == 'CX10' ~ '*osm-9<br>(ky10)*'),
    cue = case_when(
      cue == 'quinine' ~ 'Quinine',
      TRUE ~ cue
    )) %>% 
  pivot_longer(c('ci', 'ci_total'), names_to = 'measure', values_to = 'value')

# tax-4 KO response to quinine
t4 <- ci %>% 
  filter(strain == 'PR678', cue == 'Quinine') %>%
  group_nest(measure) %>% 
  mutate(t.test = map(data, ~t.test(.x$value)))

reps <- ci %>% 
  distinct(strain, genotype, treatment, conc, cue, experiment_date) %>% 
  group_by(strain, genotype, treatment, conc, cue) %>% 
  arrange(desc(experiment_date)) %>% 
  mutate(rep = row_number()) %>% 
  ungroup()

stat_layer <- ci %>% 
  mutate(strain = fct_relevel(strain, "N2")) %>%
  group_by(measure, cue) %>% 
  group_nest() %>% 
  mutate(otm = map(data, ~ PMCMRplus::kwManyOneDunnTest(.x$value ~ as.factor(.x$strain), p.adjust.method = "bonferroni"))) %>%
  mutate(summary = map(otm, ~ summary(.x))) %>%
  mutate(pval = map(summary, ~ as_tibble(pluck(summary, 1, 3), rownames = "strain"))) %>%
  mutate(estimate = map(summary, ~ as_tibble(pluck(summary, 1, 4), rownames = "strain"))) %>%
  mutate(tidy = map2(pval, estimate, ~ left_join(.x, .y, by = "strain"))) %>%
  unnest(tidy) %>%
  select(-data, -otm, -summary, -pval, -estimate, pval = `N2.x`, estimate = `N2.y`) %>%
  mutate(sig = case_when(
    pval <= 0.0001 ~ "****",
    pval <= 0.001 ~ "***",
    pval <= 0.01 ~ "**",
    pval <= 0.05 ~ "*",
    pval > 0.05 ~ ""
  )) %>% 
  left_join(., select(reps, strain, genotype)) %>% 
  mutate(value = 1.1,
         measure = case_when(
                    measure == 'ci' ~ '(Cue - Control) / (Cue + Control)',
                    measure == 'ci_total' ~ '(Cue - Control) / Total'
                  ))

(fig1d <- ci %>% 
    mutate(
      measure = case_when(
        measure == 'ci' ~ '(Cue - Control) / (Cue + Control)',
        measure == 'ci_total' ~ '(Cue - Control) / Total'
      )) %>% 
    filter(top_worms + bottom_worms > 0) %>%
    left_join(reps) %>% 
    ggplot(aes(x = genotype, y = value)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_quasirandom(aes(color = wells, shape = wells),
      size = 1) +
    geom_text(data = stat_layer,
               aes(label = sig)) +
    stat_summary(fun.data = mean_cl_normal, color = 'indianred', size = 0.3) +
    labs(x = 'Strain', y = 'Chemotaxis index', color = 'Format', shape = 'Format') +
    scale_y_continuous(limits = c(-1, 1.2)) +
    scale_x_discrete(limits = c('N2', '*tax-4<br>(p678)*', '*osm-9<br>(ky10)*', '*dyf-2<br>(m160)*', '*che-1<br>(p672)*')) +
    # MetBrewer::scale_color_met_d('Juarez') +
    scale_color_manual(values = c(alpha('black', 0.75), 'steelblue')) +
    # scale_shape_manual(values = c('circle', 'triangle')) +
    facet_grid(rows = vars(cue), cols = vars(measure)) +
    theme_minimal() +
    theme(
      strip.text = element_markdown(size = 8),
      axis.text.x = element_markdown(size = 8, angle = 0, hjust = 0.5),
      axis.text.y = element_markdown(size = 8),
      axis.title.y = element_markdown(size = 10, face = 'plain'),
      axis.line = element_line(size = 0.25, color = "black"), 
      axis.ticks = element_line(size = 0.25, color = "black"),
      legend.title = element_markdown(size = 10),
      legend.position = 'bottom'
    ) +
    NULL
)


# final -------------------------------------------------------------------

top <- image_read_pdf(here('Fig1', 'subplots', 'Fig1a.pdf'))
top <- ggdraw() + draw_image(top)


middle <- plot_grid(fig1b, fig1c, rel_widths = c(0.525, 1), 
                    align = 'h', axis = 'tb',
                    labels = c('B', 'C'), label_size = 12)

fig1 <- plot_grid(top, middle, fig1d, rel_heights = c(0.75, .8, 1),
                  labels = c('A', '', 'D'), label_size = 12,
                  nrow = 3)

save_plot(here('Fig1', 'plots', 'Fig1.pdf'), fig1,
          base_width = 7.08661, base_height = 8.26772)



# supplement --------------------------------------------------------------

(s2_fig <- ci %>% 
   mutate(
     measure = case_when(
       measure == 'ci' ~ '(Cue - Control) / (Cue + Control)',
       measure == 'ci_total' ~ '(Cue - Control) / Total'
     )) %>% 
   filter(top_worms + bottom_worms > 2) %>%
   left_join(reps) %>% 
   ggplot(aes(x = genotype, y = value)) +
   geom_hline(yintercept = 0, linetype = 'dashed') +
   geom_quasirandom(aes(color = wells, shape = wells),
                    size = 1) +
   geom_text(data = stat_layer,
             aes(label = sig)) +
   stat_summary(fun.data = mean_cl_normal, color = 'indianred', size = 0.3) +
   labs(x = 'Strain', y = 'Chemotaxis index', color = 'Format', shape = 'Format') +
   scale_y_continuous(limits = c(-1, 1.2)) +
   scale_x_discrete(limits = c('N2', '*tax-4<br>(p678)*', '*osm-9<br>(ky10)*', '*dyf-2<br>(m160)*', '*che-1<br>(p672)*')) +
   # MetBrewer::scale_color_met_d('Juarez') +
   scale_color_manual(values = c(alpha('black', 0.75), 'steelblue')) +
   # scale_shape_manual(values = c('circle', 'triangle')) +
   facet_grid(rows = vars(cue), cols = vars(measure)) +
   theme_minimal() +
   theme(
     strip.text = element_markdown(size = 8),
     axis.text.x = element_markdown(size = 8, angle = 0, hjust = 0.5),
     axis.text.y = element_markdown(size = 8),
     axis.title.y = element_markdown(size = 10, face = 'plain'),
     axis.line = element_line(size = 0.25, color = "black"), 
     axis.ticks = element_line(size = 0.25, color = "black"),
     legend.title = element_markdown(size = 10),
     legend.position = 'bottom'
   ) +
   NULL
)

save_plot(here('Fig1', 'plots', 'FigS2.pdf'), s2_fig,
          base_width = 5, base_height = 4)
