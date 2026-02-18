library(ggplot2); theme_set(theme_light())
library(tidyr)
library(dplyr)
library(tidylog)
library(ggeffects)
library(sdmTMB)
library(modelr)
library(readr)

home <- here::here()

d <- read_csv(paste0(home, "/data/HS_SKAG_matrix_biomass_abundance.csv")) |> 
  janitor::clean_names() |> 
  rename(pw = sum_biomass_not_corrected_g) |> 
  mutate(family = as.factor(family),
         year_f = as.factor(year),
         year_ct = year - min(year),
         month_f = as.factor(month))

head(d)

ggplot(d, aes(lon, lat, color = as.factor(year))) + 
  geom_point() + 
  coord_sf()

d |> 
  summarise(n = n_distinct(sample_id), .by = year) |> 
  ggplot(aes(year, n)) + 
  geom_bar(stat = "identity")

summary(d$pw)

ggplot(d, aes(pw)) + 
  geom_histogram()

ggplot(d, aes(month)) + 
  geom_histogram()

# Expand the data
d <- d |> 
  mutate(id_fam = paste(sample_id, family, sep = "_"))

cross <- expand_grid(
  sample_id = unique(d$sample_id), 
  family = unique(d$family)
  ) |> 
  mutate(id_fam = paste(sample_id, family, sep = "_")) |> 
  filter(!id_fam %in% d$id_fam) |> 
  left_join(d |>
              dplyr::select(-family, -id_fam, -pw, -sum_abundance) |> 
              distinct(sample_id, .keep_all = TRUE),
            by = "sample_id")

cross |> filter(sample_id == "A2012/05278")

d2 <- bind_rows(d, cross) |> 
  mutate(pw = replace_na(pw, 0))

d2 |> 
  arrange(sample_id)

# Fit a simple model
m1 <- sdmTMB(
  pw ~ 0 + family,
  data = d2,
  family = tweedie(),
  spatial = "off"
)

summary(m1)
sanity(m1)
tidy(m1, exponentiate = TRUE) |> arrange(desc(estimate))
plot(ggeffect(m1, terms = "family")) + coord_flip()


# Next model, family random effects?
m2 <- sdmTMB(
  pw ~ year_f + (1|family),
  data = d2,
  family = tweedie(),
  spatial = "off"
)

summary(m2)
sanity(m2)

plot(ggeffect(m2, terms = "year_f"))
tidy(m2, effects = "ran_vals") |> arrange(estimate)


# Next model, family fixed effects?
m3 <- sdmTMB(
  pw ~ year_f + family,
  data = d2,
  family = tweedie(),
  spatial = "off"
)

summary(m3)
sanity(m3)

plot(ggeffect(m3, terms = "year_f"))
plot(ggeffect(m3, terms = "family")) + coord_flip()


# Next model, month random effects?
m4 <- sdmTMB(
  pw ~ year_f + family + (1|month_f),
  data = d2,
  family = tweedie(),
  spatial = "off"
)

summary(m4)
sanity(m4)

plot(ggeffect(m4, terms = "year_f"))
plot(ggeffect(m4, terms = "family")) + coord_flip()
tidy(m4, effects = "ran_vals")

AIC(m1, m2, m3, m4) |> 
  arrange(AIC)

# 
m5 <- sdmTMB(
  pw ~ family + (1|month_f) + (1|year_f),
  data = d2,
  family = tweedie(),
  spatial = "off"
)

summary(m5)
sanity(m5)

plot(ggeffect(m5, terms = "family")) + coord_flip()
tidy(m5, effects = "ran_vals")

AIC(m1, m2, m3, m4, m5) |> 
  arrange(AIC)

# Check fit and predictions
qqnorm(residuals(m4)); qqline(residuals(m4))

# Conditional predictions
nd <- d |> 
  data_grid(year_f, family, month_f) |> 
  mutate(year = as.numeric(as.character(year_f)))

p <- predict(m4, newdata = nd, se_fit = TRUE) |> 
  mutate(est = exp(est),
         upr = exp(est + est_se*1.96),
         lwr = exp(est - est_se*1.96))

p |> 
  filter(month_f == 1) |> 
  ggplot(aes(year, est, color = family)) + 
  #geom_line() + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = NA) +
  # geom_point(data = d |> 
  #              summarise(pw = mean(pw), .by = c(year, family)),
  #            aes(year, pw)) +
  facet_wrap(~family, scales = "free") +
  guides(color = "none")

# Spatial residuals in m3?
d$res <- residuals(m3)

ggplot(d, aes(lon, lat, color = res)) + 
  geom_point() + 
  scale_color_gradient2() +
  coord_sf()

# Spatial random fields? - Nope
# What if smooth year? - Nope
# What if random walk? - Nope

# conclusions... 
d2 |> 
  summarise(mean = mean(pw), .by = family) |> 
  arrange(desc(mean))

# Let's do 9 most common groups then other?
d |> 
  summarise(mean = mean(pw), .by = family) |> 
  arrange(mean)

# TODO: also why not random year fixed family?
