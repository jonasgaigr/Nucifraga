# === Parametry souborů ======================================================
infile <- "Data/Processed/data_nuc.csv"
out_dir <- "Outputs/Analyses_nuc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# === Načtení dat ============================================================
data_raw <- readr::read_csv2(infile, show_col_types = FALSE)

# --- ošetření názvů sloupců ---
names(data_raw) <- names(data_raw) %>%
  stringr::str_replace_all("\\.", "_") %>%
  stringr::str_replace_all("[\\(\\)]", "")

# Převod sloupců na numerické
data_model <- data_raw %>%
  dplyr::mutate(
    SMRK_m2 = as.numeric(SMRK_m2),
    Pocet_nalezu = as.numeric(Pocet_nalezu),
    prob_1h_mean = as.numeric(prob_1h_mean),
    Plocha_kurovcove_kalamity_celkem_m2 = as.numeric(Plocha_kurovcove_kalamity_celkem_m2),
    Plocha_souse_m2 = as.numeric(Plocha_souse_m2),
    Plocha_tezby_m2 = as.numeric(Plocha_tezby_m2),
    Celkem_zjisteno_oresniku = as.numeric(Celkem_zjisteno_oresniku),
    Datum_1 = readr::parse_date(Datum_1, format = "%Y-%m-%d"),
    Datum_2 = readr::parse_date(Datum_2, format = "%Y-%m-%d"),
    Provokovano_1 = as.numeric(Provokovano_1),
    Doba_2 = as.numeric(Doba_2)
  )

data_model$Datum_2_num <- as.numeric(data_model$Datum_2)  # počet dní od 1970-01-01
data_model$Datum_1_num <- as.numeric(data_model$Datum_1)

# Vypočet podílů relativně k plochám
data_model <- data_model %>%
  dplyr::mutate(
    podil_kalamity = ifelse(!is.na(SMRK_m2) & SMRK_m2 > 0,
                            Plocha_kurovcove_kalamity_celkem_m2 / SMRK_m2, NA_real_),
    podil_souse = ifelse(!is.na(SMRK_m2) & SMRK_m2 > 0,
                         Plocha_souse_m2 / SMRK_m2, NA_real_),
    podil_tezby = ifelse(!is.na(SMRK_m2) & SMRK_m2 > 0,
                         Plocha_tezby_m2 / SMRK_m2, NA_real_)
  )

# === Modely pro jednotlivé hypotézy =========================================

# Hypotéza 1: SMRK_m2 -> Celkem_zjisteno_oresniku
m_h1 <- MASS::glm.nb(
  Celkem_zjisteno_oresniku ~ SMRK_m2,
  data = data_model,
  na.action = na.exclude
)

# Hypotéza 2: Pocet_nalezu -> Celkem_zjisteno_oresniku
m_h2 <- MASS::glm.nb(
  Celkem_zjisteno_oresniku ~ Pocet_nalezu,
  data = data_model,
  na.action = na.exclude
)

# Hypotéza 3: prob_1h_mean -> Celkem_zjisteno_oresniku
m_h3 <- MASS::glm.nb(
  Celkem_zjisteno_oresniku ~ prob_1h_mean,
  data = data_model,
  na.action = na.exclude
)

# Hypotéza 4: podil_kalamity -> Celkem_zjisteno_oresniku
m_h4 <- MASS::glm.nb(
  Celkem_zjisteno_oresniku ~ podil_kalamity,
  data = data_model,
  na.action = na.exclude
)

# Hypotéza 5: podil_souse -> Celkem_zjisteno_oresniku
m_h5 <- MASS::glm.nb(
  Celkem_zjisteno_oresniku ~ podil_souse,
  data = data_model,
  na.action = na.exclude
)

# Hypotéza 6: podil_tezby -> Celkem_zjisteno_oresniku
m_h6 <- MASS::glm.nb(
  Celkem_zjisteno_oresniku ~ podil_tezby,
  data = data_model,
  na.action = na.exclude
)

# Hypotéza 7: Datum_1 -> Celkem_zjisteno_oresniku
m_h7 <- MASS::glm.nb(
  Celkem_zjisteno_oresniku ~ Datum_1,
  data = data_model,
  na.action = na.exclude
)

# === Export výsledků ========================================================

# Funkce pro export CSV, DOC a TXT
export_model <- function(model, model_name){
  
  # tidy CSV
  tab <- broom::tidy(model, conf.int = TRUE, conf.level = 0.95)
  utils::write.csv(tab,
                   file = file.path(out_dir, paste0(model_name, "_tidy.csv")),
                   row.names = FALSE)
  
  # Null model pro srovnání (jen intercept)
  m_null <- MASS::glm.nb(Celkem_zjisteno_oresniku ~ 1,
                         data = data_model, na.action = na.exclude)
  
  # DOC srovnání
  sjPlot::tab_model(
    model, m_null,
    show.aic = TRUE,
    show.ci = 0.95,
    dv.labels = c(model_name, "Null model"),
    file = file.path(out_dir, paste0(model_name, "_models.doc"))
  )
  
  # TXT report
  report_txt <- report::report(model)
  base::cat(report_txt, file = file.path(out_dir, paste0(model_name, "_report.txt")))
}

# Export všech modelů
export_model(m_h1, "H1_SMRC")
export_model(m_h2, "H2_NDOP")
export_model(m_h3, "H3_Prob")
export_model(m_h4, "H4_Kalamita")
export_model(m_h5, "H5_Souse")
export_model(m_h6, "H6_Tezba")
export_model(m_h7, "H7_Datum1")

message("Analýza dokončena, výstupy uložené v: ", normalizePath(out_dir))

pred_effect_plot <- function(model, data, predictor, predictor_label = predictor, n = 100){
  
  # vybereme řádky bez NA v predictor a v odpovědi
  data_sub <- data[!is.na(data[[predictor]]) & !is.na(data$Celkem_zjisteno_oresniku), ]
  
  # rozsah prediktoru
  x_seq <- seq(min(data_sub[[predictor]], na.rm = TRUE),
               max(data_sub[[predictor]], na.rm = TRUE),
               length.out = n)
  
  # vytvoříme datový rámec pro predikci
  newdata <- data_sub[1, , drop = FALSE]  # první řádek jen pro strukturu
  newdata <- newdata[rep(1, n), ]
  newdata[[predictor]] <- x_seq
  
  # ostatní numerické prediktory nastavíme na průměr
  num_vars <- names(data_sub)[sapply(data_sub, is.numeric)]
  num_vars <- setdiff(num_vars, "Celkem_zjisteno_oresniku")
  for (v in num_vars){
    if (v != predictor) newdata[[v]] <- mean(data_sub[[v]], na.rm = TRUE)
  }
  
  # predikce s SE
  preds <- predict(model, newdata = newdata, type = "response", se.fit = TRUE)
  newdata$fit <- preds$fit
  newdata$lower <- preds$fit - 1.96*preds$se.fit
  newdata$upper <- preds$fit + 1.96*preds$se.fit
  
  # ggplot
  p <- ggplot2::ggplot(newdata, ggplot2::aes_string(x = predictor, y = "fit")) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
    ggplot2::labs(
      x = predictor_label,
      y = "Predikovaný počet ořešníků",
      title = paste("Predikce podle", predictor_label)
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}

library(ggplot2)
library(patchwork)  # pro skládání grafů

# H1
p_h1 <- pred_effect_plot(m_h1, data_model, "SMRK_m2", "SMRK_m2")

# H2
p_h2 <- pred_effect_plot(m_h2, data_model, "Pocet_nalezu", "Pocet_nalezu")

# H3
p_h3 <- pred_effect_plot(m_h3, data_model, "prob_1h_mean", "prob_1h_mean")

# H4
p_h4 <- pred_effect_plot(m_h4, data_model, "podil_kalamity", "Podíl kůrovcové kalamity")

# H5
p_h5 <- pred_effect_plot(m_h5, data_model, "podil_souse", "Podíl souše")

# H6
p_h6 <- pred_effect_plot(m_h6, data_model, "podil_tezby", "Podíl těžby")

# H7
p_h7 <- pred_effect_plot(m_h7, data_model, "Datum_1", "Datum_1")

# --- Skladání a ukládání ---
ggplot2::ggsave("Outputs/Analyses_nuc/H1_H3.png", p_h1 + p_h2 + p_h3, width = 12, height = 4)
ggplot2::ggsave("Outputs/Analyses_nuc/H4_H6.png", p_h4 + p_h5 + p_h6, width = 12, height = 4)
ggplot2::ggsave("Outputs/Analyses_nuc/H7.png", p_h7 , width = 12, height = 4)

pred_effect_plot_points <- function(model, data, predictor, predictor_label = predictor, n = 100) {
  
  # Vybereme data, která mají hodnotu prediktoru a odpovědi
  data_sub <- data[!is.na(data[[predictor]]) & !is.na(data$Celkem_zjisteno_oresniku), ]
  
  # Rozsah prediktoru pro predikci
  x_seq <- seq(min(data_sub[[predictor]], na.rm = TRUE),
               max(data_sub[[predictor]], na.rm = TRUE),
               length.out = n)
  
  # Vytvoříme datový rámec pro predikci
  newdata <- data_sub[1, , drop = FALSE]  # jen první řádek pro strukturu
  newdata <- newdata[rep(1, n), ]
  newdata[[predictor]] <- x_seq
  
  # ostatní numerické prediktory nastavíme na průměr
  num_vars <- names(data_sub)[sapply(data_sub, is.numeric)]
  num_vars <- setdiff(num_vars, "Celkem_zjisteno_oresniku")
  for (v in num_vars) {
    if (v != predictor) newdata[[v]] <- mean(data_sub[[v]], na.rm = TRUE)
  }
  
  # Predikce s SE
  preds <- predict(model, newdata = newdata, type = "response", se.fit = TRUE)
  newdata$fit <- preds$fit
  newdata$lower <- preds$fit - 1.96*preds$se.fit
  newdata$upper <- preds$fit + 1.96*preds$se.fit
  
  # ggplot s body
  p <- ggplot2::ggplot() +
    # skutečná data
    ggplot2::geom_point(data = data_sub, ggplot2::aes_string(x = predictor, y = "Celkem_zjisteno_oresniku"),
                        alpha = 0.5, color = "black") +
    # predikce
    ggplot2::geom_line(data = newdata, ggplot2::aes_string(x = predictor, y = "fit"), color = "blue") +
    ggplot2::geom_ribbon(data = newdata, ggplot2::aes_string(x = predictor, ymin = "lower", ymax = "upper"),
                         alpha = 0.2, fill = "blue") +
    ggplot2::labs(
      x = predictor_label,
      y = "Počet zjištěných ořešníků",
      title = paste("Predikce a skutečná data podle", predictor_label)
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}


# H1–H3
p_h1 <- pred_effect_plot_points(m_h1, data_model, "SMRK_m2", "SMRK_m2")
p_h2 <- pred_effect_plot_points(m_h2, data_model, "Pocet_nalezu", "Pocet_nalezu")
p_h3 <- pred_effect_plot_points(m_h3, data_model, "prob_1h_mean", "prob_1h_mean")

# H4–H6
p_h4 <- pred_effect_plot_points(m_h4, data_model, "podil_kalamity", "Podíl kůrovcové kalamity")
p_h5 <- pred_effect_plot_points(m_h5, data_model, "podil_souse", "Podíl souše")
p_h6 <- pred_effect_plot_points(m_h6, data_model, "podil_tezby", "Podíl těžby")

# H7–H8
p_h7 <- pred_effect_plot_points(m_h7, data_model, "Datum_1", "Datum_1")

# --- Skladání a ukládání ---
ggplot2::ggsave("Outputs/Analyses_nuc/H1_H3_points.png", p_h1 + p_h2 + p_h3, width = 12, height = 4)
ggplot2::ggsave("Outputs/Analyses_nuc/H4_H6_points.png", p_h4 + p_h5 + p_h6, width = 12, height = 4)
ggplot2::ggsave("Outputs/Analyses_nuc/H7_points.png", p_h7, width = 12, height = 4)

