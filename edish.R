#############################################
# LIBRARIES
#############################################
library(dplyr)
library(tidyr)
library(plotly)

options(plotly.trace_attributes = NULL)
set.seed(2026)

#############################################
# 1) ADTTE - 200 patients
#############################################

n_total <- 200

adtte <- data.frame(
  STUDYID = "STUDYX001",
  USUBJID = paste0("STUDYX001-", 1001:(1000 + n_total)),
  AGE     = sample(40:85, n_total, replace = TRUE),
  SEX     = sample(c("M", "F"), n_total, replace = TRUE),
  ECOG    = sample(0:2, n_total, replace = TRUE),
  TRTA    = sample(c("Treatment A", "Treatment B"), n_total, replace = TRUE)
)

#############################################
# 1bis) Target AE Grade ???3 (nécessaire pour ADAE)
#############################################

set.seed(999)
adtte <- adtte %>%
  mutate(
    N_G3_target = sample(
      c(0, 1, 2),
      n(),
      replace = TRUE,
      prob = c(0.55, 0.30, 0.15)
    )
  )

#############################################
# 2) MEDDRA (avec type de toxicité)
#############################################

meddra_dict <- data.frame(
  AE_TOX = c(
    "Hepatic",
    "Hepatic",
    "Non-Hepatic",
    "Non-Hepatic",
    "Non-Hepatic"
  ),
  SOC = c(
    "Hepatobiliary disorders",
    "Hepatobiliary disorders",
    "Nervous system disorders",
    "General disorders",
    "Infections and infestations"
  ),
  PT = c(
    "Drug-induced liver injury",
    "Hyperbilirubinaemia",
    "Headache",
    "Fatigue",
    "Upper respiratory tract infection"
  ),
  stringsAsFactors = FALSE
)

#############################################
# 3) ADAE - Adverse Events
#############################################

set.seed(456)

adae <- adtte %>%
  mutate(
    N_AE = pmax(N_G3_target + rbinom(n(), 2, 0.4), N_G3_target)
  ) %>%
  filter(N_AE > 0) %>%
  uncount(N_AE, .id = "AESEQ")

# Ajout MedDRA
meddra_sample <- meddra_dict %>%
  slice_sample(n = nrow(adae), replace = TRUE)

adae <- bind_cols(adae, meddra_sample) %>%
  group_by(USUBJID) %>%
  rowwise() %>%
  mutate(
    AESEV = ifelse(
      row_number() <= first(N_G3_target),
      "SEVERE",
      sample(c("MILD", "MODERATE"), 1)
    ),
    AEGR3 = ifelse(AESEV == "SEVERE", 1, 0),
    AESTDY = sample.int(1500, 1)
  ) %>%
  ungroup()

#############################################
# 4) AE hépatiques cliniques Grade ???3
#############################################

hep_ae <- adae %>%
  filter(AE_TOX == "Hepatic", AEGR3 == 1) %>%
  group_by(USUBJID) %>%
  summarise(
    HEP_AE_G3_LIST = paste(
      unique(paste0(PT, " (G3+, Day ", AESTDY, ")")),
      collapse = "<br>"
    ),
    .groups = "drop"
  )

#############################################
# 5) ADLB - Laboratoires hépatiques
#############################################

adlb <- adtte %>%
  select(USUBJID, TRTA) %>%
  mutate(N_VISIT = sample(3:7, n(), replace = TRUE)) %>%
  uncount(N_VISIT, .id = "AVISITN") %>%
  mutate(
    ADY = AVISITN * 30,
    
    ULN_ALT   = sample(c(35, 40), n(), replace = TRUE),
    ULN_TBILI = sample(c(20, 21), n(), replace = TRUE),
    ULN_ALP   = sample(c(120, 130), n(), replace = TRUE),
    
    ALT   = rlnorm(n(), log(1.5 * ULN_ALT), 0.4),
    TBILI = rlnorm(n(), log(1.2 * ULN_TBILI), 0.3),
    ALP   = rlnorm(n(), log(1.1 * ULN_ALP), 0.25)
  )


#############################################
# 5bis) Injection de cas Hy's Law
#############################################

set.seed(777)
hys_ids <- sample(unique(adlb$USUBJID), 3)

adlb <- adlb %>%
  mutate(
    ALT = ifelse(
      USUBJID %in% hys_ids & AVISITN == max(AVISITN),
      runif(1, 5, 9) * ULN_ALT,
      ALT
    ),
    TBILI = ifelse(
      USUBJID %in% hys_ids & AVISITN == max(AVISITN),
      runif(1, 2.5, 4) * ULN_TBILI,
      TBILI
    )
  )

#############################################
# 6) eDISH dérivés (lab-driven)
#############################################

edish <- adlb %>%
  group_by(USUBJID) %>%
  summarise(
    # ALT
    MAX_ALT_ULN = max(ALT / ULN_ALT, na.rm = TRUE),
    DAY_MAX_ALT = ADY[which.max(ALT / ULN_ALT)][1],
    
    # TBILI
    MAX_TBILI_ULN = max(TBILI / ULN_TBILI, na.rm = TRUE),
    DAY_MAX_TBILI = ADY[which.max(TBILI / ULN_TBILI)][1],
    
    # ALP
    MAX_ALP_ULN = max(ALP / ULN_ALP, na.rm = TRUE),
    
    # R-ratio calculé au pic d'ALT
    R_RATIO = {
      idx <- which.max(ALT / ULN_ALT)[1]
      (ALT[idx] / ULN_ALT[idx]) / (ALP[idx] / ULN_ALP[idx])
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    # Pattern hépatique standard
    LIVER_PATTERN = case_when(
      R_RATIO >= 5 ~ "Hepatocellular",
      R_RATIO > 2 & R_RATIO < 5 ~ "Mixed",
      R_RATIO <= 2 ~ "Cholestatic",
      TRUE ~ "Unknown"
    ),
    
    # Screening Hy's Law eDISH
    HYS_LAW = ifelse(
      MAX_ALT_ULN >= 3 & MAX_TBILI_ULN >= 2 & R_RATIO >= 5,
      "Yes", "No"
    )
  )



#############################################
# 7) Merge FINAL pour eDISH
#############################################

edish <- edish %>%
  left_join(
    adtte %>% select(USUBJID, TRTA, AGE, ECOG),
    by = "USUBJID"
  ) %>%
  left_join(hep_ae, by = "USUBJID")

#############################################
# 8) TOOLTIP eDISH (DMC-grade)
#############################################

edish$tooltip <- paste0(
  "<b>USUBJID:</b> ", edish$USUBJID,
  "<br><b>Treatment:</b> ", edish$TRTA,
  "<br><b>Age:</b> ", edish$AGE,
  "<br><b>ECOG:</b> ", edish$ECOG,
  
  "<br><br><b>Max ALT:</b> ",
  round(edish$MAX_ALT_ULN, 2), " × ULN",
  " (Day ", edish$DAY_MAX_ALT, ")",
  
  "<br><b>Max TBILI:</b> ",
  round(edish$MAX_TBILI_ULN, 2), " × ULN",
  " (Day ", edish$DAY_MAX_TBILI, ")",
  
  "<br><b>Max ALP:</b> ",
  round(edish$MAX_ALP_ULN, 2), " × ULN",
  
  "<br><br><b>R-ratio:</b> ",
  round(edish$R_RATIO, 2),
  
  "<br><b>Liver pattern:</b> ",
  edish$LIVER_PATTERN,
  
  "<br><br><b>Hy's Law Screen:</b> ",
  edish$HYS_LAW,
  
  ifelse(
    !is.na(edish$HEP_AE_G3_LIST),
    paste0(
      "<br><br><b>Hepatic AE ??? G3:</b><br>",
      edish$HEP_AE_G3_LIST
    ),
    ""
  )
)


#############################################
# 9) eDISH PLOT 
#############################################

fig_edish <- plot_ly(
  data = edish,
  x = ~MAX_ALT_ULN,
  y = ~MAX_TBILI_ULN,
  type = "scatter",
  mode = "markers",
  color = ~TRTA,
  colors = c(
    "Treatment A" = "#d81b60",
    "Treatment B" = "darkblue"
  ),
  marker = list(
    size = 9,
    opacity = 0.8,
    symbol = "circle"
  ),
  text = ~tooltip,
  hovertemplate = "%{text}<extra></extra>"
) %>%
  layout(
    title = paste0(
      "eDISH - Hepatic Safety Overview",
      "<br><span style='font-size:12px;'>Population: ITT</span>"
    ),
    xaxis = list(
      title = "Maximum ALT / ULN (log scale)",
      type = "log"
    ),
    yaxis = list(
      title = "Maximum Total Bilirubin / ULN (log scale)",
      type = "log"
    ),
    shapes = list(
      # ALT = 3×ULN
      list(
        type = "line",
        x0 = 3, x1 = 3,
        y0 = 0.1, y1 = max(edish$MAX_TBILI_ULN) * 1.2,
        line = list(color = "red", dash = "dash", width = 1)
      ),
      # TBILI = 2×ULN
      list(
        type = "line",
        x0 = 0.1, x1 = max(edish$MAX_ALT_ULN) * 1.2,
        y0 = 2, y1 = 2,
        line = list(color = "red", dash = "dash", width = 1)
      )
    ),
    legend = list(
      title = list(text = "Treatment Arm")
    )
  )

fig_edish