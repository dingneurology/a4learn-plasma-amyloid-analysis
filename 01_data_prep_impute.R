# Load source data
demo <- read_csv("SUBJINFO.csv") %>%
  select(SUBSTUDY, BID, AGEYR, SEX, RACE, EDCCNTU, ETHNIC,
         APOEGN, APOEGNPRSNFLG, AAPOEGNPRSNFLG) %>%
  mutate(APOE4_positive = ifelse(grepl("E4", APOEGN), 1, 0))

plasma <- read_csv("biomarker_Plasma_Roche_Results.csv") %>%
  select(BID, LBTESTCD, LABRESN)

plasma_wide <- plasma %>%
  pivot_wider(names_from = LBTESTCD, values_from = LABRESN) %>%
  select(BID, GFAP, AMYLB40, AMYLB42, TPP181)

ptau217 <- read_csv("biomarker_pTau217.csv")  %>%
  select(BID, ORRESRAW)

amyloid_diagnose <- read_csv("imaging_PET_VA.csv")

# Merge analytic dataset
all_data <- demo %>%
  inner_join(amyloid_diagnose, by = c("BID", "SUBSTUDY")) %>%
  inner_join(plasma_wide, by = "BID") %>%
  inner_join(ptau217, by = "BID")

# Data for imputation
analysis_data <- all_data %>%
  filter(!is.na(pmod_suvr), !is.na(ORRESRAW))
imp_data <- analysis_data %>%
  select(
    SUBSTUDY, BID,
    AGEYR, SEX, RACE, EDCCNTU, ETHNIC,
    APOE4_positive,  
    AMYLB40, AMYLB42, TPP181,
  ) %>%
  mutate(
    SUBSTUDY = as.factor(SUBSTUDY),
    BID = as.factor(BID),
    SEX = as.factor(SEX),
    RACE = as.factor(RACE),
    ETHNIC = as.factor(ETHNIC),
    APOE4_positive = as.factor(APOE4_positive)
  )

ini <- mice(imp_data, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

meth["BID"] <- ""
meth["SUBSTUDY"] <- ""

meth["pmod_suvr"] <- ""
meth["ORRESRAW"] <- ""

pred[, c("BID", "SUBSTUDY")] <- 0
pred[c("BID", "SUBSTUDY"), ] <- 0

meth[c("AGEYR", "EDCCNTU", "AMYLB40", "AMYLB42", "TPP181")] <- "pmm"

set.seed(123)
imp <- mice(
  imp_data,
  m = 5,
  method = meth,
  predictorMatrix = pred,
  maxit = 20,
  printFlag = TRUE
)
complete_data <- complete(imp, 1)
write_csv(complete_data, "complete_data_no_impute_pet_ptau217.csv")
