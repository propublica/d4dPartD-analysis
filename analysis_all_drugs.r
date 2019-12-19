# --------------------------------------Documentation --------------------------------------------------
#
# Author:   Hannah Fresques
# Project:  Part D/ D4D project 
# Purpose:  Data analysis
# 
#_______________________________________________________________________________________________________


# setup -------------------------------------------------------------------

library(tidyverse)
library(here)
library(broom)
library(purrr)

# load data ---------------------------------------------------------------
DOC_DRUG <- here("data/dataStore/doctor_drug_2016.csv") %>% read_csv()



# make list of drugs ------------------------------------------------------

useDrugs <- DOC_DRUG %>% select(drug_name,top50_claims,top50_cost) %>% distinct() %>% arrange(drug_name)


# run models ---------------------------------------------------------------

# covars, used in all models
covars <- c("total_claim_count_other", "generalist")

# model 1:
model1 <- as.formula(paste("partD_claim_count ~ d4d_hasPayment + ", paste(covars, collapse= "+")))

model1DF <- DOC_DRUG %>% 
  # limit to people who prescribed the drug
  filter(is.na(partD_claim_count)==F & is.na(d4d_hasPayment)==F & partD_claim_count>=11)

lm1 <- model1DF %>% 
  group_by(drug_name) %>% 
  # Run separate regressions for each drug. 
  do(tidy(lm(model1, data=.))) %>% 
  ungroup()


# sensitivity: exclude ppl who got something other than meals
lm1_meals <- model1DF %>% 
  filter(
    d4d_payment_amount_nmeals==0
  ) %>% 
  group_by(drug_name) %>% 
  do(tidy(lm(model1, data=.))) %>% 
  ungroup()



# direction 2:
# Look at the other direction - do prescribers get more $ than non-prescribers?

model2 <- as.formula(paste("d4d_payment_amount ~ partD_hasClaim + ", paste(covars, collapse= "+")))

model2DF <- DOC_DRUG %>%
  # limit to ppl with a payment for the drug
  filter(d4d_hasPayment & is.na(d4d_payment_amount)==F & is.na(partD_hasClaim)==F)

lm2 <- model2DF %>% 
  group_by(drug_name) %>% 
  do(tidy(lm(model2, data=.))) %>% 
  ungroup()




# add a few more vars and some formatting to model results ----------------

# add variables to show how many providers fell into each category, and how many claims/how much money each group had

makeRawCountsLM1 <- function(df){
  tmp <- df %>% 
    group_by(drug_name,d4d_hasPayment) %>% 
    summarize(
      n_providers = n_distinct(npi),
      n_claims = mean(partD_claim_count)
    ) 
  rawCounts <- full_join( 
    tmp %>% 
      select(-n_providers) %>%
      spread(key=d4d_hasPayment,value=n_claims) %>% 
      rename(
        n_claims_c = `FALSE`,
        n_claims_p = `TRUE`
      ),  
    tmp %>% 
      select(-n_claims) %>% 
      spread(key=d4d_hasPayment,value=n_providers) %>% 
      rename(
        n_providers_c = `FALSE`,
        n_providers_p = `TRUE`
      ),
    by="drug_name"
  )
  return(rawCounts)
}

makeRawCountsLM2 <- function(df){
  tmp <- df %>% 
    group_by(drug_name,partD_hasClaim) %>% 
    summarize(
      n_providers = n_distinct(npi),
      payment_amount = mean(d4d_payment_amount)
    ) 
  rawCounts <- full_join(
    tmp %>% 
      select(-n_providers) %>%
      spread(key=partD_hasClaim,value=payment_amount) %>% 
      rename(
        payment_amount_c = `FALSE`,
        payment_amount_p = `TRUE`
      ),
    tmp %>% 
      select(-payment_amount) %>% 
      spread(key=partD_hasClaim,value=n_providers) %>% 
      rename(
        n_providers_c = `FALSE`,
        n_providers_p = `TRUE`
      ),
    by="drug_name"
  )
  return(rawCounts)
}

lm1_formatted <- lm1 %>% 
  left_join(
    model1DF %>% makeRawCountsLM1(),
    by="drug_name"
  )  

lm1_meals_formatted <- lm1_meals %>% 
  left_join(
    model1DF %>% filter(d4d_payment_amount_nmeals==0) %>% makeRawCountsLM1(),
    by="drug_name"
  )  

lm2_formatted <- lm2 %>% 
  left_join(
    model2DF %>% makeRawCountsLM2(),
    by="drug_name"
  )  


# functions to nicely format statistical output
pvalStars <- function(pvalue){
  return(
    case_when(
      is.na(pvalue) ~ "",
      pvalue<0.01 ~ "***",
      pvalue>=0.01 & pvalue<0.05 ~ "**",
      pvalue>=0.05 & pvalue<0.1 ~ "*",
      TRUE ~ ""
    )
  )
}

formatModelTables <- function(df){
  df %>% left_join(
    useDrugs,
    by="drug_name"
  ) %>% 
    mutate(
      stars = pvalStars(p.value)
    )
}

# Function to always round 0.5 up
round2 <- function(x, digits = 0) {  
  posneg <- sign(x)
  z <- abs(x) * 10^digits
  z <- z + 0.5
  z <- trunc(z)
  z <- z / 10^digits
  z * posneg
}

lm1_formatted       <- formatModelTables(lm1_formatted)
lm1_meals_formatted <- formatModelTables(lm1_meals_formatted)
lm2_formatted       <- formatModelTables(lm2_formatted)


# table5: lm1 -------------------------------------------------------------

lm1_table <- lm1_formatted %>% 
  filter(term=="d4d_hasPaymentTRUE") %>% 
  mutate(
    diff         = n_claims_p - n_claims_c,
    estimate_pct = estimate / n_claims_c,
    conf_int_95  = paste0("(",round2((estimate - std.error*1.96),0),", ",round2((estimate+std.error*1.96),0),")")
  ) %>%
  select(drug_name,n_providers_c,n_providers_p,n_claims_c,n_claims_p,diff,estimate,conf_int_95,estimate_pct,stars,top50_claims,top50_cost) %>% 
  ungroup()



# table6: lm2 -------------------------------------------------------------


lm2_table <- lm2_formatted %>% 
  filter(term=="partD_hasClaimTRUE") %>% 
  mutate(
    diff         = payment_amount_p - payment_amount_c,
    estimate_pct = estimate / payment_amount_c,
    conf_int_95  = paste0("($",round2((estimate - std.error*1.96),0),", $",round2((estimate+std.error*1.96),0),")") %>% str_replace_all("\\$-","- $")
  ) %>%
  select(drug_name,n_providers_c,n_providers_p,payment_amount_c,payment_amount_p,diff,estimate,conf_int_95,estimate_pct,stars,top50_claims,top50_cost) %>% 
  ungroup()



