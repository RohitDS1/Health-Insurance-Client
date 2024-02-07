library(data.table)
library(DT)
library(lmtest)
library(dplyr)

adherence = fread("../Data/adherence.csv")

baseline = fread("../Data/baseline measurements.csv")

str(adherence)
str(baseline)

adherence[,id == "KQ68Ubs0yMY5ID5a",]

adherence %>%
  filter(id == "KQ68Ubs0yMY5ID5a") %>%
  select(id,t1,t2,ace,bb,statin)

unique(adherence$id)

head(adherence, 50)

unique(baseline$region)
unique(baseline$diabetes)
unique(baseline$gender)
unique(baseline$diabetes)

# Q1 - A patient's **length of follow-up** is the amount of time after diagnosis for which they 
#were under observation (with records in the data). 
# What was the median length of follow-up time?  
# What percentage of the patients had at least 1 year (365 days) of records?

# Calculate value for follow up time and group by patient ID
adherence[, max_t2 := max(t2), by = id]

# Calculate median of max_t2
median_follow_up_time <- median(adherence$max_t2)

# Display the result
print(round.numerics(median_follow_up_time,2))

# Patients with follow up time more than a year i.e. 365 days
patients_with_1_year <- adherence[max_t2 >= 365, unique(id)]

# Calculate the percentage
percentage_1_year <- length(patients_with_1_year) / length(unique(adherence$id)) * 100

# Display the result
print(round.numerics(percentage_1_year,2))

#Answer - Median Follow up time - 1335 and patients % over 365 days is 54.38%

##Q2 - For patients with at least 1 year (365 days) of follow-up, their **one-year adherence** 
# to a medication is the proportion of days in the first year after diagnosis during which the
# medication was possessed.  For each medication, what was the average one-year adherence of the
# patients?  Use only the patients with at least 1 year of follow-up records.

# Filter patients with at least 1 year of follow-up
patients_1_year <- adherence[max_t2 >= 365 & ace == 1]

# Calculate the proportion of days possessed for each medication
patients_1_year[, adherence_ace := sum(ace) / 365, by = id]
patients_1_year[, adherence_bb := sum(bb) / 365, by = id]
patients_1_year[, adherence_statin := sum(statin) / 365, by = id]

# Calculate the average one-year adherence for each medication
average_adherence_ace <- mean(patients_1_year$adherence_ace, na.rm = TRUE)
average_adherence_bb <- mean(patients_1_year$adherence_bb, na.rm = TRUE)
average_adherence_statin <- mean(patients_1_year$adherence_statin, na.rm = TRUE)

# Display the results
print(paste("Average One-Year Adherence for ACE Inhibitors:", round.numerics(average_adherence_ace, 3)))
print(paste("Average One-Year Adherence for Beta Blockers:", round.numerics(average_adherence_bb, 3)))
print(paste("Average One-Year Adherence for Statins:", round.numerics(average_adherence_statin, 3)))

#Answer - Average One-Year Adherence for ACE Inhibitors: 0.11
#         Average One-Year Adherence for Beta Blockers: 0.12
#         Average One-Year Adherence for Statins: 0.14

##Q3 - For ace inhibitors, beta blockers, and statins, we are interested in the number of these
# medications that patients take at the same time.  During the first year (365 days), 
# how frequently are the patients taking 0, 1, 2, or all 3 medications?  

# Filter patients with at least 1 year of follow-up
patients_1_year <- adherence[max_t2 >= 365]

# Calculate the number of medications taken for each patient
patients_1_year[, num_medications := ace + bb + statin]

# Count the frequency of each number of medications taken
medication_counts <- patients_1_year[, .N, by = num_medications]

# Calculate the percentage for each number of medications taken
medication_counts[, percentage := round.numerics(100 * N / sum(N), 1)]

# Format N with commas
medication_counts[, N := format(N, big.mark = ",", scientific = FALSE)]

# Sort by num_medications
medication_counts <- medication_counts[order(num_medications)]

# Display the results
print("Percentage of Patients Taking 0, 1, 2, or 3 Medications:")
print(medication_counts)

#Second Approach
# Q3 - Filter 1+ year records
q3_data <- adherence[get(t2.name) - get(t1.name) >= one.year.days]

# Q3 - Calculate concurrent meds
q3_data[, concurrent := ace * bb * statin]
str(q3_data)
# Q3 - Summarize info
q3_summary <- q3_data[, .(total_days = sum(get(t2.name) - get(t1.name)), 
                          concurrent_days = sum(concurrent, na.rm=TRUE))]
str(q3_summary)
# Q3 - Calculate percentage                       
q3_summary[, adherence_pct := round(100 * concurrent_days / total_days, 1)]

# Q3 - Print result
print(q3_summary)

##Q4 - What is the impact of diabetes, age, gender, region, and baseline condition on the
# one-year adherence to each medication? 
# Use only the patients with at least 1 year (365 days) of follow-up records. 
# Fit separate linear regression models for each medicine.  Then briefly comment on the results.

linear_adherence_model <- function(med.name) {
 
  # Calculate medicine adherence
  patients_1_year[, medicine_adherence := sum(get(med.name)) / (max(get(t2.name)) - min(get(t1.name))), by = .(id)]
  
  # Merge adherence and baseline datasets on id.name
  data_impact_adherence <- merge(baseline, patients_1_year[, .(id, medicine_adherence)], by = id.name)
  
  # Fit linear regression models for each medication
  lm_adherence_model <- lm(medicine_adherence ~ diabetes + age + gender + region + baseline.condition, data = data_impact_adherence)
  
  return(lm_adherence_model)
}

# Calling function to generate linear model for medication adherence
lm_statin_adherence_model <- linear_adherence_model(statin.name)
lm_ace_adherence_model <- linear_adherence_model(ace.name)
lm_bb_adherence_model <- linear_adherence_model(bb.name)

# Printing model summary and fetching Intercept "diabetes" coefficient
print(summary(lm_statin_adherence_model))
print(round.numerics(coef(lm_statin_adherence_model)["diabetes"],3))

print(summary(lm_ace_adherence_model))
print(round.numerics(coef(lm_ace_adherence_model)["diabetes"],3))

print(summary(lm_bb_adherence_model))
print(round.numerics(coef(lm_bb_adherence_model)["diabetes"],3))

cat("# ACE Inhibitors (ace) Beta Blockers (bb) and Statins (statin) Linear Regression model summary\n")
cat("1. Diabetes, age, gender, region, and baseline condition collectively have a statistically significant impact on adherence.\n")
cat("2. Patients with diabetes tend to have higher adherence to medications.\n")
cat("3. Older age is associated with decreased adherence.\n")

#Q5 - For each medicine, what percentage of the patients filled a prescription in the first two weeks 
# (before day t1 = 14) after their initial diagnoses? 
# Use data from all of the patients for this analysis, classifying each one as filling or not filling
# the prescription during this period.

# Approach 1 - Segregating patients by ace, bb and statin
# Calculate percentage of patients filling prescriptions in the first two weeks for ACE Inhibitors (ace)
ace_percentage <- sum(adherence$ace == 1 & adherence$t1 < 14) / nrow(adherence) * 100

# Calculate percentage of patients filling prescriptions in the first two weeks for Beta Blockers (bb)
bb_percentage <- sum(adherence$bb == 1 & adherence$t1 < 14) / nrow(adherence) * 100

# Calculate percentage of patients filling prescriptions in the first two weeks for Statins
statin_percentage <- sum(adherence$statin == 1 & adherence$t1 < 14) / nrow(adherence) * 100

# Print the percentages
cat("Percentage of patients filling prescriptions in the first two weeks for ACE Inhibitors (ace):", round.numerics(ace_percentage,2), "%\n")
cat("Percentage of patients filling prescriptions in the first two weeks for Beta Blockers (bb):", round.numerics(bb_percentage,2), "%\n")
cat("Percentage of patients filling prescriptions in the first two weeks for Statins:", round.numerics(statin_percentage,2), "%\n")

#Answer
# Percentage of patients filling prescriptions in the first two weeks for ACE Inhibitors (ace): 2.92 %
# Percentage of patients filling prescriptions in the first two weeks for Beta Blockers (bb): 3.19 %
# Percentage of patients filling prescriptions in the first two weeks for Statins: 3.7 %

# Filter to first 14 days  
early_fill_patients <- adherence[get(t1.name) < two.weeks.days]

# Get unique patients
unique_early_fill_patients <- unique(early_fill_patients$id)

# Identify bb fills
bb_fills <- early_fill_patients[get(bb.name) == 1] 

# Count unique bb patients
bb_patients <- unique(bb_fills$id)

# Calculate percentage
early_fill_percent <- round(100* length(bb_patients)/length(unique_early_fill_patients), 1)

# Print result  
print(paste0(early_fill_percent, "%"))

# Approach 2 - Keeping ace, bb, statin together

# Create a new column indicating whether any medication was filled in the first two weeks
adherence$any_medication_filled <- ifelse(adherence$t1 < 14 & (adherence$ace == 1 | adherence$bb == 1 | adherence$statin == 1), 1, 0)

# Calculate the overall percentage of patients filling prescriptions in the first two weeks
overall_percentage <- sum(adherence$any_medication_filled == 1) / nrow(adherence) * 100

# Print the overall percentage
cat("Overall percentage of patients filling prescriptions in the first two weeks:", round.numerics(overall_percentage,2), "%\n")

#Answer - Overall percentage of patients filling prescriptions in the first two weeks: 5.15 % / 5.14%

#Q6 - Now let's compare those who filled a prescription for a statin in the first two weeks 
# (before day t1 = 14) after diagnosis to those who did not.  
# Do these two groups have different baseline covariates?  Compare the groups based on their ages. 
# Then compare the distribution of baseline conditions in the two groups. 
# For continuous variables, compare their means using a t-test. 
# For the categorical variables, compare their distributions using a chi-squared test of independence.  

# Subset the data for patients who filled a prescription for a statin in the first two weeks
statin_filled <- merged_data[adherence$t1 < 14 & adherence$statin == 1, ]

# Subset the data for patients who did not fill a prescription for a statin in the first two weeks
statin_not_filled <- merged_data[adherence$t1 < 14 & adherence$statin == 0, ]

# Compare the groups based on age using a t-test
age_t_test <- t.test(statin_filled$age, statin_not_filled$age)

baseline_filled <- statin_filled$baseline.condition
baseline_not_filled <- statin_not_filled$baseline.condition

# Combine the factors to create a common level
combined_baseline <- factor(c(baseline_filled, baseline_not_filled))

# Create a contingency table
baseline_table <- table(combined_baseline)

# Perform the chi-squared test
baseline_chi_squared <- chisq.test(baseline_table)

# Print the results
cat("T-test for age:\n")
print(age_t_test)

cat("\nChi-squared test for baseline conditions:\n")
print(baseline_chi_squared)

#In summary, the t-test suggests that there is no significant difference in mean ages between the
# two groups, 
# while the chi-squared test indicates a significant association between baseline conditions.

#Q7 - How do the variables of age, gender, region, diabetes, and baseline condition impact
# the likelihood of initiating a medication within 14 days (before day t1 = 14)? 
# For each medicine, fit a logistic regression model and comment on the odds ratios. 
# Use data from all of the patients for this analysis.

glm_initiation_model <- function(data, med.name, two.weeks.days) {
  
  patients_two_weeks <- data[data[[t1.name]] < two.weeks.days]
  
  # Create ace initiated variable  
  patients_two_weeks[, initiated := get(med.name) == 1]
  
  logit_model <- glm(initiated ~ age + gender + region + diabetes + baseline.condition, data = patients_two_weeks, family = "binomial")
  
  return(logit_model)
}

glm_ace_model <- glm_initiation_model(merged_data, ace.name, two.weeks.days)
summary(glm_ace_model)

glm_bb_model <- glm_initiation_model(merged_data, bb.name, two.weeks.days)
summary(glm_bb_model)

glm_statin_model <- glm_initiation_model(merged_data, statin.name, two.weeks.days)
summary(glm_statin_model)

# Results
cat("Results:\n")

# Age
cat("Age: Older age tends to be associated with a decrease in adherence for all three medications.\n")

# Gender
cat("Gender: Being male is associated with a decrease in adherence.\n")

# Region
cat("Region: Patients in the Northeast, South, and West regions tend to have higher adherence compared to the Midwest.\n")

# Diabetes
cat("Diabetes: Presence of diabetes is associated with an increase in adherence.\n")

# Baseline Condition
cat("Baseline Condition: Having moderate symptoms or a light procedure is generally associated with a decrease in adherence.\n")

#2nd Approach
logit_ace <- glm(ace ~ age + gender + region + diabetes + baseline.condition,  
                 family = "binomial", data = merged_data)
logit_bb <- glm(bb ~ age + gender + region + diabetes + baseline.condition,   
                family = "binomial", data = merged_data)
logit_statin <- glm(statin ~ age + gender + region + diabetes + baseline.condition,    
                    family = "binomial", data = merged_data)


#### ACE Inhibitors
# Q7 - Print ACE model
print(summary(logit_ace))


#### Beta Blockers

# Q7 - Print BB model
print(summary(logit_bb))


#### Statins

# Q7 - Print statin model
print(summary(logit_statin))

#Q8 - For patients who did fill their prescriptions within 2 weeks (before day t1 = 14), 
# how long does it typically take to fill that first prescription after the initial diagnosis? 
# For each medicine, provide the mean, median, and standard deviation in units of days.

# Filter data for patients who filled prescriptions within 2 weeks for each medication
ace_filled_within_2_weeks <- adherence[adherence$ace == 1 & adherence$t1 < 14, ]
bb_filled_within_2_weeks <- adherence[adherence$bb == 1 & adherence$t1 < 14, ]
statin_filled_within_2_weeks <- adherence[adherence$statin == 1 & adherence$t1 < 14, ]
str(ace_filled_within_2_weeks)
# Calculate mean, median, and standard deviation for each medication
mean_ace <- mean(ace_filled_within_2_weeks$t2 - ace_filled_within_2_weeks$t1)
median_ace <- median(ace_filled_within_2_weeks$t2 - ace_filled_within_2_weeks$t1)
sd_ace <- sd(ace_filled_within_2_weeks$t2 - ace_filled_within_2_weeks$t1)

mean_bb <- mean(bb_filled_within_2_weeks$t2 - bb_filled_within_2_weeks$t1)
median_bb <- median(bb_filled_within_2_weeks$t2 - bb_filled_within_2_weeks$t1)
sd_bb <- sd(bb_filled_within_2_weeks$t2 - bb_filled_within_2_weeks$t1)

mean_statin <- mean(statin_filled_within_2_weeks$t2 - statin_filled_within_2_weeks$t1)
median_statin <- median(statin_filled_within_2_weeks$t2 - statin_filled_within_2_weeks$t1)
sd_statin <- sd(statin_filled_within_2_weeks$t2 - statin_filled_within_2_weeks$t1)

# Function to display results
display_results <- function(mean_vals, median_vals, sd_vals, medication) {
  result_table <- data.table(
    Medication = medication,
    Mean = mean_vals,
    Median = median_vals,
    SD = sd_vals
  )
  
  # Format the Mean, Median, and SD columns with two decimal places
  result_table[, `:=` (
    Mean = format(round.numerics(Mean, 2), nsmall = 2),
    Median = format(round.numerics(Median, 2), nsmall = 2),
    SD = format(round.numerics(SD, 2), nsmall = 2)
  )]
  
  print(result_table)
}

# Display results
display_results(mean_ace, median_ace, sd_ace, "ACE Inhibitors")
display_results(mean_bb, median_bb, sd_bb, "Beta Blockers")
display_results(mean_statin, median_statin, sd_statin, "Statins")

#Second Approach
df <- data.table(adherence)
# Create separate datasets for each medication - Removing Any patient who never filled a prescription
ace_data <- df[df[, .I[any(ace == 1)], by = id]$V1]
bb_data <- df[df[, .I[any(bb == 1)], by = id]$V1]
statin_data <- df[df[, .I[any(statin == 1)], by = id]$V1]

# Apply the function by group for "ace"
ace_data[, c("start_point_ace", "end_point_ace", "duration_ace") := calculate_start_end_duration_general(ace, t1, t2), by = id]

# Apply the function by group for "bb"
bb_data[, c("start_point_bb", "end_point_bb", "duration_bb") := calculate_start_end_duration_general(bb, t1, t2), by = id]

# Apply the function by group for "statin"
statin_data[, c("start_point_statin", "end_point_statin", "duration_statin") := calculate_start_end_duration_general(statin, t1, t2), by = id]

#Q9 - How does filling a prescription in the first two weeks impact adherence?  
# If we want to see that a medicine is working, we need to start the observation after the 
# patient has had a chance to fill the prescription. 
# To answer this question, we will follow a number of steps:

# 1.  Identify which patients filled a prescription in the first two weeks.  
#     You can call this variable **initiated** with binary values (1 for yes, 0 for no).

# 2.  Then, for each patient with at least 379 days of followup, measure the one-year adherence rate
#     (see Question 2) **starting at two weeks after the initial diagnosis**. 
#     This interval will begin at day 14 and last for 365 days.

# 3.  Fit a linear regression model of this one-year adherence including the baseline covariates
#.    (age, gender, region, diabetes, baseline condition) and an indicator of whether this patient
#.    filled a prescription for the medicine in the first two weeks.

# Perform this analysis for each medicine and comment on the results.

run_analysis <- function(medication, data) {
  # Identify patients who filled a prescription in the first two weeks
  data$initiated <- ifelse(data$t1 < 14 & data[[medication]] == 1, 1, 0)
  
  # Assuming 't2' represents the end of the follow-up period
  adherence_after_initiation <- data[data$t2 > (14 + 365) & data$initiated == 1, .(adherence = sum(get(medication))/(365 - 14)), by = id]
  
  # Merge the new adherence data with the original data
  data <- merge(data, adherence_after_initiation, by = "id", all.x = TRUE)
  
  # Fit a linear regression model
  formula <- as.formula(paste("adherence ~ age + gender + region + diabetes + baseline.condition + initiated"))
  model <- lm(formula, data = data)
  
  # Display the model summary
  print(paste("Results for", medication))
  print(summary(model))
}

# Run the analysis for ace
#run_analysis("ace", merged_data)

# Run the analysis for bb
run_analysis("bb", merged_data)

# Run the analysis for statin
run_analysis("statin", merged_data)


#Q10 - Once a patient starts a medication, how long do they continuously have a filled prescription? 
# For each patient who filled a medication, start with the first filled prescription and 
# count the duration of days until a gap occurs or follow-up ends. 
# (Hint: The first duration begins the first time ace = 1. The first duration ends at the first
# time ace = 0 after this point. Each patient will have 1 calculated duration if they filled a
# prescription and 0 otherwise.)  

#Then provide the mean, median, and standard deviation for these durations.  
# Do this separately for each medicine.

duration_filled_patients <- function(medication){
  
  filled_patients <- merged_data[get(medication) == 1, .(id, t1, t2, medication)]
  
  # Calculate the duration for each patient
  filled_patients[, duration := t2 - t1]
  
  # Aggregate durations for each medication
  medication_durations <- filled_patients[, .(mean_duration = round(mean(duration, na.rm = TRUE), 2),
                                              median_duration = round(median(duration, na.rm = TRUE), 2),
                                              sd_duration = round(sd(duration, na.rm = TRUE), 2)),
                                          by = .(medication)]
  
  # Print the result
  print(medication_durations)
}


duration_filled_patients <- function(med.name){
  
  filled_patients <- merged_data[get(med.name) == 1, .(get(id.name), get(t1.name), get(t2.name), get(med.name))]
  setnames(filled_patients, c(id.name, t1.name, t2.name, med.name))

  # Calculate the duration for each patient
  filled_patients[, duration := get(t2.name) - get(t1.name)]
 
  # Aggregate durations for each medication
  medication_durations <- filled_patients[, .(mean_duration = round(mean(duration, na.rm = TRUE), 2),
                                              median_duration = round(median(duration, na.rm = TRUE), 2),
                                              sd_duration = round(sd(duration, na.rm = TRUE), 2)),
                                          by = .(get(med.name))]
  
  # Print the result
  print(medication_durations)
}

duration_filled_patients(ace.name)
duration_filled_patients(bb.name)
duration_filled_patients(statin.name)














