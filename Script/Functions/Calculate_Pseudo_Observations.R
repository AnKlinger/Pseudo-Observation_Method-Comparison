# ========================================================================================================================================
# Calculate_Pseudo_Observations Function
# ========================================================================================================================================

source("Script/Functions/Generate_Data.R")

# Calculating pseudo observations for the restricted mean is done using the pseudo package 

# Example:

data <- Generate_Data()
pseudo::pseudomean(time = data$ObservedTime,
                   event = data$Status,
                   tmax = 10)
