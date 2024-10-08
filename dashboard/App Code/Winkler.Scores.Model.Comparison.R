#------------------------------------------------------------------------------#
#                                                                              #
#                 Calculating Winkler Scores - Model Comparison                #
#                                                                              #
#------------------------------------------------------------------------------#
# About:                                                                       #
#                                                                              #
# This function calculates the Winkler Scores for each of the forecasts        #
# produced in earlier steps of the dashboard and read into the dashboard in    #
# on the Model Comparison page. The Winkler scores provide a measure           #
# to compare prediction interval coverage for model fits and forecasts. The    #
# Winkler scores calculated in this function follow the equation given in:     # 
# https://otexts.com/fpp3/distaccuracy.html.                                   #
#------------------------------------------------------------------------------#
#                         By: Amanda Bleichrodt                                #
#------------------------------------------------------------------------------#
Winkler.Scores.Model.Comparison <- function(formatted.forecast.DASHBOARD,
                                            formatted.forecast.Other,
                                            calibrationPeriod.input,
                                            forecastHorizon.input,
                                            locations.input,
                                            date.type.input,
                                            avgWinler.input) {
  
#------------------------------------------------------------------------------#
# Creating the 'not-in' function -----------------------------------------------
#------------------------------------------------------------------------------#
# About: This section creates the 'not-in' function. Therefore, `%!in%` now    #
# can be used as the inverse of the built-in `%in%` function.                  #
#------------------------------------------------------------------------------#
  
  `%!in%` <- function(x, y) {
    
    !(x %in% y)
    
  }
  
#------------------------------------------------------------------------------#
# Reading in inputs ------------------------------------------------------------
#------------------------------------------------------------------------------#
# About: This section renames the function inputs to avoid any type of over-   #
# writing.                                                                     #
#------------------------------------------------------------------------------#
  
  #####################################################################
  # Reading in the formatted forecasts - ARIMA, GLM, GAM, and Prophet #
  #####################################################################
  formatted.forecast.input <<- formatted.forecast.DASHBOARD
  
  ##############################################
  # Reading in the formatted forecasts - Other #
  ##############################################
  formatted.forecast.other.input <<- formatted.forecast.Other 
  
  #####################################
  # Reading in the calibration period #
  #####################################
  calibrationPeriod <<- calibrationPeriod.input
  
  ##########################
  # Reading in the horizon #
  ##########################
  forecastingHorizon <<- forecastHorizon.input
  
  #####################################
  # Reading in the original locations #
  #####################################
  locationList <<- locations.input
  
  ############################
  # Reading in the date type #
  ############################
  dateType <<- date.type.input
  
  ##########################
  # Average Winkler Scores #
  ##########################
  avgWinkler <- avgWinler.input
  
  #########################################
  # Empty list to fill with renamed files #
  #########################################
  formattedForecastOtherRenamed <- list()
  
  ##########################################
  # Data frame to fill with winkler scores #
  ##########################################
  allWinklerScores <- data.frame()
  
  
#------------------------------------------------------------------------------#
# Error for running the figures without the dashboard results ------------------
#------------------------------------------------------------------------------#
# About: This section returns an error if a user trys to load other files      #
# prior to running the full dashboard.                                         #
#------------------------------------------------------------------------------#
  
  if(all(any(is.null(formatted.forecast.input) || length(formatted.forecast.input) == 0) & !is.null(formatted.forecast.other.input))){
    
    # Error to return
    return("ERROR1")
    
  }
  
#------------------------------------------------------------------------------#
# Potential errors with the loaded data ----------------------------------------
#------------------------------------------------------------------------------#
# About: This section checks for errors in the column names and file names of  #
# the loaded data.                                                             #
#------------------------------------------------------------------------------#
  
  for(i in 1:length(formatted.forecast.other.input)){
    
    # Indexed file
    data <- formatted.forecast.other.input[[i]]
    
    #############################
    # Checking the column names #
    #############################
    
    # Expected
    expectedNames <- c("Date", "data", "median", "LB", "UB")
    
    # Observed
    observedNames <- c(colnames(data))
    
    # Checking if they match each other
    if(any(expectedNames != observedNames)){
      
      # Returning an error
      return("ERROR2")
      
    }
    
    ##########################
    # Checking the file name #
    ##########################
    
    # Indexed file name
    dataName <- names(formatted.forecast.other.input)[i]
    
    # Checking for the word horizon #
    horizonModel <- qdapRegex::ex_between(dataName, "-", "-calibration")[[1]][1]
    
    # Horizon
    horizon <- qdapRegex::ex_between(horizonModel, "-", "-")[[1]][1]
    
    # Checking for the word calibration #
    calibration <- qdapRegex::ex_between(dataName, paste0(horizonModel, "-"), "-")[[1]][1]
    
    # Checking if they match what is expected
    if(any(horizon != "horizon" || calibration != "calibration")){
      
      # Returning an error
      return("ERROR3")
      
    }
    
    ###############################
    # Checking date specification #
    ###############################
    
    # Pulling the calibration period length
    caliLength <- qdapRegex::ex_between(dataName, paste0(calibration, "-"), "-")[[1]][1]
    
    # Pulling the location
    location <- qdapRegex::ex_between(dataName, paste0(calibration, "-", caliLength, "-"), "-")[[1]][1]
    
    # Pulling the date
    date <- qdapRegex::ex_between(dataName,  paste0(location, "-"), ".csv")[[1]][1]
    
    # Checking the date
    if(all(dateType == 'year' & nchar(date) != 4)){
      
      # Returning an Error
      return("ERROR4")
      
    }else if(all(dateType %in% c("week", "day") & nchar(date) != 10)){
      
      return("ERROR4")
      
    }
    
    ##########################
    # Checking the locations #
    ##########################
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  print(locationWord)
  print(locationsOrg)
    if(location %!in% c(locationList)){
      
      return("ERROR5")
      
    }
    
  }
  
#------------------------------------------------------------------------------#
# Winkler scores function ------------------------------------------------------
#------------------------------------------------------------------------------#
# About: This section creates the function to calculate Winkler scores for     #
# each formatted forecast executed in the main dashboard.                      #
#------------------------------------------------------------------------------#

  winkler_score <- function(upper.bound, lower.bound, data){
    
    ##################################################
    # Runs if the observed data is lower than the LB #
    ##################################################
    if(data < lower.bound){
      
      # Calculated Winkler Score
      score <- (upper.bound - lower.bound) + (2/0.05)*(lower.bound-data)
      
    #################################################
    # Runs if the observed data falls in the bounds #
    #################################################
    }else if(data >= lower.bound & data <= upper.bound){
      
      # Calculated Winkler score
      score <- upper.bound - lower.bound
      
    ################################################
    # Runs if the observed data falls above the UB #
    ################################################
    }else{
      
      # Calculated Winkler score 
      score <- (upper.bound - lower.bound) + (2/0.05)*(data - upper.bound)
      
    }
    
    #######################
    # Returning the score #
    #######################
    return(score)
    
  }
  

  
#------------------------------------------------------------------------------#
# Cleaning up the names of the other formatted forecasts -----------------------
#------------------------------------------------------------------------------#
# About: This section fixes the names of the other formatted forecasts to      #
# allow for easy merging and manipulation with the dashboard models.           #
#------------------------------------------------------------------------------#
  
  #######################################
  # Looping through formatted forecasts #
  #######################################
  for(i in 1:length(formatted.forecast.other.input)){
    
    # Name of forecast file
    forecastFileName <- names(formatted.forecast.other.input[i])
    
    #########################
    # Determining the model #
    #########################
    model <- qdapRegex::ex_between(forecastFileName, "", "-horizon")[[1]][1]
    
    ######################################
    # Determining the calibration period #
    ######################################
    calibration <-  qdapRegex::ex_between(forecastFileName, "calibration-", "-")[[1]][1]
    
    ############################
    # Determining the location #
    ############################
    location <- qdapRegex::ex_between(forecastFileName, paste0("calibration-", calibration, "-"), "-")[[1]][1]
    
    ####################################
    # Determining the forecast horizon #
    ####################################
    horizon <- qdapRegex::ex_between(forecastFileName, "horizon-", "-calibration")[[1]][1]
    
    ###################################
    # Determining the forecast period #
    ###################################
    
    # Forecast period 
    forecastDate <- qdapRegex::ex_between(forecastFileName, paste0(location, "-"), ".csv")[[1]][1]
    
    ################################################
    # Fixing the date format: Daily or Weekly data #
    ################################################
    if(nchar(forecastDate) > 4){
      
      # Changing the forecast date to a date YYYY-MM-DD format
      forecastDateFormat <- anytime::anydate(forecastDate)
      
      # Changing the date back to a character
      forecastDateFinal <- as.character(forecastDateFormat)
      
    #####################################
    # Setting the date: Yearly or Index #
    #####################################
    }else{ 
      
      forecastDateFinal <- as.character(forecastDate)
      
      } # End of 'if-else' for dates
    
    ######################################
    # Creating the new name for the file #
    ######################################
    fileNameNew <- paste0(model, "-", location, "-", forecastDateFinal)
    
    ##################################################
    # Adding calibration and horizon to the forecast #
    ##################################################
    forecast <- formatted.forecast.other.input[[i]] %>%
      dplyr::mutate(Calibration = as.numeric(calibration), # Calibration period 
                    Horizon = as.numeric(horizon)) # Horizon 
    
    ###############################
    # Renaming the forecast files #
    ###############################
    
    # Adding the file to the new list
    formattedForecastOtherRenamed[[i]] <- forecast
    
    # Renaming the file
    names(formattedForecastOtherRenamed)[i] <- fileNameNew
    
  }
  
#------------------------------------------------------------------------------#
# Cleaning the dashboard model forecasts ---------------------------------------
#------------------------------------------------------------------------------#
# About: This section adds the calibration and forecast horizon columns to the #
# formatted forecast files for the dashboard models. It then adds them back to #
# a list to be combined with the other dashboard models.                       #
#------------------------------------------------------------------------------#
  
  #######################################
  # Empty list to save the edited files #
  #######################################
  newDashboardForecasts <- list()
  
  ##################################
  # Looping through forecast files #
  ##################################
  for(i in 1:length(formatted.forecast.input)){
    
    ######################################
    # Pulling the forecast name and file #
    ######################################
    
    # Forecast name
    nameForecast <- names(formatted.forecast.input[i])
    
    # Forecast file
    data <- formatted.forecast.input[[i]]
    
    ####################################
    # Adding the necessary information #
    ####################################
    finaldata <- data %>%
      dplyr::mutate(Calibration = as.numeric(calibrationPeriod),
                    Horizon = as.numeric(forecastingHorizon))
    
    ####################################
    # Adding the data back to the list #
    ####################################
    
    # Adding the data 
    newDashboardForecasts[[i]] <- finaldata
    
    # Adding the name
    names(newDashboardForecasts)[i] <- nameForecast
  }
  
#------------------------------------------------------------------------------#
# Combining the other list with the dashboard list -----------------------------
#------------------------------------------------------------------------------#
# About: This section combines the newly-named other formatted forecasts and   #
# the forematted forecasts from the main dashboard into one list.              #
#------------------------------------------------------------------------------#
  
  #####################
  # Creating the list #
  #####################
  allForecasts <- c(formattedForecastOtherRenamed, newDashboardForecasts)
  
  
#------------------------------------------------------------------------------#
# Looping through formatted forecasts ------------------------------------------
#------------------------------------------------------------------------------#
# About: This section loops through the formatted forecasts to determine the   #
# forecast period (i.e., split calibration and forecasts), and to calculate    #
# the winkler scores, and saving the results for outputting.                   #
#------------------------------------------------------------------------------#
  for(w in 1:length(allForecasts)){
    
    # Indexed forecast
    indexForecast <- allForecasts[[w]]
    
    # Name indexed forecast
    nameForecast <- names(allForecasts[w])

    
    #########################################################
    # Pulling needed information from the list element name #
    #########################################################
    
    # Model name - Dashboard models
    if(grepl("ARIMA|GLM|GAM|SLR|Prophet", nameForecast)){
      
      model <- qdapRegex::ex_between(nameForecast, "", "-")[[1]][1]
      
    # Model name - Other models
    }else{
      
      model <- paste0(qdapRegex::ex_between(nameForecast, "", "-")[[1]][1], "-", qdapRegex::ex_between(nameForecast, "", "-")[[1]][3])
      
    }
    
    # Location name
    location <- qdapRegex::ex_between(nameForecast, paste0(model, "-"), "-")[[1]][1]
    
    # Adjusting for possible parenthesis in the name
    if(grepl("\\)", location) | grepl("\\(", location)) {
      
      # Adding \\ before the first instance of parenthesis 
      firstParenthesis <- gsub("\\(", "\\\\(", location)
      
      # Adding \\ before he last instance of parenthesis 
      locationForDate <- gsub("\\)", "\\\\)", firstParenthesis)
      
    }else{
      
      
      locationForDate <- location
      
    }
    
    # Forecast period
    forecastDate <- sub(paste0('.*-', locationForDate, '-'), '', nameForecast)
    
    #######################################################
    # Formatted the data frame prior to applying function #
    #######################################################
    
    formattedPreWinkler <- indexForecast %>%
      dplyr::mutate(Model = model, # Model type 
                    Location = location, # Location
                    `Forecast Date` = forecastDate, # Forecast date
                    `Winkler Score` = NA) # Empty column for later step 
    

    ########################################
    # Applying the winkler scores function #
    ########################################
    
    # Looping through rows of data set
    for(r in 1:nrow(formattedPreWinkler)){
      
      # Handling NAs in the data - Skipping that row 
      if(is.na(formattedPreWinkler[r,5]) | is.na(formattedPreWinkler[r,4]) | is.na(formattedPreWinkler[r,2])){
        
        next
        
      }
      
      # Indexed row 
      winklerScore <- winkler_score(formattedPreWinkler[r,5], formattedPreWinkler[r,4], formattedPreWinkler[r,2])
      
      # Adding the score to the data frame
      formattedPreWinkler[r,11] <- winklerScore
      
    } # End of row loop
    
    ####################################
    # Preparing the data for exporting #
    ####################################
    
    # If working with daily or weekly data
    if(nchar(formattedPreWinkler[1,1]) > 4){
      
        # Final data set to export 
        winklerData <- formattedPreWinkler %>%
          dplyr::mutate(`Forecast Date` = anytime::anydate(`Forecast Date`)) %>% # Formatted forecast date
          dplyr::mutate(CalibrationIndicator = ifelse(Date <= `Forecast Date`, 1, 0)) %>% # Indicator for calibration period
          dplyr::mutate(Date = anytime::anydate(Date)) %>% # Formatted date 
          dplyr::select(Location, Model, Date, `Forecast Date`, Calibration, Horizon, `Winkler Score`, CalibrationIndicator) %>% # Selecting needed variables 
          dplyr::group_by(CalibrationIndicator) %>% # Grouping by forecast period type 
          dplyr::mutate(`Winkler Score` = round(mean(`Winkler Score`), 2)) # Average Winkler score across forecast or calibration periods 
        
      # If working with yearly or time index data  
      }else{
        
        # Final data set to export 
        winklerData <- formattedPreWinkler %>%
          dplyr::mutate(`Forecast Date` = as.numeric(`Forecast Date`)) %>% # Formatted forecast date
          dplyr::mutate(CalibrationIndicator = ifelse(Date <= `Forecast Date`, 1, 0)) %>% # Indicator for calibration period 
          dplyr::mutate(Date = as.numeric(Date)) %>% # Formatted date 
          dplyr::select(Location, Model, Date, `Forecast Date`, Calibration, Horizon, `Winkler Score`, CalibrationIndicator) %>% # Selecting needed variables 
          dplyr::group_by(CalibrationIndicator) %>% # Grouping by forecast period type 
          dplyr::mutate(`Winkler Score` = round(mean(`Winkler Score`), 2)) # Average Winkler score across forecast or calibration periods
        
      }
    
   
    ####################################################
    # Adding the winkler score to the final data frame #
    ####################################################
    allWinklerScores <- rbind(allWinklerScores, winklerData)
    
  } # End of loop going through forecast files 
  
 
  #######################################
  # Preparing the data frame for export #
  #######################################
  
  # If non-average metrics 
  if(!avgWinkler){
    
    finalWinkler <- allWinklerScores %>%
      dplyr::ungroup() %>% # Removing any grouping from earlier 
      dplyr::mutate(`Performance Metric Type` = ifelse(CalibrationIndicator == 0, "Forecast", "Fit")) %>% # Creating the formatted variable for indicator
      dplyr::distinct(Location, Model, `Forecast Date`, CalibrationIndicator, .keep_all = T) %>% # Removing un-needed rows 
      dplyr::select(`Performance Metric Type`, Location, Model, `Forecast Date`, Calibration, Horizon,`Winkler Score`) %>% # Ordering variables 
      na.omit() # Removing NA rows 
  
  # If working with average metrics  
  }else{
    
    finalWinkler <- allWinklerScores %>%
      dplyr::ungroup() %>% # Removing any grouping from earlier 
      dplyr::mutate(`Performance Metric Type` = ifelse(CalibrationIndicator == 0, "Forecast", "Fit")) %>% # Creating the formatted variable for indicator
      dplyr::distinct(Location, Model, `Forecast Date`, CalibrationIndicator, .keep_all = T) %>% # Removing un-needed rows 
      dplyr::group_by(Location, Model, `Performance Metric Type`) %>% # Grouping variables
      dplyr::mutate(`Avg. Winkler Score` = round(mean(`Winkler Score`), 2)) %>%
      dplyr::select(`Performance Metric Type`, Location, Model, Calibration, Horizon, `Avg. Winkler Score`, `Forecast Date`) %>% # Ordering variables 
      na.omit() %>% # Removing NA rows 
      dplyr::ungroup() %>% # Ungroup
      dplyr::distinct(`Performance Metric Type`, Location, Model, Calibration, Horizon, .keep_all = T) # Removing repeat rows
    
  }
    
  ######################################
  # Returning the final Winkler Scores #
  ######################################
  return(finalWinkler)
  
}
    