#' ES_AVERAGE 
#' This function utilizes ensembles and single automatic ARIMAX models which have mean cases by AVERAGE states as exogenous variables.
#' The function fits rolling windows of N weeks for the state under analysis and rolling windows of the same size with 1 week-lag for the exogenous variables to generate forecasts.
#' It return some metrics that evaluate the performance of the models:
#' target_end_date, abs_error, cases, forecast, 'N_of_models", weighted interval score (WIS), predictive quantiles (%)  
#' The user can choose single best automatic ARIMAXs (auto=TRUE), or ensembles of 27 permutations of 0,1,2 pdq's (ES27=TRUE) or 64 permutations of 0,1,2,3 pdq's (ES64=TRUE).
#' The user also chooses the number of weeks ahead for each forecast, and the size of the rolling window.
#' 
#' @param current_state_tb A list containing Influenza Like Illness tibbles for one or many U.S. states.
#' @param auto A logical value indicating whether to use AUTO ARIMA. Default is \code{FALSE}.
#' @param ES27 A logical value indicating whether to use ensembles of 27 models. Default is \code{TRUE}.
#' @param ES64 A logical value indicating whether to use ensembles of 64 models. Default is \code{FALSE}.
#' @param n_weeks_ahead An integer specifying the number of weeks ahead for each forecast. Default is \code{1}.
#' @param week_lag An integer specifying the week lag between the exogenous variables and the ILI cases in the state under analysis. Default is \code{1}.
#' @param list_of_states A list containing Influenza Like Illness tibbles for ALL U.S. states.
#' @param window An integer specifying the size of the rolling window in number of weeks. Default is \code{104}.
#'
#' @return A list containing the forecast results and performance metrics.

ES_AVERAGE<-function(current_state_tb, auto=FALSE, n_weeks_ahead=1, week_lag=1, ES27=TRUE, ES64=FALSE, window=104, list_of_states=list_of_states){
  
  set.seed(1)
  # The model run the ES27 or the ES64.
  ES27=!ES64 
  # Empty list that will contain forecasts, predictive quantiles and number of models.
  results<-listenv()   
  if(ES27){
    pdq=c(0,1,2) # Possible ARIMA pdq's.
    my_order_params<-permutations(3,3,pdq, repeats.allowed = TRUE) # Create 27 permutations of [0,1,2]
  }
  if(ES64){
    pdq=c(0,1,2,3) # Possible ARIMA pdq's.
    my_order_params<-permutations(4,3,pdq, repeats.allowed = TRUE) # Create 64 permutations of [0,1,2,3]
  }
  # Apply the ROLLING_ARIMA function to get results. Window set to 104 weeks (2 years).
  results[[1]]<- ROLLING_ARIMAX_AVERAGE(current_state_tb, n_weeks_ahead=n_weeks_ahead, week_lag=week_lag, window = window, order_params=my_order_params, auto=auto, list_of_states=list_of_states) # %packages% "forecast" 

  # Put forecasts, prediction intervals and number of models into separate lists. 
  list_all_pred<-list() # Empty list for forecasts
  list_all_pred_quantiles<- list() # Empty list for predictive quantiles
  list_number_of_models<- list() # Empty list for number of models
  # Get forecasts and dates
  list_all_pred[[1]]<- results[[1]][[1]][[1]]
  # Get prediction intervals
  list_all_pred_quantiles[[1]]<-results[[1]][[2]][[1]]
  # Get the number of models in each ensemble by date
  list_number_of_models[[1]]<-results[[1]][[3]][[1]]
  
  #######################################################################
  # Format a tibble of predictive quantiles and point forecasts for WIS #
  predictive_quantiles_tb<- FormatForWIS(list_all_pred_quantiles=list_all_pred_quantiles, current_state_tb, model_name = "TestModel", n_weeks_ahead = n_weeks_ahead) # Format for weighted interval score calculation
  
  ########################
  # Format truth for WIS #
  truth<-current_state_tb
  truth["target_variable"]<-"cases" # rename cases as target_variable
  truth["model"]<-predictive_quantiles_tb[1,"model"]# the model name from predictive quantiles in the truth
  truth<- truth %>% rename_at("cases", ~'value') # rename the column named cases to values
  
  ########################################################################################
  # Calculate the WIS using our prediction intervals and truth values by target_end_date #
  my_forecast_scores<-score_forecasts(predictive_quantiles_tb, truth) 
  
  ################
  # Final results
  # Get the number of models in the ensembles #
  all_models<-data.frame(as.Date(unique(predictive_quantiles_tb$target_end_date)),list_number_of_models[[1]]) #Get the number of models for each forecasted date
  colnames(all_models)<-c("Dates","Number_of_models")
  # Format flu cases and dates #
  all_cases<-data.frame(current_state_tb$cases, as.Date(current_state_tb$target_end_date)) # Get flu cases from current state
  colnames(all_cases)<-c("cases","Dates")
  # Format forecasts and dates #
  all_forecasts<-data.frame(list_all_pred[[1]]$Prediction, as.Date(list_all_pred[[1]]$predicted_date)) # Get ensembled forecast from current state
  colnames(all_forecasts)<-c("forecasts","Dates")
  # Join flu cases and forecasts #
  forecasts_and_cases<-inner_join(all_forecasts,all_cases, by="Dates") # Join cases and forecasts
  colnames(forecasts_and_cases)<-c("forecasts","Dates","cases")
  # Get WIS and absolute error #
  WIS_errors<-data.frame(as.Date(c(my_forecast_scores[,"target_end_date"]$target_end_date)),c(my_forecast_scores[,"abs_error"]$abs_error),c(my_forecast_scores[,"wis"]$wis))
  colnames(WIS_errors)<-c("Dates","abs_error","WIS")  
  # Join WIS, absolute error and number of models by dates # 
  WIS_error_Nmodels<-inner_join(WIS_errors,all_models, by="Dates")
  colnames(WIS_error_Nmodels)<-c("Dates","abs_error","WIS","Number_of_models")
  # Join WIS, absolute error, number of models, forecasts and ILI cases by dates # 
  final_results<-inner_join(WIS_error_Nmodels,forecasts_and_cases, by="Dates")
  colnames(final_results)<-c("target_end_date","abs_error","WIS","Number_of_models","forecasts","cases")
  # Getting quantiles by target_end_date
  
  quantiles_by_date <- data.frame()
  for (i in 1:length(list_all_pred_quantiles[[1]])) {
    # Extract the date (name of the data frame) and the quantiles
    date <- names(list_all_pred_quantiles[[1]][i])
    quantiles <- t(list_all_pred_quantiles[[1]][[i]][2])

    # Create a temporary data frame with the quantiles and date
    temp_df <- data.frame(date = date, quantiles)
    # Bind the temporary data frame to the final data frame
    quantiles_by_date <- rbind(quantiles_by_date, temp_df)
  }
  colnames(quantiles_by_date)<-c("target_end_date","0.010", "0.025", "0.050", "0.100", "0.150", "0.200", "0.250", "0.300", "0.350", "0.400", "0.450", "0.500", "0.550",
                                 "0.600", "0.650", "0.700", "0.750", "0.800", "0.850", "0.900", "0.950", "0.975", "0.990")
  quantiles_by_date$target_end_date <- as.Date(quantiles_by_date$target_end_date)
  # Join WIS, absolute error, number of models, forecasts, ILI cases by dates and quantiles  
  results_and_quantiles<-inner_join(final_results,quantiles_by_date, by="target_end_date")
  return(results_and_quantiles)
}

#' ROLLING_ARIMAX_AVERAGE
#'
#' This function works together with the ES_AVERAGE function.
#' It fits the models, generate the forecasts and counts the number of models utilized in each ensemble.
#' 
#' @param current_state_tb A tibble containing the ILI data for the current U.S. state.
#' @param n_weeks_ahead An integer specifying the number of weeks ahead for each forecast. Default is \code{1}.
#' @param window An integer specifying the number of weeks to look back for the rolling window. Default is \code{104}.
#' @param order_params A list specifying the permutation of pdq's for the ensembles. Default is \code{NULL}.
#' @param week_lag An integer specifying the week lag between the exogenous variables and the ILI cases in the state under analysis. Default is \code{1}.
#' @param auto A logical value indicating whether to use an automatic best fitted ARIMA instead of an ensembles. Default is \code{FALSE}.
#' @param list_of_states A list containing Influenza Like Illness tibbles for ALL U.S. states.
#' 
#' @return A list containing the forecast, predictive quantiles and number of models in each ensemble.
#'

ROLLING_ARIMAX_AVERAGE <- function(current_state_tb, week_lag=week_lag, n_weeks_ahead=1, window = 104, order_params=NULL, auto=FALSE, list_of_states=list_of_states) {
  
  set.seed(1)
  ###############################
  # MEAN CASES BY AVERAGE STATE
  # A list with AVERAGE of U.S., without current state, based on our state_codes 
  
  all_my_states<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48)
  
  us_states_without_current <- list(
    Alabama = all_my_states[all_my_states != 1], 
    Arizona = all_my_states[all_my_states != 2], 
    Arkansas = all_my_states[all_my_states != 3], 
    California = all_my_states[all_my_states != 4], 
    Colorado = all_my_states[all_my_states != 5], 
    Connecticut = all_my_states[all_my_states != 6], 
    Delaware = all_my_states[all_my_states != 7], 
    Georgia = all_my_states[all_my_states != 8], 
    Idaho = all_my_states[all_my_states != 9], 
    Illinois = all_my_states[all_my_states != 10], 
    Indiana = all_my_states[all_my_states != 11], 
    Iowa = all_my_states[all_my_states != 12],
    Kansas = all_my_states[all_my_states != 13],
    Kentucky = all_my_states[all_my_states != 14],
    Louisiana = all_my_states[all_my_states != 15],
    Maine = all_my_states[all_my_states != 16],
    Maryland = all_my_states[all_my_states != 17],
    Massachusetts = all_my_states[all_my_states != 18],
    Michigan = all_my_states[all_my_states != 19],
    Minnesota = all_my_states[all_my_states != 20],
    Mississippi = all_my_states[all_my_states != 21],
    Missouri = all_my_states[all_my_states != 22],
    Montana = all_my_states[all_my_states != 23],
    Nebraska = all_my_states[all_my_states != 24],
    Nevada = all_my_states[all_my_states != 25],
    `New Hampshire` = all_my_states[all_my_states != 26],
    `New Jersey` = all_my_states[all_my_states != 27],
    `New Mexico` = all_my_states[all_my_states != 28],
    `New York` = all_my_states[all_my_states != 29],
    `North Carolina` = all_my_states[all_my_states != 30],
    `North Dakota` = all_my_states[all_my_states != 31],
    Ohio = all_my_states[all_my_states != 32],
    Oklahoma = all_my_states[all_my_states != 33],
    Oregon = all_my_states[all_my_states != 34],
    Pennsylvania = all_my_states[all_my_states != 35],
    `Rhode Island` = all_my_states[all_my_states != 36],
    `South Carolina` = all_my_states[all_my_states != 37],
    `South Dakota` = all_my_states[all_my_states != 38],
    Tennessee = all_my_states[all_my_states != 39],
    Texas = all_my_states[all_my_states != 40],
    Utah = all_my_states[all_my_states != 41],
    Vermont = all_my_states[all_my_states != 42],
    Virginia = all_my_states[all_my_states != 43],
    Washington = all_my_states[all_my_states != 44],
    `West Virginia` = all_my_states[all_my_states != 45],
    Wisconsin = all_my_states[all_my_states != 46],
    Wyoming = all_my_states[all_my_states != 47],
    Florida = all_my_states[all_my_states != 48]
  )
  
  ################################################################################### 
  # Function that calculates the mean cases in the U.S. without the current state  # 
  ##################################################################################
  
  sum_cases<-numeric(NROW(current_state_tb))
  # For the current state, compute the average using the list above

  for (given_state_ in us_states_without_current[[unique(current_state_tb$state_name)]]){ 
    cases<- list_of_states[[given_state_]]$cases
    sum_cases<-cases + sum_cases
  }
  
  mean_cases_AVERAGE_states<-sum_cases/length(us_states_without_current[[unique(current_state_tb$state_name)]]) 
###################################################################
  
  # SOME LISTS AND VARIABLES
  # All models in the ensemble
  N_of_models<-c() 
  # Final predictions list
  prediction<-list() 
  # Predictions and dates data frame 
  prediction_df<- data.frame("predicted_date"= NULL, "Prediction" = NULL) 
  # Predictive quantiles lists
  prediction_quantile<-list() 
  prediction_quantile_ls<- list() 
  
  #################################################################
  # Iterations over the dataset adapted to the number of week_lags
  for(iter in  (1+week_lag):(NROW(current_state_tb$cases)-(window))){ 
    
    # rolling window for current state ILI data  
    current_state_rolling_window<- (iter):(window+iter-1)  
    # rolling window for the exogenous data with N week_lags
    exog_rolling_window<- (iter-week_lag):(window+iter-1-week_lag)    
    
    # list that will get our ARIMA models
    fitted_models<-list()
    # list that will get the AIC scores 
    model_aic_scores<-c() 
    # Model id, utilized in the loop
    model_id<-1
    
    ##################################################
    # AVERAGE states time series for each iteration #
    # Selection of the calculated mean cases based on the exogenous window, which has 1 week lag
    AVERAGE_states_dataset= data.frame(mean_cases_AVERAGE_states[exog_rolling_window])
    #########################################
    # Exogenous variable for each iteration #
    # mean cases in the date the forecasting is being made
    exog_var<-c(mean_cases_AVERAGE_states[104+(iter-1)], mean_cases_AVERAGE_states[104+(iter-1)], mean_cases_AVERAGE_states[104+(iter-1)], mean_cases_AVERAGE_states[104+(iter-1)])

##########
# Fitting 
##########
    
    # Start if we have 104 elements in the rolling window.
    if(length(current_state_tb$cases[current_state_rolling_window])==window){ 
      
      # run 1 time for the auto.arima or run 27 or 64 times for the ensembles
      for(j in 1:nrow(order_params)){
        fit<- NULL # start with fit as NULL
        # try to fit an ARIMA model
        tryCatch(
          expr = {
            if(!auto){
              # if auto = FALSE, run the ensembles of 27 or 64
              set.seed(1)
              # fit ensembles of ARIMAs 
              fit<-Arima(current_state_tb$cases[current_state_rolling_window], xreg=AVERAGE_states_dataset[,1], order = order_params[j,], method = "CSS-ML") #
            }
            
            # if auto = TRUE, run auto.arima
            else{
              set.seed(1)
              fit<-invisible(auto.arima(current_state_tb$cases[current_state_rolling_window], xreg=AVERAGE_states_dataset[,1] ,stepwise=TRUE)) # trace my not be avaliable
            }
            
            # save each fitted ARIMA in fitted_models[[j]]
            fitted_models[[j]]<-fit
            # save the AIC of each fitted model
            model_aic_scores[model_id]<- fit$aic
          }
          # fit will be NULL if there is an error in the fitting process
          ,error = function(e){
          }
        )#end tryCatch
        
        # If fit == NULL, save the fitted model and the AIC as NAN
        if(is.null(fit) || is.null(fitted_models[[j]])){ 
          fitted_models[[j]]<-NA
          model_aic_scores[model_id]<- NA
        }
        # if auto==TRUE break the model on the first run
        # since we just need one result and not 27 or 64
        if(auto)
          break
        # model_ids are important for the ensembles
        # for the auto.arima it will be == 1
        model_id<-model_id+1 
      }
      
##############
# Forecasting 
##############
      
      # general initial variables
      predicted_value<- 0 # predicted values 
      m<- numeric(n_weeks_ahead) # mean forecast value
      s<- numeric(n_weeks_ahead) # standard deviation
      sims<-c() # simulations for the mixture of gaussians  
      
      # Ensemble weights initial variables
      model_weights<- c()
      min_aic<- min(model_aic_scores, na.rm = TRUE) # min models' aic
      total_aic<-sum(exp(-.5*(model_aic_scores-min_aic)), na.rm =TRUE ) # sum of aics without nan values
      
      # Counts the number of models utilized in each forecast  
      my_n_models<-0 
      for(my_model in fitted_models){ # for each valid model in fitted_models sum 1 to my_n_models
        if(length(my_model)>0 && !is.na(my_model[1]) && !(is.na(my_model$aic))){
          my_n_models<-my_n_models+1 
        }}
      
      ######################################################
      # Save the number of models utilized in each iteration  
      N_of_models<-append(N_of_models,my_n_models)
      
      ###############################
      # Generate the target_end_dates
      # Generates a sequence of dates based on the last date by n weeks ahead
      weekly_dates<- current_state_tb$target_end_date[current_state_rolling_window] # current data inside the 104 weeks window
      last_date <- max(weekly_dates) # last date of this window
      my_predicted_dates <- seq.Date(from = last_date + 7 , by = "week", length.out = n_weeks_ahead) 
      predicted_date<-my_predicted_dates[n_weeks_ahead]
      
      ##########################################
      # run the models stored on fitted models #
      for(my_model in fitted_models){
        ######################
        # if a model is valid
        if(length(my_model)>0 && !is.na(my_model[1]) && !(is.na(my_model$aic))){ 
          # calculate the weights based on the total AIC previously calculated
          model_weights_<- exp(-.5*(my_model$aic - min_aic))/total_aic
          
          set.seed(1)
          #######################################################################################
          # simulate values for each model based on its weights to build the mixture of Gaussians
          new.sims<-c() # list of simulations for each mixture of Gaussians                  
          fc <- forecast(my_model, h=n_weeks_ahead, xreg=exog_var[1:n_weeks_ahead], level=99)
          m <- fc$mean[n_weeks_ahead]  ## mean forecast value 
          s <- ((fc$upper[n_weeks_ahead]-fc$lower[n_weeks_ahead])/2.576/2)  # 99% confidence interval
          n<-ceiling(model_weights_*1e6) # number of simulations based on the weights
          new.sims <- rnorm(n, m=m, sd=s) # simulate values for each weighted model in as a gaussian
          sims <- c(sims, new.sims) # combine simulated values for each model
          
          #####################################
          # calculate the ensemble prediction #
          predicted_value <- model_weights_*m + predicted_value ### m = mean forecast for 99% confidence
          ###########################
          # if predicted value is NA 
          if(is.na(predicted_value)){
            print("predicted_value is na")
          }
        }
      }
      
      # get the prediction dataset
      single_prediction<- data.frame("predicted_date"= predicted_date, "Prediction" = predicted_value)# get the forecast for that date
      #rbind all forecasts and dates
      prediction_df<-rbind(prediction_df, single_prediction) 
      # Define the 23 quantiles
      probabilities <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99) 
      # get 23 predictive quantiles based on the mixture of gaussians 
      my_quantiles <- quantile(sims, probs=probabilities)
      
      ########################
      # Predictive quantiles #
      ########################
      # Creating predictive quantiles levels index
      pi_level<-c(0.01, 0.025, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
                  0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99)
      
      # reset the initial prediction_df_quantile for each iteration
      prediction_df_quantile<- data.frame("pi_level"= NULL, "quantile"= NULL, "point_forecast" = NULL)# empity dataframe of prediction intervals
      # save the 23 predictive quantiles, the upper and lower bounds for the given week ahead, and the ensemble forecasts into a data frame
      for(j in 1:23){ 
        single_df_quantile<- data.frame("pi_level"= pi_level[j],"quantile"= my_quantiles[j], "point_forecast" = predicted_value)# fill prediction intervals dataframe with correct values for each week ahead
        prediction_df_quantile<-rbind(prediction_df_quantile, single_df_quantile)
      }
      # put it into a list named with its forecast target_end_date
      prediction_quantile_ls[[toString(predicted_date)]]<-prediction_df_quantile 
    }
    # If don't have 104 weeks of data
    else
      print(paste0("Not enough values"))
  }
  # after everything is done
  print("Complete.")
  
  prediction[[1]]<-prediction_df # save the list of forecasts
  prediction_quantile[[1]]<- prediction_quantile_ls # save the list of predictive quantiles
  df_N_of_models<-data.frame(N_of_models) # save the number of models utilized in each forecast
  
  return(list("Point_ForeCast "=prediction, "Quantiles"=prediction_quantile, "Number_of_models"=df_N_of_models))
}

#' combining_states_data
#'
#' This function organizes the ILI dataset into a list of tibbles. Each tibble is a different U.S. state.
#' It also get the dates based on the year and epidemiological week.
#' 
#' @param ILI_data A ILI dataset from ILInet.
#' @param states_codes A dataset with codes and names which selects the U.S. states in a given order.
#' 
#' @return A list containing the forecast, predictive quantiles and number of models in each ensemble.
#'

combining_states_data<-function(ILI_data=NULL, state_codes=NULL){
  
  # select only the STATE, YEAR, EPI_WEEK, ILITOTAL columns from ILI data
  ILI_data = subset(ILI_data, select = c(STATE,YEAR,EPI_WEEK,ILITOTAL))
  # add a column with weekly dates
  ILI_data<-cbind(ILI_data, MMWRweek2Date(MMWRyear=ILI_data$YEAR, MMWRweek=ILI_data$EPI_WEEK))
  # select only location and location_name from state_codes
  state_codes = subset(state_codes, select = c(location,location_name))
  names(state_codes)<- c('STATE_NUMBER','STATE')
  
  # Joining datasets
  combined_data <- ILI_data %>%
    left_join(state_codes, by = "STATE")
  
  # Renaming, organizing variables types and removing NANs
  names(combined_data)<- c('state_name','MMWRyear','MMWRweek','cases','target_end_date','location')
  combined_data$location<-as.numeric(combined_data$location)
  combined_data$cases<-suppressWarnings(as.numeric(combined_data$cases))
  combined_data$target_end_date = as.Date(combined_data$target_end_date,format = "%Y/%m/%d")
  combined_data<-drop_na(combined_data)
  
  # split into different states
  states_data <-combined_data %>% group_split(location)
  return(states_data)
}

#' FormatForWIS
#'
#' This function formats the ROLLING_ARIMAX_AVERAGE results for calculating the weighted interval score.
#' 
#' @param list_all_pred_quantiles The list with predictive quantiles.
#' @param current_state_tb Tibble of the current U.S. under analysis.
#' @param model_name A name for the model you are utilizing, currently all models are named as 1. It does not afect the final results.
#' @param my_temporal_resolution Temporal resolution of the forecasts. Currently is set as "wk", which means weekly.
#' @param my_target_variable My target variable is named "cases".
#' @param n_weeks_ahead The number of weeks ahead the model is forecasting. This is important for calculating the data in which the forecast was made. I uses the same value as the ROLLING_ARIMAX_AVERAGE.
#' 
#' @return A list containing the forecast, predictive quantiles and number of models in each ensemble.
#'
FormatForWIS <- function(list_all_pred_quantiles, current_state_tb, model_name, n_weeks_ahead=1, my_temporal_resolution="wk", my_target_variable="cases") {
  my_tibble<- NULL
  
  # Create an empty tibble in the correct format
  my_tibble<-tibble(model=c(""),forecast_date=c(as.Date(c())), location=c(double()), horizon=c(double() ),
                    temporal_resolution=c(""), target_variable=c(""), target_end_date=c(as.Date(c())), type= c(""), quantile=c(double() ),
                    value =c(double()))
  
  # Loop over all the quantiles
  for(single_quantile in 1:NROW(list_all_pred_quantiles) ){
    # Get the dates of each single quantiles 
    predicted_dates_ls<- names(list_all_pred_quantiles[[single_quantile]])
    # Get the location number 
    my_location<-current_state_tb$location[1]
    
    # Get predicted dates from the list of predicted dates
    for(predicted_date_ in predicted_dates_ls){
      # Save predicted date as target_end_date
      my_target_end_date<-as.Date(predicted_date_)
      # Save the date in which the prediction was made as (predicted_date - (n_weeks_ahead*7))
      my_forecast_date<-as.Date(predicted_date_)-(7*n_weeks_ahead)
      
      # Add rows with predictions to the tibble as point_forecast
      my_tibble<- my_tibble%>%add_row(model=model_name,forecast_date=my_forecast_date, location=my_location, horizon=n_weeks_ahead,
                                      temporal_resolution=my_temporal_resolution, target_variable=my_target_variable, target_end_date=my_target_end_date, type= "point", quantile=NA,
                                      value = list_all_pred_quantiles[[1]][[predicted_date_]]$point_forecast[1]) # exponentiating the predictions back
      
      # Add rows with predictive_quantiles to the tibble as quantiles      
      for(quantile_level in list_all_pred_quantiles[[single_quantile]][predicted_date_]){
        my_quantile_value<-quantile_level$quantile # exponentiating the predictions back
        my_tibble<-my_tibble%>%add_row(model=model_name,forecast_date=my_forecast_date, location=my_location, horizon=n_weeks_ahead,
                                       temporal_resolution=my_temporal_resolution, target_variable=my_target_variable, target_end_date=my_target_end_date, type= "quantile",
                                       quantile=quantile_level$pi_level, value = my_quantile_value)
      }
    }
  }
  return(my_tibble)
}
