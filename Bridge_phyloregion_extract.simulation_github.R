no_cores <- detectCores() - 1

cl <- makeCluster(no_cores)


output.list <- parLapply(cl, rds_files, function(file) {
  
  require(dplyr)
  require(reshape2)
  
  # Use tryCatch to handle potential errors and allow continuation
  tryCatch({
    
    # Read the RDS file
    data <- readRDS(file)
    
    df <- data$output
    
    # Extract the param dataframe (one row, 8 columns)
    param_df <- as.data.frame(data$param)
    
    total_runs <- length(data$runs)
    
    successful_runs <- sum(sapply(data$runs, function(x) x$out == "Y"))
    
    # Initialize an empty list to store results for each run
    results <- list(div = list(), meta = list())
    
    
    #df$row_numbers <- row_numbers
    
    j <- 1
    # Loop over each run in the RDS file
    for (i in seq_along(data$runs)) {
      
      
      successful_array <- data$runs[[i]]$out
      
      
      if (successful_array == "Y") {
        
        # Extract the array for the current run
        array_1 <- data$runs[[i]]$all_species_array
        
        if (length(dim(array_1)) == 2) {
          
          # Perform the summing operation on array_1
          
          #sum.arr <- apply(array_1, 1:2, sum, na.rm = TRUE) %>%
          #  as.data.frame()
          sum.arr <- colSums(array_1, na.rm = TRUE) 
          
          # Store the combined result for this run
          results$div[[j]] <- sum.arr
          results$meta[[j]] <- c(sum(data$runs[[i]]$spec.id$alive),data$runs[[i]]$sim_num)
          
          j <- j + 1
          
        } else {
          # If array_1 does not have the expected dimensions, log a message or return NULL
          message("Array dimensions are incorrect for file: ", file, " run: ", i)
          results[[i]] <- NULL  # Skip this run
        }
      }
      # Ensure the array has the correct dimensions
      
    }
    
    # Return the stacked data frame for all runs in this file
    results_df <- (do.call(rbind, results$div))
    meta_df <- do.call(rbind,results$meta)
    runs <- c(successful_runs,total_runs)
    
    
    results_mean <- colMeans(results_df) %>% as_data_frame() %>% t()
    df <- cbind(df,results_mean)
    df$successful <- length(which(successful_array == "Y"))
    
    return(list(output = df, diversity_all = results_df, meta = meta_df, runs = runs))
  }, error = function(e) {
    # Handle errors by logging the error and returning NULL
    message("Error processing file: ", file, " - ", e$message)
    return(NULL)  # Continue processing other files despite the error
  })
})

stopCluster(cl)


# organize output list -------------------------------------------------------------------------



# Define the number of columns in the output
col_output <- ncol(output.list[[1]]$output)

# Use lapply to handle each element of the list
output_list2 <- lapply(output.list, function(x) {
  if (is.null(x$output)) {
    rep(NA, col_output)  # Create a vector of NAs for missing outputs
  } else {
    x$output
  }
})

# Combine the list of rows into a data frame
output_df <- do.call(rbind, output_list2)

# Convert to a data frame (if not already)
output_df <- as.data.frame(output_df)


runs <- do.call(rbind, lapply(output.list, function(x) {
  (x$runs)
}))

meta <- data_frame()

meta <- do.call(rbind, lapply(output.list, function(x) {
  c(mean(x$meta[, 1]), mean(x$meta[, 2]))
}))


# downstream analysis -------------------------------------------------------------------------

log_dif <- function(a,b){log((a+const)/(b+const))}
library(dplyr)
library(reshape2)
library(ggplot2)

# Load data bridge-------------------------------------------------------------------------

output_prepare_bridge <- function(output_df){
  
  output.df <- output_df
  
  colnames(output.df)[1:3] <- c("MPD1","MPDb","MPD2")
  colnames(output.df)[which(colnames(output.df) %in% c("1","2","3"))] <- c("div_trop","div_bridge","div_temp")
  
  #output.df <- output.df[!duplicated(output.df[,8:15]),]
  
  max_MPD1 <- max(output.df$MPD1[which(!is.infinite(output.df$MPD1))], na.rm = T)
  max_MPD2 <- max(output.df$MPD2[which(!is.infinite(output.df$MPD2))], na.rm = T)
  min_MPD1 <- min(output.df$MPD1[which(!is.infinite(output.df$MPD1))], na.rm = T)
  min_MPD2 <- min(output.df$MPD2[which(!is.infinite(output.df$MPD2))], na.rm = T)
  max_MPDb <- max(output.df$MPDb[which(!is.infinite(output.df$MPDb))], na.rm = T)
  min_MPDb <- min(output.df$MPDb[which(!is.infinite(output.df$MPDb))], na.rm = T)
  
  output.df$MPD1[which(is.infinite(output.df$MPD1) & output.df$MPD1 > 0)] <- max_MPD1
  output.df$MPD1[which(is.infinite(output.df$MPD1) & output.df$MPD1 < 0)] <- min_MPD1
  output.df$MPD2[which(is.infinite(output.df$MPD2) & output.df$MPD2 > 0)] <- max_MPD2
  output.df$MPD2[which(is.infinite(output.df$MPD2) & output.df$MPD2 < 0)] <- min_MPD2
  output.df$MPDb[which(is.infinite(output.df$MPDb) & output.df$MPDb > 0)] <- max_MPDb
  output.df$MPDb[which(is.infinite(output.df$MPDb) & output.df$MPDb < 0)] <- min_MPDb
  
  output.df$MPDdif[which(is.infinite(output.df$MPDdif))] <- output.df$MPD2[which(is.infinite(output.df$MPDdif))] - output.df$MPD1[which(is.infinite(output.df$MPDdif))]
  output.df$MPDdif[which(is.nan(output.df$MPDdif))] <- output.df$MPD2[which(is.nan(output.df$MPDdif))] - output.df$MPD1[which(is.nan(output.df$MPDdif))]
  
  
  #this line is important - incorporate the meta results for the second run
  output.df_filt <- output.df#[which(meta_results$alive_mean>450),]
  #output.df_filt <- output.df_filt[complete.cases(output.df_filt[,1:3]),]
  
  return(output.df_filt)
}
output_prepare_nb <- function(output_df){
  
  output.df <- output_df
  
  colnames(output.df)[which(colnames(output.df) %in% c("1","2"))] <- c("div_trop","div_temp")
  
  output.df <- output.df[!duplicated(output.df[,6:11]),]
  
  max_MPD1 <- max(output.df$MPD1[which(!is.infinite(output.df$MPD1))], na.rm = T)
  max_MPD2 <- max(output.df$MPD2[which(!is.infinite(output.df$MPD2))], na.rm = T)
  min_MPD1 <- min(output.df$MPD1[which(!is.infinite(output.df$MPD1))], na.rm = T)
  min_MPD2 <- min(output.df$MPD2[which(!is.infinite(output.df$MPD2))], na.rm = T)
  
  output.df$MPD1[which(is.infinite(output.df$MPD1) & output.df$MPD1 > 0)] <- max_MPD1
  output.df$MPD1[which(is.infinite(output.df$MPD1) & output.df$MPD1 < 0)] <- min_MPD1
  output.df$MPD2[which(is.infinite(output.df$MPD2) & output.df$MPD2 > 0)] <- max_MPD2
  output.df$MPD2[which(is.infinite(output.df$MPD2) & output.df$MPD2 < 0)] <- min_MPD2
  
  output.df$MPDdif[which(is.infinite(output.df$MPDdif))] <- output.df$MPD2[which(is.infinite(output.df$MPDdif))] - output.df$MPD1[which(is.infinite(output.df$MPDdif))]
  output.df$MPDdif[which(is.nan(output.df$MPDdif))] <- output.df$MPD2[which(is.nan(output.df$MPDdif))] - output.df$MPD1[which(is.nan(output.df$MPDdif))]
  
  
  #this line is important - incorporate the meta results for the second run
  output.df_filt <- output.df#[which(meta_results$alive_mean>450),]
  #output.df_filt <- output.df_filt[complete.cases(output.df_filt[,1:3]),]
  
  return(output.df_filt)
}

output.df <- output_prepare_bridge(output_df_bridge)

const <- 1e-6
output.df <- output.df %>%
  mutate(spec = log_dif(spec_trop,spec_temp),
         xtc = log_dif(xtc_trop,xtc_temp),
         disp = log_dif(disp_trop,disp_temp),
         xtc_trop_frac = xtc_trop/spec_trop,
         xtc_temp_frac = xtc_temp/spec_temp,
  )


output.df$n <- 1:nrow(output.df)
output.df_filt <- output.df[complete.cases(output.df[,1:3]),]
output.df_bridge <- output.df_filt
#filter list given the remaining rows of the output_df
output.list_bridge <- output.list_bridge[output.df_bridge_choose$n]


output.df <- output_prepare_nb(output_df)

output.df <- output.df[!duplicated(output.df[,6:11]),]
output.df <- output.df %>%
  mutate(spec = log_dif(spec_trop,spec_temp),
         xtc = log_dif(xtc_trop,xtc_temp),
         disp = log_dif(disp_trop,disp_temp),
         xtc_trop_frac = xtc_trop/spec_trop,
         xtc_temp_frac = xtc_temp/spec_temp
  )

output.df$n <- 1:nrow(output.df)
output.df_filt <- output.df[complete.cases(output.df[,1:2]),]
output.df_nobridge <- output.df_filt
output.list_nobridge <- output.list_nobridge[output.df_nobridge_choose$n]

# -------------------------------------------------------------------------

remove_param <- function(output, param, bridge) {
  param <- c(param,5e-3)
  if (bridge == TRUE) {
    output.df <- output[!(output$disp_trop %in% param | 
                            output$disp_temp %in% param | 
                            output$disp_btrop %in% param | 
                            output$disp_btemp %in% param), ]
  } else {
    output.df <- output[!(output$disp_trop %in% param | 
                            output$disp_temp %in% param), ]
  }
  
  return(output.df)
}


filt1 <- c(5e-2)
output.df_bridge_choose <- remove_param(output.df_bridge,filt1,T)
output.df_nobridge_choose <- remove_param(output.df_nobridge,filt1,F)