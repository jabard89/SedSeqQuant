#' Prepare Data with Replicates
#'
#' This function prepares count data by filtering transcripts with low counts, rounding counts to integers, and
#' dividing data into nested fractions. It then converts the data to a count matrix of transcript_IDs and Reps.
#'
#' @param count_data A data frame containing count data (output of load_counts).
#' @param min_counts Minimum count threshold for filtering low counts (default = 20).
#'
#' @return A tibble containing a column of prepared data in the form of a named list:
#' \itemize{
#'   \item{"NREP"}{Number of replicates.}
#'   \item{"NRNA"}{Number of RNA sequences.}
#'   \item{"tot_obs_counts"}{Total observed counts matrix.}
#'   \item{"sup_obs_counts"}{Supernatant observed counts matrix.}
#'   \item{"pel_obs_counts"}{Pellet observed counts matrix.}
#' }
#' @export
#' @importFrom dplyr select filter group_by arrange mutate
#' @importFrom tidyr pivot_wider pivot_longer nest
#' @importFrom purrr map
#' @importFrom tibble rownames_to_column
#' @importFrom stats setNames
#' @importFrom magrittr %>%
prepare_data_with_reps <- function(count_data,min_counts=20){
  filter_low_counts <- function(count_data, min_counts = 20) {
    wide_data <- count_data %>%
      select(-Sample_ID) %>%
      pivot_wider(names_from="Fraction",values_from="Count")
    selected_columns <- wide_data[, c('Supernatant', 'Pellet')]

    # Calculate the row sums
    row_sums <- rowSums(selected_columns, na.rm = TRUE)

    # Filter out rows where the sum is less than min_sum
    # Also filter out transcripts that don't appear in all bioreps
    count_data_filtered <- wide_data[row_sums >= min_counts, ] %>%
      dplyr::filter(Total>min_counts) %>%
      group_by(Condition) %>%
      mutate(NREPS = length(unique(Rep))) %>%
      group_by(Condition,transcript_ID) %>%
      dplyr::filter(length(transcript_ID)==NREPS) %>%
      select(-NREPS) %>%
      pivot_longer(cols=c("Total","Supernatant","Pellet"),names_to="Fraction",values_to="Count")

    return(count_data_filtered)
  }

  nested_data <- count_data %>%
    filter_low_counts(.,min_counts) %>%
    mutate(Count = round(Count)) %>% # round counts to integers
    arrange(Condition,Fraction,transcript_ID,Rep) %>%
    group_by(Condition) %>% nest %>%
    # divide up data into a list of fractions
    mutate(prepared_data=map(data,function(df){
      nested_data <- df %>% group_by(Fraction) %>% split(.$Fraction) %>%
        # convert data to a count matrix of transcript_IDs and Reps
        map(~.x%>%pivot_wider(names_from=Rep,values_from=Count) %>%
              column_to_rownames(var="transcript_ID") %>%
              select(-Fraction) %>% as.matrix(.))

      # Compute mean and sd of log-transformed tot_obs_counts
      mean_log_tot_obs <- apply(log1p(nested_data$Total), 1, mean)
      sd_log_tot_obs <- apply(log1p(nested_data$Total), 1, sd)

      # get guesses for mixing ratios
      mixing_guess <- df %>%
        group_by(transcript_ID) %>%
        mutate(Total.mean = exp(mean(log1p(Count[Fraction=="Total"])))) %>%
        ungroup %>%
        mutate(mixing_factor_guess = if_else(Fraction=="Total",Count/Total.mean,
                                             Count/(Total.mean*0.5))) %>% # guess that the pSup is 0.5
        group_by(Fraction) %>%
        summarise(mean = mean(mixing_factor_guess),
                  sd = sd(mixing_factor_guess))

      prepared_data <- list(
        "NREP"=dim(nested_data$Total)[2],
        "NRNA"=dim(nested_data$Total)[1],
        "mean_log_tot_obs" = mean_log_tot_obs,
        "sd_log_tot_obs" = sd_log_tot_obs,
        "mixing_factor_total_guess_mean" = mixing_guess %>% filter(Fraction=="Total") %>% pull(mean),
        "mixing_factor_total_guess_sd" = mixing_guess %>% filter(Fraction=="Total") %>% pull(sd),
        "mixing_factor_sup_guess_mean" = mixing_guess %>% filter(Fraction=="Supernatant") %>% pull(mean),
        "mixing_factor_sup_guess_sd" = mixing_guess %>% filter(Fraction=="Supernatant") %>% pull(sd),
        "mixing_factor_pellet_guess_mean" = mixing_guess %>% filter(Fraction=="Pellet") %>% pull(mean),
        "mixing_factor_pellet_guess_sd" = mixing_guess %>% filter(Fraction=="Pellet") %>% pull(sd),
        tot_obs_counts=nested_data$Total,
        sup_obs_counts=nested_data$Supernatant,
        pel_obs_counts=nested_data$Pellet
      )
      return(prepared_data)
    }))
  return(nested_data %>% select(Condition,prepared_data))
}

#' Model fit modified
#'
#' Fit the bayesian statistical model to the counts data.
#' @details The input data should be from prepare_data_with_reps
#' @param nested_data output of prepare_data_with_reps
#' @param chains A number, defaulting to 4
#' @param iter A number, defaulting to 1000
#' @param ... Other arguments passed to `rstan::sampling` (e.g. warmup, seed).
#' @return A dataframe with a column containing the class `stanfit` returned by `rstan::sampling`
#' @seealso [rstan::sampling()], `browseVignettes("rstan")`
#' @keywords model fit
#' @importFrom rstan stan_model sampling summary
#' @export
#'
model_fit_reps <- function(nested_data, chains = 4,
                           iter = 1000,...)
{
  nested_stanfit <- nested_data %>%
    mutate(stanfit = map(prepared_data,function(x){
      rstan::sampling(stanmodels$SedSeqQuantReps, data = x,
                      chains = chains, iter = iter,
                      ...)
    })) %>%
    select(Condition,stanfit)
  return(nested_stanfit)
}

#' parameters statistical results
#'
#' Get the estimated statistical results of parameters from a stan model
#' Extracts the names of transcript_IDs from the input data
#' @details The input data is the output of model_fit_reps
#' Output a data frame that includes a statistical result of parameters.
#' @param nested_stanfit output of model_fit_reps
#' @param nested_data input data into model_fit_reps
#' @return A data frame
#' @keywords parameter extraction
#' @import dplyr
#' @export
get_stan_summary_reps <- function(nested_stanfit,nested_data)
{
  temp_data <- full_join(nested_stanfit,nested_data,by="Condition") %>%
    mutate(stan_summary = map2(stanfit,prepared_data,function(sf,pd){
      item <- rstan::summary(sf)
      transcript_tibble <- tibble(transcript_ID=rownames(pd$tot_obs_counts),
                                  index=seq(1,pd$NRNA))
      params_summary <- tibble(Variable=rownames(item$summary),
                               mean=item$summary[,"mean"],
                               se_mean=item$summary[,"se_mean"],
                               sd=item$summary[,"sd"],
                               "x2.5" = item$summary[,"2.5%"],
                               x25 = item$summary[,"25%"],
                               x50 = item$summary[,"50%"],
                               x75 = item$summary[,"75%"],
                               "x97.5" = item$summary[,"97.5%"],
                               n_eff = item$summary[,"n_eff"],
                               Rhat = item$summary[,"Rhat"])
      # extract the scaling parameter for total (relative to Rep=1)
      long_params_mixing <- params_summary %>%
        filter(grepl("mixing_factor_",Variable)) %>%
        pivot_longer(cols=c(mean:Rhat),names_to="Term",values_to="Value")%>%
        mutate(
          Fraction = str_extract(Variable, "^[^\\[]+"),
          RepIndex = as.integer(str_extract(Variable, "\\d+")) # Get the index
        ) %>%
        rowwise() %>%
        mutate(Rep = colnames(pd$tot_obs_counts)[RepIndex]) %>% # Convert index to Rep character representation
        select(-RepIndex, -Variable)
      long_params_counts <- params_summary %>%
        filter(grepl("latent_counts",Variable)) %>%
        mutate(
          Fraction = str_extract(Variable, "^[^\\[]+"),
          index = as.integer(str_extract(Variable, "\\d+"))
        ) %>%
        select(-Variable) %>%
        left_join(transcript_tibble,by="index") %>%
        pivot_longer(cols=c(mean:Rhat),names_to="Term",values_to="Value") %>%
        select(-index)
      long_params_pSups <- params_summary %>%
        filter(grepl("^pSup",Variable)) %>%
        mutate(index = as.integer(str_extract(Variable, "\\d+")),
               Variable="pSup") %>%
        left_join(transcript_tibble,by="index") %>%
        pivot_longer(cols=c(mean:Rhat),names_to="Term",values_to="Value") %>%
        select(-index)
      long_params_lopSups <- params_summary %>%
        filter(grepl("^lopSup",Variable)) %>%
        mutate(index = as.integer(str_extract(Variable, "\\d+")),
               Variable="lopSup") %>%
        left_join(transcript_tibble,by="index") %>%
        pivot_longer(cols=c(mean:Rhat),names_to="Term",values_to="Value") %>%
        select(-index)
      long_params_other <- params_summary %>%
        filter(!grepl("mixing_factor_",Variable) &
                 !grepl("latent_counts",Variable) &
                 !grepl("pSup",Variable) &
                 !grepl("lopSup",Variable))
      stan_summary <- tibble(
        count_params = list(long_params_counts),
        mixing_params = list(long_params_mixing),
        other_params = list(long_params_other),
        pSup_est = list(long_params_pSups),
        lopSup_est = list(long_params_lopSups)
      )
      return(stan_summary)
    })) %>%
    select(Condition,stan_summary) %>% unnest(stan_summary)
}


#' write statistical summary of fit
#'
#' Write the results of get_stan_summary_reps to files
#' @details The input data is the output of the model_fit_reps function
#' Writes a data frame that includes a statistical result of parameters.
#' @param stan_summary output of get_stan_summary_reps
#' @param output_dir folder to write summaries in
#' @return NULL
#' @keywords write summary
#' @importFrom readr write_tsv
#' @import dplyr
#' @export
write_stan_summary_with_reps <- function(stan_summary,output_dir)
{
  write_tsv(stan_summary %>% select(Condition,count_params) %>%
              unnest(count_params),
            file.path(output_dir,"count_parameters.tsv"))
  write_tsv(stan_summary %>% select(Condition,mixing_params) %>%
              unnest(mixing_params),
            file.path(output_dir,"mixing_parameters.tsv"))
  write_tsv(stan_summary %>% select(Condition,other_params) %>%
              unnest(other_params),
            file.path(output_dir,"other_stan_parameters.tsv"))
  write_tsv(stan_summary %>% select(Condition,pSup_est) %>%
              unnest(pSup_est),
            file.path(output_dir,"pSup_est.tsv"))
  write_tsv(stan_summary %>% select(Condition,lopSup_est) %>%
              unnest(lopSup_est),
            file.path(output_dir,"lopSup_est.tsv"))
}
