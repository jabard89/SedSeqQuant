#' Prepare Data
#'
#' This function prepares count data by filtering transcripts with low counts, rounding counts to integers, and
#' dividing data into nested fractions. It then converts the data to a count matrix of transcript_IDs and Reps.
#'
#' @param count_data A data frame containing count data (output of load_counts).
#' @param min_counts Minimum count threshold for filtering low counts (default = 20).
#'
#' @return A tibble containing a column of prepared data in the form of a named list:
#' \itemize{
#'   \item{"NRNA"}{Number of RNA sequences.}
#'   \item{"tot_obs_counts"}{Total observed counts vector}
#'   \item{"sup_obs_counts"}{Supernatant observed counts vector}
#'   \item{"pel_obs_counts"}{Pellet observed counts vector}
#' }
#' @export
#' @importFrom dplyr select filter group_by arrange mutate
#' @importFrom tidyr pivot_wider pivot_longer nest
#' @importFrom purrr map
#' @importFrom tibble rownames_to_column
#' @importFrom stats setNames
#' @importFrom magrittr %>%
prepare_data_noreps <- function(count_data,min_counts=20){
  filter_low_counts <- function(count_data, min_counts) {
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
    group_by(Condition,Rep) %>% nest %>%
    # divide up data into a list of fractions
    mutate(prepared_data=map(data,function(df){
      nested_data <- df %>% group_by(Fraction) %>% split(.$Fraction) %>%
        # convert data to a count matrix of transcript_IDs and Reps
        map(function(x){
          counts <- x %>%
              column_to_rownames(var="transcript_ID") %>%
              select(-Fraction) %>% pull(Count)
          names(counts) <- x$transcript_ID
          return(counts)
        })
      prepared_data <- list(
        "NRNA"=length(nested_data$Total),

        tot_obs_counts=nested_data$Total,
        sup_obs_counts=nested_data$Supernatant,
        pel_obs_counts=nested_data$Pellet
      )
      return(prepared_data)
    }))
  return(nested_data %>% select(Condition,Rep,prepared_data))
}

#' Model fit
#'
#' Fit the bayesian statistical model to the counts data.
#' @details The input data should be the wide format in terms of fractions.
#'
#' @param nested_data output of prepare_data_with_noreps
#' @param chains A number, defaulting to 4
#' @param iter A number, defaulting to 1000
#' @param ... Other arguments passed to `rstan::sampling` (e.g. warmup, seed).
#' @return A dataframe with a column containing the class `stanfit` returned by `rstan::sampling`
#' @seealso [rstan::sampling()], `browseVignettes("rstan")`
#' @keywords model fit
#' @importFrom rstan stan_model sampling summary
#' @export
#'
model_fit_noreps <- function(nested_data, chains = 4,
                           iter = 1000,...)
{
  conflicted::conflicts_prefer(rstan::extract)
  conflicted::conflicts_prefer(stats::lag)
  nested_stanfit <- nested_data %>%
    mutate(stanfit = map(prepared_data,function(x){
      rstan::sampling(stanmodels$SedSeqQuantNoReps, data = x,
                      chains = chains, iter = iter,
                      ...)
    })) %>%
    select(Condition,Rep,stanfit)
  return(nested_stanfit)
}
#' parameters statistical results
#'
#' Get the estimated statistical results of parameters from a stan model
#' Extracts the names of transcript_IDs from the input data
#' @details The input data is the output of the rstan::sampling function
#' Output a data frame that includes a statistical result of parameters.
#' @param nested_stanfit output of rstan::sampling
#' @param nested_data input data into rstan::sampling
#' @return A data frame
#' @keywords parameter extraction
#' @import dplyr
#' @export
get_stan_summary_noreps <- function(nested_stanfit,nested_data)
{
  temp_data <- full_join(nested_stanfit,nested_data,by=c("Condition","Rep")) %>%
    mutate(stan_summary = map2(stanfit,prepared_data,function(sf,pd){
      item <- rstan::summary(sf)
      transcript_tibble <- tibble(transcript_ID=names(pd$tot_obs_counts),
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
                               Rhat = item$summary[,"Rhat"]) %>%
        pivot_longer(cols=c(-Variable),names_to="Term",values_to="Value")
      # extract the scaling parameter for total (relative to Rep=1)
      return(params_summary)
    })) %>%
    mutate(pSup_est = map2(stan_summary,prepared_data,function(ss,pd){
      mix_sup <- ss %>% filter(Term=="mean",
                                     Variable=="mixing_factor_sup") %>% pull(Value)
      mix_pel <- ss %>% filter(Term=="mean",
                                     Variable=="mixing_factor_pellet") %>% pull(Value)
      est <- mix_sup*pd$sup_obs_counts/(mix_sup*pd$sup_obs_counts+
                                               mix_pel*pd$pel_obs_counts)
      return(tibble(transcript_ID=names(pd$tot_obs_counts),pSup_est=est))
    })) %>%
    select(Condition,Rep,stan_summary,pSup_est)
}

#' write statistical summary of fit
#'
#' Write the results of get_stan_summary to files
#' @details The input data is the output of the rstan::sampling function
#' Output a data frame that includes a statistical result of parameters.
#' @param stan_summary output of get_stan_summary
#' @param output_dir folder to write summaries in
#' @return NULL
#' @keywords write summary
#' @importFrom readr write_tsv
#' @import dplyr
#' @export
write_stan_summary_noreps <- function(stan_summary,output_dir)
{
  write_tsv(stan_summary %>% select(Condition,Rep,stan_summary) %>%
              unnest(stan_summary),
            file.path(output_dir,"stan_summary.tsv"))
  write_tsv(stan_summary %>% select(Condition,Rep,pSup_est) %>%
              unnest(pSup_est),
            file.path(output_dir,"pSup_estimates.tsv"))
}
