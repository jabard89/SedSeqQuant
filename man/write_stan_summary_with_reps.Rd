% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fit_reps.R
\name{write_stan_summary_with_reps}
\alias{write_stan_summary_with_reps}
\title{write statistical summary of fit}
\usage{
write_stan_summary_with_reps(stan_summary, output_dir)
}
\arguments{
\item{stan_summary}{output of get_stan_summary_reps}

\item{output_dir}{folder to write summaries in}
}
\description{
Write the results of get_stan_summary_reps to files
}
\details{
The input data is the output of the model_fit_reps function
Writes a data frame that includes a statistical result of parameters.
}
\keyword{summary}
\keyword{write}
