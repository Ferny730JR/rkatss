% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{cor.pwm}
\alias{cor.pwm}
\title{Calculate the correlation coefficient between two PWM's.}
\usage{
cor.pwm(pwm1, pwm2, kmer, ret = c("cor", "pwm"))
}
\arguments{
\item{pwm1}{Position Weight Matrix to correlate}

\item{pwm2}{Position Weight Matrix to correlate}

\item{kmer}{Length of kmer to attempt to correlate from the PWM's. Should
not have a length greater than the number of rows in either PWM1/PWM2. Kmer
should be specified to be able to compare PWM's of with different number of
rows.}

\item{ret}{Return the actual correlation and p-value (cor), or the PWM which
was used for the correlation value.}
}
\value{
Correlation coefficient of the PWM's
}
\description{
Calculate the correlation coefficient between two PWM's.
}
\examples{
# Load sample data
data(rbfox2_pwms)

# Get the correlation between both PWMs
cor.pwm(rbfox2_pwms$pwm[[1]], rbfox2_pwms$pwm[[2]], kmer = 5)
}
