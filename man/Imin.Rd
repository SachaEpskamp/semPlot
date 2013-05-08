\name{Imin}
\alias{Imin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helper function to substract matrix from identity matrix and take inverse.
}
\description{
This function can be used to more easilly compute I - X or (I - X)^(-1), which are common in SEM models.
}
\usage{
Imin(x, inverse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix
}
  \item{inverse}{
Logical, should the inverse be taken?
}
}

\author{
Sacha Epskamp <mail@sachaepskamp.com>
}
