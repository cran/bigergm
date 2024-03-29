% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hergm.R
\name{install_python_dependencies}
\alias{install_python_dependencies}
\title{Install optional Python dependencies}
\usage{
install_python_dependencies(envname = "r-bigergm", method = "auto", ...)
}
\arguments{
\item{envname}{The name, or full path, of the environment in which Python packages
are to be installed. When NULL (the default), the active environment as set
by the RETICULATE_PYTHON_ENV variable will be used; if that is unset, then the
r-reticulate environment will be used.}

\item{method}{Installation method. By default,
"auto" automatically finds a method that will
work in the local environment. Change the
default to force a specific installation method.
Note that the "virtualenv" method is not available on Windows.}

\item{...}{Additional arguments, to be passed to lower-level functions}
}
\value{
No return value, called for installing the Python dependencies 'infomap' and 'numpy'
}
\description{
Install Python dependencies needed for using the Python implementation of infomap
}
