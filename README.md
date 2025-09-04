```markdown

# saemixPBPK

**Version:** 1.3.0  
**License:** GPL (>= 2)

## Overview

**saemixPBPK** is an R package that implements the Stochastic Approximation Expectation Maximization (SAEM) algorithm
for PBPK/QSP models developed with the OSP Suite. This package enables robust and efficient estimation procedures for
pharmacokinetic and quantitative systems pharmacology modeling.

## Features

- Estimation using the SAEM algorithm for PBPK/QSP models.
- Designed to work with models developed in the OSP Suite.
- Visualization and diagnostic tools.
- Support for simulation and model evaluation workflows.

## Authors

- Donato Teutonico ([Maintainer](mailto:donato.teutonico@gmail.com))
- Marc Lavielle
- Audrey Lavenu
- Emmanuelle Comets
- Belhal Karimi
- Maud Delattre
- Marilou Chanel
- Johannes Ranke
- Sofia Kaisaridi
- Lucie Fayette

## Installation

Please use the latest release file available on GitHub or
to test the development versionuse the following code to install from GitHub:

```R
# If devtools is not installed, install it first
install.packages("devtools")

# Install saemixPBPK from GitHub
devtools::install_github("Sanofi-Public/saemixPBPK")
```

## Dependencies

- graphics
- stats
- methods
- gridExtra
- ggplot2
- grid
- rlang
- abind
- progress

<!--  
## Usage

```R
library(saemixPBPK)
# Work in progress
```
-->

## Documentation

Documentation is provided within the package. Use `help(package = "saemixPBPK")` in R to see available functions and their documentation.

## Contributing

Contributions are welcome. Please open an issue or pull request for bugs, feature requests, or improvements.

## License

This package is licensed under GPL (>= 2).

```
