# msc.scrnaseq

## Contents
[Description](#description) | [Installation](#installation) |

## Description

<br>

In the interest of transparent and reproducible science, ```msc.scrnaseq``` was created to not only house all the analysis scripts of this project, but make them executable. This means a user will be able to replicate our core analysis as well as expand upon it, using the functions provided here. This is made up of the following four functions : 

* ```preprocess```
* ```homologize```
* ```integrate_samples```
* ```integrate_species```

## Installation 
### Prerequisites
Our analysis, and in consequence this package, makes use of the following pre-existing packages :
```cowplot```, ```DoubletFinder```, ```DropletQC```, ```dplyr```, ```ggplot2```, ```harmony```, ```Matrix```,
```Seurat```, ```SoupX```, ```stringr```, ```patchwork```, ```purrr```.

<br>

These must be loaded into your local R environment before installation. This can be done altogether as shown below. 

<br>

```R
packages <- c("cowplot", "DoubletFinder", "dplyr", "ggplot2", "harmony", "Matrix", "Seurat", "SoupX", "stringr", "patchwork", "purrr")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
        library(pkg)
    }
}

```

<br>

Please note that ```DropletQC``` needs to be installed from ```GitHub``` using the ```remotes``` package. 

<br>

```R
# install.packages("devtools")
library(devtools)
devtools::install_github("powellgenomicslab/DropletQC", build_vignettes = TRUE)
library(DropletQC
```

<br>

### Package installation
The ```msc.scrnaseq``` package can then be installed and loaded from ```GitHub``` as shown below. 

<br>

```R
devtools::install_github("AlicenJoyHenning/msc.scrnaseq")
library(msc.scrnaseq)
```

<br>

## Verify installation

To ensure the package has installed correctly, run the following to see if you can view the help pages. 

```R
?preprocess
?homologize
?integrate_samples
?integrate_species
```




