<br>

<p align="left">
  <img src="https://github.com/AlicenJoyHenning/msc.scrnaseq/blob/main/inst/extdata/msc.rnaseq.png" alt="limiric_logo" height="47" width="380">
</p>

🔵[Description](#-description)   ⚪[Installation](#-installation)   🟣[Quickstart](#-quickstart)

<br>

## 🔵 Description

<br>

In the interest of transparent and reproducible science, ```msc.scrnaseq``` was created to not only house the analysis scripts of this project, but make them executable in any environment. This means a user will be able to replicate our core analysis as well as expand upon it. ```msc.scrnaseq``` is made up of the following four functions : 

* ```preprocess```
* ```homologize```
* ```integrate_samples```
* ```integrate_species```

<br>
<br>


## ⚪ Installation 
### Prerequisites
Our analysis, and in consequence this package, makes use of the following pre-existing packages : <br>
```cowplot```, ```DoubletFinder```, ```DropletQC```, ```dplyr```, ```ggplot2```, ```harmony```, ```Matrix```,
```Seurat```, ```SoupX```, ```stringr```, ```patchwork```, ```purrr```.

<br>

These must be loaded into your local R environment before installation. This can be done altogether as shown below. 

```R
packages <- c("cowplot", "DoubletFinder", "dplyr", "ggplot2", "harmony", "Matrix", "Seurat", "SoupX", "stringr", "patchwork", "purrr")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
        library(pkg)
    }
}

```

> However, ```DropletQC``` needs to be installed from ```GitHub``` using the ```remotes``` package. 
>
> ```R
> # install.packages("devtools")
> library(devtools)
> devtools::install_github("powellgenomicslab/DropletQC", build_vignettes = TRUE)
> library(DropletQC)
> ```
> 

<br>

### Package installation
The ```msc.scrnaseq``` package can then be installed and loaded from ```GitHub``` as shown below. 

<br>

```R
devtools::install_github("AlicenJoyHenning/msc.scrnaseq")
library(msc.scrnaseq)
```

<br>

### Verify installation

To ensure the package has installed correctly, run the following to see if you can view the help pages. 

```R
?preprocess
?homologize
?integrate_samples
?integrate_species
```

<br>
<br>

## 🟣 Quickstart 
### ```preprocess```
This function is responsible for the heavy lifting preprocessing work. This includes ambient RNA correction with ```SoupX```, removing red blood cells, isolating immune cells, filtering damaged cells detected by ```DropletQC``` and ```limiric```, and filtering doublets detected by ```DoubletFinder```. It also provides a preliminary check of sample quality by immune cell marker visualisation (UMAP). Running this function requires the inputs mentioned below and will deposit outputs, including diagnostic plots and processed ```Seurat``` object, in the following directory structure.  

```
output_dir/
|
├── RBCQC
├── CellQC
├── IMCQC
├── SampleClustering
└── R_objects 

```

####     _Inputs_

| **Parameter**      | **Default/Required** | **Description**                                                                                             |
|--------------------|----------------------|-----------------------------------------------------------------------------------------------------------  |
| **project_name**   | Required             | String with project or sample name.                                                                         |
| **organism**       | Required             | String with one of the following Mmul, Hsap, Mmus.                                                          |
| **ptprc_threshold**| Default: 0            | Percentage CD45 over which cells will be retained. E.g: If 0, all cells must express CD45                  |
| **filtered_path**  | Required             | Directory of filtered alignment output, matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz.                    |
| **raw_path**       | Reuired             | Directory of raw alignment output, matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz.                          |
| **velocyto_path**  | Required             | Directory of Velocyto filtered alignment output, spliced.mtx.gz, unspliced.mtx.gz, ambiguous.mtx.gz, features.tsv.gz, barcodes.tsv.gz. |
| **filter_damaged** | Default: TRUE        | Should output contain no damaged cells.                                                  |
| **output_dir**     | Required             | Directory where preprocess_core output should be generated.                                                |
| **sample_list**    | Default: NULL        | Input multiple samples in list.                                                          |
> Please note that those marked as ```required``` need to be given else the function will not work. 

<br>

_Example_ 
```R
# Define a sample list 
samples <- list(

  list(project_name = "SRR1234567",
      filtered_path = "/home/username/alignment/SRR1234567/filtered",
      raw_path = "/home/username/alignment/SRR1234567/raw",
      velocyto_path = "/home/username/alignment/SRR1234567/velocyto",
      output_dir = "/home/username/postalignment/preprocess"),

  list(project_name = "SRR1234568",
      filtered_path = "/home/username/alignment/SRR1234568/filtered",
      raw_path = "/home/username/alignment/SRR1234568/raw",
      velocyto_path = "/home/username/alignment/SRR1234568/velocyto",
      output_dir = "/home/username/postalignment/preprocess"))

# Run the preprocess function
 GSE1234567 <- preprocess(sample_list = samples)

```

####     _Outputs_ 
As mentioned, the outputs include plots and a processed ```Seurat``` object. Plots are divided into directories according to the quality control process the are associated with, with ```RBCQC``` and ```IMCQC``` showing scatter plots of filtered red blood cells and isolated immune cells, respectively. Similarily, the ```CellQC``` directory houses tSNE and scatter plots where damaged cells are marked. The last plot is found in the ```SampleClustering``` directory and contains isolated visualisations of the cells in a sample. This gives a good indication of what types of immune cells are present in the sample and in what proportion. Lastly, the processed ```Seurat``` objects are housed in the ```R_objects``` output directory and will continue to the next stage where homologous genes are found. 

<br>
<br>

### ```homologize```
This function is responsible for converting genes in a ```Seurat``` object belonging to one organism, such as a mouse, to homologous genes from another organism, such as a human. Both input and output organism must be specified as one of the following: _Hsap_, _Mmus_, or _Mmul_. The inputs required are detailed in the table below and include the directory to the preprocessed ```Seurat``` object, ideally that generated from the previous ```preprocess``` function but any ```Seurat``` object will be suitable as input for this function. The output of this function is just the ```Seurat``` object with gene converted to their homologous counterparts. This is stored in a directory called ```R_objects```. 

```
output_dir/
└── R_objects 
```

####     _Inputs_
| **Parameter**      | **Default/Required** | **Description**                                                                                            |
|--------------------|----------------------|----------------------------------------------------------------------------------------------------------- |
| **project_name**   | Required             | String with project or sample name.                                                                        |
| **seurat_dir**     | Required             | Path to the Seurat object used as input.                                                                   |
| **organism_in**    | Required             | Organism of the input Seurat object (Hsap, Mmul, Mmus)                                                     |
| **organism_out**   | Required             | Organism for which gene homologs must be found                                                             |
| **check_umap**     | Default: TRUE        | Whether or not marker UMAPs should be generated.                                                           |
| **output_dir**     | Required             | Directory where homologize_core output should be generated.                                                |
| **sample_list**    | Default: NULL        | Input multiple samples in list.                                                                            |


<br>

_Example_
```R
# Define samples
samples <- list(

  list(project_name = "SRR1234567",
       seurat_dir   = "/home/username/postalignment/preprocess/R_objects/SRR1234567.rds",
       organism_in  = "human",
       organism_out = "mouse",
       output_dir   = "/home/username/postalignment/homologize"),

 list(project_name  = "SRR1234568",
       seurat_dir   = "/home/username/postalignment/preprocess/R_objects/SRR1234568.rds",
       organism_in  = "human",
       organism_out = "mouse",
       output_dir   = "/home/username/postalignment/homologize"),

# Run the homologize function
GSE1234567 <- homologize(sample_list = samples)

```

<br>

### ```integrate_samples```
This is an optional function that uses ```Seurat```'s canonical correlation method to integrate samples of the same species origin. This can be helpful to visualise what cell types are detectable in the context of an individual species with a high level of granularity. This is because integration and clustering is performed on the full set of genes in the organism's genome, rather than being isolated to those that have homologs. The output of this function is a single plot and ```Seurat``` object. This plot shows the reduced form of the integrated samples, with each sample marked in a different colour to confirm the absence of sample-specific batch effects, in a UMAP alongside the visualisation of marker genes. The ```Seurat``` object contains the integrated data from all input samples. Please note that, depending on how many samples are added as input, this is likely going to be a _large_ object. 

<br>

#### _Inputs_ 
| **Parameter**      | **Default/Required** | **Description**                                                                                           |
|--------------------|----------------------|-----------------------------------------------------------------------------------------------------------|
| **project_name**   | Required             | String with project name.                                                                                 |
| **organism**       | Required             | Which organism do all inputs belong to, either Hsap, Mmul, or Mmus.                                       |
| **sample_list**    | Required             | Required list of directories for input Seurat objects.                                                    |
| **output_dir**     | Required             | Directory where integrate_samples output should be generated.                                             |

<br>

_Examples_
```R
# Define input list
   samples <- list(
     "/home/user/postalignment/preprocess/R_objects/SRR1234567.rds",
     "/home/user/postalignment/preprocess/R_objects/SRR1234567.rds"
     )

   # Run the integrate samples function
   GSE1234567 <- integrate_samples(
     project_name = "mouse_samples",
     organism     = "mouse",
     sample_list  = samples,
     output_dir   = "/home/user/postalignment/integrate_samples/"
   )
```

<br>

### ```integrate_species```
This is the main integration function that takes the input ```Seurat``` objects from the ```homologize``` function. This means the input will contain different species, but all with a common set of genes. This is essential else integration is not possible. As with ```integrate_samples```, this requires a few default inputs but includes an additional option for selection of an integration method. The options available currently include ```cca``` and ```harmony```, of which ```harmony``` is used as default as it is best at dealing with inter-species batch effects as well as inter-sample batch effects. The output for this function, as with the above ```integrate_samples``` function, is a plot showing the cells in a reduced dimensional view coloured according to sample of origin as well as a single ```Seurat``` object. Please note, as well, that this object is going to be very large as it contains information from all samples across species. 

####     _Inputs_
| **Parameter**      | **Default/Required** | **Description**                                                                                           |
|--------------------|----------------------|-----------------------------------------------------------------------------------------------------------|
| **project_name**   | Required             | String with project name.                                                                                 |
| **organism**       | Required             | Organism that all samples gene annotations conform to.                                                    |
| **sample_list**    | Required             | Required list of directories for input Seurat objects.                                                    |
| **method**         | Default: harmony     | Integration method to use, either CCA or harmony.                                     |
| **output_dir**     | Required             | Directory where integrate_species output should be generated.                                             |

<br>

_Example_
```R
  # Define input list, here all are defined in terms of mouse genes
   samples <- list(
     "/home/user/postalignemnt/Homologize/R_objects/human_1_mouse.rds",
     "/home/user/postalignemnt/Homologize/R_objects/human_2_mouse.rds",
     "/home/user/postalignemnt/Homologize/R_objects/mouse_1_mouse.rds",
     "/home/user/postalignemnt/Homologize/R_objects/mouse_2_mouse.rds",
     "/home/user/postalignemnt/Homologize/R_objects/macaque_1_mouse.rds",
     "/home/user/postalignemnt/Homologize/R_objects/macaque_2_mouse.rds")

   # Run the integrate species function
   output <- integrate_species(
     project_name = "mouse_genes_harmony",
     organism     = "mouse",
     sample_list  = samples,
     output_dir   = "/home/user/postalignment/integrate_species/"
   )
```
