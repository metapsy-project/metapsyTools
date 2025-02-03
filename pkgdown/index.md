<h1><code style="background: white; font-size: 45px !important; font-weight: 500; color: black;">metapsyTools</code></h1>

<br>

### Overview {#overview}

The `metapsyTools` package facilitates the calculation of effect sizes and meta-analyses based on the [Metapsy](https://www.metapsy.org) database (or databases adhering to the same format).

Databases that are loaded into R using the [`metapsyData`](https://data.metapsy.org) package can automatically be analyzed using the [analysis module](https://tools.metapsy.org/articles/web/get-started#the-analysis-module) of the package.

<br></br>

### Usage Example {#usage-example}

``` r
# Load metapsyData
library(metapsyData)

# Load database, filter trials, and analyze
getData("depression-inpatients") %>% 
    filterPoolingData(
      condition_arm1 == "cbt",
      instrument == "bdi") %>% 
    runMetaAnalysis()
```

<br></br>

### Installation {#installation}

The `metapsyTools` package lives in a **GitHub repository**. It can be downloaded using the code below:

``` r
if (!require("remotes"))
  install.packages("remotes")

remotes::install_github(
    "metapsy-project/metapsyTools")
```

More details on how to install and update the package can be found in the [installation guide](articles/web/installation.html).
