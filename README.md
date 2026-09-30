# Commensurability and cross-study inference in quantitative species interaction networks

Code and data used to reproduce the analyses in this study, along with archived materials referenced in the manuscript.

## Repository structure

### Archived materials
- **Archived_taxonomic_name_corrections.pdf**  
  Modifications to the taxonomic identifiers (i.e., string names) of nodes in networks.

- **Archived_null_network_corrections.pdf**  
  (a) Percent differences between Patefield null models (Patefield, 1981) and empirical networks for specialization, nestedness, and modularity.  
  (b) Percent differences between Vaznull null models (Vazquez et al., 2007) and empirical networks for specialization, nestedness, and modularity.

- **Archived_networks_and_publication_sources.pdf**  
  A list of the open-access bipartite networks and their publication sources used in this study.

### Data
All input data live in `data/`.

- **data/general_network_information.csv**  
  Information on network type (pollination or seed-dispersal), publication grouping, and latitude and longitude for each network.

- **data/all_results_for_patefield_randomization.csv**  
  Observed network indices (weighted specialization, weighted nestedness, weighted modularity) and the corresponding results from 1000 Patefield randomizations for each network.

- **data/all_results_vaznull_randomizations.csv**  
  Observed network indices (weighted specialization, weighted nestedness, weighted modularity) and the corresponding values obtained from 1000 Vázquez null-model randomizations for each network.
  
- **data/rarefaction_sensitivity_500_iterations.csv**  
  Network indices (weighted specialization, weighted nestedness, weighted modularity) calculated after rarefying observed networks to 50, 100, and 200 interactions, with 500 iterations at each rarefaction depth.

- **data/fricke_metadata.csv**  
  Seed-dispersal network metadata from Fricke, E. C. and J.-C. Svenning (2020), *Accelerating homogenization of the global plant–frugivore meta-network*.

- **data/pollination_sampling_metadata.csv**  
  Pollination network metadata compiled by us for the analyses in this repository.
  
### Network data
- **data/networks/**  
  All 279 bipartite networks used in this study.

- **data/networks/column_names/**  
  Taxonomic identities for all animal nodes appearing in the network datasets.
  
- **data/networks/row_names/**  
  Taxonomic identities for all plant nodes appearing in the network datasets.  

### Code
All scripts live in `R/`. Each script reads from `data/` and writes its figures and tables to `output/` (created automatically if missing). Paths are resolved with the `here` package, so scripts can be run from any working directory inside the repository.

- **R/beta_diversity.R**  
  Generates Figure 3 and Table S1, comparing within- and between-publication beta diversity by guild and network type, and evaluating the effects of geographic distance and publication identity on beta diversity.

- **R/linear_models_vs_linear_mixed_models.R**  
  Generates Table 2, comparing the relationships between latitude and network indices (weighted specialization, weighted nestedness, and weighted modularity) using linear models and linear mixed models that account for publication-level variation.

- **R/taxonomic_resolution_analysis.R**  
  Generates Figures 5 and S2 and Table S2, summarizing the percentage of nodes lacking species-level identification and their finest available taxonomic resolution.

- **R/patefield_null_model_analysis.R**  
  Generates Figure S4 and Tables 4 and S6–S8, comparing observed network indices with 1,000 Patefield null-model realizations and summarizing within- and between-publication variation in null-corrected indices.
  
- **R/vaznull_model_analysis.R**  
  Generates Figure S5 and Tables 5 and S9–S11, comparing observed network indices with 1,000 Vázquez null-model realizations and summarizing within- and between-publication variation in null-corrected indices.
  
- **R/sampling_intensity_and_network_topology.R**  
  Generates Figure 6, illustrating the relationships between network size, sampling intensity, and network indices (weighted specialization, weighted nestedness, and weighted modularity).
  
- **R/interaction_recording_and_edge_weight_units.R**  
  Generates Figure 4, illustrating the publication-collapsed frequencies of sampling methods and edge-weight definitions.
  
- **R/rarefaction.R**  
  Generates Tables 3 and S3–S5, summarizing network indices and comparing within- and between-publication variation after rarefaction to 50, 100, and 200 interactions.
  
- **R/map_of_locations.R**  
  Generates Figure S1, illustrating the geographic distribution of pollination and seed-dispersal networks, with point sizes indicating the number of networks at each location.
  
## Installing dependencies

The R packages used by the scripts in `R/` are listed in the `DESCRIPTION` file. From the repository root, install them with either:

```r
# using pak
install.packages("pak")
pak::pak()

# or using remotes
install.packages("remotes")
remotes::install_deps()
```

## Running the analysis

`run_analysis.R` sources every script in `R/` in turn and reports progress for each one. From the repository root:

```sh
Rscript run_analysis.R
```

or, from an R session opened in the repository:

```r
source("run_analysis.R")
```

Scripts can also be run individually. All figures (PDF) and tables (CSV, plus one TXT file with the model summaries for Table 2) are written to `output/`.

## Contact

For questions about the data or code, please contact:

Chris Brimacombe  
- University of Guelph  
- Email: cbrimaco@uoguelph.ca
