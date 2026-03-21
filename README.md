# Commensurability and cross-study inference in quantitative species interaction networks

Code and data used to reproduce the analyses in this study.

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
- **general_network_information.csv**  
  Information on network type (pollination or seed-dispersal), publication grouping, and latitude and longitude for each network.

- **all_results_for_patefield_randomization.csv**  
  Observed network indices (weighted specialization, weighted nestedness, weighted modularity) and the corresponding results from 1000 Patefield randomizations for each network.

- **all_results_vaznull_randomizations.csv**  
  Observed network indices (weighted specialization, weighted nestedness, weighted modularity) and the corresponding values obtained from 1000 Vázquez null-model randomizations for each network.

- **fricke_metadata.csv**  
  Seed-dispersal network metadata from Fricke, E. C. and J.-C. Svenning (2020), *Accelerating homogenization of the global plant–frugivore meta-network*.

- **pollination_sampling_metadata.csv**  
  Pollination network metadata compiled by us for the analyses in this repository.

### Network data
- **networks/**  
  All 279 bipartite networks used in this study.

- **networks/column_names/**  
  Taxonomic identities for all animal nodes appearing in the network datasets.
  
- **networks/row_names/**  
  Taxonomic identities for all plant nodes appearing in the network datasets.  

### Code
- **beta_diversity.R**  
  Generates Figure 3, showing beta diversity among guilds (plants, animals) and publication groupings ("One network per publication", "Multiple networks per publication") within each network type.

- **linear_models_vs_linear_mixed_models.R**  
  Generates Table 2 by fitting linear models and linear mixed models to explain weighted specialization, weighted nestedness, and weighted modularity in pollination and seed-dispersal networks.

- **nodes_lacking_species_level_identification.R**  
  Generates Figures 5 and S2, and Table S2. Quantifies: (i) the percentage of nodes lacking species-level identification in each network, summarized by publication grouping (Figure 5), (ii) the percentage of nodes lacking species-level identification in each network across all networks (Figure S2), and (iii) the finest available taxonomic resolution for nodes lacking species-level identification (Table S2).  

- **nullModel_compared_to_empirical_patefield.R**  
  Generates Figure S4 and Tables 3, S3, S4, S5, and S6. Quantifies: (i) weighted specialization, weighted nestedness, and weighted modularity for each network and its corresponding 1000 Patefield null-model realizations, summarized by publication grouping (Figure S4); (ii) the standard deviation of delta values for these indices across 1000 Patefield null-model realizations for "One network per publication" and "Multiple networks per publication", separated by pollination and seed-dispersal (Table 3); (iii) the exact values underlying Figure S4 for the observed networks, along with the mean and standard deviation of delta values from 1000 Patefield null-model realizations for each index by publication grouping (Table S3); (iv) the standard deviation of z-scores for these indices relative to 1000 Patefield null-model realizations for each index in "One network per publication" and "Multiple networks per publication", separated by pollination and seed-dispersal (Table S4); (v) the mean and standard deviation of z-scores for these indices relative to 1000 Patefield null-model realizations for each publication grouping (Table S5); and (vi) the percent difference between each observed network index and the mean of its 1000 Patefield null-model realizations (Table S6).
  
- **nullModel_compared_to_empirical_vaznull.R**  
  Generates Figure S5 and Tables 4, S7, S8, S9, and S10. Quantifies: (i) weighted specialization, weighted nestedness, and weighted modularity for each network and its corresponding 1000 Vázquez null-model realizations, summarized by publication grouping (Figure S5); (ii) the standard deviation of delta values for these indices across 1000 Vázquez null-model realizations for "One network per publication" and "Multiple networks per publication", separated by pollination and seed-dispersal (Table 4); (iii) the exact values underlying Figure S5 for the observed networks, along with the mean and standard deviation of delta values from 1000 Vázquez null-model realizations for each index by publication grouping (Table S7); (iv) the standard deviation of z-scores for these indices, relative to 1000 Vázquez null-model realizations, for "One network per publication" and "Multiple networks per publication", separated by pollination and seed-dispersal (Table S8); (v) the mean and standard deviation of z-scores for these indices relative to 1000 Vázquez null-model realizations for each publication grouping (Table S9); and (vi) the percent difference between each observed network index and the mean of its 1000 Vázquez null-model realizations (Table S10).
  
- **sampling_intensive_and_network_size_and_topological_indices.R**  
  Generates Figure 6 (four panels), showing the relationships between (i) sampling intensity and network size, (ii) weighted specialization and sampling intensity, (iii) weighted nestedness and sampling intensity, and (iv) weighted modularity and sampling intensity.
  
- **sampling_method_and_edge_weight.R**  
  Generates Figure 4 (two panels), showing publication-collapsed frequencies of (i) sampling methods used to observe species interactions and (ii) edge-weight definitions.
  
## Contact

For questions about the data or code, please contact:

Chris Brimacombe  
- University of Guelph  
- Email: cbrimaco@uoguelph.ca