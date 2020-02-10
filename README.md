# ColombianBarcode

To reproduce the work in the manuscript entitled "Identifying departures from isolation-by-distance:  connectivity between P. falciparum populations on the Colombian-Pacific coast" complete the following four steps

1. Run in no particular order

- DataFormat.R (ouputs SNPData.RData and hmmInput.txt)
- Generate_geo_dist_info.R (outputs LonLat.RData, geo_dist_info_states.RData, and geo_dist_info_cities.RData)

2. Run in the following order

- Generate_mles_CIs.R (SLOW, outputs mles_frequencies_true.RData and mles_frequencies_unif.RData)
- Format_mle_df.R (ouputs All_results_true.RDat and All_results_unif.RData)

3. Run in no particular order

- Generate_transportdistances_time.R (outputs W_results_t.RData)
- Generate_transportdistances_city.R (outputs All_W_results_c.RData)
- Generate_proportions.R (outputs proportions_time.RData and proportions_geo.RData)
- Generate_proportions_sensitivity.R (outputs proportions_sensitivities.RData)

4. Run in no particular order

- Plot_network.R (outputs Colombia_network%d.png and Colombia_Map.png) 
- Plot_proportions_distances.R (outputs Proportions_and_W_distance.pdf)
- Plot_proportions_sensitivity.R (outputs Proportions_sensitivity.pdf)
- Plot_graph_components.R (outputs All_CCs.pdf, Graphs_timespace.pdf and Graphs.pdf)
