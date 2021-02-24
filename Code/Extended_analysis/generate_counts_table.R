#####################################################
# Function to generate yearly sample counts per city
# Takes a data frame with a "City" column and a "Year"
# column
#####################################################

generate_counts_table <- function(city_year){
  
  unique_cities <- unique(city_year$City)
  unique_years <- unique(city_year$Year)
  
  counts_table <- array(dim = c(length(unique_cities)+1, length(unique_years)+1), 
                        dimnames = list(c(unique_cities, 'Total'), c(unique_years, 'Total')))
  
  for(city in unique_cities){
    for(year in unique_years){
      counts_table[city, as.character(year)] <- sum(city_year$City == city & city_year$Year == year)
    }
  }
  
  counts_table[,'Total'] <- rowSums(counts_table, na.rm = T)
  counts_table['Total',] <- colSums(counts_table, na.rm = T)  
  
  return(counts_table)
}



