#############################################################################################
# Calculates the geodesic distance between two points specified by radian latitude/longitude 
# using the Haversine formula (hf) - assumes earth is spherical
# (more robust than the Spherical Law of Cosines for small distances)
#############################################################################################
gcd.hf <- function(A, B) {
  
  # Conver to radians
  long1 <- A[1]*pi/180
  lat1 <- A[2]*pi/180
  long2 <- B[1]*pi/180
  lat2 <- B[2]*pi/180
  
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}
