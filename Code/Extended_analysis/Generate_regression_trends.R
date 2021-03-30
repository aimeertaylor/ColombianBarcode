rm(list = ls())
load("../../RData/metadata_extended.RData")
load("../../RData/mles_CIs_extended_freqsTaylor2020_meta.RData")
eps <- 0.01 # Used to define clones in Taylor et al. 2020

# # Filter uninformative
uninformative_ind <- mle_CIs$r2.5. < eps & mle_CIs$r97.5. > 1-eps; sum(uninformative_ind)
mle_CIs <- mle_CIs[!uninformative_ind,]

# Generate a data set of samples that featured in Taylor et al. 2020
sids_Taylor2020 <- metadata$SAMPLE.CODE[metadata$PloSGen2020]
mles_CIs_Taylor2020 <- mle_CIs[mle_CIs$individual1 %in% sids_Taylor2020 & mle_CIs$individual2 %in% sids_Taylor2020, ]

# # For all estimates - need to add intra site label as in TM border
# glmfitted <- glm(highly_related ~ Coast + Port + year_diff + geo_dist, 
#                  data = mle_CIs)
# summary(glmfitted)
# 
# glmfitted <- glm(highly_related ~ Coast + Port + year_diff + geo_dist, 
#                  data = mles_CIs_Taylor2020)
# summary(glmfitted)


# For all inter-site estimates 
glmfitted <- glm(highly_related ~ Port + Coast + year_diff + geo_dist, data = mle_CIs[mle_CIs$geo_dist > 0,])
summary(glmfitted) 

glmfitted <- glm(highly_related ~ Port + Coast + year_diff + geo_dist, data = mles_CIs_Taylor2020[mles_CIs_Taylor2020$geo_dist > 0,])
summary(glmfitted)


