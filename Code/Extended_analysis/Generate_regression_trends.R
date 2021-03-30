load("../../RData/metadata_extended.RData")
load("../../RData/mles_CIs_extended_freqsTaylor2020_meta.RData")

mle_CIs$IntraCity <- NA
mle_CIs$IntraCity[mle_CIs$City1 == mle_CIs$City2] <- mle_CIs$City1[mle_CIs$City1 == mle_CIs$City2]

sids_Taylor2020 <- metadata$SAMPLE.CODE[metadata$PloSGen2020]
mles_CIs_Taylor2020 <- mle_CIs[mle_CIs$individual1 %in% sids_Taylor2020 & mle_CIs$individual2 %in% sids_Taylor2020, ]

# For all estimates - need to add intra site label as in TM border
glmfitted <- glm(highly_related ~ Coast + Port + year_diff + geo_dist, 
                 data = mle_CIs)
summary(glmfitted)

glmfitted <- glm(highly_related ~ Coast + Port + year_diff + geo_dist, 
                 data = mles_CIs_Taylor2020)
summary(glmfitted)


# For all inter-site estimates 
glmfitted <- glm(highly_related ~ Port + Coast + year_diff + geo_dist, data = mle_CIs[mle_CIs$geo_dist > 0,])
summary(glmfitted)

glmfitted <- glm(highly_related ~ Port + Coast + year_diff + geo_dist, data = mles_CIs_Taylor2020[mles_CIs_Taylor2020$geo_dist > 0,])
summary(glmfitted)


