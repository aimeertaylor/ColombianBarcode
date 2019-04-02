##############################################
# This script simply loads a very large file 
# and saves more manageable versions
##############################################
rm(list = ls())
load('../RData/mles.RData') # Few mins
mle_core = mle_df[,1:4] # Extract core

# Save a version withs rhat CIs 
# rather than bootstrapped quantities
P = c(0.025, 0.975) # CI quantiles
rhat_CIs = t(apply(mle_df, 1, function(x)quantile(x$rhats_boot, probs=P))) # Few seconds
mle_CIs = cbind(mle_core, rhat_CIs)
save(mle_CIs, file = '../RData/mle_CIs.RData')

# Save a more manageable version withs only
# 100 rhat bootstrapped quantities 
rhats_boot_100 = lapply(mle_df$rhats_boot, function(x) x[1:100])
mle_core$rhats_boot = rhats_boot_100
mle_boot = mle_core
save(mle_boot, file = '../RData/mle_boot.RData')

