######################################################
# Effect of unequal COI and transmission intensity
# 
# COI is an almost perfect correlate of cases per 1000
# people 2000-2007 [Echeverry et al. 2013]
#
# Rationale for cor.test:
# Under null hypothesis of equal connectivity, we might
# expect a negative correction between the proportion
# of samples with COI > 1 and the fraction of highly 
# related sample pairs for the reasons stated below. 
# 
# Result: non-significant negative correlation 
# -0.19 (-0.73, 0.50) driven entirely by Cauca (Guapi)
# If Cauca removed, correlation 0.32 (-0.44, 0.81)
# 
# Reason 1) unequal proportions of multiclonal samples: 
# removal of multiclonal parasite samples
# with COI > 1 is likely to remove rare parasites. 
# Since the vast majority of parasites are presumably 
# domestic versus migrant, removal of parasite samples 
# with COI > 1 is likely to remove migrant parasites. 
# In that case, if connectivity is equal between sites, 
# expect to see fewer highly related samples pairs between
# sites with greater proportions of parasite samples with
# COIs > 1 removed, i.e. a negative correlation between 
# the fraction of highly related samples pairs and the 
# weighted (to account for unequal sample sizes) average
# proportion of COIs > 1. 
#
# Unequal transmission intensities: 
# The effective population of malaria parasites is likely
# larger in high versus low transmission settings. As such, 
# the denominator used to calculate the faction of highly
# related sample pairs with be larger in high transmission 
# settings rending the fraction itself smaller and leading 
# to a negative correlation between the fraction of highly
# related samples pairs and transmission intensity. Since 
# COI is an almost perfect correlate of cases in Echeverry
# et al. 2013, we can use the weighted average proportion of 
# COIs > 1 as a correlated of cases and thus transmission 
# intensity. 
######################################################
rm(list = ls()) 
load('../RData/SNPData.RData') # For weights
load('../RData/proportions_sensitivities.RData')

# Extract COI > 1 proportions from paper (available at state level)
COI <- c('Choco' = 0.25, 'Valle' = 0.14, 'Cauca' = 0.14, 'Narino' = 0.19)

# Load highly related fractions
fract_highly_related <- proportions_states["mean",,"0.5","Unfiltered"]

# Create weights based on monoclonal sample sizes
weights <- table(SNPData$STATE)/sum(table(SNPData$STATE))

# State comparisons
state_comparisons = do.call(rbind, strsplit(names(fract_highly_related), split = '_'))

# Compute weighted average proportions of COI > 1
weighted_avCOIs <- sapply(1:length(fract_highly_related), function(i){
  states = state_comparisons[i, ]
  sum(COI[states] * weights[states]) * 0.5
})

# Indicator of same or different state (for plot colour)
same_state = state_comparisons[,1] == state_comparisons[,2]

# Plot
plot(x = weighted_avCOIs, y = fract_highly_related, 
     pch = 20, col = same_state+1, 
     ylab = "Fraction highly related",
     xlab = "(correlate of transmission intensity)")
text(x = weighted_avCOIs, y = fract_highly_related, 
     labels = names(fract_highly_related), cex = 0.5, pos = 2)
title(xlab = "Weighted proportion of removed samples with COI > 1", line = 2)

# Insignificant negative correlation 
cor.test(x = weighted_avCOIs, y = fract_highly_related)

# Correlation when 'Cauca_Cauca' is removed: insignificant positive correlation
ind <- names(fract_highly_related) != 'Cauca_Cauca'
cor.test(x = weighted_avCOIs[ind], y = fract_highly_related[ind])
