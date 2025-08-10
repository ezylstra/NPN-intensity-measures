# Compare two models based on same multi-year red maple dataset
# 13395 observations of 249 individuals over 8 years

load("output/rema-multiyear-model.RData")
summary(m_rema)
loo1 <- loo(m_rema, cores = 4)
loo1

load("output/rema-gdd-model.RData")
summary(m_rema_gdd)
loo2 <- loo(m_rema_gdd, cores = 4)
loo2

loo_compare(loo1, loo2)
# according to LOOCV, DOY*yr is a WAY better model than GDD

m_rema <- add_criterion(m_rema, "waic")
m_rema_gdd <- add_criterion(m_rema_gdd, "waic")
loo_compare(m_rema, m_rema_gdd, criterion = "waic")
# Same