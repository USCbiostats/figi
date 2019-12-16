#=============================================================================#
# Diagnostic plots
# Use with other results script
#=============================================================================#
calculate_pval <- function(data, statistic, df) {
  data$P <- pchisq(data[,statistic], df = df, lower.tail = F)
  data
}


gxe_twostep <- format_twostep_data(dat = gxe_clump_chiSqG, 'chiSqG', 5, 0.05)




# chiSqG vs chiSqGxE
## categorize chiSqG - not sig, suggestive (1e-6), sig (5e-8)
results <- gxe_clump_chiSqG %>% 
  mutate(chiSqG_logp = -log10(pchisq(chiSqG, df = 1, lower.tail = F)),
         chiSqGxE_logp = -log10(pchisq(chiSqGxE, df = 1, lower.tail = F))) %>% 
  mutate(chiSqG_cat = cut(chiSqG_logp, breaks = c(-Inf, 0.125, 0.3, 0.6, -log10(0.05), -log10(1e-6), -log10(5e-8), Inf), labels = c(0, 1, 2, 3, 4, 5, 6)))
table(results$chiSqG_cat)


plot(results$chiSqG_logp, results$chiSqGxE_logp)
abline(v = 7.30103, col = 'red')


results <- gxe_clump_chiSqG %>% 
  mutate(chiSqG_logp = -log10(pchisq(chiSqG, df = 1, lower.tail = F)),
         chiSqGxE_logp = -log10(pchisq(chiSqGxE, df = 1, lower.tail = F))) %>% 
  mutate(chiSqG_cat = cut(chiSqG_logp, breaks = c(-Inf, 2.30103, 0.3, 0.6, -log10(0.05), -log10(1e-6), -log10(5e-8), Inf), labels = c(0, 1, 2, 3, 4, 5, 6)))
table(results$chiSqG_cat)


ggplot(results, aes(chiSqG_cat, chiSqGxE_logp)) + 
  geom_boxplot()

wtf <- arrange(gxe_clump_chiSqG, desc(chiSqGxE))
head(wtf)

wtfff <- filter(results, chiSqG_cat == 6)









# chiSqG vs chiSqGxE - after LD clumping of chiSqGxE statistics
## categorize chiSqG - not sig, suggestive (1e-6), sig (5e-8)
results <- gxe_clump_chiSqGxE %>% 
  mutate(chiSqG_logp = -log10(pchisq(chiSqG, df = 1, lower.tail = F)),
         chiSqGxE_logp = -log10(pchisq(chiSqGxE, df = 1, lower.tail = F))) %>% 
  mutate(chiSqG_cat = cut(chiSqG_logp, breaks = c(-Inf, 0.125, 0.3, 0.6, -log10(0.05), -log10(1e-6), -log10(5e-8), Inf), labels = c(0, 1, 2, 3, 4, 5, 6)))
table(results$chiSqG_cat)

ggplot(results, aes(chiSqG_cat, chiSqGxE_logp)) + 
  geom_boxplot()

wtf <- arrange(gxe_clump_chiSqG, desc(chiSqGxE))
head(wtf)

wtfff <- filter(results, chiSqG_cat == 6)






# chiSqG vs chiSqGxE - after LD clumping of chiSqGxE statistics
## categorize chiSqG - not sig, suggestive (1e-6), sig (5e-8)
results <- gxe_clump_chiSqGxE %>% 
  mutate(chiSqGE_logp = -log10(pchisq(chiSqGE, df = 1, lower.tail = F)),
         chiSqGxE_logp = -log10(pchisq(chiSqGxE, df = 1, lower.tail = F))) %>% 
  mutate(chiSqGE_cat = cut(chiSqGE_logp, breaks = c(-Inf, 0.125, 0.3, 0.6, -log10(0.05), -log10(1e-6), -log10(5e-8), Inf), labels = c(0, 1, 2, 3, 4, 5, 6)))
table(results$chiSqGE_cat)

ggplot(results, aes(chiSqGE_cat, chiSqGxE_logp)) + 
  geom_boxplot()

