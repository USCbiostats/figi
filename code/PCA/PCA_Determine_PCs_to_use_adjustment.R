#=============================================================================#
#
# FIGI PCA
# Selecting PCs for adjustment
#
#=============================================================================#



eigenval <- fread("~/data/PrincipalComponents/FIGI_GxESet_181120.eigenval")

eigenval$PC <- seq(1,20)
eigenvalTotal <- sum(eigenval$V1)

pct <- round((eigenval$V1 / eigenvalTotal) * 100, 1)
cumpct <- round(cumsum(pct), 1)

eigenval$PC_pct <- pct
eigenval$PC_pct_cumsum <- cumpct

# ggplot likes long format, remember


ggplot(eigenval, aes(PC, V1)) +
  geom_point() + geom_line() + theme_bw() +
  labs(x = "Principal Component", y = "Eigenvalue")


df <- eigenval %>% 
  rename(EigenvalPct = PC_pct,
         EigenvalPctCum = PC_pct_cumsum) %>%
  dplyr::select(-V1) %>% 
  gather(key = group, value, -PC)




ggplot(df, aes(PC, value) ) +
  geom_point() +
  geom_line(aes(linetype = group)) + 
  theme_bw() + 
  theme(legend.position = 'bottom', legend.text=element_text(size=14)) + 
  labs(x = "Principal Component", y = "Eigenval Proportion") + 
  scale_linetype_discrete(name = "", labels = c("Proportion", "Cumulative Proportion")) + 
  geom_text(aes(label=ifelse(value > 40, as.character(value), ''), vjust = -3))



## broken stick??
p = 20
z <- rep(1, p) / seq(1, p) ; z
g <- cumsum(rev(z)) / p; g

df1 <- eigenval %>% 
  mutate(group = "Proportion of Var") %>% 
  dplyr::select(group, PC, PC_pct)

df2 <- data.frame(group = "Broken-Stick", 
                  PC = seq(1,20), 
                  PC_pct = rev(100*g))

df3 <- rbind(df1, df2)	


ggplot(df3, aes(PC, PC_pct) ) +
  geom_point() +
  geom_line(aes(linetype = group)) + 
  theme_bw() + 
  theme(legend.position = 'bottom', legend.text=element_text(size=14))
# labs(x = "Principal Component", y = "Eigenval Proportion") + 
# scale_linetype_discrete(name = "", labels = c("Proportion", "Cumulative Proportion")) +
# geom_text(aes(label=ifelse(value > 40, as.character(value), ''), vjust = -3))



# avg eigenvalu?
#Kaiser-Guttman test
# /* Average Root (Kaiser 1960; Guttman 1954; J. E. Jackson, p. 47) */
# 	mean = mean(Proportion);
# keepAvg = loc( Proportion >= mean )[<>];
# 
# /* Scaled Average Root (Joliffe 1972; J. E. Jackson, p. 47-48) */
# 	keepScaled = loc( Proportion >= 0.7*mean )[<>];
# print keepAvg keepScaled;
# 
# 
# 
# avg_eigen <- sum(eigenval$V1)/20

# might be as simple as the mean proportions

zz <- mean(eigenval$PC_pct)

ggplot(eigenval, aes(PC, PC_pct)) +
  geom_point() + geom_line() + theme_bw() +
  labs(x = "Principal Component", y = "Eigenvalue") + 
  geom_hline(yintercept = zz, color = 'red')


