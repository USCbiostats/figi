#Some test data
dat <- data.frame(x=runif(10),y=runif(10),
                  grp = rep(LETTERS[1:5],each = 2),stringsAsFactors = TRUE)

#Create a custom color scale
library(RColorBrewer)
myColors <- brewer.pal(5,"Set3")
# names(myColors) <- levels(dat$grp)
colScale <- scale_colour_manual(values = myColors)
colScale <- scale_fill_manual(values = c("black", "#cc0000", "#ff9933", "#00cc00", "#99ccff"))


getPlots <- function(data, group, exposure, outcome) {
  
  # quote variables
  group_quo <- enquo(group)
  exposure_quo <- enquo(exposure)
  outcome_quo <- enquo(outcome)
  
  # create a temporary grouped dataset
  cov <- data %>%
    dplyr::select(!! outcome_quo, !! group_quo, !! exposure_quo) %>%
    group_by(!! group_quo)
  
  p <- ggplot(data = cov, aes(x = !! outcome_quo)) +
    geom_bar(aes(fill = !! exposure_quo), position = 'fill') +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    colScale
    # scale_fill_manual(values = c( "black", "cyan", "red"))
}

test <- getPlots(data = gxe_set, exposure = exposure_fac, group = study_gxe, outcome = outc_fac)
test


gxe_set <- figi %>% 
  filter(drop == 0 & gxe == 1) %>%  # GxE set N = 102792 as of 3/3/2019
  mutate(outc_num = ifelse(outc == "Case", 1, 0), 
         outc_fac = factor(outc, labels = c("Case", "Control")),
         sex = factor(sex),
         exposure = asp_ref,
         exposure_fac =  fct_relevel(fct_rev(as.factor(folate_totqc2)), ""),
         exposure_num = ifelse(exposure == "", NA, ifelse(exposure == "No", 0, 1)))

table(gxe_set$asp_ref)
table(gxe_set$exposure_fac)
