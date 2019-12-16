#-----------------------------------------------------------------------------#
# Regression functions ------
#
#-----------------------------------------------------------------------------#

# plot stratified analysis Gxe

#' create_glm_stratified_plot
#'
#' Uses package 'effects' to create a stratified plot of GxE interactions. Incidentally, whenever you specify interactions in glm, ALWAYS USE G*E (otherwise, specify on 'gxe' flag)
#'
#' @param model GLM output with higher order (interaction) term in the format G*E
#' @param G Genotype dosage variable (string)
#' @param E Exposure variable (string)
#' @param gxe T/F, just to tell function if GxE or ExG (sorry)
#'
#' @return A stratified GxE plot (NOT exported as *.png)
#' @export
#'
#' @examples create_glm_stratified_plot(model1, "X5.12345678", "exposure", gxe = T)
create_glm_stratified_plot <- function(model, G, E, gxe) {

  # how you define list names right off the bat? grr
  xlevel_list <- list(c(0,1,2))
  if(gxe == T) {
    names(xlevel_list) <- G
  } else {
    names(xlevel_list) <- E # bad but let's get on with it
  }

  model_eff <- effect(paste(G,E, sep = "*"),
                      model,
                      xlevels=xlevel_list, # might have to add more here, depending on coding of E
                      se=TRUE,
                      confidence.level=.95,
                      typical=mean)

  model_eff <- as.data.frame(model_eff) # factors should be defined before fitting model

  ggplot(data=model_eff, aes(x = !! sym(G), y = fit, group = !! sym(E))) +
    # coord_cartesian(ylim = c(0.25,.5)) +
    geom_line(size=2, aes(color=!! sym(E))) +
    ylab("Model Fit")+
    xlab(paste0(G, " Dosage")) +
    ggtitle(paste0("GxE Stratified Analysis - ", G, "*", E)) +
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    scale_fill_grey()

}

# model1 <- glm(outcome ~ X5.40252294*aspirin+age_ref_imp+sex+study_gxe+PC1+PC2+PC3, data = posthoc_df, family = 'binomial')
# create_glm_stratified_plot(model1, "X5.40252294", "aspirin", gxe = T)
#
# model1 <- glm(outcome ~ aspirin*X5.40252294+age_ref_imp+sex+study_gxe+PC1+PC2+PC3, data = posthoc_df, family = 'binomial')
# create_glm_stratified_plot(model1, "aspirin", "X5.40252294", gxe = F)
#




# strat_wrapper <- function(data, snp, exposure) {
#
#   # get counts
#   ctrl <- data %>%  filter(type == 0)
#   ta <- data.frame(table(ctrl[, snp], ctrl[, exposure]))
#
#   case <- data %>% filter(type == 1)
#   tb <- data.frame(table(case[, snp], case[, exposure]))
#
#   tf <- bind_cols(ta, tb[, 3, drop = F]) %>%
#     mutate(counts = paste(Freq, Freq1, sep = "/")) %>%
#     dplyr::select(-starts_with('Freq')) %>%
#     spread(Var2, counts)
#
#   # get E estimates stratified by genotype (0,1,2)
#   covarsE <- c(covars_full, exposure)
#   form <- as.formula(paste0('type ~ ', paste(covarsE, collapse = "+")))
#
#   # just make sure ORs are correct..
#   b <-  data[, c(covarsE, 'type', snp)] %>%
#     filter(complete.cases(.)) %>%
#     group_by_(snp) %>%
#     do(confint_tidy(glm(form, data = ., family = 'binomial')))
#
#   a <- data[, c(covarsE, 'type', snp)] %>%
#     filter(complete.cases(.)) %>%
#     group_by_(snp) %>%
#     do(tidy(glm(form, data = ., family = 'binomial'))) %>%
#     bind_cols(b[, c("conf.low", "conf.high")]) %>%
#     mutate(OR = paste0(round(exp(estimate), 2), " (", round(exp(conf.low), 2), "-", round(exp(conf.high), 2), ")")) %>%
#     filter(grepl(exposure, term)) %>% dplyr::select(!!snp, OR, term) %>%
#     spread_(snp, value = 'OR')
#
#   af <- rbind(c("baseline", "1.0 (Ref)", "1.0 (Ref)", "1.0 (Ref)"), a)
#
#   # get G estimate (additive) for each E level
#   covarsG <- c(covars_full, snp)
#   form <- as.formula(paste0('type ~ ', paste(covarsG, collapse = "+")))
#
#   b <-  data[, c(covarsG, 'type', exposure)] %>%
#     filter(complete.cases(.)) %>%
#     group_by_(exposure) %>%
#     do(confint_tidy(glm(form, data = ., family = 'binomial')))
#
#   aa <- data[, c(covarsG, 'type', exposure)] %>%
#     filter(complete.cases(.)) %>%
#     group_by_(exposure) %>%
#     do(tidy(glm(form, data = ., family = 'binomial'))) %>%
#     bind_cols(b[, c("conf.low", "conf.high")]) %>%
#     mutate(OR = paste0(round(exp(estimate), 2), " (", round(exp(conf.low), 2), "-", round(exp(conf.high), 2), ")")) %>%
#     filter(grepl(snp, term)) %>% dplyr::select(!!exposure, OR, term)
#
#   final <- cbind(tf, af, aa[,2])
#
#
# }

# logit2prob <- function(logit){
#   odds <- exp(logit)
#   prob <- odds / (1 + odds)
#   return(prob)
# }
