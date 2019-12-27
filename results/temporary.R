# asp_ref
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'asp_ref',
  is_exposure_categorical = T,
  energy_adj = F
  ), output_file = "~/Dropbox/FIGI/Results/asp_ref/asp_ref_post_harmonization.html")

# aspirin
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'asp_ref',
  is_exposure_categorical = T,
  energy_adj = F
), output_file = "~/Dropbox/FIGI/Results/asp_ref/aspirin_post_harmonization.html")


# alcoholc_moderate
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'alcoholc_moderate',
  is_exposure_categorical = T,
  energy_adj = T
), output_file = "~/Dropbox/FIGI/Results/alcoholc_moderate/alcoholc_moderate_postharmonization.html")


# alcoholc_heavy
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'alcoholc_heavy',
  is_exposure_categorical = T,
  energy_adj = T
), output_file = "~/Dropbox/FIGI/Results/alcoholc_heavy/alcoholc_heavy_postharmonization.html")



# hrt_ref_pm
rmarkdown::render("~/git/FIGI_code/results/results_report_parent.Rmd", params = list(
  exposure = 'hrt_ref_pm',
  is_exposure_categorical = T,
  energy_adj = F
), output_file = "~/Dropbox/FIGI/Results/hrt_ref_pm/hrt_ref_pm_postharmonization.html")





