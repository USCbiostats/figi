#!/bin/bash

# ---- asp_ref ---- #
Rscript create_sig_hit_dataframe.R    asp_ref FIGI_v2.3_gxeset_asp_ref_basic_covars_gxescan age_ref_imp sex study_gxe pc1 pc2 pc3

# ---- aspirin ---- #
Rscript create_sig_hit_dataframe.R    aspirin FIGI_v2.3_gxeset_aspirin_basic_covars_gxescan age_ref_imp sex study_gxe pc1 pc2 pc3

# ---- nsaids ---- #
Rscript create_sig_hit_dataframe.R    nsaids FIGI_v2.3_gxeset_nsaids_basic_covars_gxescan age_ref_imp sex study_gxe pc1 pc2 pc3

# ---- diab ---- #
Rscript create_sig_hit_dataframe.R    diab FIGI_v2.3_gxeset_diab_basic_covars_gxescan age_ref_imp sex study_gxe pc1 pc2 pc3

# ---- calcium_totqc2 ---- #
Rscript create_sig_hit_dataframe.R    calcium_totqc2 FIGI_v2.3_gxeset_calcium_totqc2_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3

# ---- calcium_dietqc2 ---- #
Rscript create_sig_hit_dataframe.R    calcium_dietqc2 FIGI_v2.3_gxeset_calcium_dietqc2_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3

# ---- folate_totqc2 ---- #
Rscript create_sig_hit_dataframe.R    folate_totqc2 FIGI_v2.3_gxeset_folate_totqc2_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3

# ---- folate_dietqc2 ---- #
Rscript create_sig_hit_dataframe.R    folate_dietqc2 FIGI_v2.3_gxeset_folate_dietqc2_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3

# ---- alcoholc_moderate ---- #
Rscript create_sig_hit_dataframe.R   alcoholc_moderate FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3

# ---- alcoholc_heavy ---- #
Rscript create_sig_hit_dataframe.R   alcoholc_heavy FIGI_v2.3_gxeset_alcoholc_heavy_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3

# ---- alcoholc_heavy_by_moderate ---- #
Rscript create_sig_hit_dataframe.R    alcoholc_heavy_vs_moderate FIGI_v2.3_gxeset_alcoholc_heavy_vs_moderate_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3

# ---- hrt_ref_pm ---- #
Rscript create_sig_hit_dataframe.R   hrt_ref_pm FIGI_v2.3_gxeset_hrt_ref_pm_basic_covars_gxescan age_ref_imp sex study_gxe pc1 pc2 pc3

# ---- eo_ref_pm_gxe ---- #
Rscript create_sig_hit_dataframe.R   eo_ref_pm_gxe FIGI_v2.3_gxeset_eo_ref_pm_gxe_basic_covars_gxescan age_ref_imp sex study_gxe pc1 pc2 pc3

# ---- ep_ref_pm_gxe ---- #
Rscript create_sig_hit_dataframe.R   ep_ref_pm_gxe FIGI_v2.3_gxeset_ep_ref_pm_gxe_basic_covars_gxescan age_ref_imp sex study_gxe pc1 pc2 pc3

# ---- redmeatqc2 ---- #
Rscript create_sig_hit_dataframe.R    redmeatqc2 FIGI_v2.3_gxeset_redmeatqc2_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3

# ---- procmeatqc2 ---- #
Rscript create_sig_hit_dataframe.R    procmeatqc2 FIGI_v2.3_gxeset_procmeatqc2_basic_covars_gxescan age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
