rmarkdown::render(input = "dataset_analysis_template.Rmd",
                  output_file = "UCLA_Schizophrenia_Classification_Analysis.html",
                  params = list(data_path = "~/data/UCLA_Schizophrenia/",
                                dataset_ID = "UCLA_Schizophrenia",
                                github_dir = "~/github/fMRI_FeaturesDisorders/",
                                univariate_feature_set = "catch22",
                                pairwise_feature_set = "pyspi14",
                                ggseg_atlas = "DK",
                                brain_region_file = "Brain_Region_info.csv",
                                noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                noise_proc_for_null = "AROMA+2P+GMR"))


# rmarkdown::render(input = "univariate_analysis_template.Rmd",
#                   output_file = "ABIDE_ASD_Univariate_Analysis.html",
#                   params = list(data_path = "~/data/ABIDE_ASD/",
#                                 dataset_ID = "ABIDE_ASD",
#                                 github_dir = "~/fMRI_FeaturesDisorders/",
#                                 univariate_feature_set = "catch22",
#                                 pairwise_feature_set = "pyspi14",
#                                 ggseg_atlas = "HO",
#                                 brain_region_file = "Harvard_Oxford_cort_prob_2mm_ROI_lookup.csv",
#                                 noise_procs = c("FC1000"),
#                                 noise_proc_for_null = "FC1000"))
