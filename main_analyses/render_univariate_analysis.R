rmarkdown::render(input = "univariate_analysis_template.Rmd",
                  output_file = "UCLA_Schizophrenia_Univariate_Analysis.html",
                  params = list(data_path = "D:/Virtual_Machines/Shared_Folder/PhD_work/data/UCLA_Schizophrenia/",
                                dataset_ID = "UCLA_Schizophrenia",
                                github_dir = "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/",
                                univariate_feature_set = "catch22",
                                pairwise_feature_set = "pyspi_19",
                                ggseg_atlas = "DK",
                                noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                noise_proc_for_null = "AROMA+2P+GMR"))


rmarkdown::render(input = "univariate_analysis_template.Rmd",
                  output_file = "ABIDE_ASD_Univariate_Analysis.html",
                  params = list(data_path = "D:/Virtual_Machines/Shared_Folder/PhD_work/data/ABIDE_ASD/",
                                dataset_ID = "ABIDE_ASD",
                                github_dir = "D:/Virtual_Machines/Shared_Folder/github/fMRI_FeaturesDisorders/",
                                univariate_feature_set = "catch22",
                                pairwise_feature_set = "pyspi_19",
                                ggseg_atlas = "HO",
                                noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                noise_proc_for_null = "AROMA+2P+GMR"))