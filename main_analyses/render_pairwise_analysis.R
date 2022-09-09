rmarkdown::render(input = "pairwise_analysis_template.Rmd",
                  output_file = "UCLA_Schizophrenia_Pairwise_Analysis.html",
                  params = list(data_path = "~/data/UCLA_Schizophrenia/",
                                dataset_ID = "UCLA_Schizophrenia",
                                github_dir = "~/github/fMRI_FeaturesDisorders/",
                                univariate_feature_set = "catch22",
                                pairwise_feature_set = "pyspi14",
                                ggseg_atlas = "DK",
                                brain_region_file = "Brain_Region_info.csv",
                                noise_procs = c("AROMA+2P", "AROMA+2P+GMR", "AROMA+2P+DiCER"),
                                noise_proc_for_null = "AROMA+2P+GMR"))

