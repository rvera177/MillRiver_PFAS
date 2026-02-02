
#working with global caravan data. 
#this data comes from caravan-qual lite. Released December 11, 2025
#link to lite dataset: https://zenodo.org/records/17787066 

#"heavy" dataset can be found here: https://github.com/SustainableWaterSystems/Caravan-Qual 

#pull in global caravan data from my github. 
PFOA_Raw <- read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Caravan_PFOA.csv")
PFOS_Raw = read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Caravan_PFOS.csv")


site_info = read_csv("https://raw.githubusercontent.com/rvera177/MillRiver_PFAS/refs/heads/main/data/Caravan_wqms_site_info.csv")

PFOA <- left_join(PFOA_Raw, site_info, by = "wqms_id")
PFOS <- left_join(PFOS_Raw, site_info, by = "wqms_id")