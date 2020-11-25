## READ MAPPING



## RUN inSTRAIN
## Compare each samples to the different reference genomes
# Species 1
inStrain profile -p 6 -o "$sample_name".species_1.pid96 -c 5 -s 20 -f 0.05 -l 0.96 --min_mapq 2 "$BAM_FILE" PH2015_12U_Oscillatoriales_45_315.fa
# Species 2
inStrain profile -p 6 -o "$sample_name".species_1.pid96 -c 5 -s 20 -f 0.05 -l 0.96 --min_mapq 2 "$BAM_FILE" PH2015_13D_Oscillatoriales_45_19.fa
# Species 3
inStrain profile -p 6 -o "$sample_name".species_1.pid96 -c 5 -s 20 -f 0.05 -l 0.96 --min_mapq 2 "$BAM_FILE" PH2017_22_RUC_O_B_Oscillatoriales_46_93.fa
