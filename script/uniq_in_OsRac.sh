bcftools view -Oz 104-1-90-2_combined.flt.vcf >104-1-90-2.flt.gz
bcftools view -Oz 104-3-1-90-2_combined.flt.vcf >104-3-90-2.flt.gz
bcftools view -Oz Kin-I_combined.flt.vcf >Kin-I_combined.flt.gz
bcftools isec -p 104-1 104-1-90-2.flt.gz Kin-I_combined.flt.gz 
bcftools isec -p 104-3 104-3-90-2.flt.gz Kin-I_combined.flt.gz 

