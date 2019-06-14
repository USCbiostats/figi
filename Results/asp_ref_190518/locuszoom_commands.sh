#!/bin/bash
# locuszoom standalone
# use ASN 1000G_March2012 to calculate LD using axiom_mecc_cfr_ky 
# use EUR to actually use 1000G to calculate LD

results="examples/GxEScanR_aspref_age_sex_pc3_studygxe_N72820_RsqFilter.txt"


# 1:70992951
bin/locuszoom  --snpset NULL --metal ${results} --pop ASN --build hg19 --source 1000G_March2012  --refsnp 1:70992951 --flank 500kb --plotonly --prefix /home/rak/Dropbox/hrc

bin/locuszoom  --snpset NULL --metal ${results} --pop EUR --build hg19 --source 1000G_March2012  --refsnp 1:70992951 --flank 500kb --plotonly --prefix /home/rak/Dropbox/kgp


# 6:31312538 
bin/locuszoom  --snpset NULL --metal ${results} --pop ASN --build hg19 --source 1000G_March2012  --refsnp 6:31312538 --flank 500kb --plotonly --prefix /home/rak/Dropbox/hrc

bin/locuszoom  --snpset NULL --metal ${results} --pop EUR --build hg19 --source 1000G_March2012  --refsnp 6:31312538 --flank 500kb --plotonly --prefix /home/rak/Dropbox/kgp


# 1:183025555
bin/locuszoom  --snpset NULL --metal ${results} --pop ASN --build hg19 --source 1000G_March2012  --refsnp 1:183025555 --flank 500kb --plotonly --prefix /home/rak/Dropbox/hrc

bin/locuszoom  --snpset NULL --metal ${results} --pop EUR --build hg19 --source 1000G_March2012  --refsnp 1:183025555 --flank 500kb --cache None --plotonly --prefix /home/rak/Dropbox/kgp


# 1:55261752
bin/locuszoom  --snpset NULL --metal ${results} --pop ASN --build hg19 --source 1000G_March2012  --refsnp 1:55261752 --flank 500kb --plotonly --prefix /home/rak/Dropbox/hrc

bin/locuszoom  --snpset NULL --metal ${results} --pop EUR --build hg19 --source 1000G_March2012  --refsnp 1:55261752 --flank 500kb --plotonly --prefix /home/rak/Dropbox/kgp





# # 1:70992951
# bin/locuszoom  --snpset NULL --metal ${results} --pop ASN --build hg19 --source 1000G_March2012  --refsnp 1:70992951 --flank 500kb --cache None --plotonly --prefix /home/rak/Dropbox/hrc

# bin/locuszoom  --snpset NULL --metal ${results} --pop EUR --build hg19 --source 1000G_March2012  --refsnp 1:70992951 --flank 500kb --cache None --plotonly --prefix /home/rak/Dropbox/kgp


# # 6:31312538 
# bin/locuszoom  --snpset NULL --metal ${results} --pop ASN --build hg19 --source 1000G_March2012  --refsnp 6:31312538 --flank 500kb --cache None --plotonly --prefix /home/rak/Dropbox/hrc

# bin/locuszoom  --snpset NULL --metal ${results} --pop EUR --build hg19 --source 1000G_March2012  --refsnp 6:31312538 --flank 500kb --cache None --plotonly --prefix /home/rak/Dropbox/kgp


# # 1:183025555
# bin/locuszoom  --snpset NULL --metal ${results} --pop ASN --build hg19 --source 1000G_March2012  --refsnp 1:183025555 --flank 500kb --cache None --plotonly --prefix /home/rak/Dropbox/hrc

# bin/locuszoom  --snpset NULL --metal ${results} --pop EUR --build hg19 --source 1000G_March2012  --refsnp 1:183025555 --flank 500kb --cache None --plotonly --prefix /home/rak/Dropbox/kgp


# # 1:55261752
# bin/locuszoom  --snpset NULL --metal ${results} --pop ASN --build hg19 --source 1000G_March2012  --refsnp 1:55261752 --flank 500kb --cache None --plotonly --prefix /home/rak/Dropbox/hrc

# bin/locuszoom  --snpset NULL --metal ${results} --pop EUR --build hg19 --source 1000G_March2012  --refsnp 1:55261752 --flank 500kb --cache None --plotonly --prefix /home/rak/Dropbox/kgp




6:31312538  rs2854008
HLA-B




1:70992951  rs2651241
not in a region, but close to CTH gene



1:55261752  rs886738
TTC22
