#!/bin/bash
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --output=out.txt #output file
#SBATCH --error=err.txt #error file

#lugh version

cd /data2/smatthews/UKB/chrX/10000_chunks

#cut log into files of 10,000 individuals
awk '
    {
        for( i=c=1; i<=50; i++ ){
            o = $(c++)
            for( n=1; n<10000; n++ ) o = o OFS $(c++)
            print o > "ukb22431_cX_" i ".txt"
        }
    }
' ../ukb22431_cX_b0_v2.txt

#cut baf into files of 10,000 individuals
awk '
      {
          for( i=c=1; i<=50; i++ ){
              o = $(c++)
             for( n=1; n<10000; n++ ) o = o OFS $(c++)
              print o > "ukb22437_cX_" i ".txt"
          }
      }
  ' ../ukb22437_cX_b0_v2.txt
