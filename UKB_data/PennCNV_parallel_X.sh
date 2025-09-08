#!/bin/bash


#call CNVs from batches of 100,000
#1
perl detect_cnv.pl -test -chrx --sexfile /data/UKB/sophie/sexfile.txt -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chrX/chrX.pfb -list /data2/smatthews/UKB/chrX/xaa -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/log1.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/calls1.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
perl detect_cnv.pl -test -chrx --sexfile /data/UKB/sophie/sexfile.txt -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chrX/chrX.pfb -list /data2/smatthews/UKB/chrX/xab -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/log2.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/calls2.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
perl detect_cnv.pl -test -chrx --sexfile /data/UKB/sophie/sexfile.txt -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chrX/chrX.pfb -list /data2/smatthews/UKB/chrX/xac -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/log3.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/calls3.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
perl detect_cnv.pl -test -chrx --sexfile /data/UKB/sophie/sexfile.txt -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chrX/chrX.pfb -list /data2/smatthews/UKB/chrX/xad -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/log4.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/calls4.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
perl detect_cnv.pl -test -chrx --sexfile /data/UKB/sophie/sexfile.txt -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chrX/chrX.pfb -list /data2/smatthews/UKB/chrX/xae -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/log5.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chrX/calls5.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
