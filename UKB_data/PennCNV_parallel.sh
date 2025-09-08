#!/bin/bash


#call CNVs from batches of 100,000
#1
perl detect_cnv.pl -test -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chr2/chr2.pfb -list /data2/smatthews/UKB/chr2/xaa -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/log1.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/calls1.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
perl detect_cnv.pl -test -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chr2/chr2.pfb -list /data2/smatthews/UKB/chr2/xab -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/log2.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/calls2.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
perl detect_cnv.pl -test -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chr2/chr2.pfb -list /data2/smatthews/UKB/chr2/xac -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/log3.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/calls3.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
perl detect_cnv.pl -test -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chr2/chr2.pfb -list /data2/smatthews/UKB/chr2/xad -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/log4.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/calls4.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
perl detect_cnv.pl -test -hmm affy/libgw6/affygw6.hmm -pfb /data2/smatthews/UKB/PennCNV_outputs/chr2/chr2.pfb -list /data2/smatthews/UKB/chr2/xae -log ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/log5.log -out ~/../../data2/smatthews/UKB/PennCNV_outputs/chr2/calls5.rawcnv -gcmodel affy/libgw6/affygw6.hg19.gcmodel
