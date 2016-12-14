#!/bin/bash

chr=$1
grep -wFf $(dirname $0)/$chr.snps <(zcat $(dirname $0)/imputed/$chr.vcf.gz) > $(dirname $0)/$chr.snps.on.vcf
