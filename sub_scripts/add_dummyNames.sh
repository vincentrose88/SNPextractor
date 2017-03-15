#!/bin/bash
sed -r 's/([0-9]+)\t([0-9]+)\t\./\1\t\2\tchr\1\:\2/g' $1/genoFile.tmp > $1/genoFile.noHead
