tail -1 $1 | sed 's/#//g' | sed 's/\t/"\t"/g' | sed 's/$/"/g' | sed 's/^/"/g' > newHeader
