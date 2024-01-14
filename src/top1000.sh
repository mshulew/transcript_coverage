#!/bin/bash
# creates list of top 1000 genes by read counts

sort -n -r -k3 $1 | cut -f 2 | head -1000 > $2
