#!/usr/bin/env bash

## find every gff file
#find: -L option to include symlinks
find -L . \
  -type f \
  -name "*.vcf.gz" \
| sed 's#.vcf.gz#.filtered.vcf#' \
| xargs mk
