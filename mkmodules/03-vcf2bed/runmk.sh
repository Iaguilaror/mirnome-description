#!/usr/bin/env bash

## find every gff file
#find: -L option to include symlinks
find -L . \
  -type f \
  -name "*.vcf" \
| sed 's#.vcf#.bed#' \
| xargs mk
