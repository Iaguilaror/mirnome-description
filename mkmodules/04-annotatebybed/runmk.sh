#!/usr/bin/env bash

## find every bed file
#find: -L option to include symlinks
find -L . \
  -type f \
  -name "*.bed" \
  ! -name "*.filtered.bed" \
| sed 's#.bed#.variants.bed#' \
| xargs mk
