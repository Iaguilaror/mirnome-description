#!/usr/bin/env bash

## find every gff file
#find: -L option to include symlinks
find -L . \
  -type f \
  -name "*.gff3" \
| sed 's#.gff3#.CREATE_BED#' \
| xargs mk
