#!/bin/bash

sed $'$!N;s/\\\n/\t/' "$@" 2> /dev/null


# alternative
#awk '{printf /^>/ ? "\n"$0"\t" : $0}' "$@" | tail -n+2
#printf "\n" ""
