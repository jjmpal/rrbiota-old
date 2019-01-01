#!/bin/bash

cd "${HOME}/phd/research/articletwo/reports" || exit 1
HOSTPATH="utu:articletwo"
FILE="rrbiome.html"
BACKUP="rrbiome_$(date +%Y%m%d).html"

RET=$(rsync -aEim $HOSTPATH/$FILE $FILE)

if [ $? -eq 0 -a -n "${RET}" ]; then
	cp $FILE $BACKUP
	open $FILE
fi
