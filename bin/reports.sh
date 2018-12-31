#!/bin/bash

cd ${HOME}/phd/research/articletwo/ || exit 1

rsync -avu utu:articletwo/*.{html,pdf} reports/ #&& open reports/rrbiome.html
