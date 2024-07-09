#!/bin/bash
echo >> ./manuscript/references.bib
curl -LH "Accept: application/x-bibtex" doi.org/$1 >> ./manuscript/references.bib
echo >> ./manuscript/references.bib