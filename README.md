
# Trading genetics

This is the git repository for the paper "Trading social status for genetics 
in marriage markets: evidence from UK Biobank" by Abdel Abdellaoui, Oana
Borcan, Pierre-André Chiappori and David Hugh-Jones.

The current version of the paper is at https://github.com/hughjonesd/trading-genetics/master/trading-genetics.pdf.

Data is from UK Biobank and is not publicly available, but you can request it
from them. We can help with all other data - just get in touch.

To reproduce:

* make sure you have the required data
* clone this repository
* you'll also need the https://github.com/hughjonesd/import-ukbb-data repository
* edit paths in import-ukbb-data.R
* start R from within this project directory. This should install the R package
  "renv".
* Run `renv::restore()` to install dependencies.
* Run `Rscript -e make-plan-job.R` from the command line to build the
  paper
