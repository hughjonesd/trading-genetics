
# Trading genetics

This is the git repository for the paper "Trading social status for genetics 
in marriage markets: evidence from Great Britain and Norway" by Abdel Abdellaoui, Oana
Borcan, Pierre-André Chiappori, David Hugh-Jones, Fartein Ask Torvik and Eivind Ystrøm.

The current version of the paper is at https://github.com/hughjonesd/trading-genetics/blob/master/trading-genetics.pdf.

You can also look at previous versions, to check for specification searches,
p-hacking and the Garden of Forking Paths. (But what if we hid them from github? 
It's a wicked world out there....)

Data is from UK Biobank and is not publicly available, but you can request it
from them. We can help with all other data - just get in touch.

Most of the code for the analyses is either in `_drake.R` or in `trading-genetics.Rmd`.

To reproduce:

* Make sure you have the required data.
* Clone this repository.
* You'll also need the https://github.com/hughjonesd/import-ukbb-data repository.
* Edit paths in `import-ukbb-data.R`.
* Start R from within this project directory. This should install the R package
  "renv".
* Run `renv::restore()` to install dependencies.
* Run `Rscript -e make-plan-job.R` from the command line to build the
  paper.
