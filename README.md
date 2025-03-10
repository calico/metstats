# metstats
Functions for statistical analysis of metabolomics data

Calico Life Sciences, LLC

[Demo analysis of Metabolomics Workbench ST003519 published in quarto pub](https://delfarah.quarto.pub/metstats-demo---metabolomics-workbench-st003519-dd5a/)

# Analysis workflow

1. Convert mass spec raw files to mzML using [ProteWizard MSConvert](https://proteowizard.sourceforge.io/tools/msconvert.html) software.
e.g. Metabolomics Workbench [ST003519](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&DataMode=AllData&StudyID=ST003519&StudyType=MS&ResultType=1#DataTabs) for demo, raw files are already mzML.

2. Run dataset through [MAVEN peakdetector](https://github.com/eugenemel/maven/releases) to generate ```mzrollDB``` file for peak identification and quantification.

3. Open mzrollDB in [MAVEN 2 software](https://www.mdpi.com/2218-1989/12/8/684) and annotate peaks.

4. Create a project folder and put the following files in the folder:
- Annotated ```mzrollDB``` file
- qmd file ```metstats/vignettes/ST003519_demo.Qmd``` for the demo

5. In terminal go to the project folder and run the following to set up a new Quarto website project:
```quarto create project my-website
```
- Choose a folder name for the project/website.

6. Render the qmd to see results in localhost.
7. To publish results to Quarto pub, run the following in terminal:
```quarto publish quarto-pub ST003519_demo.qmd
```

Example result page on Quarto pub for the [ST003519 demo](https://delfarah.quarto.pub/metstats-demo---metabolomics-workbench-st003519-dd5a/)

# Network Analysis

To run Cytoscape and generate metabolic networks:

1. Download latest version of [Cytoscape](https://cytoscape.org/download.html)

2. Install [RCy3 2.5.1](https://github.com/delfarahalireza/RCy3) Cytoscape-R API
- Newer versions of RCy3 are problematic.
- Recommended method of installation:

```r
# 1. Clone the RCy3 repository (if you haven't already)
git clone git@github.com:delfarahalireza/RCy3.git

# 2. Install the package from source
install.packages("cloned directory of RCy3", repos = NULL, type = "source")

# 3. Restart R session

# 4. Check that version of RCy3 is 2.5.1 after loading library(Rcy3)
```

3. Open ```metstats/inst/extdata/Networks.cys``` file

4. In R, run the ```network plots```, ```pathways color```, ```pathways shape```, and ```pathways``` chunks similar to ```metstats/vignettes/ST003519_demo.Qmd```




