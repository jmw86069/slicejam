
# How to review slicejam output in R

This document describes two approaches to review the output
results after running the Rmarkdown analysis.
If the option `SAVE_RDATA=1` was used, there should be
a file with extension `.RData` in the output directory.

Loading data from cache files is usually much faster than
the RData file, because it only loads the promise of data,
not the data itself, see `base::lazyLoad()` help docs.
However the cache/ directory is usually substantially
larger than the RData file.


## If slicejam analysis included SAVE_RDATA=1 load the .RData file

Two ways to load the file, depending whether you want the data
loaded into your current active R workspace, or into a
specific environment of your workspace. You would load into
a specific environment of your workspace when you have data
in your R workspace already, and do not want the `.RData`
contents to overwrite the workspace data.

### Load RData into global workspace (default)

```R
load("file.RData")

# confirm the featureCounts SummarizedExperiment object is available
class(fc_se)
```

### Load RData into a specific environment

```R
# create a new environment
sliceenv <- new.env();

# load data into that environment
load("file.RData", envir=sliceenv)

# use data in the environment by using a prefix
class(sliceenv$fc_se)
# > SummarizedExperiment
```


## How to load Rmarkdown cache data

When the `RData` file is not available, you can load the cached data
from the Rmarkdown `render()` step, which will give you the exact
data used during the analysis and can be helpful for debugging,
or for further review of the contents.

* Navigate to the analysis output folder, then to the "cache/" sub-folder.
* Open R terminal.

Run these commands to load the cache files in the order they were created
during the Rmarkdown steps:

```R
# obtain the full list of cache files
fi1 <- file.info(list.files(pattern="rdx"))
# sort files in order they were created
fo1 <- rownames(fi1[order(fi1$ctime),])
# remove file extension
fo2 <- gsub(".rdx$", "", fo1);

# load files in order, printing each filename onscreen
for (fo in fo2) {
   jamba::printDebug(fo);
   lazyLoad(fo)
}
# print the objects in the working environment
print(ls())
```

At this point you should be able to re-run any R code chunk from
the Rmarkdown file, making modifications where needed.

### Alternate: Load cache files into specific environment

The cache data can be loaded into a specific R environment.

```R
# obtain the full list of cache files
fi1 <- file.info(list.files(pattern="rdx"))
# sort files in order they were created
fo1 <- rownames(fi1[order(fi1$ctime),])
# remove file extension
fo2 <- gsub(".rdx$", "", fo1);

# create new environment
sliceenv <- new.env()

# load files in order, printing each filename onscreen
for (fo in fo2) {
   jamba::printDebug(fo);
   lazyLoad(fo, envir=sliceenv)
}
# print the objects in the specific environment
print(ls(envir=sliceenv))
```
