
# How to review cached data from slicejam_analysis.Rmd

* Navigate to the analysis output folder, then the "cache/" sub-folder.
* Open R terminal.

Run these commands to load the cache files in the order they were created.

```R
fi1 <- file.info(list.files(pattern="rdx"))
fo1 <- rownames(fi1[order(fi1$ctime),])
fo2 <- gsub(".rdx$", "", fo1);
fo2;
for (fo in fo2) {
   jamba::printDebug(fo);
   lazyLoad(fo)
}
print(ls())
```

At this point you should be able to re-run any R code chunk from
the Rmarkdown file, making modifications where needed.
