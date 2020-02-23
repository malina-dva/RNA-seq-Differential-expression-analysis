## RNA-seq:Differential expression analysis

"RNA-seq:Differential expression analysis" is a continuation of another R tutorial called "RNA-seq: Data normalization and clustering" (uploaded on gitHub https://github.com/malina-dva/RNA-seq-Data-normalization-and-clustering?fbclid=IwAR0op-RofdeX91wzPSPeW2OSUKYXWv59oaJ7Btz__Noibrnzn6RBYRdIw-E)

You will have to repeat the first four sections from the "RNA-seq: Data normalization and clustering" in order to be able to perform the steps from the current tutorial. 


If you are new to R and Rstudio, please have a look at the video for how to set up your working directory and install packages.

[![Getting started with RNA-seq](http://img.youtube.com/vi/kR_iHVau8GI/0.jpg)](https://www.youtube.com/watch?v=kR_iHVau8GI "Getting started with RNA-seq" )

To maintain the flow of the current tutorial I have included the commands you need from the "RNA-seq: Data normalization and clustering" down below. 

### 1. Installing and loading the DESeq2 package

First you need to load the package DESeq2, which is one of the **MAIN** packages that will be used in this tutorial.
if you have never used DESeq2, you will have to install and load it beforehand.
The commands for installation and loading you can find here:

https://bioconductor.org/packages/release/bioc/html/DESeq2.html

By the way this is how you can install any package if you have R version >=3.6

Alternatively you can install DESEq2 by typing: 

install.packages("DESeq2", dependencies = TRUE)

Regardless of the way you decide to install DESeq2, it will take some time for all the dependences and DESeq2 itself to be installed. So be patient. 

To just load/re-load the package if you have it already installed use the line below. 



```r
require("DESeq2")
```

```
## Loading required package: DESeq2
```

```
## Loading required package: S4Vectors
```

```
## Loading required package: stats4
```



You can also use an alternative command below to load the package, once installed:

library("DESEq2")

### 2. Setting working directory

Once we have DESEq2 installed and loaded, then we can set the working directory where you will have to place all the files you will use as an input.
This is also the directory where all of the output files will be stored automatically once the directory is set.
Working directory can be any folder on your computer.

To set working directory in R studio, you use the function

setwd("")

You have to place the location of your working directory (folder) Within the double quotes of this function.

Let's create a new **folder** (working directory) on your desktop and call it 

Tutorial_RNA_seq

To retrieve the full path to the location of the newly created folder, simply copy the whole folder and paste within the double quotes.
This is what I get when I do this:
"file:///C:/Users/Malina/Desktop/Tutorial_RNA_seq"
You will have to change the forward slashes to backward double slashes and you will have to remove everything before C: in order R to be able to use the path.**NO NEED to change the type of slashes if you are using Linux and macOS machines.**

Place the full path within the double quotes in the command below. 



```r
setwd("C:\\Users\\Malina\\Desktop\\Tutorial_RNA_seq")
```

### 3. Downloading and importing the RNA-seq data

**The RNA-seq data** which is going to be analyzed here is from this paper:

https://www.nature.com/articles/ncomms12418 

DNA hydroxymethylation controls cardiomyocyte gene expression in development and hypertrophy. 

The data is from cardiac myocytes (CMs) cells isolated from mice - embryonic day 14,5, neonatal, adult and TAC (transverse aortic constriction). The first three conditions represent developmental stages and the last condition is hypertrophic (disease) state.

**Each of the four conditions has two biological replicates.**

In order to perform the analysis on this data set you will have to download the file with the raw counts (the input file) named Raw_counts_input.txt from github.
To do so go to:

https://github.com/malina-dva/RNA-seq-Differential-expression-analysis

Click on Raw_counts_input.txt file

Then click on Raw.

The file with the Raw_counts_input.txt will open in a separate page.

Right click and then click on "save page as" or "save as" to save the file, preferably on your Desktop. Make sure you save the file with its original name.

Once the file is saved on your Desktop move it to your working directory that we just created and set. (Tutorial_RNA_seq).


It is a good idea to open the input file in excel and check its content before you start to work with it in R.
**However, I dont' recommend keeping input files open in excel while you are using them in R.**

Repeating again! The file with the raw counts has to be placed in your working directory.
To read/load the input file in Rstudio use the command below.

```r
just.raw.counts = read.delim("Raw_counts_input.txt")
```

Let's check first several rows of the input file using the function 'head'.


```r
head(just.raw.counts)
```

```
##        Probe E14.5_R1 E14.5_R2 Neonatal_R1 Neonatal_R2 Adult_R1 Adult_R2 TAC_R1
## 1       Xkr4      229      363         545         417      133       96    280
## 2     Gm1992        0        3           0           0        0        0      1
## 3        Rp1        3        5           4          26       48       73      7
## 4      Sox17      206      195         285         226       51       44     53
## 5 AC129937.1        0        3           0           2        0        0      0
## 6     Mrpl15      597      599         468         480      318      316    389
##   TAC_R2
## 1    219
## 2      0
## 3     10
## 4     57
## 5      0
## 6    347
```


Now let's check the dimensions (total number of rows and columns) of your imported table. 


```r
dim(just.raw.counts)
```

```
## [1] 27195     9
```

There are 27195 genes (rows) and9 columns. The first column (e.g. Probe) holds the names of the genes, the remaining 8 columns hold the gene counts (gene expression info) for the replicates of the four conditions we have. 
We want to actually specify that the column named 'Probe' has the information for the names of each row in the table. 

We can do that when we read the input file by specifying row.names =1, meaning that names of the rows are to be found in the first column.


```r
just.raw.counts = read.delim("Raw_counts_input.txt", row.names = 1)
```

```r
head(just.raw.counts)
```

```
##            E14.5_R1 E14.5_R2 Neonatal_R1 Neonatal_R2 Adult_R1 Adult_R2 TAC_R1
## Xkr4            229      363         545         417      133       96    280
## Gm1992            0        3           0           0        0        0      1
## Rp1               3        5           4          26       48       73      7
## Sox17           206      195         285         226       51       44     53
## AC129937.1        0        3           0           2        0        0      0
## Mrpl15          597      599         468         480      318      316    389
##            TAC_R2
## Xkr4          219
## Gm1992          0
## Rp1            10
## Sox17          57
## AC129937.1      0
## Mrpl15        347
```

After defining row.names = 1, the first column is no longer named Probe as you can see from output above.

```r
dim(just.raw.counts)
```

```
## [1] 27195     8
```
Also the number of columns is 8, not 9 - exactly matching the total number of replicates in the data set.

### 4. Generating the expression set for DESeq2

To do that we need three elements

1. names of the genes, which we already defined with row.names = 1

2. raw counts of the genes for the respective condition and replicate

it seems we already have the first two elements from the 'just.raw.counts' variable

3. metadata - phenotypic data for the expression set
The main purpose of the metadata is to define the relationship between replicates and conditions

The metadata is another text file called meta_data.txt.



Again download this file from github and place it in your working directory. 


To do so go 

https://github.com/malina-dva/RNA-seq-Differential-expression-analysis

Click on meta_data.txt file

Then click on Raw

The file with the meta_data.txt will open in a separate page

Right click and then click on "save page as" or "save as" to save the file, preferably on your Desktop. Keep the original name of the file.

Move meta_data.txt to the Tutorial_RNA_seq folder.


Open it with excel to check its content.

Then let's load it in Rstudio.



```r
meta.data = read.delim(file="meta_data.txt", row.names = 1)
```
Again we specify the first column as row.names.


```r
head(meta.data)
```

```
##             condition
## E14.5_R1        E14.5
## E14.5_R2        E14.5
## Neonatal_R1  Neonatal
## Neonatal_R2  Neonatal
## Adult_R1        Adult
## Adult_R2        Adult
```
The column called 'condition' links the name of each replicate from 'just.raw.counts' to its respective condition.
For example E14.5_R1 and E14.5_R2, in the 'condition' column have the same name e.g E14.5, because they both belong to the E14.5 condition.


Now we can proceed to build DESeq object.

First we need to create DESeq matrix. 

```r
count.data.set <- DESeqDataSetFromMatrix(countData=just.raw.counts, 
                                         colData=meta.data, design= ~ condition) 
```

'countData' argument requires the table with the counts and the names of the genes
'colData' argument requires the metadata file
'design' requires the column of the metadata that holds the information for the relationships between replicates and conditions (e.g the 'condition' column from the metadata in our case)



Then we create DESeq object.

```r
count.data.set.object <- DESeq(count.data.set)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

### 5. Differential expression analysis 
For performing differential expression analysis we don't need to normalize the data beforehand.
DESeq2 uses raw counts to perform the statistical test for identifying genes that are differentially expressed between conditions.

The function that we use to obtain a list with the differentially expressed genes is called "results".

We can use this function by passing Deseq object to it. **We just created one above (count.data.set.object)**.
Here we are going to use "contrast" argument to get the differentially expressed genes (DEG) for the comparison between two **conditions** -  Neonatal compared to E14.5.



```r
res <- results(count.data.set.object, contrast=c("condition","Neonatal","E14.5"))
```

To check the first few lines of the "results" output just run res.


```r
res 
```

```
## log2 fold change (MLE): condition Neonatal vs E14.5 
## Wald test p-value: condition Neonatal vs E14.5 
## DataFrame with 27195 rows and 6 columns
##                      baseMean     log2FoldChange             lfcSE
##                     <numeric>          <numeric>         <numeric>
## Xkr4         264.493584216531  0.935393970854056  0.25066339056516
## Gm1992        0.3956712134287  -2.74785350532486  4.85964654230835
## Rp1          28.6314442544208    2.1380664991517 0.839134353214377
## Sox17        120.722246152093  0.571306786655843 0.274738774018729
## AC129937.1  0.454668627784066 -0.303968583443821  4.75128030921768
## ...                       ...                ...               ...
## AC151712.3                  0                 NA                NA
## AC151712.1 0.0836284315793158  -1.20138696387088  4.99456102608126
## AC151712.4   944.557348992265  0.753466152174321 0.237796422606015
## AC151712.5                  0                 NA                NA
## AC151712.2                  0                 NA                NA
##                           stat               pvalue                padj
##                      <numeric>            <numeric>           <numeric>
## Xkr4           3.7316736550362 0.000190211817781812 0.00126241968144792
## Gm1992      -0.565443079327251    0.571772430033651                  NA
## Rp1           2.54794299740161   0.0108360175433807  0.0378445318580725
## Sox17         2.07945452437994    0.037575593459159   0.101909160872415
## AC129937.1 -0.0639761419367553    0.948989223858412                  NA
## ...                        ...                  ...                 ...
## AC151712.3                  NA                   NA                  NA
## AC151712.1  -0.240539049897943    0.809912394261418                  NA
## AC151712.4    3.16853442922762  0.00153209601516601 0.00756077537381214
## AC151712.5                  NA                   NA                  NA
## AC151712.2                  NA                   NA                  NA
```
There are six columns in the output (in addition to the gene names). The most important for us are log2FoldChange (how many times the gene expression in the Neonatal changes relative to E14.5) and padj (p value adjusted for multiple testing) columns, which shows us the significance of the difference.
If the value of the log2FoldChange is positive, it means the gene was upregulated in Neonatal compared to E14.5, and if it is negative, the gene was downregulated in Neonatal compared to E14.5 respectively.
To get a summary of the output use:

```r
summary(res)
```

```
## 
## out of 20902 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)       : 3082, 15%
## LFC < 0 (down)     : 2526, 12%
## outliers [1]       : 0, 0%
## low counts [2]     : 5584, 27%
## (mean count < 8)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```
We can see that there are in total 3082 and 2526 genes up and downregulated respectively, when the default adjusted p-value cutoff of 0.1 is used.

Additionally DESeq2 did not identify outliers based on the cooksCutoff, because the replicates for the conditions that we are comparing are less than three. 
What cooksCutoff filtering is... 
This is simply a cutoff for excluding genes that have higher variability of counts between replicates than between the compared condition.
For example if we had at least three replicates for each of the compared conditions, this test will check if the variability within the replicates of the Neaonatal or E14.5 is not exceeding the variability between   Neonatal relative to E14.5 conditions. 
If instances like that were found, the pvalue and p.adjust column for the respective genes for which this was the case would have been labeled as NA.
If we have condition with only two replicates, even if the variability of the counts within the replicates is high, there would be no way to determine which value for the counts in the two replicates is likely to be "more true". The third replicate comes handy is such cases.
With only two replicates cooksCutoff test is not performed.

Another test that is performed is excluding genes that have too low counts relative to the mean counts of the experiment.

This test is called independentFiltering. You can see that there were 5584 genes or 27% of all genes that did not pass the minimal counts threshold.
These genes are labeled as NA in the results output. If a gene has zero counts for all replicates, baseMean column will be zero and the remaining columns will be NAs. If a gene is filtered out due to having a low mean normalized count, but is not zero in all replicates, only p.adjust will be labeled as NA.
How performing independentFiltering can improve the results we get? By removing the weakly-expressed genes,  multiple  testing adjustment of the remaining genes improves  and we are able to find  more  genes  to  be  significant  among  those which we keep. (https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf)

The default adjusted p-value cutoff (0.1) is relatively weak threshold for significance. We want to filter the data based on more stringent **adjusted** p value (padj) and to get only the genes that have padj value less than 0.05.
Let's get started.


#### 1. Filter for p adjusted <= 0.05


Let's first remove all genes that have columns with NAs. These will be the genes which didn't pass the independentFiltering. An easy way to do that is using the na.omit function and to apply it to the variable that holds our results.

```r
res = na.omit(res)
```
Let's check if there are any NAs left


```r
head(res)
```

```
## log2 fold change (MLE): condition Neonatal vs E14.5 
## Wald test p-value: condition Neonatal vs E14.5 
## DataFrame with 6 rows and 6 columns
##                baseMean     log2FoldChange             lfcSE               stat
##               <numeric>          <numeric>         <numeric>          <numeric>
## Xkr4   264.493584216531  0.935393970854056  0.25066339056516    3.7316736550362
## Rp1    28.6314442544208    2.1380664991517 0.839134353214377   2.54794299740161
## Sox17  120.722246152093  0.571306786655843 0.274738774018729   2.07945452437994
## Mrpl15   425.9268829942 -0.112065875825264 0.173950414675517 -0.644240348804624
## Lypla1 2125.25733572459 -0.212902093679033 0.138183910889204  -1.54071550232602
## Gm6104 10.7188340040208 -0.958849545023435 0.841281350766104  -1.13974896049968
##                      pvalue                padj
##                   <numeric>           <numeric>
## Xkr4   0.000190211817781812 0.00126241968144792
## Rp1      0.0108360175433807  0.0378445318580725
## Sox17     0.037575593459159   0.101909160872415
## Mrpl15    0.519419590155383   0.673193102800588
## Lypla1    0.123386042168149   0.248851533104898
## Gm6104     0.25439090305859   0.416276023186783
```

Good for us, it seems, at least when looking at the first few rows, that there are no NAs left in the output.


Now let's select only the significant genes based on padj value cutoff equals or lower than 0.05. 
We can do that by the following command:


```r
resfiltered = res[res$padj <= 0.05,]
head(resfiltered)
```

```
## log2 fold change (MLE): condition Neonatal vs E14.5 
## Wald test p-value: condition Neonatal vs E14.5 
## DataFrame with 6 rows and 6 columns
##                       baseMean    log2FoldChange             lfcSE
##                      <numeric>         <numeric>         <numeric>
## Xkr4          264.493584216531 0.935393970854056  0.25066339056516
## Rp1           28.6314442544208   2.1380664991517 0.839134353214377
## Atp6v1h       1480.46643710253  0.71491744184966 0.151692536580221
## Pcmtd1        2015.43715746467 0.868627278470412 0.144607239462168
## Adhfe1         610.25231063051  1.70080201575893 0.199941239271289
## 2610203C22Rik  204.18358523198  2.86448724741401 0.378475860286967
##                           stat               pvalue                 padj
##                      <numeric>            <numeric>            <numeric>
## Xkr4           3.7316736550362 0.000190211817781812  0.00126241968144792
## Rp1           2.54794299740161   0.0108360175433807   0.0378445318580725
## Atp6v1h       4.71293748503956 2.44171123159486e-06 2.76064092221974e-05
## Pcmtd1        6.00680354386864 1.89216555921233e-09 3.96500575048077e-08
## Adhfe1        8.50650932222747 1.79248198281982e-17 1.18350168158767e-15
## 2610203C22Rik  7.5684807090262 3.77614382875013e-14  1.6112248236433e-12
```
The commands asks to get all the rows(genes), where padj column has values equal or lower than 0.05.
Now we want to check the summary of the filtered results for padj and NAs.

```r
summary(resfiltered)
```

```
## 
## out of 4648 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)       : 2615, 56%
## LFC < 0 (down)     : 2033, 44%
## outliers [1]       : 0, 0%
## low counts [2]     : 0, 0%
## (mean count < 8)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

You can see that after all the steps of filtering we got 2615 up and 2033 down regulated genes respectively, using threshold of padj 0.05.
No NAs are present, as per the independentFiltering output.

But wait! In the summary the adjusted p-value is < 0.1. It is really? Haven't we already set a different threshold?
We can check if there are genes with padj value higher than 0.05 with the following line of code:


```r
sum(resfiltered$padj > 0.05)
```

```
## [1] 0
```

With this code we get the number of rows that have padj value > 0.05. it is clear that there are none, which means that our filtering for padj works.

Yet,the adjusted p-value will always be shown as  < 0.1 (the default) in the summary output, unless we specify a diffrent value for "alpha" (padj value) directly when we generate the res variable!
If we were to that the command will be:

res <- results(count.data.set.object, contrast=c("condition","Neonatal","E14.5"), alpha=0.05)


Finally, we want to order the results so that on the top there will be genes with the lowest (most significant padj value).

#### 2. Order the results based on increasing padj value

We use the function "order" and apply it on the values in the resfiltered$padj column and for all rows (genes) of the resfiltered variable.
The output from the reordering we keep in a new variable called res.filtered.ordered


```r
res.filtered.ordered = resfiltered[order(resfiltered$padj),] 
```

Let's check what we got now.


```r
head(res.filtered.ordered)
```

```
## log2 fold change (MLE): condition Neonatal vs E14.5 
## Wald test p-value: condition Neonatal vs E14.5 
## DataFrame with 6 rows and 6 columns
##                baseMean   log2FoldChange             lfcSE             stat
##               <numeric>        <numeric>         <numeric>        <numeric>
## Parp14 906.738982073535 3.81372984953951 0.188852056379392 20.1942722925821
## Mndal  1113.88864099369  3.7151892299422 0.187099384516676 19.8567688479545
## Rsad2  554.291380343165 3.77928939515909 0.199843963942738 18.9112011221013
## Rnf213 2769.26656489804 2.31717682379161 0.125832000934182 18.4148452427744
## Mx2    352.485286822358 4.69237762534611 0.266311395047114  17.619890521455
## Pparg  924.951149646171 3.04779703337712 0.178390873142063 17.0849381456303
##                      pvalue                 padj
##                   <numeric>            <numeric>
## Parp14  1.0993444109865e-90 1.68397576874913e-86
## Mndal  9.63109598055043e-88 7.37645641150357e-84
## Rsad2  9.22318914103121e-80 4.70936037541053e-76
## Rnf213 9.98744112309483e-76 3.82469057808917e-72
## Mx2    1.73340065573481e-69 5.31044624890916e-66
## Pparg  1.92146838252122e-65 4.90550878057667e-62
```

You see that the genes with the lowest padj values are on the top with a "leader", a gene called Parp14 (1.68397576874913e-86). Veryyy low padj value though, very highly significant difference (upregulation based on the logFoldChange, which is a positive value ~3.8) in Neonatal compared to E14.5.

Now we are going to save the filtered and ordered list with the DEG.


```r
write.table(res.filtered.ordered, sep="\t",file="Results_filtered_ordered.txt", row.names=TRUE,col.names=NA,quote=FALSE)
```

#### 3. logFoldChange based filtering and visualization
**This step of analysis is just for demonstration. Generally, cutoff based on FC is not needed, unless the is a very good reason behind it.**

But anyway, let's check what is the distribution of the logFoldChange values for all genes in the results we have.
We will need ggplot2 to do that. By now you should be a PRO in package installation and loading. :P If not google BiocManager and ggplot2 and follow the steps for installation. Then load ggplot2 using the command below:


```r
require(ggplot2)
```

```
## Loading required package: ggplot2
```

We need to convert the output with our DEG (res.filtered.ordered) into dataframe in order to be able to visualize it on ggplot2.


```r
resOrdFilt.data.frame = as.data.frame(res.filtered.ordered)
head(resOrdFilt.data.frame)
```

```
##         baseMean log2FoldChange     lfcSE     stat       pvalue         padj
## Parp14  906.7390       3.813730 0.1888521 20.19427 1.099344e-90 1.683976e-86
## Mndal  1113.8886       3.715189 0.1870994 19.85677 9.631096e-88 7.376456e-84
## Rsad2   554.2914       3.779289 0.1998440 18.91120 9.223189e-80 4.709360e-76
## Rnf213 2769.2666       2.317177 0.1258320 18.41485 9.987441e-76 3.824691e-72
## Mx2     352.4853       4.692378 0.2663114 17.61989 1.733401e-69 5.310446e-66
## Pparg   924.9511       3.047797 0.1783909 17.08494 1.921468e-65 4.905509e-62
```

Now let's plot it. We will plot a **histogram** of the distribution of values in the column with log2FoldChange (e.g. the x in the aesthetics in the ggplot function).


```r
ggplot(resOrdFilt.data.frame, aes(x=log2FoldChange)) +
  geom_histogram()
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](https://github.com/malina-dva/RNA-seq-Differential-expression-analysis/blob/master/unnamed-chunk-26-1.png)<!-- -->
 

We can zoom the x axis a bit by specifying limits=c(-5,5) to the scale_x_continuous argument (meaning the x axis). 
Additionally we can put labels to the positive and negative 5, 1 and 0.5 values of the x axis (e.g specifying that in the breaks argument) and increase the number of histogram bins from the default 30 to 100, to be able to better see where the distribution peaks.

```r
ggplot(resOrdFilt.data.frame, aes(x=log2FoldChange)) +
  geom_histogram(bins=100) +
  scale_x_continuous( name ="log2FoldChange",breaks=c(-5,-1,-0.5,0, 0.5,1,5), limits = c(-5,5))
```

```
## Warning: Removed 26 rows containing non-finite values (stat_bin).
```

![](https://github.com/malina-dva/RNA-seq-Differential-expression-analysis/blob/master/unnamed-chunk-27-1.png)<!-- -->


Based on the histogram it is obvious that most of FC are between 0.5 and 1. Just for demonstration we can set a threshold for the FC to -0.5 and +0.5 for down and up regulated genes respectively to see how many DEG will remain.

To filter for log2FoldChange threshold we use "subset". We apply subset on the resOrdFilt.data.frame. The **conditions** we want to sunset on are:  log2FoldChange column from resOrdFilt.data.frame is either >0.5 or <0.5. To achieve that, we use "|" argument in our command which means **or**. 


```r
resOrdFiltLFC = subset(resOrdFilt.data.frame,log2FoldChange >0.5 | log2FoldChange < -0.5 )
```

To check how many DEG remained after filtering for fold change > or < 0.5:

```r
dim(resOrdFiltLFC)
```

```
## [1] 3985    6
```

Well, we have reduced the number of DEG from 4648 to 3985 by applying filtering for LFC.


```r
ggplot(resOrdFiltLFC, aes(x=log2FoldChange)) +
  geom_histogram(bins=100) +
  scale_x_continuous( name ="log2FoldChange",breaks=c(-5,-1,-0.5,0, 0.5,1,5), limits = c(-5,5))
```

```
## Warning: Removed 26 rows containing non-finite values (stat_bin).
```

![](https://github.com/malina-dva/RNA-seq-Differential-expression-analysis/blob/master/unnamed-chunk-30-1.png)<!-- -->

Let's write the results with the filtered FCs.


```r
write.table(resOrdFiltLFC, sep="\t",file="Results_filtered_ordered_LFC.txt", row.names=TRUE,col.names=NA,quote=FALSE)
```
**As I mentioned in the beginning of the section, filtering on FC is not always recommended or needed. The included steps from above were only for demonstration.**

#### 4. Subsetting and visualizing the DEGs 

We will continue to work with our original list of DEG **prior to foldchange filtering**. These were the results stored in res.filtered.ordered variable. We want to get two tables with either up or down DEG and save them. 


```r
summary(res.filtered.ordered)
```

```
## 
## out of 4648 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)       : 2615, 56%
## LFC < 0 (down)     : 2033, 44%
## outliers [1]       : 0, 0%
## low counts [2]     : 0, 0%
## (mean count < 8)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

To be able to subset our up and down DEG genes we will use the resOrdFilt.data.frame variable which holds the results from res.filtered.ordered.


```r
Sign_genes_up = subset(resOrdFilt.data.frame, log2FoldChange >0)
Sign_genes_down = subset(resOrdFilt.data.frame, log2FoldChange <0)
```

Then save the up and down genes to text files.


```r
write.table(Sign_genes_up, sep="\t",file="Sign_genes_up.txt", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(Sign_genes_down, sep="\t",file="Sign_genes_down.txt", row.names=TRUE,col.names=NA,quote=FALSE)
```

We can directly visualize the up and down genes using "EnhancedVolcano" package.

To download the package please copy and execute the commands from:
http://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html


Now... let's create a volcano plot.


```r
require("EnhancedVolcano")
```

```
## Loading required package: EnhancedVolcano
```

```
## Loading required package: ggrepel
```

For labels on the volcano plot, we are going to use the rownames(resOrdFilt.data.frame) which are in fact the names of the genes in the resOrdFilt.data.frame.
Not a bad idea to have a look at the first few lines of the resOrdFilt.data.frame, right?

```r
head(resOrdFilt.data.frame)
```

```
##         baseMean log2FoldChange     lfcSE     stat       pvalue         padj
## Parp14  906.7390       3.813730 0.1888521 20.19427 1.099344e-90 1.683976e-86
## Mndal  1113.8886       3.715189 0.1870994 19.85677 9.631096e-88 7.376456e-84
## Rsad2   554.2914       3.779289 0.1998440 18.91120 9.223189e-80 4.709360e-76
## Rnf213 2769.2666       2.317177 0.1258320 18.41485 9.987441e-76 3.824691e-72
## Mx2     352.4853       4.692378 0.2663114 17.61989 1.733401e-69 5.310446e-66
## Pparg   924.9511       3.047797 0.1783909 17.08494 1.921468e-65 4.905509e-62
```
The x axis on the volcano plot is going to be log2FoldChange column from resOrdFilt.data.frame. The y axis will be padj column respectively. For the FC cutoff (FoldChange cutoff) let's specify 0.3 (that will refer to +/- 3), and pCutoff (which in our case is actually padj) let's specify 0.05.
Ideally, one would like to change the labels of the y axis and legend to padj instead of p value. You can try to do that as a **homework** by following the steps in this tutorial: https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html (point 4.6)


```r
EnhancedVolcano(resOrdFilt.data.frame,

        lab = rownames(resOrdFilt.data.frame),

        x = "log2FoldChange",
        y = "padj",
        FCcutoff = 0.3,
         pCutoff = 0.05)
```

![](https://github.com/malina-dva/RNA-seq-Differential-expression-analysis/blob/master/unnamed-chunk-37-1.png)<!-- -->

You see that all of our genes from resOrdFilt.data.frame have a difference in FC greater than 0.3 and the padj value is < 0.05, which is what we would expect anyway!

### 6. GO enrichment analysis and visualization

Let's check what biological processes the Up and Down genes are involved in by performing GO enrichment analysis. 
Of course, will have to load the packages that will help us find this out.
Again, please install them as per usual if you don't have them already (Bioconductor, BiocManager and the name of the package). To load them use the commands below.

```r
require("clusterProfiler")
```

```
## Loading required package: clusterProfiler
```

```
## 
```

```
## Registered S3 method overwritten by 'enrichplot':
##   method               from
##   fortify.enrichResult DOSE
```

```
## clusterProfiler v3.12.0  For help: https://guangchuangyu.github.io/software/clusterProfiler
## 
## If you use clusterProfiler in published research, please cite:
## Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
```

```
## 
## Attaching package: 'clusterProfiler'
```

```
## The following object is masked from 'package:DelayedArray':
## 
##     simplify
```

```r
require("org.Mm.eg.db")
```

```
## Loading required package: org.Mm.eg.db
```

```
## Loading required package: AnnotationDbi
```

```
## 
```
Now we are going to create a variable with all unique mouse genes.
We will need that as a gene background for our GO analysis.

```r
all_genes_mouse = unique(just.raw.counts$Probe)
```

To perform GO enrichment we will need a few arguments:

1. The names of the genes we want to perform GO enrichment test on (e.g rownames(Sign_genes_up) or rownames(Sign_genes_down) respectively)

2. OrgDb - org.Mm.eg.db, Genome wide annotation data base for mouse. The enrichGO function will perform the tests using this data base

3. keyType - the type of the annotation of our gene set

4. ont - type of Ontology - can be BP - Biological process, MF - Molecular function,  and CC - cell component subontologies

5. pAdjustMethod - method for the FDR correction, "BH" is Benjamini and Hochberg test

6. universe - background genes (e,g all_genes_mouse variable)

7. pvalueCutoff - cutoff for the pvalue of the GO test

8. qvalueCutoff - cutoff for the corrected for mutiple testing pvalue (e.g the qvalue) of the GO test


```r
GO_up <- enrichGO(gene         = rownames(Sign_genes_up),
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 universe      = all_genes_mouse,
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
```

enrichGO function takes a little bit of time to complete, so please be patient. 

To have a look at the results of the GO enrichment, we can create a dataframe from the GO_up varaible. Then let's check only the first three rows using the head function. You can also opne the file in excel to quickly check its content, but then don't forget to close it.


```r
GO_up_dataframe=as.data.frame(GO_up)
head(GO_up_dataframe, n=3)
```

```
##                    ID                           Description GeneRatio   BgRatio
## GO:1901342 GO:1901342 regulation of vasculature development  102/2237 340/23174
## GO:0003013 GO:0003013            circulatory system process  126/2237 489/23174
## GO:0001667 GO:0001667         ameboidal-type cell migration  109/2237 406/23174
##                  pvalue     p.adjust       qvalue
## GO:1901342 2.614672e-26 1.572987e-22 8.518327e-23
## GO:0003013 2.145316e-25 6.453110e-22 3.494606e-22
## GO:0001667 1.177161e-23 1.912059e-20 1.035453e-20
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  geneID
## GO:1901342                                                                                                                                        Pparg/Xdh/Cd36/Tcf4/Adamts9/Stat1/Ptprm/C3/Sash1/Rhoj/Kit/Stim1/Tlr3/Sema6a/Hspb6/Angpt2/Fgf1/Tek/Ceacam1/Jak1/Agtr1a/Pml/Flt1/Abcc8/Cxcl10/Aqp1/Gpr4/Pgf/Col4a2/Cybb/Meox2/Cflar/Map2k5/Smad1/Sema5a/Serpine1/Cdh5/Dcn/Adam12/Akt3/Tmem100/Itgb3/Nfe2l2/Prkd2/Efnb2/Hdac9/Sulf1/S100a1/Itgb2/Cma1/Tspan12/Prkca/Lrg1/Agt/Pik3r6/Adm/Flcn/Nfatc1/Cd40/Tnfrsf1a/Hgf/Itga5/Ptgis/Gab1/Fgf9/Cyp1b1/Rap1a/Fbxw7/Rapgef3/Hipk2/Dll4/Pdgfd/Dll1/Thbs2/Creb3l1/Card10/Ccl2/Srpx2/Enpp2/Cd34/Kdr/Wt1/Pxn/Notch1/Mecp2/Col4a3/Pdgfb/Ppp1r16b/C3ar1/Aplnr/Hspb1/Fgf2/Smoc2/Jup/Rock2/Cd59a/Ets1/Stard13/Sema4a/Sod2/Cysltr1/Nf1
## GO:0003013 Pparg/Cd36/Pde3a/Epas1/Thrb/Myh6/Ephx2/Sgcg/Ace/Stat1/Mgll/Ptprm/Ppara/Pde2a/Adra1b/Dock4/Nav2/Wwtr1/Ceacam1/Dmpk/Pde4b/Agtr1a/Ryr2/Cacna1g/Fli1/Smad6/Olfr78/Arhgap42/Atp2a2/Heg1/Rgs4/Uty/Kcnq1/Fyn/Ada/Tgfb1/Icam1/Tmem38a/Ramp3/Mef2a/Kcnj5/Egfr/Tnni3/P2ry2/Dmd/Gch1/Ank2/Kcnj2/Trpv4/Bmpr2/Ptp4a3/Coro2b/Wnk1/Akap6/Tgfbr3/Ednrb/Ppard/Kcna5/S100a1/Sgcd/Tmem65/Mrvi1/Cav1/Prkca/Agt/Tnf/Adm/Akap13/Camk2d/Akap12/Kcne4/Pln/Ar/Prcp/Ccl4/Adra1a/Lrp5/Nr3c1/Stk39/Nedd4l/Dll1/Ncald/Pik3r1/Ptk2/Kcne1/Aoc3/Cd34/Hrc/Nox1/Kdr/Nts/Mcpt4/Prkg1/Corin/Smtn/Kcnn2/Grip2/Comp/Ptafr/Npr1/Slc1a1/Mecp2/Adora2b/Slc6a4/Rnf207/Snta1/Jak2/Tcap/Pdgfb/Cacna2d1/C3ar1/Nampt/Nos2/Jup/Rock2/Chrm2/Tnni3k/Cacna1c/Smad7/Cacnb2/Ednra/Dbh/Rps6ka2/Casq2/Sod2/Cysltr1
## GO:0001667                                                                                            Pparg/Sema7a/Adamts9/Ptprm/Sash1/Prkce/Nrp2/Rhoj/Kit/Itga2/Gab2/Sema6a/Kitl/Angpt2/Fgf1/Tek/Ceacam1/Ptprg/Cdh13/Pml/Rreb1/Bcas3/Sema3f/Aqp1/Cygb/Tns1/Meox2/Srgap2/Arid5b/Rffl/Map2k5/Sema5a/Fgf16/Tgfb1/Serpine1/Dock1/Sema6b/Zeb2/Dcn/Bmpr2/Stat5a/Fgfr1/Akt3/Ptp4a3/Smurf2/Itgb3/Nfe2l2/Appl2/Mcc/Tgfbr3/Prkd2/Ednrb/Ppard/Efnb2/Hdac9/Itgb2/Prkca/Agt/Evl/Fap/Jun/Ddr2/Akap12/Map4k4/Cd40/Mapre2/Dnaja4/Prcp/Dusp10/Pecam1/Lrp5/Cyp1b1/Fbxw7/Dll4/Fer/Ptk2/Card10/Gna12/Srpx2/Enpp2/Sema4d/Kdr/Ret/Acvr1/Pik3r3/Clasp1/Sgpl1/Nus1/Itgav/Prox1/Pxn/Pmp22/Notch1/Mecp2/Adora2b/Pdgfb/Pten/Hspb1/Fgf2/Smoc2/Clasp2/Jup/Rock2/Ets1/Has2/Stard13/Sema4a/Sema3g/Nf1
##            Count
## GO:1901342   102
## GO:0003013   126
## GO:0001667   109
```

The output is a table with a first column holding the GO term unique identifier, then second column, the respective description of the GO term, the third - "Gene Ratio" is the ratio between the genes found in that GO category relative to all queried genes, the fourth - "BgRatio" is the ratio between all the genes in the GO category vs all genes in the background gene list, the fifth is pvalue, the sixth padj value (corrected for multuple testing p value using BH method, the seventh is the qvalue - direct estimate of the FDR associated with padj, then "gene ID" - long list ey :), the ninth - "Count" is the number of genes found to be enriched for the respective GO category.

Once we have performed the GO enrichment analysis we can plot the results using the dotplot function from the clusterProfiler package.

1.GO_up -  compareClusterResult object

2.x - the x variable 

3.showCategory - the number of GO categories to be shown

4.color - dots on the plot are going to be colored based on padj value

5.title - title of the plot

GeneRatio in the legend is from the GeneRatio column in the GO_up dataframe above.


```r
dotplot(GO_up, x="count", showCategory=8, color = 'p.adjust', title  = "DEG_UP P1 vs E14.5")
```

```
## wrong orderBy parameter; set to default `orderBy = "x"`
```

![](https://github.com/malina-dva/RNA-seq-Differential-expression-analysis/blob/master/unnamed-chunk-41-1.png)<!-- -->
        
        
Let's now check the GO of down-regulated genes        

```r
GO_down <- enrichGO(gene         = rownames(Sign_genes_down),
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 universe      = all_genes_mouse,
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
```


```r
dotplot(GO_down, x="count", showCategory=8, color = 'p.adjust', title  = "DEG_DOWN P1 vs E14.5")
```

```
## wrong orderBy parameter; set to default `orderBy = "x"`
```

![](https://github.com/malina-dva/RNA-seq-Differential-expression-analysis/blob/master/unnamed-chunk-44-1.png)<!-- -->

Save tables for the GO results (up and Down DEG respectively)


```r
GO_up_dataframe=as.data.frame(GO_up) 

write.table(GO_up_dataframe,sep="\t",file="GO_up_dataframe.txt", row.names=TRUE,col.names=NA,quote=FALSE)


GO_down_dataframe=as.data.frame(GO_down)

write.table(GO_down_dataframe,sep="\t",file="GO_down_dataframe.txt", row.names=TRUE,col.names=NA,quote=FALSE)
```

GO resuls look good! Now you can try to perform the analysis by yourself but comparing Adult vs E14.5 condition for example.

Let me know what top GO terms you are getting. Do they make sense?
    
        
