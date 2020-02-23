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


```{r}
require("DESeq2")

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


```{r}

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

https://github.com/malina-dva/RNA-seq-Data-normalization-and-clustering

Click on Raw_counts_input.txt file

Then click on Raw.

The file with the Raw_counts_input.txt will open in a separate page.

Right click and then click on "save page as" or "save as" to save the file, preferably on your Desktop. Make sure you save the file with its original name.

Once the file is saved on your Desktop move it to your working directory that we just created and set. (Tutorial_RNA_seq).


It is a good idea to open the input file in excel and check its content before you start to work with it in R.
**However, I dont' recommend keeping input files open in excel while you are using them in R.**

Repeating again! The file with the raw counts has to be placed in your working directory.
To read/load the input file in Rstudio use the command below.
```{r}

just.raw.counts = read.delim("Raw_counts_input.txt")
```

Let's check first several rows of the input file using the function 'head'.

```{r}
head(just.raw.counts)
```


Now let's check the dimensions (total number of rows and columns) of your imported table. 

```{r}
dim(just.raw.counts)
```

There are 27195 genes (rows) and9 columns. The first column (e.g. Probe) holds the names of the genes, the remaining 8 columns hold the gene counts (gene expression info) for the replicates of the four conditions we have. 
We want to actually specify that the column named 'Probe' has the information for the names of each row in the table. 

We can do that when we read the input file by specifying row.names =1, meaning that names of the rows are to be found in the first column.

```{r}
just.raw.counts = read.delim("Raw_counts_input.txt", row.names = 1)
```
```{r}
head(just.raw.counts)
```

After defining row.names = 1, the first column is no longer named Probe as you can see from output above.
```{r}
dim(just.raw.counts)
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

https://github.com/malina-dva/RNA-seq-Data-normalization-and-clustering

Click on meta_data.txt file

Then click on Raw

The file with the meta_data.txt will open in a separate page

Right click and then click on "save page as" or "save as" to save the file, preferably on your Desktop. Keep the original name of the file.

Move meta_data.txt to the Tutorial_RNA_seq folder.


Open it with excel to check its content.

Then let's load it in Rstudio.


```{r}
meta.data = read.delim(file="meta_data.txt", row.names = 1)
```
Again we specify the first column as row.names.

```{r}
head(meta.data)
```
The column called 'condition' links the name of each replicate from 'just.raw.counts' to its respective condition.
For example E14.5_R1 and E14.5_R2, in the 'condition' column have the same name e.g E14.5, because they both belong to the E14.5 condition.


Now we can proceed to build DESeq object.

First we need to create DESeq matrix. 
```{r}
count.data.set <- DESeqDataSetFromMatrix(countData=just.raw.counts, 
                                         colData=meta.data, design= ~ condition) 
```

'countData' argument requires the table with the counts and the names of the genes
'colData' argument requires the metadata file
'design' requires the column of the metadata that holds the information for the relationships between replicates and conditions (e.g the 'condition' column from the metadata in our case)



Then we create DESeq object.
```{r}
count.data.set.object <- DESeq(count.data.set)
```

### 5. Differential expression analysis 
For performing differential expression analysis we don't need to normalize the data beforehand.
DESeq2 uses raw counts to perform the statistical test for identifying genes that are differentially expressed between conditions.

The function that we use to obtain a list with the differentially expressed genes is called "results".

We can use this function by passing Deseq object to it. **We just created one above (count.data.set.object)**.
Here we are going to use "contrast" argument to get the differentially expressed genes (DEG) for the comparison between two **conditions** -  Neonatal compared to E14.5.


```{r}

res <- results(count.data.set.object, contrast=c("condition","Neonatal","E14.5"))
```

To check the first few lines of the "results" output just run res.

```{r}

res 
```
There are six columns in the output (in addition to the gene names). The most important for us are log2FoldChange (how many times the gene expression in the Neonatal changes relative to E14.5) and padj (p value adjusted for multiple testing) columns, which shows us the significance of the difference.
If the value of the log2FoldChange is positive, it means the gene was upregulated in Neonatal compared to E14.5, and if it is negative, the gene was downregulated in Neonatal compared to E14.5 respectively.
To get a summary of the output use:
```{r}
summary(res)
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
```{r}
res = na.omit(res)
```
Let's check if there are any NAs left

```{r}
head(res)
```

Good for us, it seems, at least when looking at the first few rows, that there are no NAs left in the output.


Now let's select only the significant genes based on padj value cutoff equals or lower than 0.05. 
We can do that by the following command:

```{r}
resfiltered = res[res$padj <= 0.05,]
head(resfiltered)
```
The commands asks to get all the rows(genes), where padj column has values equal or lower than 0.05.
Now we want to check the summary of the filtered results for padj and NAs.
```{r}
summary(resfiltered)
```

You can see that after all the steps of filtering we got 2615 up and 2033 down regulated genes respectively, using threshold of padj 0.05.
No NAs are present, as per the independentFiltering output.

But wait! In the summary the adjusted p-value is < 0.1. It is really? Haven't we already set a different threshold?
We can check if there are genes with padj value higher than 0.05 with the following line of code:

```{r}
sum(resfiltered$padj > 0.05)
```

With this code we get the number of rows that have padj value > 0.05. it is clear that there are none, which means that our filtering for padj works.

Yet,the adjusted p-value will always be shown as  < 0.1 (the default) in the summary output, unless we specify a diffrent value for "alpha" (padj value) directly when we generate the res variable!
If we were to that the command will be:

res <- results(count.data.set.object, contrast=c("condition","Neonatal","E14.5"), alpha=0.05)


Finally, we want to order the results so that on the top there will be genes with the lowest (most significant padj value).

#### 2. Order the results based on increasing padj value

We use the function "order" and apply it on the values in the resfiltered$padj column and for all rows (genes) of the resfiltered variable.
The output from the reordering we keep in a new variable called res.filtered.ordered

```{r}
res.filtered.ordered = resfiltered[order(resfiltered$padj),] 
```

Let's check what we got now.

```{r}
head(res.filtered.ordered)

```

You see that the genes with the lowest padj values are on the top with a "leader", a gene called Parp14 (1.68397576874913e-86). Veryyy low padj value though, very highly significant difference (upregulation based on the logFoldChange, which is a positive value ~3.8) in Neonatal compared to E14.5.

Now we are going to save the filtered and ordered list with the DEG.

```{r}
write.table(res.filtered.ordered, sep="\t",file="Results_filtered_ordered.txt", row.names=TRUE,col.names=NA,quote=FALSE)
```

#### 3. logFoldChange based filtering and visualization
**This step of analysis is just for demonstration. Generally, cutoff based on FC is not needed, unless the is a very good reason behind it.**

But anyway, let's check what is the distribution of the logFoldChange values for all genes in the results we have.
We will need ggplot2 to do that. By now you should be a PRO in package installation and loading. :P If not google BiocManager and ggplot2 and follow the steps for installation. Then load ggplot2 using the command below:

```{r}
require(ggplot2)

```

We need to convert the output with our DEG (res.filtered.ordered) into dataframe in order to be able to visualize it on ggplot2.

```{r}

resOrdFilt.data.frame = as.data.frame(res.filtered.ordered)
head(resOrdFilt.data.frame)
```

Now let's plot it. We will plot a **histogram** of the distribution of values in the column with log2FoldChange (e.g. the x in the aesthetics in the ggplot function).

```{r}
ggplot(resOrdFilt.data.frame, aes(x=log2FoldChange)) +
  geom_histogram()

```
 

We can zoom the x axis a bit by specifying limits=c(-5,5) to the scale_x_continuous argument (meaning the x axis). 
Additionally we can put labels to the positive and negative 5, 1 and 0.5 values of the x axis (e.g specifying that in the breaks argument) and increase the number of histogram bins from the default 30 to 100, to be able to better see where the distribution peaks.
```{r}
ggplot(resOrdFilt.data.frame, aes(x=log2FoldChange)) +
  geom_histogram(bins=100) +
  scale_x_continuous( name ="log2FoldChange",breaks=c(-5,-1,-0.5,0, 0.5,1,5), limits = c(-5,5))

```


Based on the histogram it is obvious that most of FC are between 0.5 and 1. Just for demonstration we can set a threshold for the FC to -0.5 and +0.5 for down and up regulated genes respectively to see how many DEG will remain.

To filter for log2FoldChange threshold we use "subset". We apply subset on the resOrdFilt.data.frame. The **conditions** we want to sunset on are:  log2FoldChange column from resOrdFilt.data.frame is either >0.5 or <0.5. To achieve that, we use "|" argument in our command which means **or**. 

```{r}
resOrdFiltLFC = subset(resOrdFilt.data.frame,log2FoldChange >0.5 | log2FoldChange < -0.5 )

```

To check how many DEG remained after filtering for fold change > or < 0.5:
```{r}
dim(resOrdFiltLFC)
```

Well, we have reduced the number of DEG from 4648 to 3985 by applying filtering for LFC.

```{r}
ggplot(resOrdFiltLFC, aes(x=log2FoldChange)) +
  geom_histogram(bins=100) +
  scale_x_continuous( name ="log2FoldChange",breaks=c(-5,-1,-0.5,0, 0.5,1,5), limits = c(-5,5))

```

Let's write the results with the filtered FCs.

```{r}
write.table(resOrdFiltLFC, sep="\t",file="Results_filtered_ordered_LFC.txt", row.names=TRUE,col.names=NA,quote=FALSE)
```
**As I mentioned in the beginning of the section, filtering on FC is not always recommended or needed. The included steps from above were only for demonstration.**

#### 4. Subsetting and visualizing the DEGs 

We will continue to work with our original list of DEG **prior to foldchange filtering**. These were the results stored in res.filtered.ordered variable. We want to get two tables with either up or down DEG and save them. 

```{r}
summary(res.filtered.ordered)
```

To be able to subset our up and down DEG genes we will use the resOrdFilt.data.frame variable which holds the results from res.filtered.ordered.

```{r}
Sign_genes_up = subset(resOrdFilt.data.frame, log2FoldChange >0)
Sign_genes_down = subset(resOrdFilt.data.frame, log2FoldChange <0)
```

Then save the up and down genes to text files.

```{r}
write.table(Sign_genes_up, sep="\t",file="Sign_genes_up.txt", row.names=TRUE,col.names=NA,quote=FALSE)

write.table(Sign_genes_down, sep="\t",file="Sign_genes_down.txt", row.names=TRUE,col.names=NA,quote=FALSE)
```

We can directly visualize the up and down genes using "EnhancedVolcano" package.

To download the package please copy and execute the commands from:
http://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html


Now... let's create a volcano plot.

```{r}
require("EnhancedVolcano")
```

For labels on the volcano plot, we are going to use the rownames(resOrdFilt.data.frame) which are in fact the names of the genes in the resOrdFilt.data.frame.
Not a bad idea to have a look at the first few lines of the resOrdFilt.data.frame, right?
```{r}
head(resOrdFilt.data.frame)
```
The x axis on the volcano plot is going to be log2FoldChange column from resOrdFilt.data.frame. The y axis will be padj column respectively. For the FC cutoff (FoldChange cutoff) let's specify 0.3 (that will refer to +/- 3), and pCutoff (which in our case is actually padj) let's specify 0.05.
Ideally, one would like to change the labels of the y axis and legend to padj instead of p value. You can try to do that as a **homework** by following the steps in this tutorial: https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html (point 4.6)

```{r}
EnhancedVolcano(resOrdFilt.data.frame,

        lab = rownames(resOrdFilt.data.frame),

        x = "log2FoldChange",
        y = "padj",
        FCcutoff = 0.3,
         pCutoff = 0.05)
```      

You see that all of our genes from resOrdFilt.data.frame have a difference in FC greater than 0.3 and the padj value is < 0.05, which is what we would expect anyway!

###6. GO enrichment analysis and visualization

Let's check what biological processes the Up and Down genes are involved in by performing GO enrichment analysis. 
Of course, will have to load the packages that will help us find this out.
Again, please install them as per usual if you don't have them already (Bioconductor, BiocManager and the name of the package). To load them use the commands below.
```{r}
require("clusterProfiler")
require("org.Mm.eg.db")
```
Now we are going to create a variable with all unique mouse genes.
We will need that as a gene background for our GO analysis.
```{r}
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

```{r}
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

```{r}
GO_up_dataframe=as.data.frame(GO_up)
head(GO_up_dataframe, n=3)
```

The output is a table with a first column holding the GO term unique identifier, then second column, the respective description of the GO term, the third - "Gene Ratio" is the ratio between the genes found in that GO category relative to all queried genes, the fourth - "BgRatio" is the ratio between all the genes in the GO category vs all genes in the background gene list, the fifth is pvalue, the sixth padj value (corrected for multuple testing p value using BH method, the seventh is the qvalue - direct estimate of the FDR associated with padj, then "gene ID" - long list ey :), the ninth - "Count" is the number of genes found to be enriched for the respective GO category.

Once we have performed the GO enrichment analysis we can plot the results using the dotplot function from the clusterProfiler package.

1.GO_up -  compareClusterResult object

2.x - the x variable 

3.showCategory - the number of GO categories to be shown

4.color - dots on the plot are going to be colored based on padj value

5.title - title of the plot

GeneRatio in the legend is from the GeneRatio column in the GO_up dataframe above.

```{r}
dotplot(GO_up, x="count", showCategory=8, color = 'p.adjust', title  = "DEG_UP P1 vs E14.5")
```
        
        
Let's now check the GO of down-regulated genes        
```{r}
GO_down <- enrichGO(gene         = rownames(Sign_genes_down),
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 universe      = all_genes_mouse,
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

```

```{r}
dotplot(GO_down, x="count", showCategory=8, color = 'p.adjust', title  = "DEG_DOWN P1 vs E14.5")
```  

Save tables for the GO results (up and Down DEG respectively)

```{r}
GO_up_dataframe=as.data.frame(GO_up) 

write.table(GO_up_dataframe,sep="\t",file="GO_up_dataframe.txt", row.names=TRUE,col.names=NA,quote=FALSE)


GO_down_dataframe=as.data.frame(GO_down)

write.table(GO_down_dataframe,sep="\t",file="GO_down_dataframe.txt", row.names=TRUE,col.names=NA,quote=FALSE)
```

GO resuls look good! Now you can try to perform the analysis by yourself but comparing Adult vs E14.5 condition for example.

Let me know what top GO terms you are getting. Do they make sense?
    
        
