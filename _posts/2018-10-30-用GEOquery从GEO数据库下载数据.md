---
layout:     post
title:      DNA甲基化入门
subtitle:   高通量测序领域常用名词解释
date:       2018-10-15
author:     DL
header-img: img/home-bg-o.jpg
catalog: true
tags:
    - NGS
    - Bioinformatics
---

# 用GEOquery从GEO数据库下载数据
===

>[Gene Expression Omnibus database (GEO)](http://www.ncbi.nlm.nih.gov/geo/)是由NCBI负责维护的一个数据库，设计初衷是为了收集整理各种表达芯片数据，但是后来也加入了甲基化芯片，甚至高通量测序数据！

##[GEO数据库基础知识](http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/geo/)

  * GEO Platform (GPL) 芯片平台
  
  * GEO Sample (GSM) 样本ID号
  
  * GEO Series (GSE)  study的ID号
 
  * GEO Dataset (GDS) 数据集的ID号
  
##GEOquery包的用法
>只需要记住三个函数和两个对象,以及它们把数据下载到哪里了~

>getGEO/getGEOfile/getGEOSuppFiles/GDS对象/expression set 对象

>这三个函数根据上面的四种ID号下载数据时候，返回的对象不一样！


###首先是下载和加载包：
      source("http://www.bioconductor.org/biocLite.R")
      
      biocLite("GEOquery")
      
      library(GEOquery)
  
##然后是使用它！

首先，我们介绍getGEO函数，它可以针对三种ID号来下载数据
    
      gds858 <- getGEO('GDS858', destdir=".") ##根据GDS号来下载数据，下载soft文件
      
      gpl96 <- getGEO('GPL96', destdir=".") ##根据GPL号下载的是芯片设计的信息！
      
      gse1009 <- getGEO('GSE1009', destdir=".")##根据GSE号下载数据，下载matrix.txt.gz
      
  
下载的文件都会保存在本地，返回的对象不一样！
  
gds858返回的对象很复杂:
用Table(gds858)可以得到表达矩阵！
用Meta(gds858)可以得到描述信息
```R
names(Meta(gds858))
Table(gds858)[1:5,1:5]
```
然后还可以用 **GDS2eSet**函数把它转变为*expression set 对象*

eset <- GDS2eSet(gds858, do.log2=TRUE)

也就是直接根据GSE号返回的对象：gse1009  ##expression set 对象非常重要，我会单独讲解它

我们的处理函数有：geneNames/sampleNames/pData/

但是根据GPL号下载返回的对象跟GDS一样，也是用Table/Meta处理！

```R
names(Meta(gpl96))
Table(gpl96)[1:10,1:4]
##下面这个就是芯片ID的基因注释信息
Table(gpl96)[1:10,c("ID","GB_LIST","Gene.Title","Gene.Symbol","Entrez.Gene")]
```

getGEO除了可以下载数据，还可以打开本地数据！

gds858 <- getGEO(filename='GDS858.soft.gz')


```R
tmp=getGEOSuppFiles(GSE1009)
if (is.null(tmp)) {
    warning("Supplementary data files not provided!\n")
}
```
把代码都自己运行一遍， 就明白了！