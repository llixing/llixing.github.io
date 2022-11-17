---
layout:     post
title:      2020-10-28-Skin_traits_preanalysis
subtitle:   R代码整理
date:       2020-10-28
author:     DL
header-img: img/home-bg-o.jpg
catalog: true
tags:
    - R
---

### 1.Venn图

[Draw Venn Diagram](http://bioinformatics.psb.ugent.be/webtools/Venn/)

### 2.Summary_description

```
skin_data <- read.table(file="your file",sep="\t",header=T,row.names=1)
summary_result <- summary(skin_data)
write.csv(summary_result,file="summary_result.csv")

```

### 3. Boxplot and density_distribution

```
boxplot(skin_data)

install.package("ggplot2")
library(ggplot2)
p<-ggplot(data, aes(x = your column))
p + geom_density(color = "red", fill = "red")

```

### 4. Correlation_analysis

```
library(psych)
library(corrplot)
library(RColorBrewer)

rm(list=ls())
a <- read.table(file="clipboard",sep="\t",header=T,row.names=1)
ct <- corr.test(a)

corrplot(as.matrix(ct$r),method = "square",tl.cex = 0.8,cl.cex = 0.8,
         col = colorRampPalette(c("royalblue","white","firebrick3"))(10),tl.col = "black",cl.pos = "b")

write.table(ct$r,file="Pain_skin_cor.txt",sep="\t",quote=F)
write.table(ct$p,file="Pain_skin_pvalue.txt",sep="\t",quote=F)
```

### 5. Linear_regression_analysis

```
data_548 <- read.table(file="covariate.txt",sep="\t",header=T,row.names = 1)
data_548$Sex <-as.factor(data_548$Sex)


result_lm <- c()

for (i in c(4:117)) {
    fit <- lm(get(colnames(data_548)[i]) ~ + Age + Sex + BMI,data = data_548)
    
    b <- coef(summary(fit))
    rownames(b)[1] <- colnames(data_548)[i]
    rownames(b)[1] <-paste(rownames(b)[1],"-Intercept",sep="")
    result_lm <- rbind(result_lm,b)
}


write.csv(result_lm,file="result_lm.csv")

```

