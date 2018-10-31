---
layout:     post
title:      Ubuntu下conda的安装及使用
subtitle:   用conda安装RNA-seq所需要的工具
date:       2018-10-31
author:     DL
header-img: img/post-bg-universe.jpg
catalog: true
tags:
    - Conda
    - Ubuntu
---

# 一、安装Anaconda
直接从官网下载对应于自己python版本的conda即可：

	~ bash Anaconda2-5.3.0-Linux-x86_64.sh
## 配置conda环境

先启动conda环境

	source ~/miniconda3/bin/activate
添加频道

	conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/ 
	conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/ 
	conda config --set show_channel_urls yes
添加第三方频道	

	conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/

安装bioconda频道

	conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ 
# 二、学习conda基本操作的管理环境

搜索bwa进行安装（注意是在conda环境下（base））

	conda search bwa
结果可以发现有很多bwa version可以安装，我们用以下命令安装(y代表yes）

	conda install bwa -y 
下一步安装软件可以先去bioconda官方网站https://bioconda.github.io/recipes.html#recipes查看相应版本等
可以进行搜索，比如samtools，可以看到其很多版本，目前最高1.9（7/24/2018 4:23:35 PM ）

安装samtools

	conda install samtools=1.9 -y 
## 体验下环境的差别
	echo $PATH 
	/home/dingli/anaconda2/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games
当前的环境是anaconda2环境，下面这个命令就是启动这个环境

	source ~/anaconda2/bin/activate
现在退出miniconda3环境看，还能不能运行bwa和samtools

	source deactivate
其实还是可以，因为bwa已经自动添加进环境变量了（ubuntu16.04）

## 总结conda安装小技巧

-1 根据软件所用的编程语言确定安装策略

-2 安装conda不要添加到环境变量中，用source activate启动

-3 官方的channel靠后，避免channel之间依赖关系混乱

-4 新建一个或多个安装环境安装生信软件

-5 国内用户利用好清华源镜像

-6 搜索生信软件用https://bioconda.github.io/

# 三、用conda安装转录组分析软件

-hisat2 samtools sratoolkit

-htseq-count

-fastqc trimmomatics

	python2环境（base）  
 
	conda install fastqc trimmomatic（conda可以同时指定两个软件安装） 
注意要在python2环境下安装htseq，先启动python2环境

	source activate python2
	conda install htseq(/一定注意不能在python3环境安装，要启动python2环境！我不小心按了y，在python3安装完成了，然后用conda install htseq卸载） 
	conda install htseq -y
	启动htseq-count
	htseq-count

继续搜索hisat2
	source deactivate
	conda install hisat2 
	conda install hisat2 sra-tools -y 

From Y大宽  生信猿