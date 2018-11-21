---
layout:     post
title:      2018-10-31-Linux软件安装
subtitle:   用conda安装RNA-seq所需要的工具
date:       2018-10-31
author:     DL
header-img: img/home-bg-o.jpg
catalog: true
tags:
    - Linux
---

> 正所谓前人栽树，后人乘凉。
> 
> 感谢[Y大宽](https://www.jianshu.com/u/51a71446d509)、[徐洲更](https://www.jianshu.com/u/9ea40b5f607a)等生信前辈。

# 1. conda软件介绍
对于生信初学者而言，最困难的事情某过于安装各种生信软件，如果一切所有软件都能像`sudo apt-get intall` 或者是`sudo yum install`那样多好。本文就介绍了目前大家认为的最强的非root软件管理器-**conda**。

## 1.1 什么是conda

想要了解什么是conda，需要先要了解什么是[Anaconda](https://anaconda.org/)。

Anaconda是Python的科学发行版，它将各种科学计算工具整合到一个安装包之中，从而使得Python变得无比的强大，就像Linux本身也只是内核，通过整合不同的软件之后才会变得如何的实用。

![](https://s1.ax1x.com/2018/10/30/i26rSx.png)

Anaconda为了避免Python原生`pip`安装软件会出现的问题，比如说Windows下安装科学计算必备的`numpy`和`pandas`时就非常的麻烦，于是它就自己编译了好一些安装包，仅仅使用`conda install`就能下载编译好的二进制包。

因此，conda最开始是Anaconda提供的Python包安装管理工具哦。

## 1.2 为什么用conda 

**conda**最开始只是Anaconda用于管理Python包的工具，但由于它为了避免Python包安装时出现的依赖库不全的问题，相当于又安装了一个虚拟系统，于是乎它能够管理的软件越来越多。

- 官方频道：https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
- conda软件管理依赖环境解决频道：https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
- 生物信息软件官方频道： bioconda
- 生物信息软件清华镜像频道： https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/


基本上，大部分你能想到的软件都能用`conda`安装，如果这些软件还不能的话，还可以基础conda环境进行编译。当然，像docker,mysql这类系统级软件，无法使用conda管理。

因此，使用conda的第一个好处就是**安装方便**。

第二个优点叫做，**环境容易管理**。 当你担心Python2会和Python3冲突的时候，使用conda专门建立一个虚拟环境（下面教程会说），相当于重新开了一台电脑工作。再也不担心Python版本冲突了。而且当你想试用最新版本的工具的时候，完全可以新建一个环境，这样子就不用担心软件不好用无法返回原先版本了。

还有一个优点就是不需要root权限，当管理员没空搭理你，或者处于系统安全考虑不能安装某一个软件的时候，**conda**这类不需要root权限的软件包管理器就是你最好的选择。
PS: 你当然可以选择自己编译，然后解决不断出现的依赖包缺失问题。

# 2. conda的安装及配置

## 2.1 安装Anaconda
直接从官网下载对应于自己python版本的conda即可：

	~ bash Anaconda2-5.3.0-Linux-x86_64.sh

**注**1：miniconda是anaconda的简化版，包括最核心的一些功能，如conda

**注**2：选择Miniconda2和Miniconda3任一都可，因为可以通过虚拟环境创建另一个版本Python环境。

## 2.2 配置conda环境

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

## 2.3 学习conda基本操作的管理环境

搜索bwa进行安装（注意是在conda环境下（base））

	conda search bwa
结果可以发现有很多bwa version可以安装，我们用以下命令安装(y代表yes）

	conda install bwa -y 
下一步安装软件可以先去bioconda官方网站https://bioconda.github.io/recipes.html#recipes查看相应版本等
可以进行搜索，比如samtools，可以看到其很多版本，目前最高1.9（7/24/2018 4:23:35 PM ）

安装samtools

	conda install samtools=1.9 -y 

## 2.4 体验下环境的差别
	echo $PATH 
	/home/dingli/anaconda2/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games
当前的环境是anaconda2环境，下面这个命令就是启动这个环境

	source ~/anaconda2/bin/activate
现在退出miniconda3环境看，还能不能运行bwa和samtools

	source deactivate
其实还是可以，因为bwa已经自动添加进环境变量了（ubuntu16.04）

## 2.5 总结conda安装小技巧

-1 根据软件所用的编程语言确定安装策略

-2 安装conda不要添加到环境变量中，用source activate启动

-3 官方的channel靠后，避免channel之间依赖关系混乱

-4 新建一个或多个安装环境安装生信软件

-5 国内用户利用好清华源镜像

-6 搜索生信软件用https://bioconda.github.io/

# 3. conda的使用
## 3.1 conda的基本操作

安装conda之后，我们需要学习一点最基本的conda使用方法，当然哪里不懂可以到[https://docs.anaconda.com/docs_oss/conda/get-started](https://docs.anaconda.com/docs_oss/conda/get-started) 找到解决方法。不过我相信，下面的已经够用了。

## 3.2 用conda安装转录组分析软件

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
	conda install hisat2 
	conda install sra-tools -y 

