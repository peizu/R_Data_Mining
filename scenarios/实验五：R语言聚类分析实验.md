# 实验目的

了解基本聚类算法的基本原理，并掌握R语言中实现聚类算法的函数方法

# 实验原理

聚类就是将相似的事物聚集在一起，而将不相似的事物划分到不同的类别的过程，是数据分析之中十分重要的一种手段。在数据分析的术语之中，聚类和分类是两种技术。分类是指我们已经知道了事物的类别，需要从样品中学习分类的规则，是一种监督式学习；而聚类则是由我们来给定简单的规则，从而得到类别，是一种无监督学习。

在K均值算法中，质心是定义聚类原型（也就是机器学习获得的结果）的核心。在介绍算法实施的具体过程中，我们将演示质心的计算方法。而且你将看到除了第一次的质心是被指定的以外，此后的质心都是经由计算均值而获得的。

首先，选择K个初始质心（这K个质心并不要求来自于样本数据集），其中K是用户指定的参数，也就是所期望的簇的个数。每个数据点都被收归到距其最近之质心的分类中，而同一个质心所收归的点集为一个簇。然后，根据本次分类的结果，更新每个簇的质心。重复上述数据点分类与质心变更步骤，直到簇内数据点不再改变，或者等价地说，直到质心不再改变。

# 算法源码
```
> kmeans
function (x, centers, iter.max = 10L, nstart = 1L, algorithm = c("Hartigan-Wong", 
    "Lloyd", "Forgy", "MacQueen"), trace = FALSE) 
{
    .Mimax <- .Machine$integer.max
    do_one <- function(nmeth) {
        switch(nmeth, {
            isteps.Qtran <- as.integer(min(.Mimax, 50 * m))
            iTran <- c(isteps.Qtran, integer(max(0, k - 1)))
            Z <- .Fortran(C_kmns, x, m, p, centers = centers, 
                as.integer(k), c1 = integer(m), c2 = integer(m), 
                nc = integer(k), double(k), double(k), ncp = integer(k), 
                D = double(m), iTran = iTran, live = integer(k), 
                iter = iter.max, wss = double(k), ifault = as.integer(trace))
            switch(Z$ifault, stop("empty cluster: try a better set of initial centers", 
                call. = FALSE), Z$iter <- max(Z$iter, iter.max + 
                1L), stop("number of cluster centres must lie between 1 and nrow(x)", 
                call. = FALSE), warning(gettextf("Quick-TRANSfer stage steps exceeded maximum (= %d)", 
                isteps.Qtran), call. = FALSE))
        }, {
            Z <- .C(C_kmeans_Lloyd, x, m, p, centers = centers, 
                k, c1 = integer(m), iter = iter.max, nc = integer(k), 
                wss = double(k))
        }, {
            Z <- .C(C_kmeans_MacQueen, x, m, p, centers = as.double(centers), 
                k, c1 = integer(m), iter = iter.max, nc = integer(k), 
                wss = double(k))
        })
        if (m23 <- any(nmeth == c(2L, 3L))) {
            if (any(Z$nc == 0)) 
                warning("empty cluster: try a better set of initial centers", 
                  call. = FALSE)
        }
        if (Z$iter > iter.max) {
            warning(sprintf(ngettext(iter.max, "did not converge in %d iteration", 
                "did not converge in %d iterations"), iter.max), 
                call. = FALSE, domain = NA)
            if (m23) 
                Z$ifault <- 2L
        }
        if (nmeth %in% c(2L, 3L)) {
            if (any(Z$nc == 0)) 
                warning("empty cluster: try a better set of initial centers", 
                  call. = FALSE)
        }
        Z
    }
    x <- as.matrix(x)
    m <- as.integer(nrow(x))
    if (is.na(m)) 
        stop("invalid nrow(x)")
    p <- as.integer(ncol(x))
    if (is.na(p)) 
        stop("invalid ncol(x)")
    if (missing(centers)) 
        stop("'centers' must be a number or a matrix")
    nmeth <- switch(match.arg(algorithm), `Hartigan-Wong` = 1L, 
        Lloyd = 2L, Forgy = 2L, MacQueen = 3L)
    storage.mode(x) <- "double"
    if (length(centers) == 1L) {
        k <- centers
        if (nstart == 1L) 
            centers <- x[sample.int(m, k), , drop = FALSE]
        if (nstart >= 2L || any(duplicated(centers))) {
            cn <- unique(x)
            mm <- nrow(cn)
            if (mm < k) 
                stop("more cluster centers than distinct data points.")
            centers <- cn[sample.int(mm, k), , drop = FALSE]
        }
    }
    else {
        centers <- as.matrix(centers)
        if (any(duplicated(centers))) 
            stop("initial centers are not distinct")
        cn <- NULL
        k <- nrow(centers)
        if (m < k) 
            stop("more cluster centers than data points")
    }
    k <- as.integer(k)
    if (is.na(k)) 
        stop(gettextf("invalid value of %s", "'k'"), domain = NA)
    if (k == 1L) 
        nmeth <- 3L
    iter.max <- as.integer(iter.max)
    if (is.na(iter.max) || iter.max < 1L) 
        stop("'iter.max' must be positive")
    if (ncol(x) != ncol(centers)) 
        stop("must have same number of columns in 'x' and 'centers'")
    storage.mode(centers) <- "double"
    Z <- do_one(nmeth)
    best <- sum(Z$wss)
    if (nstart >= 2L && !is.null(cn)) 
        for (i in 2:nstart) {
            centers <- cn[sample.int(mm, k), , drop = FALSE]
            ZZ <- do_one(nmeth)
            if ((z <- sum(ZZ$wss)) < best) {
                Z <- ZZ
                best <- z
            }
        }
    centers <- matrix(Z$centers, k)
    dimnames(centers) <- list(1L:k, dimnames(x)[[2L]])
    cluster <- Z$c1
    if (!is.null(rn <- rownames(x))) 
        names(cluster) <- rn
    totss <- sum(scale(x, scale = FALSE)^2)
    structure(list(cluster = cluster, centers = centers, totss = totss, 
        withinss = Z$wss, tot.withinss = best, betweenss = totss - 
            best, size = Z$nc, iter = Z$iter, ifault = Z$ifault), 
        class = "kmeans")
}
```

# 实验步骤

我们使用r语言中的数据集iris实现k-means聚类算法。  
首先载入iris数据集并命名为iris1，移除Species类别属性：

```
> iris1 <- iris
> iris1$Species <- NULL
```

对iris1数据集调用函数kmeans\(\)，并将结果存储在变量kmeans.result中，在命令外面加小括号可以直接查看函数调用结果。**每次操作输出结果困难不同，相应的，后面画出的散点图也会不同**

```
> (kmeans.result <- kmeans(iris1,3))

#输出结果为：
K-means clustering with 3 clusters of sizes 33, 21, 96
Cluster means:
  Sepal.Length Sepal.Width Petal.Length Petal.Width
1     5.175758    3.624242     1.472727   0.2727273
2     4.738095    2.904762     1.790476   0.3523810
3     6.314583    2.895833     4.973958   1.7031250
Clustering vector:
[1] 1 2 2 2 1 1 1 1 2 2 1 1 2 2 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1 2 2 1 1 1 2 1 1 1 2 1 1 2 2 1 1 2 1
[48] 2 1 1 3 3 3 3 3 3 3 2 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2
[95] 3 3 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
[142] 3 3 3 3 3 3 3 3 3
Within cluster sum of squares by cluster:
[1]   6.432121  17.669524 118.651875
(between_SS / total_SS =  79.0 %)
Available components:
[1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss" "betweenss"    "size"         "iter"         "ifault"
```

绘制散点图并标记所有的簇以及簇中心，需要注意的是数据集有四个维度但是绘图只使用前两个维度

```
> plot(iris1[c("Sepal.Length","Sepal.Width")],col=kmeans.result$cluster)
> points(kmeans.result$centers[,c("Sepal.Length","Sepal.Width")],col=1:3,pch=8,cex=2)
```

![](/images/2-2-5-1_20171107091014.014.jpeg)
