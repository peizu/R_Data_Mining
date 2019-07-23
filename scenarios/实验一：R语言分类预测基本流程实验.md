# 实验目的

了解分类预测的基本流程，并掌握R语言实现分类算法的基本方法

# 实验原理

分类是数据挖掘中应用领域极其广泛的重要技术之一,至今已经提出很多算法。分类是根据数据集的特点构造一个分类器,利用分类器对未知类别的样本赋予类别的一种技术。构造分类器的过程一般分为训练和测试两个步骤。在训练阶段,分析训练数据集的特点,为每个类别产生一个对相应数据集的准确描述或模型。在测试阶段,利用类别的描述或模型对测试进行分类,测试其分类准确度。

本节教程使用k-近邻\(knn\)算法作为例子解释R语言中实现分类算法的基本流程。

k-近邻\(kNN，k-Nearest Neighbors\)算法是一种基于实例的分类方法。该方法就是找出与未知样本x距离最近的k个训练样本，看这k个样本中多数属于哪一类，就把x归为那一类。k-近邻方法是一种懒惰学习方法，它存放样本，直到需要分类时才进行分类，如果样本集比较复杂，可能会导致很大的计算开销，因此无法应用到实时性很强的场合。

K值大小的选取会对K近邻算法的结果会产生重大影响。如果选择较小的K值，就相当于用较小的领域中的训练实例进行预测，“学习”近似误差会减小，只有与输入实例较近或相似的训练实例才会对预测结果起作用，与此同时带来的问题是“学习”的估计误差会增大，换句话说，K值的减小就意味着整体模型变得复杂，容易发生过拟合；如果选择较大的K值，就相当于用较大领域中的训练实例进行预测，其优点是可以减少学习的估计误差，但缺点是学习的近似误差会增大。这时候，与输入实例较远（不相似的）训练实例也会对预测器作用，使预测发生错误。在实际应用中，K值一般取一个比较小的数值，例如采用交叉验证法（简单来说，就是一部分样本做训练集，一部分做测试集）来选择最优的K值。

# 算法源码
```
> kknn
function (formula = formula(train), train, test, na.action = na.omit(), 
    k = 7, distance = 2, kernel = "optimal", ykernel = NULL, 
    scale = TRUE, contrasts = c(unordered = "contr.dummy", ordered = "contr.ordinal")) 
{
    if (is.null(ykernel)) 
        ykernel = 0
    weight.y = function(l = 1, diff = 0) {
        k = diff + 1
        result = matrix(0, l, l)
        diag(result) = k
        for (i in 1:(k - 1)) {
            for (j in 1:(l - i)) {
                result[j, j + i] = k - i
                result[j + i, j] = k - i
            }
        }
        result
    }
    kernel <- match.arg(kernel, c("rectangular", "triangular", 
        "epanechnikov", "biweight", "triweight", "cos", "inv", 
        "gaussian", "rank", "optimal"), FALSE)
    ca <- match.call()
    response = NULL
    old.contrasts <- getOption("contrasts")
    options(contrasts = contrasts)
    formula = as.formula(formula)
    mf <- model.frame(formula, data = train)
    mt <- attr(mf, "terms")
    mt2 <- delete.response(mt)
    cl <- model.response(mf)
    d <- sum(attr(mt, "order"))
    if (is.ordered(cl)) {
        response <- "ordinal"
        lev <- levels(cl)
    }
    if (is.numeric(cl)) 
        response <- "continuous"
    if (is.factor(cl) & !is.ordered(cl)) {
        response <- "nominal"
        lev <- levels(cl)
    }
    if (distance <= 0) 
        stop("distance must >0")
    if (k <= 0) 
        stop("k must >0")
    learn <- model.matrix(mt, mf)
    valid <- model.matrix(mt2, test)
    m <- dim(learn)[1]
    p <- dim(valid)[1]
    q <- dim(learn)[2]
    ind <- attributes(learn)$assign
    d.sd <- numeric(length(ind)) + 1
    we <- numeric(length(ind)) + 1
    d.sd = apply(learn, 2, stats::var)
    for (i in unique(ind)) {
        d.sd[ind == i] = sqrt(mean(d.sd[ind == i]))
        we[ind == i] = 1/sum(ind == i)
    }
    we[d.sd == 0] = 0
    d.sd[d.sd == 0] = 1
    if (scale) {
        learn <- sweep(learn, 2L, d.sd, "/", check.margin = FALSE)
        valid <- sweep(valid, 2L, d.sd, "/", check.margin = FALSE)
    }
    ord = order(we * apply(learn, 2, sd), decreasing = TRUE)
    we = we[ord]
    learn = learn[, ord, drop = FALSE]
    valid = valid[, ord, drop = FALSE]
    Euclid <- FALSE
    if (distance == 2) 
        Euclid <- TRUE
    if (Euclid) 
        dmtmp <- .C("dmEuclid", as.double(learn), as.double(valid), 
            as.integer(m), as.integer(p), as.integer(q), dm = double((k + 
                1L) * p), cl = integer((k + 1L) * p), k = as.integer(k + 
                1), as.double(distance), as.double(we), PACKAGE = "kknn")
    else dmtmp <- .C("dm", as.double(learn), as.double(valid), 
        as.integer(m), as.integer(p), as.integer(q), dm = double((k + 
            1L) * p), cl = integer((k + 1L) * p), k = as.integer(k + 
            1), as.double(distance), as.double(we), PACKAGE = "kknn")
    D <- matrix(dmtmp$dm, nrow = p, ncol = k + 1)
    C <- matrix(dmtmp$cl, nrow = p, ncol = k + 1)
    maxdist <- D[, k + 1]
    maxdist[maxdist < 1e-06] <- 1e-06
    D <- D[, 1:k]
    C <- C[, 1:k] + 1
    CL <- matrix(cl[C], nrow = p, ncol = k)
    if (response != "continuous") {
        l <- length(lev)
        weightClass <- matrix(0, p, l)
    }
    if (response == "continuous") {
        weightClass <- NULL
    }
    W <- D/maxdist
    W <- pmin(W, 1 - (1e-06))
    W <- pmax(W, 1e-06)
    if (kernel == "rank") 
        W <- (k + 1) - t(apply(as.matrix(D), 1, rank))
    if (kernel == "inv") 
        W <- 1/W
    if (kernel == "rectangular") 
        W <- matrix(1, nrow = p, ncol = k)
    if (kernel == "triangular") 
        W <- 1 - W
    if (kernel == "epanechnikov") 
        W <- 0.75 * (1 - W^2)
    if (kernel == "biweight") 
        W <- dbeta((W + 1)/2, 3, 3)
    if (kernel == "triweight") 
        W <- dbeta((W + 1)/2, 4, 4)
    if (kernel == "cos") 
        W <- cos(W * pi/2)
    if (kernel == "triweights") 
        W <- 1
    if (kernel == "gaussian") {
        alpha = 1/(2 * (k + 1))
        qua = abs(qnorm(alpha))
        W = W * qua
        W = dnorm(W, sd = 1)
    }
    if (kernel == "optimal") {
        W = rep(optKernel(k, d = d), each = p)
    }
    W <- matrix(W, p, k)
    if (response != "continuous") {
        for (i in 1:l) {
            weightClass[, i] <- rowSums(W * (CL == lev[i]))
        }
        weightClass <- weightClass/rowSums(weightClass)
        colnames(weightClass) <- lev
    }
    if (response == "ordinal") {
        blub = length(lev)
        weightClass = weightClass %*% weight.y(blub, ykernel)
        weightClass <- weightClass/rowSums(weightClass)
        weightClass <- t(apply(weightClass, 1, cumsum))
        colnames(weightClass) <- lev
        fit <- numeric(p)
        for (i in 1:p) fit[i] <- min((1:l)[weightClass[i, ] >= 
            0.5])
        fit <- ordered(fit, levels = 1:l, labels = lev)
    }
    if (response == "nominal") {
        fit <- apply(weightClass, 1, order, decreasing = TRUE)[1, 
            ]
        fit <- factor(fit, levels = 1:l, labels = lev)
        if (kernel == "rectangular" && k > 1) {
            blub <- apply(weightClass, 1, rank, ties.method = "max")
            indices = (1:p)[colSums(blub == l) > 1]
            blub = t(blub)
            nM = matrix(0, p, l)
            colnames(nM) = lev
            for (i in 1:l) nM[, i] = apply((CL == lev[i]) %*% 
                diag(1:k), 1, max)
            nM = (blub == l) * nM
            nM[nM == 0] <- k + 1
            fitv = numeric(p)
            for (i in indices) fitv[i] = which(nM[i, ] == min(nM[i, 
                ]))
            fit[indices] <- factor(fitv[indices], levels = 1:l, 
                labels = lev)
        }
    }
    if (response == "continuous") 
        fit <- rowSums(W * CL)/pmax(rowSums(W), 1e-06)
    options(contrasts = old.contrasts)
    result <- list(fitted.values = fit, CL = CL, W = W, D = D, 
        C = C, prob = weightClass, response = response, distance = distance, 
        call = ca, terms = mt)
    class(result) = "kknn"
    result
}
```

# 实验步骤

我们使用R语言中的kknn\(\)函数对iris数据集进行knn算法的实现。  
首先加载kknn包

```
> library(kknn)
```

加载iris数据集，并生成采样数据，将iris数据集分为两部分，分别为训练数据集iris.learn和测试数据集iris.valid

```
> m <- dim(iris)[1]
> val <- sample(1:m, size = round(m/3), replace = FALSE, prob = rep(1/m, m))
> iris.learn <- iris[-val,]
> iris.valid <- iris[val,]
```

调用kknn\(\)函数

```
> iris.kknn <- kknn(Species~., iris.learn, iris.valid, distance = 1, kernel = "triangular")
> summary(iris.kknn)

#输出结果为：
#输出结果每次不同
Call:
kknn(formula = Species ~ ., train = iris.learn, test = iris.valid,     distance = 1, kernel = "triangular")

Response: "nominal"
          fit prob.setosa prob.versicolor prob.virginica
1   virginica           0      0.00000000     1.00000000
2   virginica           0      0.01624687     0.98375313
3  versicolor           0      1.00000000     0.00000000
4   virginica           0      0.10041431     0.89958569
5  versicolor           0      1.00000000     0.00000000
6   virginica           0      0.00000000     1.00000000
7  versicolor           0      1.00000000     0.00000000
8  versicolor           0      1.00000000     0.00000000
9      setosa           1      0.00000000     0.00000000
10 versicolor           0      0.83072576     0.16927424
11     setosa           1      0.00000000     0.00000000
12     setosa           1      0.00000000     0.00000000
13  virginica           0      0.00000000     1.00000000
14  virginica           0      0.07954460     0.92045540
15 versicolor           0      1.00000000     0.00000000
16 versicolor           0      1.00000000     0.00000000
17  virginica           0      0.00000000     1.00000000
18     setosa           1      0.00000000     0.00000000
19  virginica           0      0.10028773     0.89971227
20  virginica           0      0.00000000     1.00000000
21 versicolor           0      1.00000000     0.00000000
22 versicolor           0      1.00000000     0.00000000
23  virginica           0      0.00000000     1.00000000
24     setosa           1      0.00000000     0.00000000
25 versicolor           0      0.68610975     0.31389025
26 versicolor           0      0.77887494     0.22112506
27 versicolor           0      1.00000000     0.00000000
28     setosa           1      0.00000000     0.00000000
29     setosa           1      0.00000000     0.00000000
30  virginica           0      0.24798470     0.75201530
31     setosa           1      0.00000000     0.00000000
32 versicolor           0      0.96768116     0.03231884
33 versicolor           0      1.00000000     0.00000000
34     setosa           1      0.00000000     0.00000000
35     setosa           1      0.00000000     0.00000000
36     setosa           1      0.00000000     0.00000000
37  virginica           0      0.00000000     1.00000000
38     setosa           1      0.00000000     0.00000000
39  virginica           0      0.00000000     1.00000000
40     setosa           1      0.00000000     0.00000000
41     setosa           1      0.00000000     0.00000000
42  virginica           0      0.10639202     0.89360798
43 versicolor           0      1.00000000     0.00000000
44     setosa           1      0.00000000     0.00000000
45 versicolor           0      0.96070841     0.03929159
46     setosa           1      0.00000000     0.00000000
47 versicolor           0      1.00000000     0.00000000
48  virginica           0      0.13230826     0.86769174
49  virginica           0      0.00000000     1.00000000
50 versicolor           0      1.00000000     0.00000000
```

使用交叉验证法查看测试数据集分类情况

```
> fit <- fitted(iris.kknn)
> table(iris.valid$Species, fit)

#输出结果为：
             fit
             setosa versicolor virginica
setosa         14          0         0
versicolor      0         15         1
virginica       0          2        18
```

绘制测试数据集的散点矩阵图，其中错误分类样本用红色标注

```
> pcol <- as.character(as.numeric(iris.valid$Species))
> pairs(iris.valid[1:4], pch = pcol, col = c("green3", "red")[(iris.valid$Species != fit)+1])
```

![](/images/2-2-1-1_20171107090129.029.jpeg)
