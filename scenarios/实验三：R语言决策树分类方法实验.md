# 实验目的

了解决策树分类算法的基本原理，并掌握R语言中实现决策树算法的函数方法

# 实验原理

决策树\(Decision Tree\)是一种十分常用的分类算法，是在已知各种情况发生概率的基础上，通过构成决策树来求取净现值的期望值大于等于零的概率，评价项目风险，判断其可行性的决策分析方法，是直观运用概率分析的一种图解法。由于这种决策分支画成图形很像一棵树的枝干，故称决策树。

# 算法源码
```
> ctree
function (formula, data = list(), subset = NULL, weights = NULL, 
    controls = ctree_control(), xtrafo = ptrafo, ytrafo = ptrafo, 
    scores = NULL) 
{
    ls <- dpp(conditionalTree, formula, data, subset, xtrafo = xtrafo, 
        ytrafo = ytrafo, scores = scores)
    fit(conditionalTree, ls, controls = controls, weights = weights)
}

> conditionalTree
An object of class "StatModel"
Slot "name":
[1] "unbiased conditional recursive partitioning"

Slot "dpp":
function (formula, data = list(), subset = NULL, na.action = NULL, 
    xtrafo = ptrafo, ytrafo = ptrafo, scores = NULL, ...) 
{
    dat <- ModelEnvFormula(formula = formula, data = data, subset = subset, 
        designMatrix = FALSE, responseMatrix = FALSE, ...)
    inp <- initVariableFrame(dat@get("input"), trafo = xtrafo, 
        scores = scores)
    response <- dat@get("response")
    if (any(is.na(response))) 
        stop("missing values in response variable not allowed")
    resp <- initVariableFrame(response, trafo = ytrafo, response = TRUE, 
        scores = scores)
    RET <- new("LearningSampleFormula", inputs = inp, responses = resp, 
        weights = rep(1, inp@nobs), nobs = inp@nobs, ninputs = inp@ninputs, 
        menv = dat)
    RET
}
<environment: namespace:party>

Slot "fit":
function (object, controls, weights = NULL, ...) 
{
    if (!extends(class(object), "LearningSample")) 
        stop(sQuote("object"), " is not of class ", sQuote("LearningSample"))
    if (!extends(class(controls), "TreeControl")) 
        stop(sQuote("controls"), " is not of class ", sQuote("TreeControl"))
    if (is.null(weights)) 
        weights <- object@weights
    storage.mode(weights) <- "double"
    if (length(weights) != object@nobs || storage.mode(weights) != 
        "double") 
        stop(sQuote("weights"), " are not a double vector of ", 
            object@nobs, " elements")
    if (max(abs(floor(weights) - weights)) > sqrt(.Machine$double.eps)) 
        stop(sQuote("weights"), " contains real valued elements; currently\n             only integer values are allowed")
    tree <- .Call(R_TreeGrow, object, weights, controls)
    where <- tree[[1]]
    tree <- tree[[2]]
    tree <- prettytree(tree, names(object@inputs@variables), 
        object@inputs@levels)
    RET <- new("BinaryTree")
    RET@tree <- tree
    RET@where <- where
    RET@weights <- weights
    RET@responses <- object@responses
    if (inherits(object, "LearningSampleFormula")) 
        RET@data <- object@menv
    RET@update <- function(weights = NULL) {
        ctreefit(object = object, controls = controls, weights = weights, 
            ...)
    }
    RET@get_where <- function(newdata = NULL, mincriterion = 0, 
        ...) {
        if (is.null(newdata) && mincriterion == 0) {
            if (all(where > 0)) 
                return(where)
        }
        newinp <- newinputs(object, newdata)
        .R_get_nodeID(tree, newinp, mincriterion)
    }
    RET@cond_distr_response <- function(newdata = NULL, mincriterion = 0, 
        ...) {
        wh <- RET@get_where(newdata = newdata, mincriterion = mincriterion)
        response <- object@responses
        if (any(response@is_censored)) {
            swh <- sort(unique(wh))
            RET <- vector(mode = "list", length = length(wh))
            resp <- response@variables[[1]]
            for (i in 1:length(swh)) {
                w <- weights * (where == swh[i])
                RET[wh == swh[i]] <- list(mysurvfit(resp, weights = w))
            }
            return(RET)
        }
        RET <- .Call(R_getpredictions, tree, wh)
        return(RET)
    }
    RET@predict_response <- function(newdata = NULL, mincriterion = 0, 
        type = c("response", "node", "prob"), ...) {
        type <- match.arg(type)
        if (type == "node") 
            return(RET@get_where(newdata = newdata, mincriterion = mincriterion, 
                ...))
        cdresp <- RET@cond_distr_response(newdata = newdata, 
            mincriterion = mincriterion, ...)
        if (type == "prob") 
            return(cdresp)
        if (object@responses@ninputs > 1) 
            return(cdresp)
        response <- object@responses
        if (all(response@is_nominal || response@is_ordinal)) {
            lev <- levels(response@variables[[1]])
            RET <- factor(lev[unlist(lapply(cdresp, which.max))], 
                levels = levels(response@variables[[1]]))
            return(RET)
        }
        if (any(response@is_censored)) {
            RET <- sapply(cdresp, mst)
            return(RET)
        }
        RET <- unlist(cdresp)
        RET <- matrix(unlist(RET), nrow = length(RET), byrow = TRUE)
        if (response@ninputs == 1) 
            colnames(RET) <- names(response@variables)
        return(RET)
    }
    RET@prediction_weights <- function(newdata = NULL, mincriterion = 0, 
        ...) {
        wh <- RET@get_where(newdata = newdata, mincriterion = mincriterion)
        swh <- sort(unique(wh))
        RET <- vector(mode = "list", length = length(wh))
        for (i in 1:length(swh)) RET[wh == swh[i]] <- list(weights * 
            (where == swh[i]))
        return(RET)
    }
    return(RET)
}
```

# 实验步骤

决策树是一种树形结构，其中每个内部节点表示一个属性上的测试，每个分支代表一个测试输出，每个叶节点代表一种类别。下面我们在iris数据集上演示如何使用party包中的ctree\(\)函数来建立一棵决策树。  
调用决策树函数需要加载party包，因此首先安装party包并加载

```
> library(party)
```

加载iris数据集，并生成采样数据，将iris数据集分为两部分，分别为训练数据集trainData和测试数据集testData

```
> data(iris)
> set.seed(1234)
> ind <- sample(2,nrow(iris),replace=TRUE,prob=c(0.7,0.3))
> testData <- iris[ind==2,]
> trainData <- iris[ind==1,]
```

调用ctree\(\)函数，在训练集上建立决策树模型

```
> myFormular <- Species~Sepal.Length+Sepal.Width+Petal.Length+Petal.Width
> iris_ctree <- ctree(myFormular,data=trainData)
```

查看决策树模型在训练数据集上的拟合效果

```
> table(predict(iris_ctree,trainData),trainData$Species)

#执行结果为：
           setosa versicolor virginica
setosa         40          0         0
versicolor      0         37         3
virginica       0          1        31
```

由上表可以看出决策树模型在训练集上的拟合效果，其中错误分类的个数为4，分类正确率为96.43%

绘制生成的决策树

```
> plot(iris_ctree,type="simple")
```

![](/images/Rplot_20180409022530.030.png)
接下来查看决策树模型在测试集上的预测效果

```
> testPred <- predict(iris_ctree,newdata=testData)
> table(testPred,testData$Species)

#执行结果为：
testPred     setosa versicolor virginica
setosa         10          0         0
versicolor      0         12         2
virginica       0          0        14
```

由上表可以看出，决策树模型在测试集上错误分类个数为2，分类正确率为94.74%
