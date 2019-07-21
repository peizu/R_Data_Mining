# 实验目的

了解关联分析的基本原理，并掌握R语言中实现关联分析的函数方法

# 实验原理

关联规则是表示两个项集之间的关联度或相关性的规则，关联规则的形式为下式：

![](https://kfcoding-static.oss-cn-hangzhou.aliyuncs.com/gitcourse-bigdata/2-2-6-1_20171107091219.019.png)

其中A 与 B是两个互斥的项集，分别位于规则的左侧和右侧。

3个常用于选择有趣关联规则的度量是支持度、置信度和提升度。支持度是数据集中既包含项集A又包含项集B的实例所占的百分比，置信度是包含项集A的实例中也包含项集B所占的百分比，提升度是置信度与包含项集B的实例所占百分比的比值。

本节教程使用R语言中的数据集Titanic作为例子，Titanic数据集是包含在datasets包中的一个4维的数据表，其中包含了泰坦尼克上的乘客的综合信息，包括社会地位、性别、年龄以及是否幸存。为了在该数据集上进行关联规则挖掘，首先对原始数据进行重构，使得数据集中的每一行记录代表一个乘客。

# 算法源码
```
> apriori
function (data, parameter = NULL, appearance = NULL, control = NULL) 
{
    data <- as(data, "transactions")
    items <- data@data
    if (is(appearance, "list")) 
        appearance <- as(c(appearance, list(labels = itemLabels(data))), 
            "APappearance")
    appearance <- as(appearance, "APappearance")
    control <- as(control, "APcontrol")
    parameter <- as(parameter, "APparameter")
    if (control@verbose) {
        cat("Apriori\n")
        cat("\nParameter specification:\n")
        print(parameter)
        cat("\nAlgorithmic control:\n")
        print(control)
    }
    abs_supp <- as.integer(parameter@support * length(data))
    if (control@verbose) {
        cat("\nAbsolute minimum support count:", abs_supp, "\n\n")
    }
    result <- .Call(R_rapriori, items@p, items@i, items@Dim, 
        parameter, control, appearance, data@itemInfo)
    call <- match.call()
    result@info <- list(data = call$data, ntransactions = length(data), 
        support = parameter@support, confidence = parameter@confidence)
    quality(result)$count <- quality(result)$support * length(data)
    if (is(result, "rules")) {
        validObject(result@lhs@data)
        validObject(result@rhs@data)
    }
    else {
        validObject(result@items@data)
    }
    result
}
```

# 实验步骤

载入Titanic数据集，并另命名为df

```
> data("Titanic")
> df <- as.data.frame(Titanic)
```

对数据集df进行重构并保存为titanic.raw数据集，使得该数据集中每一行记录代表一个乘客

```
> titanic.raw <- NULL
> for(i in 1:4){
       titanic.raw <- cbind(titanic.raw,rep(as.character(df[,i]),df$Freq))
  }
> titanic.raw <- as.data.frame(titanic.raw)
> names(titanic.raw) <- names(df)[1:4]
```

调用str\(\)和summary\(\)函数查看数据集titanic.raw结构

```
> dim(titanic.raw)
[1] 2201    4
> summary(titanic.raw)
    Class           Sex          Age   Survived  
 1st :325   Female: 470   Adult:2092   No :1490  
 2nd :285   Male  :1731   Child: 109   Yes: 711  
 3rd :706                                        
 Crew:885
```

安装并加载arules包，并进行关联规则分析

```
> library(arules)
> rules<-apriori(titanic.raw,control=list(verbose=F),parameter=list(minlen=2,supp=0.005,conf=0.8),appearance=list(rhs=c("Survived=No","Survived=Yes"),default="lhs"))
> quality(rules) <- round(quality(rules),digits=3)
> rules.sorted <- sort(rules,by="lift")
> inspect(rules.sorted)
     lhs                                  rhs            support confidence lift 
[1]  {Class=2nd,Age=Child}             => {Survived=Yes} 0.011   1.000      3.096
[2]  {Class=2nd,Sex=Female,Age=Child}  => {Survived=Yes} 0.006   1.000      3.096
[3]  {Class=1st,Sex=Female}            => {Survived=Yes} 0.064   0.972      3.010
[4]  {Class=1st,Sex=Female,Age=Adult}  => {Survived=Yes} 0.064   0.972      3.010
[5]  {Class=2nd,Sex=Female}            => {Survived=Yes} 0.042   0.877      2.716
[6]  {Class=Crew,Sex=Female}           => {Survived=Yes} 0.009   0.870      2.692
[7]  {Class=Crew,Sex=Female,Age=Adult} => {Survived=Yes} 0.009   0.870      2.692
[8]  {Class=2nd,Sex=Female,Age=Adult}  => {Survived=Yes} 0.036   0.860      2.663
[9]  {Class=2nd,Sex=Male,Age=Adult}    => {Survived=No}  0.070   0.917      1.354
[10] {Class=2nd,Sex=Male}              => {Survived=No}  0.070   0.860      1.271
[11] {Class=3rd,Sex=Male,Age=Adult}    => {Survived=No}  0.176   0.838      1.237
[12] {Class=3rd,Sex=Male}              => {Survived=No}  0.192   0.827      1.222
```

可视化查看关联规则结果

```
> library(arulesViz)
> plot(rules,method="graph")
```

![](https://kfcoding-static.oss-cn-hangzhou.aliyuncs.com/gitcourse-bigdata/2-2-6-2_20171107091747.047.jpeg)