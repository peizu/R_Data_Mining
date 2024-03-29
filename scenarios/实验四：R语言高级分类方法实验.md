# 实验目的

了解支持向量机\(svm\)分类算法的基本原理，并掌握R语言中实现svm算法的函数方法

# 实验原理

支持向量机首先是一种二分类模型，它的基本思想是找到线性空间中的一个超平面，可以将两类别正确分离并使几何间隔最大。当然能线性正确分隔的数据集毕竟很少，如果数据集线性不可分，即存在特异点，除去这些特异点后，数据集是线性可分的，则可以引入松弛变量和惩罚因子，使得分离间隔尽量大的同时误分类的点数目尽量少。对于非线性的分类问题，可以利用核技巧的方法，将数据集映射到高维空间中，在高维空间中训练数据集，找到能使分离间隔较大并且误分类点较少的超平面，这种方法叫做非线性支持向量机。支持向量机中对分离超平面的求解可以形式化为凸优化问题的求解，下面通过调用R语言中e1071包中的函数svm\(\)展示支持向量机。

# 算法源码
```

C_CLASSFICATION <- 0;
EPSILON_SVR <- 3;

is.bigmatrix.refer<-function(x)
{
    return( inherits(x, "BigMatrix.refer") );
}


svm <- function (x, ...)
    UseMethod ("svm")

svm.formula <- function (formula, data = NULL, ..., subset, na.action = na.omit, scale = TRUE)
{
    call <- match.call()
    if (!inherits(formula, "formula"))
        stop("method is only for formula objects")

    m <- match.call(expand.dots = FALSE)
    if (identical(class(eval.parent(m$data)), "matrix"))
        m$data <- as.data.frame(eval.parent(m$data))

    m$... <- NULL
    m$scale <- NULL
    m[[1]] <- as.name("model.frame")
    m$na.action <- na.action
    m <- eval(m, parent.frame())

    Terms <- attr(m, "terms")

    attr(Terms, "intercept") <- 0

    x <- model.matrix(Terms, m)

    y <- model.extract(m, "response")

    attr(x, "na.action") <- attr(y, "na.action") <- attr(m, "na.action")

    if (length(scale) == 1)
        scale <- rep(scale, ncol(x))

    if (any(scale))
    {
        remove <- unique(c(which(labels(Terms) %in%
                                 names(attr(x, "contrasts"))),
                           which(!scale) ) )
        scale <- !attr(x, "assign") %in% remove
    }

    ret <- svm.default (x, y, scale = scale, ..., na.action = na.action)
    ret$call <- call
    ret$call[[1]] <- as.name("svm")
    ret$terms <- Terms
    if (!is.null(attr(m, "na.action")))
        ret$na.action <- attr(m, "na.action")

    class(ret) <- c("svm.formula", class(ret))

    return (ret)
}

check_var_info <- function(x, y, sparse, type, scale, subset, na.action )
{
    x.scale <- y.scale <- NULL
    formula <- inherits(x, "svm.formula")

    nac <- attr(x, "na.action")

    ## scaling, and NA handling
    if (sparse)
    {
        scale <- rep(FALSE, ncol(x))
        if(!is.null(y)) na.fail(y)
        x <- SparseM::t(SparseM::t(x)) ## make shure that col-indices are sorted
    }
    else
    {
        ## na-handling for matrices
        if (!formula)
        {
            if (!missing(subset))
            {
                bigm.subset(x, subset);
                y <- y[subset];
            }

            if (is.null(y))
                nac <- bigm.naction(x, na.action)
            else
            {
                #df <- na.action(data.frame(y, x))
                #y <- df[,1]
                #x <- as.matrix(df[,-1])
                #nac <-
                #    attr(x, "na.action") <-
                #        attr(y, "na.action") <-
                #            attr(df, "na.action")

                sum.row <- rowSums(x) + as.numeric(y);
                sum.row <- na.action(sum.row);
                nac <- attr( y, "na.action") <- attr( sum.row, "na.action");
                bigm.naction( x, na.action, nac );
            }
        }

        ## scaling
        if (length(scale) == 1)
            scale <- rep(scale, ncol(x))

        if (any(scale))
        {
            #co <- !apply( x[, scale, drop = FALSE], 2, var )
            co <- !(unlist(lapply( 1:NCOL(x), function(i) {var(x[,i]);} ))[scale] )
            if (any(co))
            {
                warning(paste("Cannot scale data. Variable(s)", paste(sQuote(colnames(x[,scale,  drop = FALSE])[co]),   sep="", collapse=" and "), "constant. "))
                scale <- rep(FALSE, ncol(x))
            }
            else
            {
                x.scale <- bigm.scale( x, scale);
                # scale Y for regression
                if (is.numeric(y) && (type>2))
                {
                    y <- scale(y)
                    y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
                    y <- as.vector(y)
                }
            }
        }
    }

    ## further parameter checks
    nr <- nrow(x);
    if ( length(y) != nr )
        stop("x and y don't match.")

    if (!is.vector(y) && !is.factor (y) )
        stop("y must be a vector or a factor.")

    #### C/C++ part of rgtTrain requires the Y is sorted by -1 and 1
    #### y.idx will be used later!

    y.org <- y;
    y.idx <- c();
    for( y0 in sort(unique(y) ) )
        y.idx <- c( y.idx, which( y == y0 ) );
    y <- y[y.idx];

    ### x <- x[y.idx,];
    if (sparse)
        x <- x[y.idx,]
    else
        bigm.subset( x, rows=y.idx, cols=NULL );

    lev <- NULL;
    ## in case of classification: transform factors into integers

    if( type < EPSILON_SVR )
    {
        if (is.factor(y))
        {
            lev <- levels(y)
        }
        else
        {
            if(any(as.integer(y) != y))
               stop("dependent variable has to be of factor or integer type for classification mode.")
            y <- as.factor(y)
            lev <- levels(y)
        }
    }

    nclass <- length(lev);

    return( list( nclass = nclass,
                lev   = lev,
                scale = scale,
                x.scale = x.scale,
                y.scale = y.scale,
                y.index = y.idx,
                y.orignal = y.org,
                nr      = nrow(x),
                nac     = nac,
                x       = x,
                y       = y ) );


}

svm.default <- function (x,
          y           = NULL,
          scale       = TRUE,
          type        = "C-classification",
          kernel      = "radial",
          degree      = 3,
          gamma       = if (is.vector(x)) 1 else 1 / ncol(x),
          coef0       = 0,
          cost        = 1,
          class.weights= NULL,
          tolerance   = 0.001,
          epsilon     = 0.1,
          shrinking   = TRUE,
          cross       = 0,
          probability = FALSE,
          fitted      = TRUE,
          rough.cross = 0,
          no.change.x = TRUE,
          verbose     = FALSE,
          ...,
          subset,
          na.action = na.omit)
{
    if (verbose) cat("cost=", cost, " gamma=", gamma, " epsilon=", epsilon, "coef0=", coef0, "degree=", degree, "\n");

    if ((is.vector(x) && is.atomic(x)))
        x <- t(t(x));

    if(inherits(x, "Matrix"))
    {
        requireNamespace("SparseM");
        requireNamespace("Matrix");
        x <- as(x, "matrix.csr");
    }

    if(inherits(x, "simple_triplet_matrix"))
    {
        requireNamespace("SparseM")
        ind <- order(x$i, x$j)
        x <- new("matrix.csr",
                 ra = x$v[ind],
                 ja = x$j[ind],
                 ia = as.integer(cumsum(c(1, tabulate(x$i[ind])))),
                 dimension = c(x$nrow, x$ncol))
    }

    if (sparse <- inherits(x, "matrix.csr") ||  inherits(x, "dgCMatrix" ))
    {
        requireNamespace("SparseM")
    }

    ## NULL parameters?
    if(is.null(degree)) stop(sQuote("degree"), " must not be NULL!")
    if(is.null(gamma)) stop(sQuote("gamma"), " must not be NULL!")
    if(is.null(coef0)) stop(sQuote("coef0"), " must not be NULL!")
    if(is.null(cost)) stop(sQuote("cost"), " must not be NULL!")
    if(is.null(epsilon)) stop(sQuote("epsilon"), " must not be NULL!")
    if(is.null(tolerance)) stop(sQuote("tolerance"), " must not be NULL!")
    if(is.null(shrinking)) stop("shrinking argument must not be NULL!")
    if(is.null(cross)) stop("cross argument must not be NULL!")
    if(is.null(probability)) stop("probability argument must not be NULL!")
    if(is.null(verbose)) stop("verbose argument must not be NULL!")

    ## only support C-classification
    if (is.null(type))
        type <- ifelse ( is.factor(y), "C-classification", "eps-regression");

    type.name <- type;
    type <- C_CLASSFICATION;
    if (type.name != "C-classification" && type.name != "eps-regression")
        stop("Rgtsvm only support C-classification and eps-regression!")
    else
        if(type.name == "eps-regression") type <- EPSILON_SVR;

    if (type != C_CLASSFICATION && length(class.weights) > 0)
    {
        class.weights <- NULL;
        warning(sQuote("class.weights"), " are set to NULL for regression mode. For classification, use a _factor_ for ", sQuote("y"),", or specify the correct ", sQuote("type"), " argument.");
    }

    ## kernel type
    kernel <- pmatch(kernel, c("linear",
                               "polynomial",
                               "radial",
                               "sigmoid"), 99) - 1;

    if (kernel > 10) stop("wrong kernel specification!");
    if (is.null(kernel)) stop("kernel argument must not be NULL!");
    if (is.null(sparse)) stop("sparse argument must not be NULL!");

    if( !sparse && ( class(x) %in% c("matrix", "data.frame") ) )
    {
        if( class(x) == "data.frame" ) x <- as.matrix(x);

        x.backup <- x;
        x <- attach.bigmatrix(x.backup);
    }

    # After this push, x$data will be scaled, so we have to save it to a RDS file(rds.save=TRUE).
    if ( inherits(x, "BigMatrix.refer") && no.change.x)
        bigm.push(x, rds.save=TRUE);

    var.info <- check_var_info( x, y, sparse, type, scale, subset, na.action )

    if ( cross > var.info$nr )
        stop(sQuote("cross"), " cannot exceed the number of observations!")

    biased <- ifelse(var.info$nclass<=2, TRUE, FALSE);

    if ( probability && !fitted) fitted <- TRUE;

    param <- list(type=type, type.name = type.name, kernel=kernel, degree=degree, gamma=gamma,
             coef0=coef0, cost=cost, tolerance=tolerance, epsilon=epsilon,
             shrinking=shrinking, cross=cross, rough.cross=rough.cross,
             sparse=sparse, probability=probability,
             biased = biased, fitted=fitted, nclass = var.info$nclass, class.weights= class.weights,
             verbose = verbose);

    x <- var.info$x;
    y <- var.info$y;
    var.info$y <- NULL;
    x.scale <- var.info$x.scale;
    y.scale <- var.info$y.scale;

    if( type == C_CLASSFICATION)
        cret <- gtsvmtrain.classfication.call( y, x, param, verbose=verbose );
    if( type == EPSILON_SVR)
        cret <- gtsvmtrain.regression.call( y, x, param, verbose=verbose );

    gtsvm.class <- ifelse( cret$nclasses==2, 1, cret$nclasses );
    if (missing(subset)) subset <- NULL;

    ret <- list (
                 call      = match.call(),
                 type.name = type.name,
                 type      = type,
                 kernel    = kernel,
                 cost      = cost,
                 degree    = degree,
                 gamma     = gamma,
                 coef0     = coef0,
                 class.weights = class.weights,

                 tolerance = tolerance,
                 epsilon   = epsilon,
                 sparse    = sparse,
                 scaled    = var.info$scale,
                 x.scale   = var.info$x.scale,
                 y.scale   = var.info$y.scale,
                 biased    = biased,
                 subset    = subset,

                 #number of classes
                 nclasses  = cret$nclasses,
                 levels    = var.info$lev,
                 # total number of sv
                 tot.nSV   = cret$nr,

                 # number of SV in diff. classes
                 nSV       = cret$nSV,

                 # labels of the SVs.
                 labels    = cret$label,

                 # SV matrix
                 SV       = cret$SV,

                 # indexes of sv in x
                 index     = cret$index,

                 ##constants in decision functions
                 rho       = cret$rho,

                 ## coefficiants of sv
                 coefs    = cret$coefs,

                 ## probability
                 compprob  = probability,
                 probA     = NULL,
                 probB     = NULL,

                 totalIter = cret$totalIter,
                 t.elapsed = cret$t.elapsed,
                 na.action = var.info$nac );

    if (cross > 0 )
    {
        if (inherits(x, "BigMatrix.refer") ) bigm.push(x);
        cross.ret <- cross_validation( y, x, param );
        if (inherits(x, "BigMatrix.refer") ) bigm.pop(x);

        ret$cross <- param$cross;
        ret$rough.cross <- param$rough.cross ;

        if ( type > 2)
        {
            scale.factor     <- if (any(scale)) crossprod(y.scale$"scaled:scale") else 1;
            ret$MSE          <- cross.ret$cresults * scale.factor;
            ret$tot.MSE      <- cross.ret$ctotal1  * scale.factor;
            ret$scorrcoeff   <- cross.ret$ctotal2;
        }
        else
        {
            ret$accuracies   <- cross.ret$cresults;
            ret$tot.accuracy <- cross.ret$ctotal1;
        }
    }

    ret$host <-  try(system("hostname", intern = TRUE),silent=T);
    class (ret) <- "gtsvm";

    if (fitted) {
        if(type == C_CLASSFICATION)
        {
            org.idx <- sort.int( var.info$y.index, index.return=T )$ix;

            ret$decision.values <- cret$decision[org.idx]
            ret$fitted <- as.factor(cret$predict[org.idx]);
            levels( ret$fitted ) <- var.info$lev;
            attr(ret$fitted, "decision.values") <- NULL;

            ret$fitted.accuracy <- length( which( ret$fitted == var.info$y.orignal ) )/length(y);

            if(param$nclass==2)
            {
                if(probability)
                {
                    prob <- svc_binary_train_prob( var.info$y.orignal, ret$decision.values );
                    ret$probA = prob$A;
                    ret$probB = prob$B;
                }
            }
            else
            {
                if(probability)
                {
                    org.idx <- sort.int( var.info$y.index, index.return=T )$ix;
                    ret$decision.values <- matrix( cret$decision,ncol=cret$nclasses )[org.idx, ];
                    ret$probA <- svc_one_again_all_train_prob( as.numeric(var.info$y.orignal), ret$decision.values );
                }
            }
        }
        else
        {
            org.idx <- sort.int( var.info$y.index, index.return=T )$ix;

            ret$decision.values <- matrix( cret$predict[ org.idx ], ncol=1) ;
            ret$fitted <- cret$predict[ org.idx ];
            y1 <- y[ org.idx ];

            if(!is.null(y.scale))
            {
                ret$fitted <- ret$fitted * y.scale$"scaled:scale" + y.scale$"scaled:center";
                y1 <- y1 * y.scale$"scaled:scale" + y.scale$"scaled:center";
            }

            y1.v <- ret$fitted;
            ret$fitted.MSE <- sum(( y1 - y1.v )^2, na.rm=T)/length(y1);
            ret$fitted.r2  <- ( length(y1)* sum(y1.v*y1, na.rm=T) - sum(y1.v, na.rm=T)*sum(y1, na.rm=T) )^2 / ( length(y1)*sum(y1.v^2, na.rm=T) - (sum(y1.v, na.rm=T))^2) / ( length(y1)*sum(y1^2, na.rm=T)- (sum(y1, na.rm=T))^2);
            ret$residuals <- (y1 - y1.v);

            if(probability)
               ret$probA <- eps_train_prob( y1, y1.v );

        }
    }

    ### restore x$data
    if(inherits(x, "BigMatrix.refer") && no.change.x) bigm.pop(x);

    ret
}

# Cross-Validation-routine from svm-train in R code
cross_validation <- function ( y, x, param )
{
    y.pre <- rep(NA, length(y));

    idx.shuffle <- sample(1:length(y));
    y <- y[ idx.shuffle ];
    if ( inherits(x, "BigMatrix.refer") ) bigm.subset( x, rows=idx.shuffle, cols=NULL ) else x <- x[idx.shuffle, ];

    breaks <- round(seq(1, length(y), length.out=param$cross+1));
    cresults <- c();

    for(i in 1:(length(breaks)-1))
    {
        idx.cross <- c(breaks[i]:breaks[i+1]);

        if( param$type == C_CLASSFICATION)
        {
            if ( inherits(x, "BigMatrix.refer") )
            {
                bigm.push(x);
                bigm.subset(x, rows = -idx.cross );
                x0 <- x;
            }
            else
                x0 <- x[-idx.cross,];

            sret <- gtsvmtrain.classfication.call( y[ -idx.cross ], x0, param, final.result=TRUE, verbose=FALSE, ignoreNoProgress=TRUE );
            if ( inherits(x, "BigMatrix.refer") ) bigm.pop(x);

            if ( inherits(x, "BigMatrix.refer") )
            {
                bigm.push(x);
                bigm.subset(x, rows = idx.cross );
                x0 <- x;
            }
            else
                x0 <- x[idx.cross,];

            pret <- gtsvmpredict.classfication.call( x0, param$sparse, sret, verbose=FALSE );
            if ( inherits(x, "BigMatrix.refer") ) bigm.pop(x);

            if ( sret$nclasses==2 )
                y.pre [idx.cross] <- levels(y)[as.factor(pret$ret)]
            else
            {
                ret2 <- matrix( pret$dec[ 1:(length(idx.cross)*param$nclass) ],
                                nrow = length(idx.cross), ncol= param$nclass  );
                y.pre [idx.cross] <- apply(ret2, 1, which.max);
            }

            cresults[i] = 100.0 * sum( as.integer(y.pre [idx.cross] == y[idx.cross]) )/length(idx.cross);
        }


        if( param$type == EPSILON_SVR)
        {
            if ( inherits(x, "BigMatrix.refer") )
            {
                bigm.push(x);
                bigm.subset(x, rows = -idx.cross );
                x0 <- x;
            }
            else
                x0 <- x[-idx.cross,];

            sret <- gtsvmtrain.regression.call( y[ -idx.cross ], x0, param, final.result=TRUE, verbose=FALSE, ignoreNoProgress=TRUE );
            if ( inherits(x, "BigMatrix.refer") ) bigm.pop(x);

            if (inherits(x, "BigMatrix.refer") )
            {
                bigm.push(x);
                bigm.subset(x, rows = idx.cross );
                x0 <- x;
            }
            else
                x0 <- x[idx.cross,];

            pret <- gtsvmpredict.regression.call( x0, param$sparse, sret, verbose=FALSE);
            if (inherits(x, "BigMatrix.refer") ) bigm.pop(x);

            y.pre [idx.cross] <- pret$ret;
            cresults[i] = sum(( y.pre[idx.cross] - y[idx.cross])^2 )/length(idx.cross)
        }

        if( param$rough.cross >0 && param$rough.cross==i)
            break;

    }

    if( any( which( is.na(y.pre) | is.na(y) ) ) )
    {
        y.sel <- !is.na(y.pre) & !is.na(y);
        y.pre <- y.pre[y.sel];
        y     <- y[y.sel];
    }

    if(param$type == EPSILON_SVR )
    {
        # MSE
        ctotal1 <- sum((y.pre-y)^2)/length(y);
        # R2
        ctotal2 <- (length(y)* sum(y.pre*y) - sum(y.pre)*sum(y) )^2 / ( length(y)*sum(y.pre^2) - (sum(y.pre))^2) / ( length(y)*sum(y^2)- (sum(y))^2);
    }
    else
    {
        ctotal1 <- 100.0 * sum(y.pre==y)/length(y);
        ctotal2 <- NA;
    }

    return(list(ctotal1=ctotal1, ctotal2=ctotal2, cresults=cresults));
}


predict.gtsvm <- function (object, newdata,
          decision.values = FALSE,
          probability = FALSE,
          verbose = FALSE,
          ...,
          na.action = na.omit)
{
    if (missing(newdata))
        return(fitted(object))

    if (object$tot.nSV < 1)
        stop("Model is empty!");

    if(inherits(newdata, "Matrix"))
    {
        requireNamespace("SparseM");
        requireNamespace("Matrix");
        newdata <- as(newdata, "matrix.csr");
    }

    if(inherits(newdata, "simple_triplet_matrix"))
    {
       requireNamespace("SparseM");
       ind <- order(newdata$i, newdata$j);
       newdata <- new("matrix.csr",
                      ra = newdata$v[ind],
                      ja = newdata$j[ind],
                      ia = as.integer(cumsum(c(1, tabulate(newdata$i[ind])))),
                      dimension = c(newdata$nrow, newdata$ncol));
    }

    sparse <- inherits(newdata, "matrix.csr");
    if (object$sparse || sparse)
        requireNamespace("SparseM");

    act <- NULL;
    if ((is.vector(newdata) && is.atomic(newdata)))
        newdata <- t(t(newdata));

    if (sparse)
        newdata <- SparseM::t(SparseM::t(newdata));

    preprocessed <- !is.null(attr(newdata, "na.action"))
    rowns <- if (!is.null(rownames(newdata)))
                rownames(newdata)
            else
                1:nrow(newdata);

    if( inherits(newdata, "BigMatrix.refer") ) bigm.push(newdata);

    if (!object$sparse)
    {
        if (inherits(object, "svm.formula"))
        {
            if(is.null(colnames(newdata)))
                colnames(newdata) <- colnames(object$SV);

            newdata <- na.action(newdata);
            act <- attr(newdata, "na.action");
            newdata <- model.matrix(delete.response(terms(object)),
                                    as.data.frame(newdata));
        }
        else
        {
            if( !inherits(newdata, "BigMatrix.refer") )
            {
                newdata <- na.action(as.matrix(newdata));
                act <- attr(newdata, "na.action");
            }
            else
            {
                act <- bigm.naction( newdata, na.action );
            }
        }
    }

    if (!is.null(act) && !preprocessed)
        rowns <- rowns[-act];

    if (any(object$scaled))
    {
        if( !inherits( newdata, "BigMatrix.refer") )
             newdata[,object$scaled] <-
                scale( newdata[, object$scaled, drop = FALSE], center = object$x.scale$"scaled:center", scale  = object$x.scale$"scaled:scale")
        else
            bigm.scale( newdata, object$scaled, center = object$x.scale$"scaled:center", scale  = object$x.scale$"scaled:scale" );
    }

    if (ncol(object$SV) != ncol(newdata))
        stop ("test data does not match model !");

    param <- list( decision.values = decision.values, probability = probability );

    # Call C/C++ interface to do predict
    if(object$type == C_CLASSFICATION)
        ret <- gtsvmpredict.classfication.call( newdata, sparse, object, param, verbose=verbose)
    else if(object$type == EPSILON_SVR)
        ret <- gtsvmpredict.regression.call( newdata, sparse, object, param, verbose=verbose)
    else
        stop("only 'C-classification' and 'eps-regression' are implemented in this package!");

    ret2 <- ret$ret;
    if (is.character(object$levels)) # classification: return factors
    {
        ret2 <- as.factor(ret$ret);
        levels( ret2 ) <- object$levels;
    }
    else if (any(object$scaled) && !is.null(object$y.scale)) # return raw values, possibly scaled back
        ret2 <- ret$ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center";

    names(ret2) <- rowns;
    ret2 <- napredict(act, ret2);

    if (decision.values)
    {
        colns = c();
        #for (i in 1:(object$nclasses - 1))
        #    for (j in (i + 1):object$nclasses)
        #        colns <- c(colns,
        #                   paste(object$levels[object$labels[i]],
        #                         "/", object$levels[object$labels[j]],
        #                         sep = ""));

        if(object$nclasses==2)
            attr(ret2, "decision.values") <- napredict(act, matrix(ret$dec, nrow = nrow(newdata), ncol=1 ) )
        else
        {
            colns <- object$labels;
            attr(ret2, "decision.values") <-
                napredict(act, matrix(ret$dec, nrow = nrow(newdata), ncol=length(colns), byrow = FALSE, dimnames = list(rowns, colns) ) );
        }

    }

    if (probability && object$type < 2) {
        if (!object$compprob)
            warning("SVM has not been trained using `probability = TRUE`, probabilities not available for predictions.")
        else
        {
            if(object$nclasses==2)
                attr(ret2, "probabilities") <- napredict(act, svc_binary_predict_prob( ret$dec, object$probA, object$probB ) )
            else
            {
                colns <- object$labels;
                ydeci <- matrix(ret$dec, nrow = nrow(newdata), ncol=length(colns), byrow = FALSE, dimnames = list(rowns, colns) );
                ret$prob <- svc_one_again_all_predict_prob( ydeci, object$probA );
                colnames(ret$prob) <- object$labels;
                attr(ret2, "probabilities") <- napredict(act,ret$prob );
             }
        }
    }

    if( inherits(newdata, "BigMatrix.refer") ) bigm.pop(newdata);

    ret2
}

print.gtsvm <- function (x, ...)
{
    cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep="\n");
    cat("Parameters:\n");
    cat("   SVM-Type: ", x$type.name, "\n");
    cat(" SVM-Kernel: ", c("linear",
                           "polynomial",
                           "radial",
                           "sigmoid")[x$kernel+1], "\n");

    cat("       cost: ", x$cost, "\n");
    if (x$kernel==1)
        cat("     degree: ", x$degree, "\n");
    cat("      gamma: ", x$gamma, "\n");
    if (x$kernel==1 || x$kernel==3)
        cat("     coef.0: ", x$coef0, "\n");

    cat("    tolerance: ", x$tolerance, "\n");
    cat(" time elapsed: ", x$t.elapsed[3], "\n\n");

    cat("\nNumber of Support Vectors: ", x$tot.nSV);
    cat("\n\n");

}

summary.gtsvm <- function(object, ...)
    structure(object, class="summary.gtsvm")

print.summary.gtsvm <- function (x, ...)
{
    print.gtsvm(x);

    cat(" (", x$nSV, ")\n\n");
    cat("\nNumber of Classes: ", x$nclasses, "\n\n");
    cat("Levels:", if(is.numeric(x$levels)) "(as integer)", "\n", x$levels);

    cat("\n\n");
}

#scale.data.frame <- function(x, center = TRUE, scale = TRUE)
#{
#    i <- sapply(x, is.numeric)
#    if (ncol(x[, i, drop = FALSE])) {
#        x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], na.omit(center), na.omit(scale))
#        if(center || !is.logical(center))
#            attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
#        if(scale || !is.logical(scale))
#            attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
#    }
#    x
#}

plot.gtsvm <- function(x, data, formula = NULL, fill = TRUE,
         grid = 50, slice = list(), symbolPalette = palette(),
         svSymbol = "x", dataSymbol = "o", ...)
{
    if (is.null(formula) && ncol(data) == 3)
    {
        formula <- formula(delete.response(terms(x)));
        formula[2:3] <- formula[[2]][2:3];
    }

    if (is.null(formula))
        stop("missing formula.");

    if (fill)
    {
        sub <- model.frame(formula, data);
        xr <- seq(min(sub[, 2]), max(sub[, 2]), length = grid);
        yr <- seq(min(sub[, 1]), max(sub[, 1]), length = grid);
        l <- length(slice);

        if (l < ncol(data) - 3)
        {
            slnames <- names(slice);
            slice <- c(slice, rep(list(0), ncol(data) - 3 - l));
            names <- labels(delete.response(terms(x)));
            names(slice) <- c(slnames, names[!names %in% c(colnames(sub), slnames)]);
        }

        for (i in names(which(sapply(data, is.factor))))
            if (!is.factor(slice[[i]])) {
                levs <- levels(data[[i]]);
                lev <- if (is.character(slice[[i]])) slice[[i]] else levs[1];
                fac <- factor(lev, levels = levs);
                if (is.na(fac))
                   stop(paste("Level", dQuote(lev), "could not be found in factor", sQuote(i)));
                slice[[i]] <- fac;
            }

        lis <- c(list(yr), list(xr), slice);
        names(lis)[1:2] <- colnames(sub);
        new <- expand.grid(lis)[, labels(terms(x))];
        preds <- predict(x, new);

        filled.contour(xr, yr,
               matrix(as.numeric(preds),
               nrow = length(xr), byrow = TRUE),
               plot.axes = {
                               axis(1)
                               axis(2)
                               colind <- as.numeric(model.response(model.frame(x, data)))
                               dat1 <- data[-x$index,]
                               dat2 <- data[x$index,]
                               coltmp1 <- symbolPalette[colind[-x$index]]
                               coltmp2 <- symbolPalette[colind[x$index]]
                               points(formula, data = dat1, pch = dataSymbol, col = coltmp1)
                               points(formula, data = dat2, pch = svSymbol, col = coltmp2)
                           },
               levels = 1:(length(levels(preds)) + 1),
               key.axes = axis(4, 1:(length(levels(preds))) + 0.5,
               labels = levels(preds),
               las = 3),
               plot.title = title(main = "SVM classification plot",
               xlab = names(lis)[2], ylab = names(lis)[1]), ...);
    }
    else {
        plot(formula, data = data, type = "n", ...);
        colind <- as.numeric(model.response(model.frame(x, data)));

        dat1 <- data[-x$index,];
        dat2 <- data[x$index,];
        coltmp1 <- symbolPalette[colind[-x$index]];
        coltmp2 <- symbolPalette[colind[x$index]];
        points(formula, data = dat1, pch = dataSymbol, col = coltmp1);
        points(formula, data = dat2, pch = svSymbol, col = coltmp2);

        invisible();
    }
}

predict.batch <- function (object, file.rds, decision.values = TRUE, probability = FALSE, verbose = FALSE, ..., na.action = na.omit)
{
    if (missing(file.rds))
        stop("No RDS files are specified.\n");

    if (object$tot.nSV < 1)
        stop("Model is empty!");

    if (object$sparse)
        requireNamespace("SparseM");

    x.count <- 0;
    rowns <- c();
    for(rds in file.rds)
    {
        newdata <- readRDS(rds);
        if (NCOL(newdata) != NCOL(object$SV))
            stop(paste("X data in RDS file doesn't have same number of feature vectors. File=", rds, sep=""));
        x.count <- x.count + NROW(newdata);
        rowns <- c( rowns, rownames(newdata));

        newdata <- na.action( as.matrix(newdata) );
        act <- attr(newdata, "na.action");
        if(!is.null(act))
            stop( paste( "Missing values in RDS file, file=", rds ) );
    }

    param <- list( decision.values = decision.values, probability = probability );

    # Call C/C++ interface to do predict
    if(object$type == C_CLASSFICATION)
        ret <- gtsvmpredict.classfication.batch.call( file.rds, x.count, object, param, verbose=verbose)
    else if(object$type == EPSILON_SVR)
        ret <- gtsvmpredict.regression.batch.call( file.rds, x.count, object, param, verbose=verbose)
    else
        stop("only 'C-classification' and 'eps-regression' are implemented in this package!");

    ret2 <- ret$ret;
    if (is.character(object$levels)) # classification: return factors
    {
        ret2 <- as.factor(ret$ret) ;
        levels( ret2 ) <- object$levels;
    }
    else if (any(object$scaled) && !is.null(object$y.scale)) # return raw values, possibly scaled back
        ret2 <- ret$ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center";

    names(ret2) <- rowns;
    act <- NULL;

    if (decision.values)
    {
        colns = c();
        for (i in 1:(object$nclasses - 1))
            for (j in (i + 1):object$nclasses)
                colns <- c(colns,
                           paste(object$levels[object$labels[i]],
                                 "/", object$levels[object$labels[j]],
                                 sep = ""));

        attr(ret2, "decision.values") <-
            napredict(act, matrix(ret$dec, nrow = x.count, ncol=length(colns), byrow = TRUE, dimnames = list(rowns, colns) ) );
    }

    if (probability && object$type < 2)
    {
        if (!object$compprob)
            warning("SVM has not been trained using `probability = TRUE`, probabilities not available for predictions.")
        else
        {

            if(object$nclasses==2)
                attr(ret2, "probabilities") <- napredict(act, svc_binary_predict_prob( ret$dec, object$probA, object$probB ) )
            else
            {
                colns <- object$labels;
                ydeci <- matrix(ret$dec, nrow = nrow(newdata), ncol=length(colns), byrow = FALSE, dimnames = list(rowns, colns) );
                ret$prob <- svc_one_again_all_predict_prob( ydeci, object$probA );
                colnames(ret$prob) <- object$labels;
                attr(ret2, "probabilities") <- napredict(act,ret$prob );
             }
        }
    }

    ret2
}
```

# 实验步骤

首先安装e1071包并加载

```
> library(e1071)
```

导入iris数据集

```
> data(iris)
```

首先通过散点图观察数据集的大致分布情况，此步需要加载lattice包，因此先安装

```
> library(lattice)
> xyplot(Petal.Length ~ Petal.Width, data = iris, groups = Species,  auto.key=list(corner=c(1,0)))
```

![](/images/2-2-4-1_20171107090711.011.jpeg)

调用svm函数并观察分类后的超平面

```
> svm_model <- svm(Species~Petal.Length+Petal.Width,data=iris)
> plot(svm_model,iris,Petal.Length~Petal.Width)
```

![](/images/2-2-4-2_20171107090804.004.jpeg)
