biplot.default2 <-
  function(x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
           xlabs = NULL, ylabs = NULL, expand=1, xlim = NULL, ylim = NULL,
           arrow.len = 0.1, pt.cex=1, axis.text.cex=1, draw.box=T,
           main = NULL, sub = NULL, xlab = NULL, ylab = NULL, textaxis=T,textaxisbg="#b2f4f4c0",textaxisborder="black",...) {
    n <- nrow(x)
    p <- nrow(y)
    if(missing(xlabs)) {
      xlabs <- dimnames(x)[[1L]]
      if(is.null(xlabs)) xlabs <- 1L:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
    if(missing(ylabs)) {
      ylabs <- dimnames(y)[[1L]]
      if(is.null(ylabs)) ylabs <- paste("Var", 1L:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2L]])
    
    if(length(cex) == 1L) cex <- c(cex, cex)
    if(missing(col)) {
      col <- par("col")
      if (!is.numeric(col)) col <- match(col, palette(), nomatch=1L)
      col <- c(col, col + 1L)
    }
    else if(length(col) == 1L) col <- c(col, col)
    
    unsigned.range <- function(x)
      c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
    rangx1 <- unsigned.range(x[, 1L])
    rangx2 <- unsigned.range(x[, 2L])
    rangy1 <- unsigned.range(y[, 1L])
    rangy2 <- unsigned.range(y[, 2L])
    
    if(missing(xlim) && missing(ylim))
      xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if(missing(xlim)) xlim <- rangx1
    else if(missing(ylim)) ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")
    if(!is.null(main))
      op <- c(op, par(mar = par("mar")+c(0,0,1,0)))
    plot(x, type = "n", xlim = xlim, ylim = ylim, col = col[1L],
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    # text(x, xlabs, cex = cex[1L], col = col[1L], ...)
    points(x, col = col[1L], cex=pt.cex,...)
    par(new = TRUE)
    dev.hold(); on.exit(dev.flush(), add = TRUE)
    plot(y, axes = FALSE, type = "n", xlim = xlim*ratio, ylim = ylim*ratio,
         xlab = "", ylab = "", col = col[1L], ...)
    # axis(3, col = col[2L], ...)
    # axis(4, col = col[2L], ...)
    if(draw.box) box(col = col[1L])
    if(textaxis){
      if(!is.null(textaxisbg)){
        boxtext(y[,1],y[,2], labels=ylabs, cex = axis.text.cex, col = col[2L], col.bg = textaxisbg, border.bg = textaxisborder, ...)
      }else{
        text(y, labels=ylabs, cex = axis.text.cex, col = col[2L], ...)
      }
    }
    if(var.axes)
      arrows(0, 0, y[,1L] * 0.8, y[,2L] * 0.8, col = col[2L], length=arrow.len)
    invisible()
  }

boxtext <- function(x, y, labels = NA, col.text = NULL, col.bg = NA, 
                    border.bg = NA, adj = NULL, pos = NULL, offset = 0.5, 
                    padding = c(0.5, 0.5), cex = 1, font = graphics::par('font'),...){
  
  ## The Character expansion factro to be used:
  theCex <- graphics::par('cex')*cex
  
  ## Is y provided:
  if (missing(y)) y <- x
  
  ## Recycle coords if necessary:    
  if (length(x) != length(y)){
    lx <- length(x)
    ly <- length(y)
    if (lx > ly){
      y <- rep(y, ceiling(lx/ly))[1:lx]           
    } else {
      x <- rep(x, ceiling(ly/lx))[1:ly]
    }       
  }
  
  ## Width and height of text
  textHeight <- graphics::strheight(labels, cex = theCex, font = font)
  textWidth <- graphics::strwidth(labels, cex = theCex, font = font)
  
  ## Width of one character:
  charWidth <- graphics::strwidth("e", cex = theCex, font = font)
  
  ## Is 'adj' of length 1 or 2?
  if (!is.null(adj)){
    if (length(adj == 1)){
      adj <- c(adj[1], 0.5)            
    }        
  } else {
    adj <- c(0.5, 0.5)
  }
  
  ## Is 'pos' specified?
  if (!is.null(pos)){
    if (pos == 1){
      adj <- c(0.5, 1)
      offsetVec <- c(0, -offset*charWidth)
    } else if (pos == 2){
      adj <- c(1, 0.5)
      offsetVec <- c(-offset*charWidth, 0)
    } else if (pos == 3){
      adj <- c(0.5, 0)
      offsetVec <- c(0, offset*charWidth)
    } else if (pos == 4){
      adj <- c(0, 0.5)
      offsetVec <- c(offset*charWidth, 0)
    } else {
      stop('Invalid argument pos')
    }       
  } else {
    offsetVec <- c(0, 0)
  }
  
  ## Padding for boxes:
  if (length(padding) == 1){
    padding <- c(padding[1], padding[1])
  }
  
  ## Midpoints for text:
  xMid <- x + (-adj[1] + 1/2)*textWidth + offsetVec[1]
  yMid <- y + (-adj[2] + 1/2)*textHeight + offsetVec[2]
  
  ## Draw rectangles:
  rectWidth <- textWidth + 2*padding[1]*charWidth
  rectHeight <- textHeight + 2*padding[2]*charWidth    
  graphics::rect(xleft = xMid - rectWidth/2, 
                 ybottom = yMid - rectHeight/2, 
                 xright = xMid + rectWidth/2, 
                 ytop = yMid + rectHeight/2,
                 col = col.bg, border = border.bg, lwd=0.7)
  
  ## Place the text:
  graphics::text(xMid, yMid, labels, col = col.text, cex = theCex, font = font, 
                 adj = c(0.5, 0.5))    
  
  ## Return value:
  if (length(xMid) == 1){
    invisible(c(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                yMid + rectHeight/2))
  } else {
    invisible(cbind(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                    yMid + rectHeight/2))
  }    
}


biplot.prcomp2 <- function(x, choices = 1L:2L, scale = 1, pc.biplot=FALSE, ...)
{
  if(length(choices) != 2L) stop("length of choices must be 2")
  if(!length(scores <- x$x))
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
         domain = NA)
  if(is.complex(scores))
    stop("biplots are not defined for complex PCA")
  lam <- x$sdev[choices]
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  biplot.default2(t(t(scores[, choices]) / lam),
                 t(t(x$rotation[, choices]) * lam),...)
  invisible()
}
