library(shiny)
library(org.Hs.eg.db)
library(reactome.db)
library(doParallel)
library(foreach)
registerDoParallel(cores=2)
suppressPackageStartupMessages(library(GGtools))

buildSMLS = function(pop, egid, chr, annodb, numpc, lbmaf) {
          if (pop == "CEU") pack2use = "GGdata"
          else if (pop == "YRI") pack2use = "hmyriB36"
          library(annodb, character.only=TRUE)
          annopref = gsub(".db", "", annodb)
          entrezmap = get(paste(annopref, "ENTREZID", sep=""))
          gsets = split(egid, chr)
          chrinplay = names(gsets)
          get1 = function(x) sapply(x, "[", 1)
          psets = lapply(gsets, function(x) unlist(get1(mget(x, 
              revmap(entrezmap)))))
          cvec = rep(chrinplay, lapply(gsets,length))
          if (numpc > 0)
             smls = lapply(1:length(gsets), function(chr) 
                     clipPCs(getSS(pack2use, chrinplay[chr]), 1:as.integer(numpc))[probeId(psets[[chr]]), ])
          else smls = lapply(1:length(gsets), function(chr) 
                     getSS(pack2use, chrinplay[chr])[probeId(psets[[chr]]), ])
          if (lbmaf > 0)
               smls = lapply(smls, function(x)
                        MAFfilter(x, lower=input$lbMAF))
          smls
          }

getTopHits = function(smls) {
  ns = length(smls)
  foreach(i= 1:ns, .combine=rbind) %dopar% {
    fn = featureNames(smls[[i]])
    ann = annotation(smls[[1]])
    annpref = gsub(".db", "", ann)
    syms = unlist(mget(fn, get(paste(annpref, "SYMBOL", sep="")), ifnotfound=NA))
    tes = eqtlTests(smls[[i]], targdir=tempfile())
    ans = topFeats(probeId(fn), tes, n=1)
    data.frame(probe=fn, chr=names(smList(smls[[i]]))[1], rsid=names(ans), scores=as.numeric(ans), 
     syms=syms,
     ml10p=-log10(1-pchisq(as.numeric(ans),1)), stringsAsFactors=FALSE)
  }
}

setupPlots = function(df, smls) {
  rs = df$rsid
  ac = as.character
  an = as.numeric
  genocodes = lapply(1:length(rs), function(x) ac(as(smList(smls[[x]])[[1]][,rs[x]], "character")))
  pid = df$probe
  expvals = lapply(1:length(rs), function(x) an(exprs(smls[[x]])[pid[x],]))
  splits = lapply(1:length(rs), function(x) split(expvals[[x]], genocodes[[x]]))
  curans = list(splits=splits, genes=pid, rsids=rs, syms=df$syms, genocodes=genocodes, 
      expvals=expvals)
  class(curans) = "plotset"
  curans
}

plot.plotset = function(x, ind){
   boxplot( x$splits[[ind]], xlab=x$rsids[ind], ylab="expr", main=x$syms[ind] )
   points(jitter(as.numeric(factor(x$genocodes[[ind]])), .4), x$expvals[[ind]], col="gray", pch=19)
   }

shinyServer( function(input, output) {
 options(error=recover)

 output$parms = reactiveTable(function() {
   pway <<- " "
# need something to clean these out when pathway not selected
   A = input$pway2useA
   H = input$pway2useH
   R = input$pway2useR
   pwchoices = c(A,H,R)
   okind = which(pwchoices != " ")
   if (length(okind) == 1) pway <<- pwchoices[okind]
   data.frame(population=input$population, 
     type=input$featType,
     gene=input$geneName,
     path = pway,
     cisradius=input$radius, 
     numPC=input$numPC,
     lbMAF=input$lbMAF,
     stringsAsFactors=FALSE)})  # end parms
  
  output$egids = reactiveTable(function() {
    if (input$featType == "gene") {
        egid=get(input$geneName, revmap(org.Hs.egSYMBOL))
        nm=get(egid, org.Hs.egGENENAME)
        chr=get(egid, org.Hs.egCHR)
        sym=get(egid, org.Hs.egSYMBOL)
        ans = data.frame(chr=chr,egid=egid, name=nm, sym=sym, stringsAsFactors=FALSE)
    }
    else if (input$featType == "pathway") {
       pway <<- " "
# need something to clean these out when pathway not selected
       A = input$pway2useA
       H = input$pway2useH
       R = input$pway2useR
       pwchoices = c(A,H,R)
       okind = which(pwchoices != " ")
       if (length(okind) == 1) pway <<- pwchoices[okind]
       syms = chrs = nms = egids = " "
       if (pway != " ") {
        pway = paste("Homo sapiens: ", pway, sep="")
        pwid = get(pway, reactomePATHNAME2ID)
        egids = get(pwid, reactomePATHID2EXTID)
        nms = unlist(mget(egids, org.Hs.egGENENAME,ifnotfound=NA))
        chrs = unlist(mget(egids, org.Hs.egCHR,ifnotfound=NA))
        syms = unlist(mget(egids, org.Hs.egSYMBOL,ifnotfound=NA))
       }
       ans = na.omit(data.frame(chr=chrs, egid=egids, name=nms, sym=syms, stringsAsFactors=FALSE))
       if (input$ready == "yes") {
          smls <- buildSMLS(input$population, ans$egid, ans$chr, "illuminaHumanv1.db", input$numPC, input$lbMAF)
          ans <- getTopHits(smls)
          pl <- setupPlots(ans , smls)
          output$activeControls = reactiveUI(function() {
               selectInput("options", "Gene to show",ans$syms)
               })
          output$pickedGene2 = reactivePlot(function(){ plot(pl, input$plotSelector) })
          }
       ans
       }
  })   # end egids


})
