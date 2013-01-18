
reactomeList = function() {
library(reactome.db)
mk = mappedkeys(reactomePATHNAME2ID)
mkhuman = mk[grep("^Homo sapiens", mk)]
mkhuman[1] = sub("  ", " ", mkhuman[1]) # bugged NS1 mediated...
mkpath = gsub("Homo sapiens: ", "", mkhuman)
class = toupper( substr(mkpath,1,1) )
grps = list(c(as.character(0:9), LETTERS[1:7]), LETTERS[8:17],
 LETTERS[18:26])
pclasses <<- lapply(grps, function(g) mkpath[which(class %in% g)])
names(pclasses) <<- c("0-9..A-G", "H-Q", "R-Z")
pclasses
}

