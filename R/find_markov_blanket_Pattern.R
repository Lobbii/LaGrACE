#' Find the Markov blanket around a target variable.
#'
#' This function returns the Markov blanket around a target vertex from a
#' graph
#' @param target A string containing the name of the target vertex
#' @param graph A list containing graph information
#' @return A list of the names of vertices in the Markov blanket
#' @export


find_markov_blanket_Pattern<-function(target,graph){
  edges <- as.data.frame(graph$edges,stringsAsFactors=F)

  parent.e<-edges[grep(paste0('^',target,'$'),edges[,2]),]
  if (nrow(parent.e) ==0) {
    parent.v<-c()
  } else {
    parent.v<-parent.e[,1]
    # print(c('parent node:'))
    # print(parent.e)
  }

  children.e<-edges[grep(paste0('^',target,'$'),edges[,1]),]
  if (nrow(children.e) ==0) {
    children.v<-c()
    co.parent.v<-c()
  } else {
    children.v<-children.e[,2]
    # print(c('children node:'))
    # print(children.e)
    co.parent.v<-c()
    for(i in 1:length(children.v)){
      if(children.e[i,3] == '-->'){
        # directed edge
        co.parent.e<-edges[grep(paste0('^',children.v[i],'$'),edges[,2]),]
        co.parent.e<-co.parent.e[co.parent.e[,3]=='-->',]
        co.parent.v<-c(co.parent.v,co.parent.e[,1])
        # print(c('coparent node',children.v[i]))
        # print(co.parent.e)
      }
    }
  }

  MB<-unique(c(parent.v,children.v,co.parent.v))
  MB<-setdiff(MB,target)
  return(MB)}


