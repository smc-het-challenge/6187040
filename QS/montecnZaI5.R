library(gtools)

##################################################################

library(data.table)
#library(ggplot2)
#library(emdbook)

############################## genere W test ###############################

GW<-function(C){
  a<-array(c(1), dim = (C+1))
  W<-rdirichlet(1,a)
  #W[1]<-W[1]*rbeta(1,1,100)
  
  Pw<-log(ddirichlet(W,a))
  W<-t(W)
  
  l<-list("W"=W,"Pw"=Pw)
  return(l)
}

############################################################################

################################# genere C #################################

GC<-function(esp){
  
  C<-rgeom(1,1/esp)+1
  #C<-rpois(1,esp)
  return(C)
  
}

GCt<-function(C,esp){
  
  Ct<-GC(esp)
  
  Pif <- dgeom(Ct-1,1/esp)
  Pfi <- dgeom(C-1,1/esp)
  
  # Pif <- dpois(Ct,esp)
  # Pfi <- dpois(C,esp)
  
  l<-list("C"=Ct,"Pif"=Pif, "Pfi"=Pfi)
  return(l)
}

############################################################################

GetA<-function(data,lZ,lW,datat,lCt,lZt,lWt){
  
  n<-dim(data$data)[1]
  C<-dim(lZ$Z)[2]-1
  vect<-matrix(c(1),nrow=C+1,ncol=1)
  a<-lZ$Z%*%vect
  vect<-matrix(c(1),nrow=dim(lZt$Z)[2],ncol=1)
  b<-lZt$Z%*%vect
  P1<-lCt$Pif+sum(log(dbinom(data$data[a[,1]!=0,1],data$data[a[,1]!=0,2],lZ$Z[a[,1]!=0,]%*%lW$W/(lZ$Zal[a[,1]!=0,]%*%lW$W))))+lZ$Pz+lW$Pw+data$prob
  P2<-lCt$Pfi+sum(log(dbinom(datat$data[b[,1]!=0,1],datat$data[b[,1]!=0,2],lZt$Z[b[,1]!=0,]%*%lWt$W/(lZt$Z[b[,1]!=0,]%*%lWt$W))))+lZt$Pz+lWt$Pw+datat$prob
  
  A<-P2-P1
  
  return(A)
  
}

allele<-function(data,ligne){
  
  a<-matrix(c(0,2),nrow=1,ncol=2)
  
  n<-dim(data)[2]
  
  for (i in seq(3,n,2)){
    
    if (data[ligne,i]!=0 & data[ligne,i+1]!=0 & length(a[a[,1]==data[ligne,i] & a[,2]==data[ligne,i+1],])==0){
      
      a<-rbind(a,data[ligne,c(i,i+1)])
      
      if (length(a[a[,1]==1 & a[,2]==data[ligne,i+1],])==0){
        
        a<-rbind(a,c(1,data[ligne,i+1]))
        
      }
      
      if (length(a[a[,1]==0 & a[,2]==data[ligne,i+1],])==0){
        
        a<-rbind(a,c(0,data[ligne,i+1]))
        
      }
      
    }
    
  }
  
  if(dim(a)[1]>1){
    a<-a[order(a[,2],decreasing = FALSE),]
    a<-a[order(a[,1],decreasing = TRUE),]
    a<-a[order(a[,1]/a[,2],decreasing = TRUE),]
  }
  
  return(a)
  
}

PostZ<-function(C,data,lW){
  
  n<-dim(data)[1]
  Z<-matrix(c(1),nrow=n,ncol=(C+1))
  Zal<-matrix(c(0),nrow=n,ncol=(C+1))
  
  for (i in 1:n){
    
    a<-allele(data,i)
    Z[i,]<-a[1,1]
    Zal[i,]<-a[1,2]
    
  }
  
  Zal[,1]<-2
  Z[,1]<-0
  #Z[,2]<-0
  
  E<-abs(Z%*%lW$W/(Zal%*%lW$W)-data[,1]/data[,2])
  
  ordre<-1+order(lW$W[-1],decreasing = TRUE)
  
  for (i in 1:n){
    
    a<-NULL
    a<-allele(data,i)
    
    for (j in 1:C){
      
      modif<-Z[i,ordre[j]]
      modifal<-Zal[i,ordre[j]]
      
      for (k in 1:dim(a)[1]){
        
        Z[i,ordre[j]]<-a[k,1]
        Zal[i,ordre[j]]<-a[k,2]
        
        e<-abs((Z[i,]%*%lW$W/(Zal[i,]%*%lW$W)-data[i,1]/data[i,2]))
        
        if (e<E[i]){
          
          E[i]<-e
          modif<-Z[i,ordre[j]]
          modifal<-Zal[i,ordre[j]]
          
        }
        
      }
      
      Z[i,ordre[j]]<-modif
      Zal[i,ordre[j]]<-modifal
      
    }
    
  }
  
  # for (i in 1:n){
  #   
  #   if (sum(Z[i,])==0){Z[i,order(lW$W[-1],decreasing = FALSE)+1]<-1}
  #   
  # }
  
  l<-list("Z"=Z,"Zal"=Zal,"Pz"=0)
  return(l)
}

Monte<-function(data,esp,np,nb){
  
  l<-bruit(data)
  C<-GC(esp)
  lW<-GW(C)
  lZ<-PostZ(C,l$data,lW)
  
  h<-array(dim=np)
  v<-array(dim=np)
  
  for (i in 1:nb){
    
    lt<-bruit(data)
    lCt<-GCt(C,esp)
    lWt<-GW(lCt$C)
    lZt<-PostZ(lCt$C,lt$data,lWt)
    
    A<- GetA(l,lZ,lW,lt,lCt,lZt,lWt)
    print(A)
    
    if (is.nan(A)){
      
      print(lZt$Z)
      print(lWt$W)
      
    }
    
    if (A>0){
      
      C<-lCt$C
      lZ<-lZt
      lW<-lWt
      l<-lt
      
    }
    
    else{
      
      if (runif(1,exp(A))==1){
        
        C<-lCt$C
        lZ<-lZt
        lW<-lWt
        l<-lt
        
      }
      
    }
    
  }
  
  for (i in 1:np){
    
    lt<-bruit(data)
    
    lCt<-GCt(C,esp)
    lWt<-GW(lCt$C)
    lZt<-PostZ(lCt$C,lt$data,lWt)
    
    A<- GetA(l,lZ,lW,lt,lCt,lZt,lWt)
    print(A)
    
    if (is.nan(A)){
      
      print(lZt$Z)
      print(lWt$W)
      
    }
    
    if (A>0){
      
      C<-lCt$C
      lZ<-lZt
      lW<-lWt
      l<-lt
      
    }
    
    else{
      
      if (runif(1,exp(A))==1){
        
        C<-lCt$C
        lZ<-lZt
        lW<-lWt
        l<-lt
        
      }
      
    }
    
    h[i]<-C
    v[i]<-lCt$C
    
  }
  
  l<-bruit(data)
  
  C<-strtoi(names(sort(table(h),decreasing = TRUE)))[1]
  lW<-GW(C)
  lZ<-PostZ(C,l$data,lW)
  
  for (i in 1:np){
    
    lt<-bruit(data)
    
    lWt<-GW(C)
    lZt<-PostZ(C,lt$data,lWt)
    a<-lZ$Z%*%matrix(c(1),nrow=dim(lZ$Z)[2],ncol=1)
    b<-lZt$Z%*%matrix(c(1),nrow=dim(lZ$Z)[2],ncol=1)
    A<-sum(log(dbinom(lt$data[b[,1]!=0,1],lt$data[b[,1]!=0,2],(lZt$Z[b[,1]!=0,]%*%lWt$W/(lZt$Zal[b[,1]!=0,]%*%lWt$W)))))-sum(log(dbinom(l$data[a[,1]!=0,1],l$data[a[,1]!=0,2],lZ$Z[a[,1]!=0,]%*%lW$W/(lZ$Zal[a[,1]!=0,]%*%lW$W))))+lt$prob-l$prob
    print(A)
    
    if (is.nan(A)){
      
      print(lZt$Z)
      print(lWt$W)
      
    }
    
    if (A>0){
      
      lZ<-lZt
      lW<-lWt
      
    }
    
  }
  
  ret<-list("Z"=lZ$Z,"Zal"=lZ$Zal,"W"=lW$W,"data"=l$data)
  return(ret)
  
}

dataprod<-function(table){
  
  data<-matrix(nrow=dim(table)[1],ncol=26)
  
  data[,1]<-table$tum_alt_allele_depth
  data[,2]<-table$tum_depth
  data[,3]<-table$nMaj1_A
  data[,4]<-table$nMaj1_A+table$nMin1_A
  data[,5]<-table$nMaj2_A
  data[,6]<-table$nMaj2_A+table$nMin2_A
  data[,7]<-table$nMaj1_B
  data[,8]<-table$nMaj1_B+table$nMin1_B
  data[,9]<-table$nMaj2_B
  data[,10]<-table$nMaj2_B+table$nMin2_B
  data[,11]<-table$nMaj1_C
  data[,12]<-table$nMaj1_C+table$nMin1_C
  data[,13]<-table$nMaj2_C
  data[,14]<-table$nMaj2_C+table$nMin2_C
  data[,15]<-table$nMaj1_D
  data[,16]<-table$nMaj1_D+table$nMin1_D
  data[,17]<-table$nMaj2_D
  data[,18]<-table$nMaj2_D+table$nMin2_D
  data[,19]<-table$nMaj1_E
  data[,20]<-table$nMaj1_E+table$nMin1_E
  data[,21]<-table$nMaj2_E
  data[,22]<-table$nMaj2_E+table$nMin2_E
  data[,23]<-table$nMaj1_F
  data[,24]<-table$nMaj1_F+table$nMin1_F
  data[,25]<-table$nMaj2_F
  data[,26]<-table$nMaj2_F+table$nMin2_F
  
  return(data)
  
}

bruit<-function(data){
  
  rdata<-data
  b<-rbeta(dim(data)[1],1,100)
  rdata[,1]<-rdata[,1]-b
  rdata[,2]<-rdata[,2]-b
  
  prob<-sum(log(dbeta(b,1,100)))
  
  l<-list("data"=data,"prob"=prob)
  return(l)
  
}

reponses<-function(vect){
  
  purete<-1-vect$W[1]
  
  n<-dim(vect$Z)[1]
  c<-dim(vect$Z)[2]
  
  mut<-array(dim=n)
  
  prop<-array(dim=n)
  
  for (i in 1:n){
    
    vect$Z[i,vect$Z[i,]!=0]<-vect$Z[i,vect$Z[i,]!=0]/vect$Z[i,vect$Z[i,]!=0]
    
  }
  
  prop<-vect$Z%*%vect$W
  
  a<-as.data.frame(table(prop))
  b<-as.array(as.numeric(names(table(prop))))
  
  for (i in 1:dim(a)[1]){
    
    mut[prop==a$prop[i]]<-rownames(a)[i]
    
  }
  
  nbclones<-dim(a)[1]
  
  arbre<-matrix(ncol=2,nrow=nbclones)
  
  for (i in 1:nbclones){
    
    arbre[i,1]<-i
    
  }
  
  arbre[nbclones,2]<-0
  
  while(length(arbre[is.na(arbre[,2]),])!=0){
    
    for (i in (nbclones-1):1){
      
      if (is.na(arbre[i,2])){
        
        for(j in nbclones:(i+1)){
          
          if (is.na(sum(b[arbre[,2]==j & is.na(arbre[,2])==FALSE]))){
            
            if (b[i]<=b[j]){
              
              arbre[i,2]<-j
              break
              
            }
            
          }
          
          else{
            
            if (b[i]+sum(b[arbre[,2]==j & is.na(arbre[,2])==FALSE])<=b[j]){
              
              arbre[i,2]<-j
              break
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  a2 = a
  a2[1] = a[2]
  a2[2] = a[1]
  
  l<-list("1A"=purete,"1B"=nbclones,"1C"=a2,"2A"=mut,"3A"=arbre)
  return(l)
  
}

test<-function(data){
  
  n<-dim(data)[1]
  a<-matrix(nrow=n,ncol=2)
  a[,1]<-0
  a[,2]<-2
  
  Z<-a[,1]
  Zal<-a[,2]
  
  b<-data[,c(1,2)]%/%a
  b[1,1]<-b[order(b[,1],decreasing = FALSE)[1],1]
  b[1,2]<-b[order(b[,2],decreasing = FALSE)[1],2]
  w<-b[1,order(b[1,],decreasing = FALSE)[1]]
  data[,c(1,2)]<-data[,c(1,2)]-w*a
  
  c<-matrix(nrow=n,ncol=2)
  d<-a
  
  for (i in 1:n){
    
    a[i,]<-allele(data,i)[1,]
    c[i,2]<-dim(allele(data,i))[1]
    
  }
  
  c[,1]<-1
  j<-0
  f<-0
  
  while(f!=1){
    print("a")
    j<-j+1
    
    for (i in 1:n){
      
      a[i,]<-allele(data,i)[1,]
      
    }
    
    c[,1]<-1
    
    b<-data[,c(1,2)]%/%a
    
    b[b[,1]<0 | b[,2]<0,]<-0
    b[is.na(b[,1]),1]<-Inf
    b[is.na(b[,2]),2]<-Inf
    
    while (length(b[b[,1]==0 | b[,2]==0,])!=0){
      print("b")
      c[b[,1]==0 | b[,2]==0,1]<-c[b[,1]==0 | b[,2]==0,1]+1
      a[c[,1]>c[,2],]<-c(0,0)
      
      for (k in 1:n){
        
        if (c[k,1]<=c[k,2]){
          
          a[k,]<-allele(data,k)[c[k,1],]
          
        }
        
      }
      
      a[data[,2]==0,]<-c(0,0)
      
      b<-data[,c(1,2)]%/%a
      b[is.na(b[,1]) | is.infinite(b[,1]),1]<-Inf
      b[is.na(b[,2]) | is.infinite(b[,2]),2]<-Inf
      b[b[,1]<0 | b[,2]<0,]<-0
      #print(data[b[,1]==0 | b[,2]==0,])
      
    }
    
    #print(a[is.na(a[,1]) | is.na(a[,2]),])
    
    b[1,1]<-b[order(b[,1],decreasing = FALSE)[1],1]
    b[1,2]<-b[order(b[,2],decreasing = FALSE)[1],2]
    w<-rbind(w,b[1,order(b[1,],decreasing = FALSE)[1]])
    data[,c(1,2)]<-data[,c(1,2)]-w[j]*a
    
    Z<-as.matrix(cbind(Z,a[,1]))
    Zal<-as.matrix(cbind(Zal,a[,2]))
    
    if (sum(abs(a-d))==0){f<-1}
    
    d<-a
    
  }
  
  w<-as.vector(w)
  
  l<-list("Z"=Z,"Zal"=Zal,"W"=w)
  return(l)
  
}

##################################################################

args <- commandArgs(trailingOnly = TRUE)
table<-suppressWarnings(fread(args[1]))


data<-dataprod(table)

h1<-Monte(data,2,20000,10000)

l<-reponses(h1)

plot(density(data[,1]/data[,2]))
lines(density(h1$Z%*%h1$W/(h1$Zal%*%h1$W)),col='red')

# sortie<-matrix(nrow=(dim(h1$Z)[1]+1),ncol=dim(h1$Z)[2])
# 
# for (i in 1:(dim(data)[1])){
#   sortie[i,]<-h1$Z[i,]
# }
# 
# sortie[(dim(data)[1]+1),]<-h1$W



write.table(x=l$`1A`, quote = FALSE, file=args[2], col.names = FALSE, row.names = FALSE)
write.table(x=l$`1B`, quote = FALSE, file=args[3], col.names = FALSE, row.names = FALSE)
write.table(x=l$`1C`, quote = FALSE, file=args[4], row.names = TRUE, col.names = FALSE,sep="\t")
write.table(x=l$`2A`, quote = FALSE, file=args[5], col.names = FALSE, row.names = FALSE,sep="\t")
write.table(x=l$`3A`, quote = FALSE, file=args[6], col.names = FALSE, row.names = FALSE,sep="\t")
