####for linear, ok, runing with old delta1 sequence, try also new one (more scripts!)

delta0=c(0.2)
#delta1=seq(from=-1,to=1,by=0.1) #for LM it is symmetric, so better use seq(from=0,to=1,by=0.1) 
delta1=seq(from=0.2,to=1,by=0.4)
task=c("medium")   
var.sigma =c(3) 
nn<-c(100,250,500) 


txt1<-"
source(\"source_sim.R\")\n
seeds<-1:20000\n



model<-\"linear\"\n
phi.beta<-NULL \n
intercept<-1 \n
mean_y<-NULL \n
B<-1000 \n



no<-1\n



for (ii in no){\n
if (ii==1) alter=FALSE else alter<-TRUE \n

  
    for (kk in alter){\n
m=mm=0\n
while(m<B){\n

  m=m+1\n
mm=mm+1\n
if (mm>20000)  set.seed(max(seeds)+mm-20000)  else      set.seed(seeds[mm])\n

  mydata<-mysim_data(n="

txt1.1<-",dif=\"medium\",no=ii, delta0=0.2,delta1="
txt1.2<-",alter=kk,model=model,intercept=intercept,var.sigma=3,phi.beta = phi.beta)


 
  dd<-try(my_sim(i=m,delta0="


  txt2<-"     )   ,silent=TRUE)\n
  if (class(dd)==\"try-error\") m<-m-1\n
}\n



}}\n
"
mm=0
for (ii in delta0){
  for (jj in delta1){
    for (zz in task){
      for (kk in var.sigma){
for (nnn in nn){
        mm=mm+1
    cat(txt1,nnn,txt1.1,jj,txt1.2,ii,
        ",delta1=",jj,
        ",task=\"",zz,"\"
        ,alter=kk,
        no=ii,
        n=",nnn,
        ",model=model,
        var.sigma=",kk,"\n,phi.beta=phi.beta,
        intercept=intercept,mean_y=mean_y",txt2,sep="",file=paste("script",mm,".R",sep=""))    
        
      }}}}}






txt1<-"
source(\"source_sim.R\")\n
seeds<-1:20000\n



model<-\"linear\"\n
phi.beta<-NULL \n
intercept<-1 \n
mean_y<-NULL \n
B<-1000 \n



no<-0\n



for (ii in no){\n
 alter=FALSE  \n

  
    for (kk in alter){\n
m=mm=0\n
while(m<B){\n

  m=m+1\n
mm=mm+1\n
if (mm>20000)  set.seed(max(seeds)+mm-20000)  else      set.seed(seeds[mm])\n

  mydata<-mysim_data(n="




for (ii in delta0){
  for (jj in 0){
    for (zz in task){
      for (kk in var.sigma){
for (nnn in nn){
        mm=mm+1
    cat(txt1,nnn,txt1.1,jj,txt1.2,ii,
        ",delta1=",jj,
        ",task=\"",zz,"\"
        ,alter=kk,
        no=0,
        n=",nnn,
        ",model=model,
        var.sigma=",kk,"\n,phi.beta=phi.beta,
        intercept=intercept,mean_y=mean_y",txt2,sep="",file=paste("script",mm,".R",sep=""))    
        
      }}}}}


