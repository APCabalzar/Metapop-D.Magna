###############################################################################################
##Script to iteratively solve the Euler-Lotka equation for the intrinsic growth rate, small r##
###############################################################################################

mydata<-read.table("PATH", header=T) #provide a data frame with the following columns: 1:replicate; 2:Names ; columns with ages at reproduction (x1-xn); columns with corresponding clutch sizes (g1-gn)
mydata<-na.omit(mydata)#remove NA's in the data set

nom_vector<-(mydata[,2]) #Vectorize the names (here column 2), to add the names back to lists later on

#Vectorize ages at reproduction (Vx) and clutch sizes (Vg) for each sample

Vx=NULL ##define a list in which the vectors will be stored
for(n in 1:n) #Where n is the number of rows
{Vx[[n]]=as.numeric(mydata[n,3:5])} ###loop to vectorize x1:x3 (here columns 3 to 5)
names(Vx)<-nom_vector #Add names back

Vg=NULL
for(n in 1:n)
{Vg[[n]]=as.numeric(mydata[n,6:8])} ###loop to vectorize g1:g3
names(Vg)<-nom_vector

Res=NULL ##define a list in which the results will be stored
eulerlotka <- function(r) sum (Vg[[n]] * exp(-r * Vx[[n]])) - 1  #set the Euler-Lotka equation to zero in order to use "uniroot"
r.range <- c(-10,10) #define a range in which the root should be searched

for(n in 1:n)
{Res[[n]]=uniroot(f = eulerlotka, interval = r.range, tol = 1e-8)}  ### loop to find root for each sample

#Extract the roots from the generated list and store them in a vector
r_table=NULL
for (n in 1:n)
{r_table[[n]]=Res[[n]]$root}
names(r_table)<-nom_vector # put names back
