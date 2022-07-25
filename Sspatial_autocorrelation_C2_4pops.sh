### Purpose: spatial autocorrelation of C2 values
# Author: Carole Smadja
# Date : 29-07-2019


#### window based analysis ##############

# sort chrID (in the form of 1, 2, ...X) so that order is : 1, 10, 11,.,19,2,3,..9,X
sort -k1,1 -k2,2n res.ana4pops_allchr_C2.bed > res.ana4pops_allchr_C2_sorted.bed


########## example given for analysis of 50kb windows

bedtools map -c 4 -o mean -a Mus_50KB_0KB_overlap.bed -b res.ana4pops_allchr_C2_sorted.bed > res.ana4pops_allchr_C2_50kbw.bed

### split data per chromosome

for i in {1..19} X; do
	awk -v chr=$i '{if($1==chr) {print}}' res.ana4pops_allchr_C2_50kbw.bed > res.ana4pops_chr${i}_C2_50kbw.bed;
done


### calculate autocorrelation coefficient for each window

#In R:
## example for chromosome 1 
data1 <- read.table("res.ana4pops_chr1_C2_50kbw.bed", header=FALSE, sep="\t")
x <- as.numeric(data1$V4)

print(data.frame(result_acf$lag,result_acf$acf)[1:30,])

# autocorrelation is very strong at lag 1
result_acf_k1=result_acf$acf[which(result_acf$lag==1)]
print(result_acf_k1)

# autocorrelation is still strong at lag 3
result_acf_k3=result_acf$acf[which(result_acf$lag==3)]
print(result_acf_k3)


# check if series is autocorrelated with itself for a lag of 3

plot  ( 1:length(x),   x,type="l")
points((1:length(x))-3,x,type="l",col="red")

# test of significance of the autocorrelation
n=3910
mu=mu_estim
sigma=sigma_estim
gamma=rep(NA,1000)
for (i in 1:1000){
  x1=rnorm(n,mu,sigma)
  x2=rnorm(n,mu,sigma)
  gamma[i]=mean((x1-mu)*(x2-mu))/(sigma^2)
}
# quantiles d'ordre 5% et 95% de gamma (x et y sont indépendants):
quantile(gamma,c(0.05,0.95))

         5%         95% 
-0.02665532  0.02634829 


print(1.96/sqrt(n))


# partial autocorrelation

result_pacf=pacf(x)

print(data.frame(result_pacf$lag,result_pacf$acf)[1:10,])








