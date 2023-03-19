# Factor Analysis

library(psych)

raisin_data <- read.csv("/Users/parthvi/Desktop/Desktop/Parthvi/MITA/SEM2/MVA/week 7/Raisin_Dataset.csv", row.names=1)

attach(raisin_data)
raisin_data[1:6]
raisin_data
fit.pc <- principal(raisin_data[1:6], nfactors=4, rotate="varimax")

fit.pc
round(fit.pc$values, 3)
fit.pc$loadings
# Loadings with more digits
for (i in c(1,3,2,4)) { print(fit.pc$loadings[[1,i]])}
# Communalities
fit.pc$communality
# Rotated factor scores, Notice the columns ordering: RC1, RC3, RC2 and RC4
fit.pc$scores
# Play with FA utilities

fa.parallel(raisin_data[1:6]) # See factor recommendation
fa.plot(fit.pc) # See Correlations within Factors
fa.diagram(fit.pc) # Visualize the relationship
vss(raisin_data[1:6]) # See Factor recommendations for a simple structure




# Computing Correlation Matrix
corrm.emp <- cor(raisin_data[1:6])
corrm.emp
plot(corrm.emp)
raisin_data_pca <- prcomp(raisin_data[1:6], scale=TRUE)
summary(raisin_data_pca)
plot(raisin_data_pca)
# A table containing eigenvalues and %'s accounted, follows. Eigenvalues are the sdev^2
(eigen_raisin_data <- round(raisin_data_pca$sdev^2,3))
round(fit.pc$values, 3)
names(eigen_raisin_data) <- paste("PC",1:6,sep="")
eigen_raisin_data
sumlambdas <- sum(eigen_raisin_data)
sumlambdas
propvar <- round(eigen_raisin_data/sumlambdas,2)
propvar
cumvar_raisin_data <- cumsum(propvar)
cumvar_raisin_data
matlambdas <- rbind(eigen_raisin_data,propvar,cumvar_raisin_data)
matlambdas
rownames(matlambdas) <- c("Eigenvalues","Prop. variance","Cum. prop. variance")
rownames(matlambdas)
eigvec.emp <- raisin_data_pca$rotation
print(raisin_data_pca)
# Taking the first four PCs to generate linear combinations for all the variables with four factors
pcafactors.emp <- eigvec.emp[,1:4]
pcafactors.emp
# Multiplying each column of the eigenvector’s matrix by the square-root of the corresponding eigenvalue in order to get the factor loadings
unrot.fact.emp <- sweep(pcafactors.emp,MARGIN=2,raisin_data_pca$sdev[1:4],`*`)
unrot.fact.emp
# Computing communalities
communalities.emp <- rowSums(unrot.fact.emp^2)
communalities.emp
# Performing the varimax rotation. The default in the varimax function is norm=TRUE thus, Kaiser normalization is carried out
rot.fact.emp <- varimax(unrot.fact.emp)
#View(unrot.fact.emp)
rot.fact.emp
# The print method of varimax omits loadings less than abs(0.1). In order to display all the loadings, it is necessary to ask explicitly the contents of the object $loadings
fact.load.emp <- rot.fact.emp$loadings[1:6,1:4]
fact.load.emp
# Computing the rotated factor scores for the 30 European Countries. Notice that signs are reversed for factors F2 (PC2), F3 (PC3) and F4 (PC4)
scale.emp <- scale(raisin_data[1:6])
scale.emp
as.matrix(scale.emp)%*%fact.load.emp%*%solve(t(fact.load.emp)%*%fact.load.emp)


