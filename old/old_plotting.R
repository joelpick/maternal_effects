
load(file="/Users/joelpick/Dropbox/0_fitness/fitness_timing/Data/Intermediate/gaussian_sims.Rdata")

list2array <- function(x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x)), dimnames=list(NULL,c("A","Me","Mg","cov_AMg","E"),NULL))


sim_sum <- function(x,models=4){
	s_array <- list2array(x)[1:models,,]
	s_means<-apply(s_array,c(1,2),mean)
	s_ses<-apply(s_array,c(1,2),function(x) sd(x)/sqrt(length(x)) ) 
	list(mean=s_means,se=s_ses)
}

plot_model <- function(hs, fs, sim_values, models=4, mod_names){
	
	h_sum <- sim_sum(hs, models=models)
	f_sum <- sim_sum(fs, models=models)
	n_traits <- 5

	plot(NA,ylim=c(0,1),xlim=c(1,n_traits * models + models-1), xaxt="n", xlab="", ylab="Variance")
	# axis(1, (1:models)*n_traits+(1:models)-mean(1:n_traits), mod_names)
	axis(1, 1:(models*n_traits + models), rep(c("Va","Vme","Vmg","COVa,mg","Ve",""),models))
	axis(3, (1:models)*n_traits+(1:models)-mean(1:n_traits), mod_names)

	# plot(NA,ylim=c(0,1),xlim=c(1,models+1))
	for(i in 1:models){

		x<- rep(i*n_traits+i,n_traits)+(1:n_traits)-mean(1:n_traits)*2

		# sim values
		arrows(x-0.5, sim_values,x+0.5, sim_values, code=0, col="grey")

		arrows(x,h_sum$mean[i,]+h_sum$se[i,],x,h_sum$mean[i,]-h_sum$se[i,],code=3, angle=90, col=1:5, length=0.05)
		points(x,h_sum$mean[i,], pch=19, col=1:5)

		arrows(x,f_sum$mean[i,]+f_sum$se[i,],x,f_sum$mean[i,]-f_sum$se[i,],code=3, angle=90, col=1:5, length=0.05)
		points(x,f_sum$mean[i,], pch=18, col=1:5)
	}
	abline(v=c((1:models)*(n_traits+1)))

	# legend("topleft", c("") col=1:5)

}


plot_trait <- function(hs, fs, sim_values, models=5){
	
	h_sum <- sim_sum(hs, models=models)
	f_sum <- sim_sum(fs, models=models)
	traits <- c("V_A","V_Me","V_Mg","cov_A,Mg","V_E")
	n_traits <- 5

	plot(NA,ylim=c(0,1),xlim=c(1,n_traits * models + n_traits), xaxt="n")
	axis(1, (1:n_traits)*models+(1:n_traits)-mean(1:models), traits)
	for(i in 1:n_traits){

		x<- rep(i*models+i,models)+(1:models)-mean(1:models)*2

		# sim values
		arrows(x[1]-0.5, sim_values[i],x[models]+0.5, sim_values[i], code=0, col="grey")

		arrows(x,h_sum$mean[,i]+h_sum$se[,i],x,h_sum$mean[,i]-h_sum$se[,i],code=3, angle=90, col=1:models, length=0.1)
		points(x,h_sum$mean[,i], pch=19, col=1:models)

		arrows(x,f_sum$mean[,i]+f_sum$se[,i],x,f_sum$mean[,i]-f_sum$se[,i],code=3, angle=90, col=1:models, length=0.1)
		points(x,f_sum$mean[,i], pch=18, col=1:models)
	}
	abline(v=c((1:models)*(n_traits+1)))

}



# plot_trait(s1_hs,s1_fs, c(Va=0.5, Vme=0, Vmg=0.000001, r_amg=0, Ve=0.5), models=4)
plot_model(s1_hs,s1_fs, c(Va=0.5, Vme=0, Vmg=0.000001, r_amg=0, Ve=0.5), models=4, mod_names=paste0("M",c(1,2,3,5)))
