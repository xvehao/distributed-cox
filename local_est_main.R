rm(list=ls())
library(glmnet)
library(splines)
library(MASS)
library(dplyr)
library(survival)
library(statmod)
library(ggplot2)
source('function_local.R')
set.seed(123)
### The commented code is for simulating data only
# n_clinical_feature = 15
# num_bh_knot = 2
# d_phi = 3
# T0=300
# lambda=0.05
# Xcov=1*ar1_cor(n_clinical_feature,0.2)
# rateC=0.1
# missing_rate = c(rep(.1,5),rep(.3,5),rep(.5,5))
# betas = data.frame(X1=-.8,X2=.8,X3=.1,X4=-.1,
#                    M1=.8,M2=.8,M3=-.1,M4=.1,
#                    'fw1.X1'=.6,'fw1.X2'=-.6,'fw1.M1'=.6,'fw1.M2'=-.6,fw=1)
# fw_coef=1
# fw_intercept=1
# dist_family = 'gompertz'
# dat_tab=local_est_tab=sum_stat_tab=params_tab=list()
# 
# shape_noise = rnorm(2,0,0.1) # generate purterbation for coefficients and intercept of shape function
# local_fw_coef = shape_noise[1]+fw_coef
# local_fw_intercept = shape_noise[2]+fw_intercept
# rho=runif(1,0.3,0.6)
# n_patient = sample(c(3000,1000),1)
# beta_noise = rnorm(ncol(betas),0,0.1)
# local_betas = beta_noise+betas
# ### simulate data for each site given parameters
# dat_obj=generate_data(n_patient,n_clinical_feature,local_betas,T0,missing_rate,
#                       lambda,rho,rateC,Xcov,local_fw_coef,local_fw_intercept,d_phi,dist_family=dist_family)
# sample_data = dat_obj$dat.obs %>% select(id:X15,calendar_time)
# colnames(sample_data) = c('id','y','failed',
#                           'a','b','c','d','e',
#                           'f','g','h','i','j',
#                           'k','l','m','n','o','calendar_time')
# write.csv(sample_data,'sample_data.csv',row.names=F)


### The simulated data is stored in the sample_data.csv file
sample_data <- read.csv("sample_data.csv")

### user have to specify 
### x: clinical feature; 
### time: survival time; 
### event: censorship (1 for not censored); 
### calendar_time: calendar_time 
x = sample_data %>% select(a:o)
time = sample_data$y
event = sample_data$failed
calendar_time = sample_data$calendar_time
missing_col=1:15
# see surv_local in function_local.R for details
sum_stat = surv_local(time=time, event=event, 
                x=x, calendar_time=calendar_time, 
                missing_col=missing_col,
                num_bh_knot=2, d_phi=3)
# to access the 
sum_stat$grad 
sum_stat$Hessian
### the commented code is used to plot predicted values and true values
eigen(sum_stat$Hessian)$value
# sum_stat$eta
# scale_factor = local_fw_coef+local_fw_coef
# dat = dat_obj$dat.true
# hist(dat$y,breaks = 100)
# mean(dat$failed)
# xdata = dat_obj$dat.obs
# alasso_model = sum_stat$alasso
# breslow_cbh = sum_stat$cbh
# pe_cbh = data.frame(time=xdata$y, hazard=sum_stat$cbhhat.obj$int)
# actual_cbh = get_true_cbh(breslow_cbh$time,lambda,rho,dist_family=dist_family)
# ## the arugument scale is used to adjust the hazard of second input into hazard*exp(scale),
# ## you can use the defalut value 0 if you do not want to adjust the hazard
# cbh.plot = const_plot(actual_cbh,
#            breslow_cbh[breslow_cbh$time<quantile(breslow_cbh$time,1),],
#            pe_cbh[pe_cbh$time<quantile(pe_cbh$time,1),],
#            scale=-scale_factor,
#            label1='actual',label2='breslow',label3='piecewise exp',
#            xlab='time',ylab='cumulative baseline hazard',
#            title='cumulative baseline hazard approximated by piecewise exp')
# cbh.plot
# bh_true = get_true_bh(sum_stat$psi.obj$cbh$time,
#                   lambda,rho,dist_family=dist_family)
# bh_approx = data.frame(time=sum_stat$psi.obj$cbh$time,
#                         hazard=as.array(exp(t(sum_stat$psi)%*%(sum_stat$eta))))
# colnames(bh_approx)=c('time','hazard')
# bh.plot=const_plot(bh_approx[bh_approx$time<quantile(bh_approx$time,1),],
#            bh_true,scale=scale_factor,label1='piecewise exp approx',label2='true value',
#            title='baseline hazard')
# bh.plot
# 
# fct_actual = data.frame(x=1:T0,y=fw(1:T0,local_fw_coef,local_fw_intercept),group='actual')
# fct_lasso_pred = data.frame(x=dat$calendar_time,
#                             y=as.matrix(xdata[,paste0('phi',1:d_phi)])%*%sum_stat$xi,
#                             group='pred')
# fct_lasso = rbind(fct_lasso_pred,fct_actual)
# ggplot()+
#   geom_point(data=fct_lasso,aes(x=x,y=y,col=group))+
#   xlab('time')+
#   ylab('shape function')+
#   ggtitle('shape function using lasso')
# 
