#analysis of soft tissue versus disequilibrium pump O2 and D14C stoichiometries
#uses CM2Mc model output and DIC decomposition from egglestone and galbraith 2018
#model output and code available here:
#https://earthsystemdynamics.org/models/cm2mc-simulation-library/ 
#input file is vol weighted annual means of deep ocean (>1.5 km) 
#regression of deep ocean means to get pump sensitivities
#inversion with LGM proxy data 
#with mcmc error propagation 

#read in deep ocean tracer means from CM2Mc library 
dat<- read.csv('cm2mc_o2_c14_dic_soft_dis_deep_avg_export_Fe.csv')

#define experiments
exp_co2<- dat$exp_co2
exp_ice<- dat$exp_ice
exp_Fe<- dat$exp_Fe

#symbols for plotting based on experients
cex<- exp_co2/200
#cex<- log(exp_co2)-4.5
pch<- ifelse(exp_ice == TRUE, 21, 1)
lwd<- ifelse(exp_Fe == TRUE, 2.5, 1) 
bg_soft<-  ifelse(exp_ice == TRUE, 'cadetblue', NA)
bg_dis<-  ifelse(exp_ice == TRUE, 'coral', NA)

#define variables used for regression 
o2_soft_mean<- dat$o2_soft #o2 soft umol/kg
o2_dis_mean<- dat$o2_dis #o2 dis umol/kg
dic_soft_mean<- dat$dic_soft #dic soft umol/kg
dic_dis_mean<- dat$dic_dis #dic dis umol/kg
c14_pre_mean<- dat$c14_pre #c14_pre (=c14_dis as atm=0) per mil
c14_age_mean<- dat$c14_age #c14_age per mil 
export_total<- dat$export_total #global export produciton at 100m Pg C


#O2 fits
############
#o2 soft versus dic soft 
fit_o2_dic_soft<- lm(o2_soft_mean~dic_soft_mean)
summary(fit_o2_dic_soft)
#Call:
#lm(formula = o2_soft_mean ~ dic_soft_mean)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-1.34754 -0.38374 -0.03639  0.40509  1.83222 
#
#Coefficients:
#               Estimate Std. Error  t value Pr(>|t|)    
#(Intercept)   -1.162902   0.434562   -2.676    0.012 *  
#dic_soft_mean -1.392224   0.005646 -246.572   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.7098 on 30 degrees of freedom
#Multiple R-squared:  0.9995,	Adjusted R-squared:  0.9995 
#F-statistic: 6.08e+04 on 1 and 30 DF,  p-value: < 2.2e-16

#o2 dis versus dic dis 
fit_o2_dic_dis<- lm(o2_dis_mean~dic_dis_mean)
summary(fit_o2_dic_dis)
#Call:
#lm(formula = o2_dis_mean ~ dic_dis_mean)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-12.8655  -4.8622  -0.8994   4.2983  15.1741 
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  -36.96962    3.54806  -10.42 1.74e-11 ***
#dic_dis_mean  -0.78172    0.07128  -10.97 5.11e-12 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 7.447 on 30 degrees of freedom
#Multiple R-squared:  0.8004,	Adjusted R-squared:  0.7937 
#F-statistic: 120.3 on 1 and 30 DF,  p-value: 5.109e-12


####c14 fits 
#c14 dis versus dic dis
fit_c14_dic_dis<- lm(c14_pre_mean~dic_dis_mean)
summary(fit_c14_dic_dis)
#Call:
#lm(formula = c14_pre_mean ~ dic_dis_mean)
#
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-16.7452  -2.5609  -0.2158   2.0394  11.9749 
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  -67.78008    2.84662  -23.81   <2e-16 ***
#dic_dis_mean  -0.99962    0.05719  -17.48   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1#

#Residual standard error: 5.975 on 30 degrees of freedom
#Multiple R-squared:  0.9106,	Adjusted R-squared:  0.9076 
#F-statistic: 305.5 on 1 and 30 DF,  p-value: < 2.2e-16

#c14 soft versus dic soft 
fit_c14_dic_soft<- lm(c14_age_mean~dic_soft_mean)
summary(fit_c14_dic_soft)
#Call:
#lm(formula = c14_age_mean ~ dic_soft_mean)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-8.2838 -2.3208 -0.1342  2.3684  7.6038 
#
#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   -51.12218    2.28151  -22.41   <2e-16 ***
#dic_soft_mean  -0.48173    0.02964  -16.25   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 3.726 on 30 degrees of freedom
#Multiple R-squared:  0.898,	Adjusted R-squared:  0.8946 
#F-statistic: 264.1 on 1 and 30 DF,  p-value: < 2.2e-16

#residual standard error in dic 
summary(fit_c14_dic_soft)$sigma*1/coef(summary(fit_c14_dic_soft))[2]
#-7.735604

#hist(export_total)
#account for influence of export production changes 
export_total_mean<- mean(export_total)
export_total_anomaly<- export_total - mean(export_total)
#plot(export_total_anomaly, residuals(fit_c14_dic_soft))

#take slope from multiple regression to account for export effect  
fit_c14_dic_soft_export<- lm(c14_age_mean~dic_soft_mean+export_total_anomaly)
summary(fit_c14_dic_soft_export)
#Call:
#lm(formula = c14_age_mean ~ dic_soft_mean + export_total_anomaly)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-3.8911 -1.4289 -0.1669  1.2620  4.5076 
#
#Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          -48.73903    1.25267  -38.91  < 2e-16 ***
#dic_soft_mean         -0.51407    0.01631  -31.52  < 2e-16 ***
#export_total_anomaly  13.36459    1.53784    8.69 1.44e-09 ***
#---
#ignif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1#
#
#Residual standard error: 1.996 on 29 degrees of freedom
#Multiple R-squared:  0.9717,	Adjusted R-squared:  0.9697 
#F-statistic: 497.8 on 2 and 29 DF,  p-value: < 2.2e-16

#residual standard error in dic 
summary(fit_c14_dic_soft_export)$sigma*1/coef(summary(fit_c14_dic_soft_export))[2]
#-3.883514

#exported corrected DIC_soft fit for plotting (this has the same slope as multiple regression above)
c14_age_mean_export_cor<- c14_age_mean - (export_total_anomaly*coef(fit_c14_dic_soft_export)[3])
fit_c14_cor_dic_soft<- lm(c14_age_mean_export_cor~dic_soft_mean)
summary(fit_c14_cor_dic_soft)
#lm(formula = c14_age_mean_export_cor ~ dic_soft_mean)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-3.8911 -1.4289 -0.1669  1.2620  4.5076 
#
#Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   -48.73903    1.20174  -40.56   <2e-16 ***
#dic_soft_mean  -0.51407    0.01561  -32.92   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 1.963 on 30 degrees of freedom
#Multiple R-squared:  0.9731,	Adjusted R-squared:  0.9722 
#F-statistic:  1084 on 1 and 30 DF,  p-value: < 2.2e-16



############################################
####plot dis and soft O2 and 14C regressions 
pdf('14_o2_dic_dis_soft.pdf',width=7, height=3.5, useDingbats=FALSE, encoding="WinAnsi") 
#dev.new(width=7, height=3.5)
layout(matrix(c(1,2),ncol=2), TRUE) 
par(mar=c(3.5,3.5,1,1)); par(ps = 10, cex = 1, cex.main = 1); par(mgp=c(2.25,0.75,0)); par(las=1);par(tck=-0.01)

plot(-999,-999, xlim=c(0,150), ylim=c(-25,-175), xlab=expression(paste(DIC," (",mu*mol/kg,")")), ylab=expression(paste(Delta^14*C," (\u2030)")))
grid (NULL,NULL, lty = 6, col = adjustcolor("cornsilk2", alpha=0.7)) 
points(dic_dis_mean, c14_pre_mean,  cex=cex, pch=pch, lwd=lwd, col=adjustcolor('coral3', alpha=0.9), bg=adjustcolor(bg_dis, alpha=0.5))
abline(fit_c14_dic_dis, lwd=0.9, col=adjustcolor('coral3', alpha=0.9))
points(dic_soft_mean, c14_age_mean_export_cor,  cex=cex, pch=pch, lwd=lwd, col=adjustcolor('cadetblue4', alpha=0.15), bg=adjustcolor(bg_soft, alpha=0.1))
#abline(fit_c14_dic_soft, lwd=0.9, col=adjustcolor('cadetblue4', alpha=0.9))
#abline(-54.75166, -0.48895, lwd=0.9, col=adjustcolor('cadetblue4', alpha=0.9), lty=3)
#abline(-54.75166, -0.4040909, lwd=0.9, col=adjustcolor('cadetblue4', alpha=0.9), lty=3)
points(dic_soft_mean, c14_age_mean,  cex=cex, pch=pch, lwd=lwd, col=adjustcolor('cadetblue4', alpha=0.9), bg=adjustcolor(bg_soft, alpha=0.5))
abline(fit_c14_cor_dic_soft, lwd=0.9, col=adjustcolor('cadetblue4', alpha=0.9))
abline(coef(fit_c14_cor_dic_soft)[1], coef(fit_c14_cor_dic_soft)[2]+0.1*coef(fit_c14_cor_dic_soft)[2], lwd=0.9, col=adjustcolor('cadetblue4', alpha=0.9), lty=3)
abline(coef(fit_c14_cor_dic_soft)[1], coef(fit_c14_cor_dic_soft)[2]-0.1*coef(fit_c14_cor_dic_soft)[2], lwd=0.9, col=adjustcolor('cadetblue4', alpha=0.9), lty=3)
#lines(p, p1_min, lwd=0.9, col=adjustcolor('cadetblue4', alpha=0.9), lty=3)
#lines(p, p1_max, lwd=0.9, col=adjustcolor('cadetblue4', alpha=0.9), lty=3)

text(25,c(-160,-170), c('soft tissue','disequilibrium'), col=c('cadetblue4', 'coral3'))

text(115,-150, expression(paste('~-1 \u2030/',mu*mol/kg)), col='coral3', cex=0.9)
text(120,-80, expression(paste('~-0.5 \u2030/',mu*mol/kg)), col='cadetblue4', cex=0.9)
text(147,-110, '+10%', col='cadetblue4', cex=0.7)
text(140,-137, '-10%', col='cadetblue4', cex=0.7)

plot(-999,-999, xlim=c(0,150), ylim=c(0,-200), xlab=expression(paste(DIC," (",mu*mol/kg,")")), ylab=expression(paste(O[2]," (",mu*mol/kg,")")))
grid (NULL,NULL, lty = 6, col = adjustcolor("cornsilk2", alpha=0.7)) 
points(dic_dis_mean, o2_dis_mean,  cex=cex, pch=pch, lwd=lwd, col=adjustcolor('coral3', alpha=0.9), bg=adjustcolor(bg_dis, alpha=0.5))
abline(fit_o2_dic_dis, lwd=0.9, col=adjustcolor('coral3', alpha=0.9))
abline(0, 1/-0.7072, lwd=0.9, lty=4, col=adjustcolor('cadetblue4', alpha=0.6)) 
points(dic_soft_mean, o2_soft_mean,  cex=cex, pch=pch, lwd=lwd, col=adjustcolor('cadetblue4', alpha=0.9), bg=adjustcolor(bg_soft, alpha=0.5))
abline(fit_o2_dic_soft, lwd=0.9, col=adjustcolor('cadetblue4', alpha=0.9))

#text(25,c(-160,-170), c('soft tissue','disequilibrium'), col=c('cadetblue4', 'coral3'))

text(130,-160, '~-1.4', col='cadetblue4', cex=0.9)
text(90,-80, '~-0.8', col='coral3', cex=0.9)
legend('topleft', pch=c(1,1,1,1,1,1,21,1), pt.cex=c(unique(cex), 1,1,1), c(as.character(unique(exp_co2)), 'no ice sheet' , 'ice sheet','Fe'), col='grey37', pt.bg=c(NA, NA,NA,NA,NA,NA,'grey87',NA), pt.lwd=c(1,1,1,1,1,1,1,2.5), cex=0.8, bty='n')
legend('bottomright', lty=4, col=adjustcolor('cadetblue4', alpha=0.9),'anderson94',text.col=adjustcolor('cadetblue4', alpha=0.9), cex=0.8, bty='n')

dev.off()
#end plot 



################################################
####LGM data import and inversion

#define D14C O2 vectors for plot
delta_dic<- seq(-25, 250) #pred range
delta_c14_dis<- delta_dic*coef(fit_c14_dic_dis)[2]
delta_c14_soft<- delta_dic*coef(fit_c14_cor_dic_soft)[2]
delta_c14_soft_export_plus<- delta_dic*(coef(fit_c14_cor_dic_soft)[2]-0.1*coef(fit_c14_cor_dic_soft)[2]) #10% increase in export - slope change calculated offline using Egglestone and Galbraith Eq.1 with age converted into C14 per mil 
delta_c14_soft_export_minus<- delta_dic*(coef(fit_c14_cor_dic_soft)[2]+0.1*coef(fit_c14_cor_dic_soft)[2]) #10% decrease in export - slope change calculated offline using Egglestone and Galbraith Eq.1 with age converted into C14 per mil 
delta_o2_dis<- delta_dic*coef(fit_o2_dic_dis)[2]
delta_o2_dis_low<- delta_dic*-0.5801711 #5th percentile of within simulation o2_dis/dic_dis slope - calculated offline
delta_o2_dis_high<- delta_dic*-0.8771248 #5th percentile of within simulation o2_dis/dic_dis slope - calculated offline
delta_o2_soft<- delta_dic*coef(fit_o2_dic_soft)[2]


#load error bar functions for plot 
load('xerrbar.Rdata')
load('yerrbar.Rdata')

#read in deep ocean proxy values and uncertainties
#c14 data from Rafter 2022
#o2 data from hoogakker 2015;2018; umling+thunell 2018; anderson 2019; gotschalk 2016
lgm_proxy_data<- read.csv('proxy_o2_d14c_so2_14c_rafter_err.csv') #13ka diff SO

#o2 solubility change for -2.3 oC + 1 PSU LGM 
DO2_sol<- 16.79792
lgm_proxy_data$o2_anom_sol<- lgm_proxy_data$o2_anom - DO2_sol #correct o2 for solubility change


#################
#plot LGM proxy data
pdf('dd14c_do2_lgm.pdf',width=3.5, height=3.5, useDingbats=FALSE, encoding="WinAnsi") 
#dev.new(width=3.5, height=3.5)
par(mar=c(3.5,3.5,1,1)); par(ps = 10, cex = 1, cex.main = 1); par(mgp=c(2.25,0.75,0)); par(las=1);par(tck=-0.01)
plot(-999,-999, xlim=c(-10,-200), ylim=c(10,-200), xlab=expression(paste(Delta*Delta^14*C[0-adj]," (\u2030)")), ylab=expression(paste(Delta*O[paste(2,",",sol-corr)]," (",mu*mol/kg,")")))
grid (NULL,NULL, lty = 6, col = adjustcolor("cornsilk2", alpha=0.7)) 
lines(delta_c14_dis, delta_o2_dis, lwd=1.2, col=adjustcolor('coral3', alpha=0.8))
lines(delta_c14_dis, delta_o2_dis_low, lwd=1, lty=3, col=adjustcolor('coral1', alpha=0.8))
lines(delta_c14_dis, delta_o2_dis_high, lwd=1, lty=3, col=adjustcolor('coral1', alpha=0.8))
lines(delta_c14_soft, delta_o2_soft, lwd=1.2, col=adjustcolor('cadetblue4', alpha=0.8))
lines(delta_c14_soft_export_plus, delta_o2_soft, lwd=1.2, lty=3, col=adjustcolor('cadetblue4', alpha=0.8))
lines(delta_c14_soft_export_minus, delta_o2_soft, lwd=1.2, lty=3, col=adjustcolor('cadetblue4', alpha=0.8))

text(-30,-145, 'soft', col='cadetblue4', cex=0.9, font=2)
text(-160,-90, 'disequilibrium', col='coral3', cex=0.9, font=2)
text(-40,-162, '+10%', col='cadetblue4', cex=0.7)
text(-80,-162, '-10%', col='cadetblue4', cex=0.7)
text(-195,-187, '-0.9', col=adjustcolor('coral1', alpha=0.8), cex=0.7)
text(-195,-105, '-0.6', col=adjustcolor('coral1', alpha=0.8), cex=0.7)

xerrbar(x=lgm_proxy_data$D14C_anom, y=lgm_proxy_data$o2_anom_sol, y_err=lgm_proxy_data$o2_err, errbar_lwd=0.8, errbar_lty=1, errbar_col=adjustcolor('grey57', alpha=0.8))
xerrbar(x=lgm_proxy_data$D14C_anom, y=lgm_proxy_data$o2_anom_sol, y_err=1.95*lgm_proxy_data$o2_err, errbar_lwd=0.6, errbar_lty=1, errbar_col=adjustcolor('grey57', alpha=0.6))
yerrbar(x=lgm_proxy_data$D14C_anom, y=lgm_proxy_data$o2_anom_sol, x_err=lgm_proxy_data$D14C_err, errbar_lwd=0.8, errbar_lty=1, errbar_col=adjustcolor('grey57', alpha=0.8))
yerrbar(x=lgm_proxy_data$D14C_anom, y=lgm_proxy_data$o2_anom_sol, x_err=1.95*lgm_proxy_data$D14C_err, errbar_lwd=0.6, errbar_lty=1, errbar_col=adjustcolor('grey57', alpha=0.6))

points(lgm_proxy_data$D14C_anom, lgm_proxy_data$o2_anom_sol, pch=c(25,24,22,21), col=c('grey17','grey17','grey17', 'grey37'), bg=adjustcolor(c('grey57','grey57','grey57', 'grey87'), alpha=0.6), cex=1.6)

points(lgm_proxy_data$D14C_anom[4], -168-DO2_sol, pch=1, col=adjustcolor('grey37', alpha=0.8), cex=1.6) #SO O2 anomaly from clim

legend('bottomright', pch=c(24,25,22,21), col=c('grey17','grey17','grey17','grey37'),pt.bg=adjustcolor(c('grey57','grey57','grey57','grey87'), alpha=0.8), c('Pacific mid','Pacific deep','Atlantic deep','SO deep'), cex=0.9, bty='n')


dev.off()
#end plot 


#################
####LGM inversion
#a<- coef(fit_c14_dic_soft)[2]
a<- coef(fit_c14_cor_dic_soft)[2]
b<- coef(fit_c14_dic_dis)[2]
c<- coef(fit_o2_dic_soft)[2]
d<- coef(fit_o2_dic_dis)[2]

DDIC_soft_LGM<- (d*lgm_proxy_data$D14C_anom - b*lgm_proxy_data$o2_anom_sol)/(a*d-b*c)
DDIC_dis_LGM<- (-c*lgm_proxy_data$D14C_anom + a*lgm_proxy_data$o2_anom_sol)/(a*d-b*c)

DDIC_total_LGM<- DDIC_soft_LGM+DDIC_dis_LGM

frac_soft_lgm<- DDIC_soft_LGM/DDIC_total_LGM

lgm_out<- data.frame(lgm_proxy_data, DDIC_soft_LGM, DDIC_dis_LGM,DDIC_total_LGM, frac_soft_lgm)
print(lgm_out)

write.csv(lgm_out, 'lgm_dis_soft_out.csv', row.names=FALSE)
#######

#deep ocean means
lgm_ddic_soft_mean<- mean(DDIC_soft_LGM[c(1,3,4)])
print(lgm_ddic_soft_mean)
#63.94212
lgm_ddic_dis_mean<- mean(DDIC_dis_LGM[c(1,3,4)])
print(lgm_ddic_dis_mean)
#73.05612

#total deep ocean DICsoft+dis
lgm_ddic_tot<- mean(DDIC_total_LGM[c(1,3,4)])
print(lgm_ddic_tot)
#136.9982

####substitutions 
########
#approx age component of lgm c14 in years
lgm_ddic_soft_mean*a/(-10/83)
# 272.8261  yrs
#approx dis component of lgm c14 in years
lgm_ddic_dis_mean*b/(-10/83)
#606.1357 years 

#o2
#soft component of lgm o2
lgm_ddic_soft_mean*c
#-89.02175 
#dis component of lgm o2
lgm_ddic_dis_mean*d
#-57.10951 


#####################################
####inversion export sensitivity test
#a<- coef(fit_c14_dic_soft)[2]
a<- coef(fit_c14_cor_dic_soft)[2]-0.1*coef(fit_c14_cor_dic_soft)[2]
b<- coef(fit_c14_dic_dis)[2]
c<- coef(fit_o2_dic_soft)[2]
d<- coef(fit_o2_dic_dis)[2]

DDIC_soft_LGM<- (d*lgm_proxy_data$D14C_anom - b*lgm_proxy_data$o2_anom_sol)/(a*d-b*c)
DDIC_dis_LGM<- (-c*lgm_proxy_data$D14C_anom + a*lgm_proxy_data$o2_anom_sol)/(a*d-b*c)

DDIC_total_LGM<- DDIC_soft_LGM+DDIC_dis_LGM

frac_soft_lgm<- DDIC_soft_LGM/DDIC_total_LGM

lgm_out_exportplus<- data.frame(lgm_proxy_data, DDIC_soft_LGM, DDIC_dis_LGM,DDIC_total_LGM, frac_soft_lgm)
print(lgm_out_exportplus)

write.csv(lgm_out_exportplus, 'lgm_dis_soft_out_Fe.csv', row.names=FALSE)
#######################



##################################
##################################
####mcmc uncertainty propogation 
a<- coef(fit_c14_cor_dic_soft)[2]
b<- coef(fit_c14_dic_dis)[2]
c<- coef(fit_o2_dic_soft)[2]
d<- coef(fit_o2_dic_dis)[2]

a_sigma<- coef(summary(fit_c14_cor_dic_soft))[2,2]
b_sigma<- coef(summary(fit_c14_dic_dis))[2,2]
c_sigma<- coef(summary(fit_o2_dic_soft))[2,2]
d_sigma<- coef(summary(fit_o2_dic_dis))[2,2]

mcmc_iterations<-9999

r_DDIC_soft_LGM<- matrix(, nrow = nrow(lgm_proxy_data), ncol = mcmc_iterations)
r_DDIC_dis_LGM<- matrix(, nrow = nrow(lgm_proxy_data), ncol = mcmc_iterations)
r_DDIC_soft_LGM_d_mean<- matrix(, nrow = 1, ncol = mcmc_iterations)
r_DDIC_dis_LGM_d_mean<- matrix(, nrow = 1, ncol = mcmc_iterations)

for(i in 1:mcmc_iterations){
	D14C_anom_i<- rnorm(nrow(lgm_proxy_data), lgm_proxy_data$D14C_anom, lgm_proxy_data$D14C_err)
	O2_anom_i<- rnorm(nrow(lgm_proxy_data), lgm_proxy_data$o2_anom_sol, lgm_proxy_data$o2_err)
	a_i<- rnorm(1, a, a_sigma)
	b_i<- rnorm(1, b, b_sigma)
	c_i<- rnorm(1, c, c_sigma)
	d_i<- rnorm(1, d, d_sigma)

	DDIC_soft_LGM_i<- (d_i*D14C_anom_i - b_i*O2_anom_i)/(a_i*d_i-b_i*c_i)
	DDIC_dis_LGM_i<- (-c_i*D14C_anom_i + a_i*O2_anom_i)/(a_i*d_i-b_i*c_i)

	r_DDIC_soft_LGM[,i]<- DDIC_soft_LGM_i
	r_DDIC_dis_LGM[,i]<- DDIC_dis_LGM_i
	
	r_DDIC_soft_LGM_d_mean[i]<- mean(DDIC_soft_LGM_i[c(1,3,4)])
	r_DDIC_dis_LGM_d_mean[i]<- mean(DDIC_dis_LGM_i[c(1,3,4)])
}

r_DDIC_total_LGM<- r_DDIC_soft_LGM+r_DDIC_dis_LGM
r_frac_soft_lgm<- r_DDIC_soft_LGM/r_DDIC_total_LGM

#means
DDIC_soft_LGM<- apply(r_DDIC_soft_LGM,1,mean)
DDIC_soft_LGM_sigma<- apply(r_DDIC_soft_LGM,1,sd)

DDIC_dis_LGM<- apply(r_DDIC_dis_LGM,1,mean)
DDIC_dis_LGM_sigma<- apply(r_DDIC_dis_LGM,1,sd)

DDIC_total_LGM<- apply(r_DDIC_total_LGM,1,mean)
DDIC_total_LGM_sigma<- apply(r_DDIC_total_LGM,1,sd)

frac_soft_lgm<- apply(r_frac_soft_lgm,1,mean)
frac_soft_lgm_sigma<- apply(r_frac_soft_lgm,1,sd)

lgm_out_err<- data.frame(lgm_proxy_data, DDIC_soft_LGM, DDIC_soft_LGM_sigma, DDIC_dis_LGM,DDIC_dis_LGM_sigma,DDIC_total_LGM,DDIC_total_LGM_sigma, frac_soft_lgm,frac_soft_lgm_sigma)
print(lgm_out_err)

write.csv(lgm_out_err, 'lgm_dis_soft_out_err.csv', row.names=FALSE)
#the end!