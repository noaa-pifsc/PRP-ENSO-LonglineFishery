#BRT for Bigeye tuna in deep-set fishery for the Climate Variability Project

# Load libraries
library(matrixStats)
library(fmsb)
library(here)

# Set working directory
mainDir <- here()
setwd(file.path(mainDir, 'BRTs/'))

# Load Justin's BRT code
source("ENSO_BRT_Eval_Function_JJS.R")

# Load dataset
# This dataset is created in the code CombineDataForBRTs.R and is 1 degree gridded data across the deep-set fishing grounds 1995-2024
# Here we add some presence/absence, and log transform the CPUE to make it more normally distributed
dfAll <- readRDS(file.path(mainDir, 'FisheryData/BRT data/ENSO_BRTdata.rds')) %>% 
  filter(Lat >= 10, Lat <= 40.125, Lon >= 180, Lon <= 230) %>% 
  mutate(BET_PA = if_else(BIGEYE_TUNA>0, 1, 0, missing=NA), YFT_PA=if_else(YELLOWFIN >0, 1, 0, missing=NA), DOL_PA=if_else(MAHIMAHI>0, 1, 0, missing=NA), BRZ_PA=if_else(POMFRET>0, 1, 0, missing=NA), SWO_PA=if_else(SWORDFISH>0, 1, 0, missing=NA),
         Log_BET_CPUE=log(BET_CPUE+0.5), Log_YFT_CPUE=log(YFT_CPUE+0.5), Log_BRZ_CPUE=log(BRZ_CPUE+0.5), Log_DOL_CPUE=log(DOL_CPUE+0.5), Log_SWO_CPUE=log(SWO_CPUE+0.5), Random=runif(54909,0,10)) %>% 
  as.data.frame()
dfNE <- dfAll %>% filter(Longitude >= (209.875) & Latitude >= 20.125)
dfCW <- dfAll %>% filter(Longitude <= (209.875),  Latitude >= 20.125, Latitude <= 26.125) 
dfNW <- dfAll %>% filter(Longitude <= (209.875),  Latitude >= 26.125)
dfSW <- dfAll %>% filter(Longitude <= (209.875),  Latitude <= 20.125) 
dfSE <- dfAll %>% filter(Longitude >= (209.875),  Latitude <= 20.125)

df <- dfAll[idx <- which(!is.na(df$Log_BET_CPUE)),]

# for (i in 1:5) {
#   spp <- c('YFT', 'BRZ', 'DOL', 'SWO','BET')[i]
#   #for (j in 1:5) {
#     region <- c('NE', 'CW', 'NW', 'SW', 'SE')[j]
#     df <- list(dfNE, dfCW, dfNW, dfSW, dfSE)[j][[1]]
# 
# 
#     Predictors <- which(colnames(df) %in% c("Oxy_1ml", "Oxy_2ml", "Catchability250", "Catchability400", "ONI", "Phase", "PDO", "NPGO", "Effort", "Random"))
#     #Predictors <- which(colnames(df) %in% c("Lon", "Lat", "Oxy_1ml", "Catchability250", "Catchability400", "Year", "Month", "ONI", "Phase", "PDO", "NPGO", "Random"))
#     Response <- which(colnames(df) %in% paste("Log",spp, 'CPUE', sep="_"))
#     #Response <- which(colnames(df) %in% paste(spp, 'PA', sep="_"))
#     
#     Abund_Model_Step<-fit.brt.n_eval_Balanced(df, gbm.x=Predictors, gbm.y= c(Response), lr=0.001, tc=3, family = "gaussian", bag.fraction=0.75, n.folds=5, 5)
#     #PA_Model_Step<-fit.brt.n_eval_Balanced(df, gbm.x=Predictors, gbm.y= c(Response), lr=0.001, tc=3, family = "bernoulli",bag.fraction=0.75, n.folds=5, 5)
#     saveRDS(Abund_Model_Step, paste0('ENSO_', spp, '_testBRT_logCPUE', region, '_wEffort.rds'))
#     #saveRDS(Abund_Model_Step, paste0('ENSO_', spp, '_testBRT_logCPUE_fullRegion_wEffort.rds'))
#     #saveRDS(PA_Model_Step, paste0('ENSO_', spp, '_testBRT_PA_', region, '_wEffort.rds'))
#     #saveRDS(PA_Model_Step, paste0('ENSO_', spp, '_testBRT_PA_fullRegion_wEffort.rds'))
#     
#   #}
# }
  
#df <- dfNE
#region <- 'NE'
#spp <- 'YFT'

# Plot the distribution to make sure a gaussian distribution is appropriate
# ggplot(df, aes(x=Log_CPUE)) + 
#   geom_histogram(color='black', fill='lightgray', binwidth=0.1)

# # Set your predictors
# #Predictors <- which(colnames(df) %in% c("Lon", "Lat", "Oxy_1ml", "Catchability250", "Catchability400", "Year", "Month", "ONI", "Phase", "PDO", "NPGO", "Random"))
# Predictors <- which(colnames(df) %in% c("Oxy_1ml", "Catchability250", "Catchability400", "ONI", "Phase", "PDO", "NPGO", "Random"))
# 
# # And response variable to use
# Response <- which(colnames(df) %in% paste("Log",spp, 'CPUE', sep="_"))
# 
# # Fit model to all predictors
# # Here we are running test models so we are only running 5 models with a higher learning rate. For the final model I'd use 50-100 models and a learning rate of 0.001
# Abund_Model_Step<-fit.brt.n_eval_Balanced(dfCW, gbm.x=Predictors, gbm.y= c(Response), lr=0.001, tc=3, family = "gaussian",bag.fraction=0.75, n.folds=5, 5)
# saveRDS(Abund_Model_Step, paste0('ENSO_', spp, '_testBRT_logCPUE_', region, '.rds'))
# 
# Abund_Model<-Abund_Model_Step[[1]]
# 
# # Check model fit for R2 and RMSE 
# Model_Evals_Abund<- data.frame(matrix(unlist(Abund_Model_Step[[2]]), nrow=length(Abund_Model_Step[[2]]), byrow=TRUE))
# colnames(Model_Evals_Abund)<-c("R2","RMSE")
# 
# print(summary(Model_Evals_Abund[,1]))
# print(summary(Model_Evals_Abund[,2]))
# 
# 
# # Now reduce to 'non-random' predictors
# # Remove everything scoring lower than the Random variable. 
# var_tested<-names(df[,Predictors])
# 
# iters=length(Abund_Model)
# percent_contrib<-NULL#list()
# for(q in 1:iters){                               
#   sum1<-summary(Abund_Model[q][[1]]  , plot=F )
#   sum2<-sum1[order(sum1[,1], levels = var_tested),]
#   percent_contrib<-cbind(percent_contrib, sum2[,2])
#   rownames(percent_contrib)<-sum1[order(sum1[,1], levels = var_tested),1]
# }
# 
# # Calculate the mean contributions of each predictor variable
# Mean_PA_Contributions<-as.data.frame(t(rowMeans(percent_contrib)))
# 
# # List the predictors that score higher than Random
# Predictors_to_Keep_Index<-which(Mean_PA_Contributions>Mean_PA_Contributions$Random)
# 
# Predictors_to_Keep<-Mean_PA_Contributions[,Predictors_to_Keep_Index]
# Reduced_Predictors<-which(colnames(df) %in% colnames(Predictors_to_Keep))
# 
# # Refit model without the <= random predictors
# Abund_Model_Reduced<-fit.brt.n_eval_Balanced(dfSW, gbm.x=Reduced_Predictors, gbm.y= c(Response), lr=0.001, tc=3, family = "gaussian",bag.fraction=0.75, n.folds=5, 5)
# saveRDS(Abund_Model_Reduced, 'ENSO_testBRT_noYrMoLatLon_logCPUE_SW_reduced.rds')

#rm(list=ls())
spp <- 'YFT'
region <- 'NW'

for (j in 1:5) {
  spp <- c('YFT', 'DOL', 'BRZ', 'SWO','BET')[j]
  for (i in 1:5) {
    region <- c('CW', 'NW', 'SW', 'NE', 'SE')[i]
    
    Abund_Model_Reduced <- readRDS(paste0('ENSO_', spp, '_testBRT_logCPUE', region,'.rds'))
    Reduced_Predictors <- which(colnames(df) %in% c("Oxy_1ml", "Oxy_2ml", "Catchability250", "Catchability400", "ONI", "Phase", "PDO", "NPGO", "Random"))
    
    
    # Re-evaluate model fit
    Abund_Model<-Abund_Model_Reduced[[1]]
    
    Model_Evals_Abund<- data.frame(matrix(unlist(Abund_Model_Reduced[[2]]), nrow=length(Abund_Model_Reduced[[2]]), byrow=TRUE))
    colnames(Model_Evals_Abund)<-c("R2","RMSE")
    
    print(summary(Model_Evals_Abund[,1]))
    print(summary(Model_Evals_Abund[,2]))
    # summary(Model_Evals_Abund)[4,]
    
    # Pot variable importance and partial dependence plots.
    var_tested<-names(df[,Reduced_Predictors])
    
    percent_contrib<-NULL
    iters=length(Abund_Model)
    part_plot<-list()
    part_plot<-list()
    percent_contrib<-NULL#list()
    Cont_Preds<-names(Filter(is.numeric,df[,Reduced_Predictors]))
    Num_Preds<-which(var_tested %in% Cont_Preds)
    
    for(q in 1:iters){                                #this was 50 
      mod<-Abund_Model[q][[1]] 
      ###
      part_plot1<-data.frame(row.names=1:100)
      for(x in Num_Preds){ ###
        pp<-plot(mod ,var_tested[x],return.grid=T) ###
        part_plot1<-cbind(part_plot1, pp) ###
      }###
      
      ###
      part_plot[[q]]<-part_plot1 ###
      
      sum1<-summary(Abund_Model[q][[1]]  , plot=F )
      sum2<-sum1[order(sum1[,1], levels = var_tested),]
      percent_contrib<-cbind(percent_contrib, sum2[,2])
      rownames(percent_contrib)<-sum1[order(sum1[,1], levels = var_tested),1]
    }
    All_percent_contribution<-cbind(rownames(percent_contrib), paste(round(rowMeans(percent_contrib),2), round(rowSds(percent_contrib),2), sep=" ± "))
    Combined_All_percent_contribution<-All_percent_contribution
    
    # Make spider plot of variable importance
    Mean_Abund_Contributions<-as.data.frame(t(rowMeans(percent_contrib)))
    Abund_Predictors_Plot<- rbind(rep(max(Mean_Abund_Contributions),length(var_tested)) , rep(0,length(var_tested)) , Mean_Abund_Contributions)
    Abund_Predictors_Plot[]<-sapply(Abund_Predictors_Plot, as.numeric)
    
    # png('ENSO_radar_noYrMoLatLon_logCPUE_SW.png')
    # par(mfrow=c(1,1))
    # radarchart(Abund_Predictors_Plot,  pfcol=rgb(0.0,0.3,0.5,0.5), pcol=rgb(0.0,0.3,0.5,0.5), title="Bigeye All" )
    # dev.off()
    
    Variable_List<-as.data.frame(t(Mean_Abund_Contributions))
    Variable_List$Variables<-rownames(Variable_List)
    Variable_List$SD <- rowSds(percent_contrib)
    Variable_List<-Variable_List[order(-Variable_List$V1),]
    
    Num_Preds<-which(rownames(Variable_List) %in% Cont_Preds)
    
    # And a point plot which I prefer
    p <- ggplot(Variable_List, aes(y=reorder(Variables, V1), x=V1)) + 
      geom_point(shape=21, size=5, fill='#03959f', color='black') + 
      geom_errorbar(aes(xmin=V1-SD, xmax=V1+SD), width=0.5) +
      theme_bw() + 
      xlab('Variable Importance') + 
      ylab('Predictor Variable') + 
      ggtitle(paste(spp, ',', region))
    ggsave(paste0('ENSO_', spp, '_varImp_logCPUE_', region, '.png'), plot = p, height = 7, width=7)
    
    
    # Make the partial dependence plots
    ncols <- 3
    nrows <- ceiling(length(Num_Preds)/ncols)
    png(paste0('ENSO_', spp, '_partDep_logCPUE_', region, '.png'), width=10.5, height=nrows*2.5, units='in', res=300)
    par(mfrow=c(nrows,ncols), mar=c(5, 3, 0.5, 0.5))
    mn_part_plot<-list()  
    for(y in Num_Preds){
      id<-which(colnames(part_plot[[1]])==Variable_List$Variables[y])
      all1<-NULL
      all2<-NULL
      for(z in 1:iters){											 #this was 50 
        all1<-rbind(all1, cbind(c(part_plot[[z]][,id])))
        all2<-rbind(all2, cbind(c(part_plot[[z]][,id+1])))
      }
      all3<-cbind(all1, all2)
      all1<-all3[order(all3[,1]),]
      
      plot(all1, xlab=Variable_List$Variables[y], col="white", ylab=paste("f(",Variable_List$Variables[y], ")", sep=""),cex.axis=1.2, cex.lab=1.2) #, ylim=c(-8,2))
      plx<-predict(loess(all1[,2] ~ all1[,1], span = 0.3), se=T)
      mn_part_plot[[y]]<- cbind(all1[,1], plx$fit)      
      lines(all1[,1],plx$fit)
      lines(all1[,1],plx$fit - qt(0.975,plx$df)*plx$se, lty=2)#0.975
      lines(all1[,1],plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
      rug(na.omit(unlist(df[Variable_List$Variables[y]])))
      legend("bottomright", paste(All_percent_contribution[which(All_percent_contribution[,1]==Variable_List$Variables[y]),2],"%", sep=" "), bty="n", cex=1.4)
    }
    dev.off()
  }
}


#-----------------------------------------------------------
#---PA-----------
#-----------------------------------------------------------
for (j in 1:5) {
  spp <- c('YFT', 'DOL', 'BRZ', 'SWO','BET')[j]
  for (i in 1:5) {
    region <- c('CW', 'NW', 'SW', 'NE', 'SE')[i]
    #region <- 'NW'
    PA_Model_Reduced <- readRDS(paste0('ENSO_', spp, '_testBRT_PA_', region,'.rds'))
    Reduced_Predictors <- which(colnames(df) %in% c("Oxy_1ml", "Oxy_2ml", "Catchability250", "Catchability400", "ONI", "Phase", "PDO", "NPGO", "Random"))
    
    PA_Model<-PA_Model_Reduced[[1]]
    
    Model_Evals_PA<-unlist(unlist(PA_Model_Reduced[[2]]))
    Model_PA_Eval<-matrix(,length(PA_Model),2)
    
    for (i in 1:length(PA_Model)){
      Model_PA_Eval[i,1]<-Model_Evals_PA[[i]]@auc
      Model_PA_Eval[i,2]<-max(Model_Evals_PA[[i]]@TPR+Model_Evals_PA[[i]]@TNR-1)
    }
    
    print(summary(Model_PA_Eval[,1]))
    print(summary(Model_PA_Eval[,2]))
    # summary(Model_PA_Eval)[4,]
    
    
    #recalculate variable importance for the reduced model
    #
    var_tested<-names(df[,Reduced_Predictors])
    
    percent_contrib<-NULL
    iters=length(PA_Model)
    part_plot<-list()
    part_plot<-list()
    percent_contrib<-NULL#list()
    Cont_Preds<-names(Filter(is.numeric,df[,Reduced_Predictors]))
    Num_Preds<-which(var_tested %in% Cont_Preds)
    
    for(q in 1:iters){                                #this was 50 
      mod<-PA_Model[q][[1]] 
      ###
      part_plot1<-data.frame(row.names=1:100)
      for(x in Num_Preds){ ###
        pp<-plot(mod ,var_tested[x],return.grid=T) ###
        part_plot1<-cbind(part_plot1, pp) ###
      }###
      
      ###
      part_plot[[q]]<-part_plot1 ###
      
      sum1<-summary(PA_Model[q][[1]]  , plot=F )
      sum2<-sum1[order(sum1[,1], levels = var_tested),]
      percent_contrib<-cbind(percent_contrib, sum2[,2])
      rownames(percent_contrib)<-sum1[order(sum1[,1], levels = var_tested),1]
    }
    All_percent_contribution<-cbind(rownames(percent_contrib), paste(round(rowMeans(percent_contrib),2), round(rowSds(percent_contrib),2), sep=" ± "))
    Combined_All_percent_contribution<-All_percent_contribution
    
    
    Mean_PA_Contributions<-as.data.frame(t(rowMeans(percent_contrib)))
    PA_Predictors_Plot<- rbind(rep(max(Mean_PA_Contributions),length(var_tested)) , rep(0,length(var_tested)) , Mean_PA_Contributions)
    PA_Predictors_Plot[]<-sapply(PA_Predictors_Plot, as.numeric)
    
    # par(mfrow=c(1,1))
    # radarchart(PA_Predictors_Plot,  pfcol=rgb(0.0,0.3,0.5,0.5), pcol=rgb(0.0,0.3,0.5,0.5), title="BET P/A" )
    
    Variable_List<-as.data.frame(t(Mean_PA_Contributions))
    Variable_List$Variables<-rownames(Variable_List)
    Variable_List<-Variable_List[order(-Variable_List$V1),]
    Variable_List$SD <- rowSds(percent_contrib)
    
    
    Num_Preds<-which(rownames(Variable_List) %in% Cont_Preds)
    
    # And a point plot which I prefer
    p <- ggplot(Variable_List, aes(y=reorder(Variables, V1), x=V1)) + 
      geom_point(shape=21, size=5, fill='#03959f', color='black') + 
      geom_errorbar(aes(xmin=V1-SD, xmax=V1+SD), width=0.5) +
      theme_bw() + 
      xlab('Variable Importance') + 
      ylab('Predictor Variable') + 
      ggtitle(paste(spp, region, 'PA', sep=', '))
    ggsave(paste0('ENSO_', spp, '_varImp_PA_', region, '.png'), plot = p, height = 7, width=7)
    
    
    # Make the partial dependence plots
    ncols <- 3
    nrows <- ceiling(length(Num_Preds)/ncols)
    png(paste0('ENSO_', spp, '_partDep_PA_', region, '.png'), width=10.5, height=nrows*2.5, units='in', res=300)
    par(mfrow=c(nrows,ncols), mar=c(5, 3, 0.5, 0.5))
    mn_part_plot<-list()  
    for(y in Num_Preds){
      id<-which(colnames(part_plot[[1]])==Variable_List$Variables[y])
      all1<-NULL
      all2<-NULL
      for(z in 1:iters){											 #this was 50 
        all1<-rbind(all1, cbind(c(part_plot[[z]][,id])))
        all2<-rbind(all2, cbind(c(part_plot[[z]][,id+1])))
      }
      all3<-cbind(all1, all2)
      all1<-all3[order(all3[,1]),]
      
      plot(all1, xlab=Variable_List$Variables[y], col="white", ylab=paste("f(",Variable_List$Variables[y], ")", sep=""),cex.axis=1.2, cex.lab=1.2) #, ylim=c(-8,2))
      plx<-predict(loess(all1[,2] ~ all1[,1], span = 0.3), se=T)
      mn_part_plot[[y]]<- cbind(all1[,1], plx$fit)      
      lines(all1[,1],plx$fit)
      lines(all1[,1],plx$fit - qt(0.975,plx$df)*plx$se, lty=2)#0.975
      lines(all1[,1],plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
      rug(na.omit(unlist(df[Variable_List$Variables[y]])))
      legend("bottomright", paste(All_percent_contribution[which(All_percent_contribution[,1]==Variable_List$Variables[y]),2],"%", sep=" "), bty="n", cex=1.4)
    }
    dev.off()
  }
}



ggplot(df_pred, aes(x=Lon, y=Lat, z=BET_CPUE)) + 
  stat_summary_2d(binwidth=c(1,1)) + 
  #scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev=T) +
  scale_fill_viridis_c() +
  theme_bw() + 
  borders('world2', fill='black', xlim=c(180,230), ylim=c(10,40)) + 
  ggtitle('Model prediction anomaly', subtitle = 'No year or location, cropped area') + 
  labs(x = "Longitude (°)", y = "Latitude (°)") + 
  coord_map("ortho", orientation = c(0, median(df_pred$Lon), 0))
