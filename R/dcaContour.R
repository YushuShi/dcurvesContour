#' Plot DCA for survival data in contour plot
#'
#' @param data The dataframe contains patients' information.
#' @param time A sequence of time considered for contour plot.
#' @param threshold A sequence of threshold probabilities considered in DCA.
#' @param model The model used, a coxph object.
#' @param option The plot option, can be "Treat All", "Treat None", "Compare with All", "Compare two models", default is to plot the effect of the model.
#' @param model2 Another model to compare
#'
#' @return a plotly object
#'
#' @examples
#' library(dplyr)
#' library(dcurves)
#' library(survival)
#' library(plotly)
#' df_time_to_cancer_dx <-
#'   readr::read_csv(
#'     file = "https://raw.githubusercontent.com/ddsjoberg/dca-tutorial/main/data/df_time_to_cancer_dx.csv"
#'   ) %>%
#'   labelled::set_variable_labels(
#'     patientid = "Patient ID",
#'     cancer = "Cancer Diagnosis",
#'     ttcancer = "Years to Diagnosis/Censor",
#'     risk_group = "Risk Group",
#'     age = "Patient Age",
#'     famhistory = "Family History",
#'     marker = "Marker",
#'     cancerpredmarker = "Prediction Model",
#'     cancer_cr = "Cancer Diagnosis Status"
#'   )
#'
#'
#' coxmod <- coxph(Surv(ttcancer, cancer) ~ age + famhistory + marker,
#' data = df_time_to_cancer_dx)
#'
#' time<-seq(min(df_time_to_cancer_dx$ttcancer),
#' quantile(df_time_to_cancer_dx$ttcancer,0.9),length.out=21)
#' thresSeq<-seq(0, 0.5, 0.01)
#'
#' tempPlot<-dcaContour(df_time_to_cancer_dx,time,thresSeq,
#' coxmod,"Compare with All")
#' tempPlot
#' coxmod2 <- coxph(Surv(ttcancer, cancer) ~ age,
#' data = df_time_to_cancer_dx)
#' tempPlot2<-dcaContour(df_time_to_cancer_dx,time,thresSeq,
#' coxmod,"Compare two models",coxmod2)
#' tempPlot2
#'
#' @import survival
#' @importFrom plotly plot_ly layout
#' @importFrom dcurves dca
#'
#' @export

dcaContour<-function(data,time,threshold,model,option=NULL,model2=NULL){
  treatAllMat<-matrix(NA,ncol=length(time),nrow=length(threshold))
  treatNoneMat<-matrix(NA,ncol=length(time),nrow=length(threshold))
  treatModelMat<-matrix(NA,ncol=length(time),nrow=length(threshold))
  for(i in 1:length(time)){
    data$tempVar<-1 - summary(survfit(model, newdata = data), times = time[i])$surv[1, ]

    eval(parse(text = paste0("temp<-dca(",sub("(.*) ~.*", "\\1", as.character(model$formula))[2],
                             "~ tempVar,data=data,
               time = time[i],
               thresholds = threshold)")))
    nb<-temp$dca$net_benefit
    if(sum(is.na(nb[(1+2*length(threshold)):(3*length(threshold))]))>0){
      tp_backup<-rep(0,length(temp$dca$tp_rate))
      fp_backup<-temp$dca$pos_rate
      nb_backup<-tp_backup-fp_backup*rep(threshold/(1-threshold),3)
      nb_backup[nb_backup<0]<-0
      print(paste0("At time ",round(time[i],3),
                   ", with threshold value ",
                   threshold[is.na(nb[(1+2*length(threshold)):(3*length(threshold))])],
                   ", there are not enough observations to calculate the survival probability, and we impute NA value with ",
                   nb_backup[is.na(nb[(1+2*length(threshold)):(3*length(threshold))])],
                   "."))
      nb[is.na(nb)]<-nb_backup[is.na(nb)]
      tp<-temp$dca$tp_rate
      tp[is.na(tp)]<-tp_backup[is.na(tp)]
      fp<-temp$dca$fp_rate
      fp[is.na(fp)]<-fp_backup[is.na(fp)]
    }
    treatAllMat[,i]<-nb[1:length(threshold)]
    treatNoneMat[,i]<-nb[(1+length(threshold)):(2*length(threshold))]
    treatModelMat[,i]<-nb[(1+2*length(threshold)):(3*length(threshold))]
  }

  if(!is.null(model2)){
    treatModel2Mat<-matrix(NA,ncol=length(time),nrow=length(threshold))
    for(i in 1:length(time)){
      data$tempVar<-1 - summary(survfit(model2, newdata = data), times = time[i])$surv[1, ]
      eval(parse(text = paste0("temp<-dca(",sub("(.*) ~.*", "\\1", as.character(model$formula))[2],
                               "~ tempVar,data=data,
               time = time[i],
               thresholds = threshold)")))
      nb<-temp$dca$net_benefit
      if(sum(is.na(nb[(1+2*length(threshold)):(3*length(threshold))]))>0){
        tp_backup<-rep(0,length(temp$dca$tp_rate))
        fp_backup<-temp$dca$pos_rate
        nb_backup<-tp_backup-fp_backup*rep(threshold/(1-threshold),3)
        nb_backup[nb_backup<0]<-0
        print(paste0("At time ",round(time[i],3),
                     ", with threshold value ",
                     threshold[is.na(nb[(1+2*length(threshold)):(3*length(threshold))])],
                     ", there are not enough observations to calculate the survival probability, and we impute NA value with ",
                     nb_backup[is.na(nb[(1+2*length(threshold)):(3*length(threshold))])],
                     " for model 2."))
        nb[is.na(nb)]<-nb_backup[is.na(nb)]
        tp<-temp$dca$tp_rate
        tp[is.na(tp)]<-tp_backup[is.na(tp)]
        fp<-temp$dca$fp_rate
        fp[is.na(fp)]<-fp_backup[is.na(fp)]
      }
      treatModel2Mat[,i]<-nb[(1+2*length(threshold)):(3*length(threshold))]
    }
  }
  outplot<-plot_ly(x=time,y=threshold,z=treatModelMat,type="contour",
                   hovertemplate = paste('At time %{x:.2f} <br>with',
                                         'threshold %{y:.2f},<br>the net',
                                         'benefit is %{z:.2f}<extra></extra>'))
  outplot<-outplot %>% layout(title=list(text="Net benefit"))
  if(option=="Treat All"){
    outplot<-plot_ly(x=time,y=threshold,z=treatAllMat,type="contour",
                     hovertemplate = paste('At time %{x:.2f} <br>with',
                                           'threshold %{y:.2f},<br>the net',
                                           'benefit is %{z:.2f}<extra></extra>'))
    outplot<-outplot %>% layout(title=list(text="Net benefit by treating all"))
  }
  if(option=="Treat None"){
    outplot<-plot_ly(x=time,y=threshold,z=treatNoneMat,type="contour",
                     hovertemplate = paste('At time %{x:.2f} <br>with',
                                           'threshold %{y:.2f},<br>the net',
                                           'benefit is %{z:.2f}<extra></extra>'))
    outplot<-outplot %>% layout(title=list(text="Net benefit by treating none"))
  }
  if(option=="Compare with All"){
    outplot<-plot_ly(x=time,y=threshold,z=treatModelMat-treatAllMat,type="contour",
                     hovertemplate = paste('At time %{x:.2f} <br>with',
                                           'threshold %{y:.2f},<br>the net',
                                           'benefit is %{z:.2f}<extra></extra>'))
    outplot<-outplot %>% layout(title=list(text="Net benefit compared with treating all"))
  }
  if(option=="Compare two models"){
    outplot<-plot_ly(x=time,y=threshold,z=treatModelMat-treatModel2Mat,type="contour",
                     hovertemplate = paste('At time %{x:.2f} <br>with',
                                           'threshold %{y:.2f},<br>the difference in',
                                           'net benefit is %{z:.2f}<extra></extra>'))
    outplot<-outplot %>% layout(title=list(text="Net benefit difference of two models"))
  }
  outplot<- outplot %>% layout(xaxis=list(title="Time"),
                               yaxis=list(title="Threshold"))
  outplot
}
