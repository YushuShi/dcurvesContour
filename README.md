# dcurvesContour

## General info
Plot DCA for survival data in contour plot.

## Install

```
library(devtools)
install_github("YushuShi/dcurvesContour")
```
	
## Usage

```
dcaContour(data, time, threshold, model, option = NULL, model2 = NULL)
```
## Arguments

* **data** The dataframe contains patients' information.
* **time** A sequence of time considered for contour plot.
* **threshold** A sequence of threshold probabilities considered in DCA.
* **model** The model used, a coxph object.
* **option** The plot option, can be "Treat All", "Treat None", "Compare with All", "Compare two models", default is to plot the effect of the model.
* **model2** Another model to compare (optional).

## Output
a plotly object

## Examples

```
library(dplyr)
library(dcurves)
library(survival)
library(plotly)
library(dcurvesContour)
df_time_to_cancer_dx <-
  readr::read_csv(
    file = "https://raw.githubusercontent.com/ddsjoberg/dca-tutorial/main/data/df_time_to_cancer_dx.csv"
  ) %>%
  labelled::set_variable_labels(
    patientid = "Patient ID",
    cancer = "Cancer Diagnosis",
    ttcancer = "Years to Diagnosis/Censor",
    risk_group = "Risk Group",
    age = "Patient Age",
    famhistory = "Family History",
    marker = "Marker",
    cancerpredmarker = "Prediction Model",
    cancer_cr = "Cancer Diagnosis Status"
  )


coxmod <- coxph(Surv(ttcancer, cancer) ~ age + famhistory + marker,
data = df_time_to_cancer_dx)

time<-seq(min(df_time_to_cancer_dx$ttcancer),
quantile(df_time_to_cancer_dx$ttcancer,0.9),length.out=21)
thresSeq<-seq(0, 0.5, 0.01)

tempPlot<-DCAcontour(df_time_to_cancer_dx,time,thresSeq,
coxmod,"Compare with All")
tempPlot
coxmod2 <- coxph(Surv(ttcancer, cancer) ~ age,
data = df_time_to_cancer_dx)
tempPlot2<-DCAcontour(df_time_to_cancer_dx,time,thresSeq,
coxmod,"Compare two models",coxmod2)
tempPlot2

```
