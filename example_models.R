#Package requirements
require(data.table)
require(sandwich)
require(lmtest)

#Read in input file
analysis <- fread(" ")

#List of outcome phenotypes
phenotypes <- c("Cancer")

#Prepare output data.frame: WF = within-family (within-sibship) estimates.
output <- data.frame(PHEN=phenotypes, BETA_MODEL1_0=NA, BETA_MODEL2_0=NA, BETA_TOTAL=NA, BETA_BF=NA, BETA_WF=NA,
                     SE_BETA_MODEL1_0=NA, SE_BETA_MODEL2_0=NA, SE_BETA_TOTAL=NA, SE_BETA_BF=NA, SE_BETA_WF=NA,
                     P_BETA_MODEL1_0=NA, P_BETA_MODEL2_0=NA, P_BETA_TOTAL=NA, P_BETA_BF=NA, P_BETA_WF=NA,
                     VCV_MODEL1_0=NA, VCV_MODEL1_0_TOTAL=NA, VCV_MODEL1_TOTAL=NA, 
                     VCV_MODEL2_0=NA, VCV_MODEL2_0_BF=NA, VCV_MODEL2_0_WF=NA, VCV_MODEL2_BF=NA, VCV_MODEL2_BF_WF=NA, VCV_MODEL2_WF=NA)



#Extract exposure (height) family-mean exposure and outcome data
merge <- data.table(FID=analysis$FID, OUTCOME=analysis$Cancer), EXPOSURE=analysis$Height_std, AGE=analysis$Age, SEX=analysis$Sex, FAM_MEAN=ave(analysis$Height_std, analysis$FID, FUN=mean))
merge2 <- na.omit(merge[,CENTREDEXPOSURE:=EXPOSURE-FAM_MEAN])

# Run unified regression for observed height on cancer using population and within-sibship models

fit_pop <- glm(formula = OUTCOME ~ EXPOSURE, family = "binomial", data = merge2)
fit_ws <- glm(formula = OUTCOME ~ FAM_MEAN + CENTREDEXPOSURE, family = "binomial", data = merge2)

#Population model
output$BETA_MODEL1_0[i] <- fit_pop$coefficients[1]
output$BETA_TOTAL[i] <- fit_pop$coefficients[2]

# save the variance covariance matrix
vcv_matrix <- vcovCL(fit_pop, cluster=merge2$FID)
	
if(  is.na(output$BETA_MODEL1_0[i]) | is.na(output$BETA_TOTAL[i])) {
        output$VCV_MODEL1_0[i] <-NA
        output$VCV_MODEL1_0_TOTAL[i] <-NA
        output$VCV_MODEL1_TOTAL[i] <-NA
} else {
        output$VCV_MODEL1_0[i] <- vcv_matrix[1,1]
        output$VCV_MODEL1_0_TOTAL[i] <- vcv_matrix[1,2]
        output$VCV_MODEL1_TOTAL[i] <- vcv_matrix[2,2]
 
}

#Derive the clustered SEs for the total effect and P-values


test_matrix <- coeftest(fit_pop, vcov.=vcv_matrix)
	
if(  is.na(output$BETA_MODEL1_0[i]) | is.na(output$BETA_TOTAL[i])) {
        output$SE_BETA_MODEL1_0[i] <- NA
        output$SE_BETA_TOTAL[i] <- NA
        output$P_BETA_MODEL1_0[i] <- NA
        output$P_BETA_TOTAL[i] <- NA
} else {
        output$SE_BETA_MODEL1_0[i] <- test_matrix[1,2] 
        output$SE_BETA_TOTAL[i] <- test_matrix[2,2]
        output$P_BETA_MODEL1_0[i] <- test_matrix[1,4] 
        output$P_BETA_TOTAL[i] <- test_matrix[2,4] 

}

#Within-sibship model

output$BETA_MODEL2_0 <- fit_ws$coefficients[1]
output$BETA_BF[i] <- fit_ws$coefficients[2]
output$BETA_WF[i] <- fit_ws$coefficients[3]
 

# save the variance covariance matrix
vcv_matrix = vcovCL(fit_ws, cluster=merge2$FID)
    if(  is.na(output$BETA_MODEL2_0[i]) | is.na(output$BETA_BF[i]) | is.na(output$BETA_WF[i]) ) {
        output$VCV_MODEL2_0[i] <-NA
        output$VCV_MODEL2_0_BF[i] <-NA
        output$VCV_MODEL2_0_WF[i] <-NA
        output$VCV_MODEL2_BF[i] <-NA
        output$VCV_MODEL2_BF_WF[i] <-NA
        output$VCV_MODEL2_WF[i] <-NA
    } else {
        output$VCV_MODEL2_0[i] <- vcv_matrix[1,1]
        output$VCV_MODEL2_0_BF[i] <- vcv_matrix[1,2]
        output$VCV_MODEL2_0_WF[i] <- vcv_matrix[1,3]
        output$VCV_MODEL2_BF[i] <- vcv_matrix[2,2]
        output$VCV_MODEL2_BF_WF[i] <- vcv_matrix[2,3]
        output$VCV_MODEL2_WF[i] <- vcv_matrix[3,3]
    }

    # save the clustered SE's and corresponding p-values
    test_matrix <- coeftest(fit_ws, vcov.=vcv_matrix)
    if(  is.na(output$BETA_MODEL2_0[i]) | is.na(output$BETA_BF[i]) | is.na(output$BETA_WF[i]) ) {
        output$SE_BETA_MODEL2_0[i] <- NA
        output$SE_BETA_BF[i] <- NA
        output$SE_BETA_WF[i] <- NA
        output$P_BETA_MODEL2_0[i] <- NA
        output$P_BETA_BF[i] <- NA
        output$P_BETA_WF[i] <- NA
    } else {
        output$SE_BETA_MODEL2_0[i] <- test_matrix[1,2] 
        output$SE_BETA_BF[i] <- test_matrix[2,2] 
        output$SE_BETA_WF[i] <- test_matrix[3,2] 
        output$P_BETA_MODEL2_0[i] <- test_matrix[1,4] 
        output$P_BETA_BF[i] <- test_matrix[2,4] 
        output$P_BETA_WF[i] <- test_matrix[3,4] 
    }
   



