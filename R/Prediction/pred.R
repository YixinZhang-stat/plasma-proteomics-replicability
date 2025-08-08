
### final model: DL1 + P1 +log(sample size)

## load data
rm(list = ls())

load("data/train.rda")
load("data/test.rda")


#######
# Scatter plots: quantile normalized rho v.s. data metrics
######

library(ggplot2)

DL_names <- paste0("DL", 1:3)
P_names  <- paste0("PLOD", 1:4)
name.AL7 <- c(DL_names, P_names)

term7 <- c(name.AL7, "sample.l")
par(mfrow = c(2, 4), 
    mar = c(4.5, 4.5, 2, 1)
    )
for (var in name.AL7) {
  plot( train[[var]],
        train$z,
        xlab = var,
        ylab = expression("Quantile-normalized "*rho),
        pch = 16)
}
plot( train$sample.l,
      train$z,
      xlab = expression(log("sample size")),
      ylab = expression("Quantile-normalized "*rho),
      pch = 16)
par(mfrow = c(1, 1))

#############
# -------- fit model: DL1 + P1 --------
############
formula_str <- as.formula("class~ DL1+PLOD1+sample.l")

response <- "class"

#------------weighted option---------
case_w <- nrow(train) / sum(train[[response]] == "cases")
ctrl_w <- nrow(train) / sum(train[[response]] == "controls")
train_w <- ifelse(train[[response]] == "cases", case_w, ctrl_w)
test_w <- ifelse(test[[response]] == "cases", case_w, ctrl_w)

fit_w <- glm(formula_str, data = train, 
           family = quasibinomial(link = "logit"), 
           weights = train_w)

summary(fit_w)

test_prob_w <- predict(fit_w, newdata = test, type = "response")

# ---------- AUC ----------
library(pROC)
par(bty = "l")  # keep axis only

test_auc_w <- roc(test$b ~ test_prob_w,smooth=FALSE,ci=T)

roc_w <- data.frame(
  specificity = test_auc_w$specificities,
  sensitivity = test_auc_w$sensitivities
)
auc_text <- paste0("AUC = ", round(auc(test_auc_w), 3))
plot_roc_w <- ggplot(roc_w, aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(color = "black", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  annotate("text", x = 0.9, y = 0.1, label = auc_text, size = 4) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  coord_equal() +  ##used to fix the ratio of height and width
  theme_classic(
     base_size = 13
    )
plot_roc_w

#---------calibration curve----------
test$pred_w <- test_prob_w
test$bin_w <- cut(test_prob_w, 
                  breaks = quantile(test_prob_w, probs = seq(0, 1, 0.05)), ##quantile equal sample size (1000+)
                  include.lowest = TRUE, labels = FALSE)

###----- calculate phat and predicted mean
bins <- sort(unique(na.omit(test$bin_w)))
calib_result_w <- data.frame()

for (bi in bins) {
  idx <- which(test$bin_w == bi)
  #----- calculate phat via glm ----
  try <- glm(class~1, data = test[idx,],
             family = quasibinomial(link = "logit"),
             weights = test_w[idx])
  phat <-unique( predict(try, newdata = test[idx,],
                         type = "response"))
  se <- as.numeric(sqrt(vcov(try)))  
  logit_phat <- unique(predict(try, newdata = test[idx,],type = "link") ) 
  lCI <- plogis(logit_phat - 1.96 * se)
  uCI <- plogis(logit_phat + 1.96 * se)
  
  
  pred_mean = mean(test$pred_w[idx])
  n <- length(idx)
  
  
  calib_result_w <- rbind(calib_result_w, data.frame(
    bin = bi,
    pred_mean = pred_mean,
    obs_mean = phat,
    lower = lCI,
    upper = uCI,
    sample = n
  ))
}

##########---- Calibration curve------
plot_cali_w <- ggplot(calib_result_w, aes(x = pred_mean, y = obs_mean)) +
                     geom_smooth(method = "loess", se = TRUE, 
                                span = 0.75, # default span
                                color = "black") +
                     geom_abline(slope = 1, 
                           intercept = 0, 
                           linetype = "dashed", 
                           color = "red") +
                     labs(x = "Predicted probability", 
                          y = "Observed proportion") +
                     coord_equal() +  ##used to fix the ratio of height and width
                     theme_classic(base_size = 13)
plot_cali_w


## #--------Simulation: Possible cases----
case1_dl <- c(0, 0, 1); case1_p <- c(1, 0, 0, 0)
case2_dl <- c(1, 0, 0); case2_p <- c(0, 0, 0, 1)
panel1_dl <- c(831, 353, 277)/1461; panel1_p <- c(922, 228, 201, 110)/1461
panel2_dl <- c(1110, 171, 178)/1459; panel2_p <- c(90, 479, 343, 547)/1459

sample_sizes <- seq(10000, 100000, by = 10000)

generate_structure <- function(dl, p, s, label) {
  cbind(
    class = NA,
    DL1 = dl[1],
    DL2 = dl[2],
    DL3 = dl[3],
    PLOD1  = p[1],
    PLOD2  = p[2],
    PLOD3  = p[3],
    PLOD4  = p[4],
    sample = s,
    sample.l = log(s),
    label = label
     )
}

# 
sim_list <- list()
labels <- c("Scenario1", "Scenario2", "panel1", "panel2")


for (s in sample_sizes) {
  sim_list[[paste0("Scenario1_", s)]]  <- generate_structure(case1_dl, case1_p, s, "Scenario1")
  sim_list[[paste0("Scenario2_", s)]]  <- generate_structure(case2_dl, case2_p, s, "Scenario2")
  sim_list[[paste0("Panel1_", s)]] <- generate_structure(panel1_dl, panel1_p, s, "Panel1")
  sim_list[[paste0("Panel2_", s)]] <- generate_structure(panel2_dl, panel2_p, s, "Panel2")
}

name_parts <- do.call(rbind, strsplit(names(sim_list), "_"))

sim <- as.data.frame(do.call(rbind, sim_list))
sim$sample_group <- as.numeric(name_parts[, 2])
sim$label   <- factor(sim$label)
num_vars <- setdiff(names(sim), c("class", "label"))
sim[num_vars] <- lapply(sim[num_vars], as.numeric)



## ----- Get predicted 95%CI via logit link \beta ----
# get the linear prediction and se

pred_link_w <- predict(fit_w, newdata = sim, type = "link", se.fit = TRUE)

# logit 95%CI
z_val <- qnorm(0.975)  # 95% CI
link_lower <- pred_link_w$fit - z_val * pred_link_w$se.fit
link_upper <- pred_link_w$fit + z_val * pred_link_w$se.fit

# Prob 95%CI
sim$prob_fit_w <- plogis(pred_link_w$fit)
sim$prob_lower_w <- plogis(link_lower)
sim$prob_upper_w <- plogis(link_upper)


plot_pred_w<- ggplot(sim, aes(x = sample_group, y = prob_fit_w, color = label, group = label)) 
plot_pred_w<- plot_pred_w + geom_line(linewidth = 1.2) + 
  geom_point(size = 2) + 
  geom_vline(xintercept = 16886, ## current sample size
             linetype = "dashed", 
             color = "brown", size = 0.8) +
  labs(
    x = "Sample Size", 
    y = "Predicted Overall Replicability Probability") +
  theme_classic(
     base_size = 13
    ) +
  theme(
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = c(0.95, 0.15),  # 
    legend.justification = c("right", "bottom"),
    legend.background = element_blank(),     #  legend 
    legend.key = element_blank(),            #  legend 
    legend.title = element_blank(),          # legend title
   )

print(plot_pred_w)


#------------unweighted option---------
test_w<- rep(1, nrow(test))

fit <- glm(formula_str, data = train, 
             family = quasibinomial(link = "logit"))

test_prob <- predict(fit, newdata = test, type = "response")
summary(fit)  
  # ---------- AUC ----------
  
test_auc <- roc(test$class, test_prob)
roc <- data.frame(
    specificity = test_auc$specificities,
    sensitivity = test_auc$sensitivities
  )
  
  auc_text <- paste0("AUC = ", round(auc(test_auc), 3))
  plot_roc <- ggplot(roc, aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(color = "black", linewidth = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    annotate("text", x = 0.9, y = 0.1, label = auc_text, size = 4) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    coord_equal() +  
    theme_classic(
      base_size = 13
    )
  plot_roc
  
  
  # ---------- Calibration plot ----------
  test$pred <- test_prob
  test$bin <- cut(test_prob, 
                  breaks = quantile(test_prob, probs = seq(0, 1, 0.05)), ##quantile equal sample size (1000+)
                  include.lowest = TRUE, labels = FALSE)
  
  
  ###----- calculate phat and predicted mean
  bins <- sort(unique(na.omit(test$bin)))
  calib_result <- data.frame()
  
  for (bi in bins) {
    idx <- which(test$bin == bi)
    #----- calculate phat via glm ----
    try <- glm(class~1, data = test[idx,],
               family = quasibinomial(link = "logit"),
               weights = test_w[idx])
    phat <-unique( predict(try, type = "response"))
    se <- as.numeric(sqrt(vcov(try)))  
    logit_phat <- unique(predict(try, type = "link") ) 
    lCI <- plogis(logit_phat - 1.96 * se)
    uCI <- plogis(logit_phat + 1.96 * se)
    
    
    pred_mean = mean(test$pred[idx])
    n <- length(idx)
    
    
    calib_result <- rbind(calib_result, data.frame(
      bin = bi,
      pred_mean = pred_mean,
      obs_mean = phat,
      lower = lCI,
      upper = uCI,
      sample = n
    ))
  }
  
  ## local smooth:3 times loess,mean-upper-lower 

  plot_cali <- ggplot(calib_result, aes(x = pred_mean, y = obs_mean)) +
    geom_smooth(method = "loess", se = TRUE, 
                span = 0.75, # default span
                color = "black") +
     geom_abline(slope = 1, 
                intercept = 0, 
                linetype = "dashed", 
                color = "red") +
    labs(
    x = "Predicted probability", y = "Observed proportion") +
    coord_equal() +  ##used to fix the ratio of height and width
    theme_classic(
      base_size = 13
    )
  plot_cali

  ###  -----Simulated cases-----
  
  ## ----- Get predicted 95%CI via logit link \beta ----
  pred_link <- predict(fit, newdata = sim, type = "link", se.fit = TRUE)
  
  # logit 95%CI
  z_val <- qnorm(0.975)  # 95% CI
  link_lower <- pred_link$fit - z_val * pred_link$se.fit
  link_upper <- pred_link$fit + z_val * pred_link$se.fit
  
  # Prob 95%CI
  sim$prob_fit <- plogis(pred_link$fit)
  sim$prob_lower <- plogis(link_lower)
  sim$prob_upper <- plogis(link_upper)
 
  
  plot_pred<- ggplot(sim, aes(x = sample_group, y = prob_fit, color = label, group = label)) 
  plot_pred<- plot_pred + geom_line(linewidth = 1.2) + 
    geom_point(size = 2) + 
    geom_vline(xintercept = 16886, ## current sample size
               linetype = "dashed", 
               color = "brown", size = 0.8) +
    labs(
      x = "Sample Size", 
      y = "Predicted Overall Replicability Probability") +
    theme_classic(base_size = 13) +
    theme(
      panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = c(0.95, 0.15),  # 
      legend.justification = c("right", "bottom"),
      legend.background = element_blank(),     #  legend 
      legend.key = element_blank(),            #  legend 
      legend.title = element_blank(),          # legend title
        )
  
  print(plot_pred)
  
  

  
  
 