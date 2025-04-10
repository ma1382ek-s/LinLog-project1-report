#### PROJECT 1 ####


library(tidyverse)
library(readxl)
library(plyr)
library(rstatix)
library(patchwork)
library(car)
library(GGally)

carotene_data <- read_excel("Data/carotene.xlsx")

# 1a #

# fitting models
model1a <- lm(betaplasma ~ bmi, data = carotene_data)
model1b <- lm(log(betaplasma) ~ bmi, data = carotene_data)

# creating data set with predictions, residuals etc.
model_predictions <- mutate(carotene_data, 
                            yhat_model1a = predict(model1a),
                            yhat_model1b = predict(model1b),
                            e_model1a = model1a$residuals,
                            e_model1b = model1b$residuals)

# residual plot of model 1a
ggplot(data = model_predictions, aes(x = yhat_model1a, y = e_model1a)) + 
  geom_point(size = 2) +
  geom_hline(yintercept = 0) +
  xlab("Predicted betaplasma levels in model 1a") + 
  ylab("Residual") +
  ggtitle("Residuals in model 1a")

# residual plot of model 1b
ggplot(data = model_predictions, aes(x = yhat_model1b, y = e_model1b)) + 
  geom_point(size = 2) +
  geom_hline(yintercept = 0) +
  xlab("Predicted betaplasma levels in model 1b") + 
  ylab("Residual") +
  ggtitle("Residuals in model 1b")


# Norm QQ plot of model 1a
ggplot(data = model_predictions, aes(sample = e_model1a)) +
  geom_qq(size = 2) + geom_qq_line() +
  ggtitle("Normal Q-Q-plot of residuals in model 1a")

# Norm QQ plot of model 1b
ggplot(data = model_predictions, aes(sample = e_model1b)) +
  geom_qq(size = 2) + geom_qq_line() +
  ggtitle("Normal Q-Q-plot of residuals in model 1b")


# 1 b #

# estimated parameters and confidence intervals in model 1b
model1b$coefficients
confint(model1b)

# plotting real values vs. model and confidence intervals
bmi_seq <- data.frame(bmi = seq(10, 60))
bmi_seq |> mutate(
  fit = predict(model1b, newdata = bmi_seq),
  conf = predict(model1b, newdata = bmi_seq, interval = "confidence"),
  pred = predict(model1b, newdata = bmi_seq, interval = "prediction")
) -> carotene_ints

ggplot(carotene_ints, aes(x = bmi)) + 
  geom_point(data = carotene_data, aes(y = log(betaplasma)), size = 1) +
  geom_line(aes(y = fit), color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = conf[, "lwr"], ymax = conf[, "upr"]), alpha = 0.3) +
  geom_line(aes(y = pred[, "lwr"]), color = "red", linetype = "dashed", linewidth = 0.7) +
  geom_line(aes(y = pred[, "upr"]), color = "red", linetype = "dashed", linewidth = 0.7)

# plotting in the original domain, i.e. not in log(betaplasma)

carotene_ints <- mutate(carotene_ints, 
                  fit_og = exp(fit),
                  conf_og = exp(conf),
                  pred_og = exp(pred))


ggplot(carotene_ints, aes(x = bmi)) + 
  geom_point(data = carotene_data, aes(y = betaplasma), size = 1) +
  geom_line(aes(y = fit_og), color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = conf_og[, "lwr"], ymax = conf_og[, "upr"]), alpha = 0.3) +
  geom_line(aes(y = pred_og[, "lwr"]), color = "red", linetype = "dashed", linewidth = 0.7) +
  geom_line(aes(y = pred_og[, "upr"]), color = "red", linetype = "dashed", linewidth = 0.7) 


# 1c #

# changes in beta-carotene with different changes of bmi

# dubbelkolla med Mattis på onsdag med vilka interval som är rätt
exp(model1b$coefficients[2])
exp(confint(model1b))

exp(-model1b$coefficients[2])
exp(-confint(model1b))

exp(-10*model1b$coefficients[2])
exp(-10*confint(model1b))


# 1d #

# statistical test whether there is a linear relationship between bmi and log plasma

confint(model1b) # bmi interval doesnt include 0 --> reject H_0
summary(model1b) # gives t value, degrees of freedom and p-value


# 2a #

# making smokstat into a categorical variable
mutate(model_predictions,
       smokstat = factor(smokstat,
                          levels = c(1, 2, 3),
                          labels = c("never", "former", "current"))) -> model_predictions


# simple frequency table
table(model_predictions$smokstat)

# getting mean and stddev of (log)betaplasma for each category 
ddply(model_predictions, .(category = smokstat), summarise, avg = mean(betaplasma))
ddply(model_predictions, .(category = smokstat), summarise, stddev = sd(betaplasma))
ddply(model_predictions, .(category = smokstat), summarise, logstddev = mean(log(betaplasma)))
ddply(model_predictions, .(category = smokstat), summarise, logmean = sd(log(betaplasma)))

# box plots
ggplot(model_predictions) +
  geom_boxplot(data = model_predictions , aes(x = smokstat , y = betaplasma))

ggplot(model_predictions) +
  geom_boxplot(data = model_predictions , aes(x = smokstat , y = log(betaplasma)))


# 2b #

model2b <- lm(log(betaplasma) ~ smokstat, data = model_predictions) # "never" is ref
model2a <- lm(log(betaplasma) ~ relevel(smokstat, "current"), data = model_predictions) # "current" is ref
# we see that model 2b is better, perhaps because "never" has more obs.


# 3a #
mutate(model_predictions,
       sex = factor(sex,
                    levels = c(1, 2),
                    labels = c("male", "female")),
       vituse = factor(vituse, 
                       levels = c(1, 2, 3),
                       labels = c("often", "seldom", "no"))) -> model_predictions



# have to motivate which category is best suited as reference category 
# -> 2 things: #obs of ref category should be large, and should also be domain relevant
# as a baseline if applicable

# for sex (male, female) -> choose largest #obs
# for vituse -> choose "no" as baseline unless large diff in #obs
# for smokestat -> choose "never" as baseline unless large diff in #obs

table(model_predictions$sex) # "female" many more #obs
table(model_predictions$vituse) # "often" most #obs, no would make sense in domain


# 3b #

# pairwise correlations:
cor_results <- model_predictions |> select(bmi, age, calories, fat, 
                            cholesterol, fiber, alcohol, betadiet) |>
  cor_test() |>  filter(var1 < var2)

# we get high correlations between calories/fat, calories/cholesterol and cholesterol/fat
filter(cor_results, cor > abs(0.6))


# need other comments on problem with data that we find
plot1 <- ggplot(model_predictions) + geom_point(data = model_predictions, aes(y = calories, x = fat))
plot2 <- ggplot(model_predictions) + geom_point(data = model_predictions, aes(y = calories, x = cholesterol))
plot3 <- ggplot(model_predictions) + geom_point(data = model_predictions, aes(y = fat, x = cholesterol))

plot1 + plot2 + plot3


# person drinking 200 alcoholic drinks a week:
filter(model_predictions, alcohol >= 200) # he consumes 27.9 MJ energy per day 
# how to know if he's extreme in other cases?


# 3c #
model3a <- lm(log(betaplasma) ~ bmi + age + calories + fat + cholesterol + 
                                fiber + alcohol + betadiet + smokstat + sex + vituse,
                                data = model_predictions)

# we will look at VGIF^(1/2f) > 2.24 for the 50 % in R_j 
vif(model3a)
# we see once again that it's fat calories  that are above the threshold

model3c <- lm(log(betaplasma) ~ bmi + age + fat + cholesterol + 
                                fiber + alcohol + betadiet + smokstat + sex + vituse,
                                data = model_predictions)
vif(model3c) # we see that the GVIF values are reasonable, still the highest for fat/cholesterol but ok


# 3d #
model3c$coefficients # beta coefficients
exp(model3c$coefficients) # coefficients in original domain
exp(confint(model3c)) # confidence intervals of coeffs. in original domain

summary(model3c) # using for statistical tests

confint(model3c)


# for test ii)
Ftest1 <- anova(model3c, model1b)
Ftest2 <- anova(model2b, model3c)



# 3e #
model_predictions <- mutate(model_predictions, 
                            yhat_model3c = predict(model3c), 
                            e_model3c = model3c$residuals, 
                            r_model3c = rstudent(model3c),
                            )

highlightcolors <- c("|r*|>3" = "red")


# plotting studentized residuals against predictions
ggplot(model_predictions, aes(x = yhat_model3c, y = r_model3c)) +
  geom_point() +
  geom_hline(yintercept = c(-2, 0, 2)) +
  geom_hline(yintercept = c(-3, 3), linetype = 2) +
  geom_point(data = filter(model_predictions, abs(r_model3c) > 3), 
             aes(color = "|r*|>3"), size = 3) +
  labs(title = "Studentized residuals vs linear predictor",
       subtitle = "model3c",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors) +
  theme(legend.position = "bottom")

# plotting sqrt(abs(studentized residuals)) against predictions
ggplot(model_predictions, aes(x = yhat_model3c, y = sqrt(abs(r_model3c)))) +
  geom_point() +
  geom_hline(yintercept = c(0, sqrt(qnorm(0.75)), sqrt(2))) +
  geom_hline(yintercept = sqrt(3), linetype = 2) +
  geom_point(data = filter(model_predictions, abs(r_model3c) > 3), 
             aes(color = "|r*|>3"), size = 3) +
  labs(title = "Sqrt absolute studentized residuals vs linear predictor",
       subtitle = "model3c",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors) +
  theme(legend.position = "bottom")

# from the two plots above we can conclude that there are not too many extreme outliers,
# and variance seems constant. Also no non-linear patterns can be seen in residuals


# Norm QQ plot of studentized residuals. We see heavy tails, indicating they are not norma
ggplot(data = model_predictions, aes(sample = r_model3c)) +
  geom_qq(size = 2) + geom_qq_line(color = "red") +
  labs(title = "Normal Q-Q-plot of studentized residuals in model3c")


# 3f #

#calculating leverage
model_predictions <- mutate(model_predictions, 
                            v_model3c = hatvalues(model3c))


#identifying data point with highest leverage
highest_r <- filter(model_predictions, r_model3c == max(abs(r_model3c)))

pplus1 <- length(model3c$coefficients)
n <- nobs(model3c)

ggplot(cbind(model_predictions), aes(x = yhat_model3c, y = v_model3c)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1/n) +
  geom_hline(yintercept = 2*pplus1/n, color = "red") +
  labs(title = "Model3c: leverage vs linear predictor",
       caption = "y = 1/n (black) and 2(p+1)/n (red)",
       color = "Highlight") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)


# we see the observation which is the guy drinking 12 vodkas is far from the gravity of X
ggplot(cbind(model_predictions), aes(y = v_model3c, x = alcohol)) +
  geom_point(size = 2) +
  geom_point(data = filter(model_predictions, alcohol > 200), 
             aes(color = "length>200"), size = 3) +
  geom_hline(yintercept = 1/n) +
  geom_hline(yintercept = 2*pplus1/n, color = "red") +
  labs(title = "Model3c: leverage vs linear predictor",
       caption = "y = 1/n (black) and 2(p+1)/n (red)",
       color = "Highlight") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)


# 3g #
model_predictions <- mutate(model_predictions, 
                            D_model3c = cooks.distance(model3c))

# preparing for Cook's distance
f1.model3c <- pplus1
f2.model3c <- model3c$df.residual
cook.limit.model3c <- qf(0.5, f1.model3c, f2.model3c)

# Cook's distance plot
ggplot(model_predictions, aes(yhat_model3c, D_model3c)) + 
  geom_point(size = 3) +
  geom_point(data = filter(model_predictions, alcohol > 200),
             aes(color = "blue"), size = 3) +
  geom_hline(yintercept = cook.limit.model3c, color = "purple") +
  geom_hline(yintercept = 4/n, linetype = 2, color = "red") +
  xlab("Fitted values (log(betaplasma") +
  ylab("D_i") +
  labs(title = "model3c: Cook's D",
       caption = "4/n (dashed), F_0.5, p+1, n-(p+1) (solid)",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors)


# DFBETAS
head(dfbetas(model3c))


model_predictions <- mutate(
  model_predictions, 
  df0 = dfbetas(model3c)[, "(Intercept)"],
  df1 = dfbetas(model3c)[, "bmi"],
  df2 = dfbetas(model3c)[, "age"],
  df3 = dfbetas(model3c)[, "fat"],
  df4 = dfbetas(model3c)[, "cholesterol"],
  df5 = dfbetas(model3c)[, "fiber"],
  df6 = dfbetas(model3c)[, "alcohol"],     # perhaps most interested in this one
  df7 = dfbetas(model3c)[, "betadiet"],
  df8 = dfbetas(model3c)[, "smokstatformer"],
  df9 = dfbetas(model3c)[, "smokstatcurrent"],
  df10 = dfbetas(model3c)[, "sexfemale"],
  df11 = dfbetas(model3c)[, "vituseseldom"],
  df12 = dfbetas(model3c)[, "vituseno"],
)

# impact on intercept
ggplot(model_predictions, aes(x = yhat_model3c, y = df0)) +
  geom_point(size = 2) +
  geom_point(data = filter(model_predictions, abs(r_model3c) > 3),
             aes(color = "|r*|>3"), size = 3) +
  geom_point(data = filter(model_predictions, D_model3c > 0.1),
             aes(shape = "Cook's D>0.1"),
             size = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = sqrt(cook.limit.model3c)*c(-1, 1),
             color = "red") +
  geom_hline(yintercept = 2/sqrt(n)*c(-1, 1),
             color = "red", linetype = "dashed") +
  ylab("DFBETAS_0(i)") +
  xlab("Fitted values") +
  labs(title = "model3c: DFBETAS_0: impact on the intercept",
       caption = "y = sqrt(F_0.5) and 2/sqrt(n)") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)


# looking at data point with largest Cook's distance and plotting in DFBETAS plot
outlier_data <-  (filter(model_predictions, D_model3c == max(model_predictions$D_model3c)))

# DFBETAS plot: The only one that's very affected are: df4 (cholesterol)
ggplot(model_predictions, aes(x = yhat_model3c, y = df4)) +
  geom_point(size = 2) +
  geom_point(data = filter(model_predictions, abs(r_model3c) > 3),
             aes(color = "|r*|>3"), size = 3) +
  geom_point(data = filter(model_predictions, D_model3c > 0.1),
             aes(shape = "Cook's D>0.1"),
             size = 3) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = sqrt(cook.limit.model3c)*c(-1, 1),
             color = "red") +
  geom_hline(yintercept = 2/sqrt(n)*c(-1, 1),
             color = "red", linetype = "dashed") +
  ylab("DFBETAS_4(i)") +
  xlab("Fitted values") +
  labs(title = "model3c: DFBETAS_4: cholesterol",
       caption = "y = sqrt(F_0.5) and 2/sqrt(n)") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = highlightcolors)


# 4a #

no_outlier <- filter(carotene_data, cholesterol < max(carotene_data$cholesterol))

# (just repetitive code to see which factors to choose as reference)
table(no_outlier$sex) # "female" many more #obs
table(no_outlier$vituse) # "often" most #obs
table(no_outlier$smokstat)

mutate(no_outlier,
       sex = factor(sex,
                    levels = c(2, 1),
                    labels = c("female", "male")),
       vituse = factor(vituse, 
                       levels = c(1, 2, 3),
                       labels = c("often", "seldom", "no")),
       smokstat = factor(smokstat,
                         levels = c(1, 2, 3),
                         labels = c("never", "former", "current"))
       ) -> no_outlier


new_model3c <- update(model3c, data = no_outlier)

no_outlier <- mutate(no_outlier,
                    D_new = cooks.distance(new_model3c),
                    yhat_new = predict(new_model3c))


f1.new_model <- pplus1
f2.new_model <- new_model$df.residual
cook.limit.new_model <- qf(0.5, f1.new_model, f2.new_model)
n_new = nobs(new_model)

# new Cook's D plot, we see a difference between the old dataset and this one in that there
# is no significan outlier value anymore. Good!
ggplot(no_outlier, aes(yhat, D)) + 
  geom_point(size = 3) +
  geom_hline(yintercept = cook.limit.new_model, color = "red") +
  geom_hline(yintercept = 4/n_new, linetype = 2, color = "red") +
  xlab("Fitted values (log(betaplasma") +
  ylab("D_i") +
  labs(title = "new model: Cook's D",
       caption = "4/n (dashed), F_0.5, p+1, n-(p+1) (solid)",
       color = "Highlight") +
  scale_color_manual(values = highlightcolors)

# looking at pairwise correlation
cor_results_new <- no_outlier |> select(bmi, age, calories, fat, 
                                           cholesterol, fiber, alcohol, betadiet) |>
  cor_test() |>  filter(var1 < var2)


# there are no apparent differences and we wouldn't expect that when removing one data point
filter(cor_results_new, cor > abs(0.6))
filter(cor_results, cor > abs(0.6))

# no significant difference and we don't expect that either?
vif(new_model3c)
vif(model3c)


# 4b #
stepwise_bic <- step(model3a, direction = "both", k = log(nrow(no_outlier)))
stepwise_bic <- step(model3a, direction = "both", k = log(nrow(model_predictions)))

