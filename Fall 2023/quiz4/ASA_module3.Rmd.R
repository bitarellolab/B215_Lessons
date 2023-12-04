---
  title: 'Module 3: Bayesian Regressions'
author: "[Kevin Ross (Cal Poly)](https://statistics.calpoly.edu/Kevin-Ross) and [Jingchen (Monika) Hu (Vassar)](https://pages.vassar.edu/jihu/)"
format:
  html:
  toc: true
number-sections: true
embed-resources: true
---
  
  
 ```{r}
#| warning: false
#| message: false
#| echo: false

library(knitr)

knitr::opts_chunk$set(warning = FALSE, message = FALSE, echo = TRUE)
```

Before we start...

Modules 1 and 2 introduced the principles of Bayesian analysis in simple settings encountered in introductory statistics courses. **Modules 3 and 4 apply Bayesian principles to fit more sophisticated models encountered in intermediate or advanced undergraduate statistics courses.**
  
  In the interest of time, Modules 3 and 4 assume familiarity with the Bayesian paradigm, and possibly also with commonly used frequentist counterparts (e.g., simple linear regression) if you want to compare the two paradigms. The focus of Modules 3 and 4 is on implementation of Bayesian models: how to set up these models in a Bayesian framework, how to implement model estimation in R and specifically the ```brms``` R package, what to look for in the output and how to interpret results.

# Introduction {#sec-introduction}

One commonly used class of Bayesian models in practice is Bayesian regressions. Just like classical inference with regression models, the choice of regression models depends on the nature of the outcome variable. For example,

1. Continuous outcome variable: **Bayesian linear regression**.

2. Binary outcome variable: **Bayesian logistic regression**.

3. Unordered categorical outcome variable: **Bayesian multinomial logistic regression**.

4. Count outcome variable: **Bayesian Poisson regression**.

In Bayesian inference, model parameters are considered as random and therefore one can use **prior** distributions to quantify the degree of uncertainty in these parameters. With appropriate **Markov chain Monte Carlo (MCMC)** estimation tools, one could arrive at the **posterior** distributions of these parameters, based on which one can answer relevant research questions.

In this Module, we will introduce a sample of the Consumer Expenditure Surveys (CE) which contains 5 variables of 4 different data types. We will spend most of the time in understanding how a Bayesian simple linear regression can be estimated using the ```brms``` package, including how to choose a prior, how to perform MCMC estimation, how to perform MCMC diagnostics, and present a selection of Bayesian inference analyses, with highlights in teaching these topics. Sample scripts for multiple linear regression, logistic regression, multinomial logistic regression, and Poisson regression are included for interested participants for further self-directed exploration.

# The Consumer Expenditure Surveys (CE) (~5 minutes) {#sec-CEdata}

::: {.panel-tabset}

## The CE sample

The [Consumer Expenditure Survey](https://www.bls.gov/cex/) is conducted by the U.S. Census Bureau for the U.S. Bureau of Labor Statistics. It contains data on expenditures, income, and tax statistics about consumer units (CU) across the country and provides information on the buying habits of U.S. consumers. The CE program releases data in two ways: tabular data (aggregated) and micro-level data, also called public-use microdata (PUMD).

We will be working with a CE sample comes from the 2019 Q1 CE PUMD: 5 variables, 5133 consumer units (CU).

Variable \& Description:
  
  - **UrbanRural**: Binary; the urban / rural status of CU: 1 = Urban, 2 = Rural.

- **Income**: Continuous; the amount of CU income before taxes in past 12 months (in USD).

- **Race**: Categorical; the race category of the reference person: 1 = White, 2 = Black, 3 = Native American, 4 = Asian, 5 = Pacific Islander, 6 = Multi-race.

- **Expenditure**: Continuous; CU's total expenditures in last quarter (in USD).

- **KidsCount**: Count; the number of CU members under age 16. 

## Code

Here are R code to load and preview the CE sample.

```{r}
#| label: load-data
#| echo: true

library(tidyverse)

CEdata <- readr::read_csv(file = "CEdata.csv")
CEdata[1:3, ]

```

:::


# Simple linear regression background (~15 minutes) {#sec-SLRbackground}

We start this Module considering a simple linear regression model, where we aim at using ```Income``` to learn and predict ```Expenditure```. Since both ```Income``` and ```Expenditure``` are highly right-skewed, we will first take log-transformation of the two and use ```LogIncome``` and ```LogExpenditure``` as the predictor and outcome variables, respectively.

::: {.panel-tabset}

## Variable transformation

- ```LogIncome```: log of ```Income```, predictor variable

- ```LogExpenditure```: log of ```Expenditure```, outcome variable

## Code

```{r}
CEdata <- CEdata %>%
  mutate(LogIncome = log(Income)) %>%
  mutate(LogExpenditure = log(Expenditure))
```

:::

In simple linear regression models, the normal data model is assumed for the outcome variable, and the mean of the normal data model is linked to the predictor variable through a linear combination. 



## Review of normal model

When you have continuous outcomes, you can use a normal model:
$$
Y_i \mid \mu, \sigma \overset{i.i.d.}{\sim} \textrm{Normal}(\mu, \sigma), \,\,\, i = 1, \cdots, n.
$$ {#eq-normal-model}

This model assumes each observation follows the same Normal density with mean $\mu$ and standard deviation $\sigma$.

Suppose now you have another continuous variable available, $x_i$. And you want to use the information in $x_i$ to learn about $Y_i$.

- $Y_i$ is the log of expenditure of CU's
- $x_i$ is the log of total income of CU's

::: {.panel-tabset}

### Questions
Is the model in @eq-normal-model flexible to include $x_i$?

### Short solution
Since $Y_i$ is record-indexed, we can consider record-indexed $\mu_i$ as well. What remains is how to include $x_i$ as part of $\mu_i$.

:::

## An observation specific mean

We can adjust the model in @eq-normal-model to @eq-SLR-model, where the common mean $\mu$ is replaced by an observation specific mean $\mu_i$:
$$
Y_i \mid \mu_i, \sigma \overset{ind}{\sim} \textrm{Normal}(\mu_i, \sigma), \,\,\, i = 1, \cdots, n.
$$ {#eq-SLR-model}

::: {.panel-tabset}

### Questions
How to link $\mu_i$ and $x_i$?

### Short solution
In simple linear regressions, the linking between $\mu_i$ and $x_i$ is through a linear combination, i.e., $\mu_i$ is a linear function of $x_i$. Later we will introduce other regressions, where the link might not be linear and the outcome variable might not follow a normal model.

:::

## Linear relationship between the mean and the predictor

In simple linear regressions, we use a linear relationship between $\mu_i$ and $x_i$ as follows.

$$
\mu_i = \beta_0 + \beta_1 x_i, \,\,\, i = 1, \cdots, n.
$$ {#eq-mu}

- $x_i$'s are known constants.

- $\beta_0$ and $\beta_1$ are unknown parameters.

Let's first focus on the interpretation of this linear relationship and some terminology.

1. The linear function $\beta_0 + \beta_1 x_i$ is the **expected outcome** with $x_i$.
2. $\beta_0$ is the **intercept**: the **expected outcome** when $x_i = 0$.
3. $\beta_1$ is the **slope**: the increase in the **expected outcome** when $x_i$ increases by 1 unit.
    

In the Bayesian paradigm, we have three general steps in our inference procedure.

1. Assign a prior distribution to $(\beta_0, \beta_1, \sigma)$.
2. Perform posterior simulation through MCMC.
3. Summarize posterior distribution of these parameters.

## The simple linear regression model

To put everything together, we have a linear regression model as follows.

\begin{eqnarray}
Y_i \mid x_i, \beta_0, \beta_1, \sigma \overset{ind}{\sim} \textrm{Normal}(\beta_0 + \beta_1 x_i, \sigma), \,\,\, i = 1, \cdots, n.
\end{eqnarray}

Alternatively, we can write:

\begin{eqnarray}
Y_i &=& \mu_i + \epsilon_i,  \\
\mu_i &=& \beta_0 + \beta_1 x_i, \\
\epsilon_i &\overset{i.i.d.}{\sim}& \textrm{Normal}(0, \sigma)\,\,\, i = 1, \cdots, n.
\end{eqnarray}

::: {.panel-tabset}

### Questions
What assumptions does this model make?

### Short solution
$Y_i$'s are independently following their own normal model where the mean $\mu_i = \beta_0 + \beta_1 x_i$. Each normal model shares the same standard deviation $\sigma$. 

### A figure illustration

![Display of linear regression model. The line represents the unknown regression line $\beta_0 + \beta_1 x$ and the normal curves represent the distribution of the response $Y$ about the line.](figures/Regression_View.png)

### Scatterplot of CE sample

```{r fig.align = "center", size = "footnotesize"}
ggplot(CEdata, aes(x = LogIncome, y = LogExpenditure)) +
  geom_point(size=1) + 
  labs(x = "LogIncome", y = "LogExpenditure") +
  theme_bw(base_size = 15, base_family = "") 
```

:::
  
  # A simple linear regression on the CE sample (~30 minutes) {#sec-SLR-CE}
  
  We can now specify a simple linear regression on the CE sample with ```LogIncome``` as the predictor ($x_i$) and ```LogExpenditure``` as the outcome ($y_i$).

\begin{eqnarray}
Y_i \mid \mu_i, \sigma &\overset{ind}{\sim}& \textrm{Normal}(\mu_i, \sigma), \\
\mu_i &=& \beta_0 + \beta_1 x_i.
\end{eqnarray}


Note that:
  
  1. $Y_i$ is ```LogExpenditure```, and $x_i$ is ```LogIncome```.
2. The intercept $\beta_0$: the **expected** LogExpenditure $\mu_i$ for a CU $i$ that has zero LogIncome (i.e., $x_i = 0$).
3. The slope $\beta_1$: the change in the **expected** LogExpenditure $\mu_i$ when the LogIncome of CU $i$ increases by 1 unit.

## Bayesian inference with the brms package

We will be using the [```brms``` R package](https://cran.r-project.org/web/packages/brms/index.html) to perform MCMC estimation for simple linear regression and other regressions in this Module. The ```brms``` R package includes a good number of regression models and prior choices and it is popular among researchers and practitioners.

```{r}
# make sure to install and load the library
## install.packages("brms")
library(brms)
```

## The model and Bayesian inference framework

\begin{eqnarray}
Y_i \mid \mu_i, \sigma &\overset{ind}{\sim}& \textrm{Normal}(\mu_i, \sigma), \\
\mu_i &=& \beta_0 + \beta_1 x_i.
\end{eqnarray}

- Model parameters: $\{\beta_0, \beta_1, \sigma\}$
  
  - Bayesian inference: 
  - Prior for $\{\beta_0, \beta_1, \sigma\}$
  - Sampling model for $Y_1, \cdots, Y_n \mid \mu_i, \sigma$
  - Posterior for $\{\beta_0, \beta_1, \sigma\}$ using MCMC estimation

::: {.panel-tabset}

### Questions
How many parameters and what are they in this simple linear regression model? What types of prior distributions you would like to give to them and why?
  
  ### Short solution
  There are three parameters in this model: $\beta_0, \beta_1, \sigma$. Unless I have useful prior information, I might want to use priors that would have little impact and let the data drive the inference.

:::
  
  
  
  ## Choosing a prior
  
  The ```brms``` package provides default priors if we do not have specific priors we want to use. To check the default prior provided by the ```brms``` package.

```{r}
brms::get_prior(data = CEdata,
                family = gaussian,
                formula = LogExpenditure ~ 1 + LogIncome)
```

We can also use our own priors, e.g., using the following in the ```brm()``` function, which gives a Normal(0, 10) prior for $\beta_0$, a Normal(0, 10) prior for $\beta_1$, and a Caucy(0, 1) prior for $\sigma$.

```{r, eval = FALSE}
prior = c(prior(normal(0, 10), class = Intercept),
          prior(normal(0, 10), class = b),
          prior(cauchy(0, 1), class = sigma))
```

## MCMC estimation

Important concepts in MCMC estimation:
  
  - ```iter```: the number of MCMC iterations.
- ```warmup```: the number of burn-in iterations (i.e., the number of MCMC iterations to be discarded from the beginning of the MCMC chain).
- ```thin```: thinning rate (i.e., every ```thin```-th draw is kept; e.g., when ```thin = 10```, every 10-th draw is saved.)
- ```chains```: the number of MCMC chains.
- ```seed```: the seed to obtain reproducible results.

```{r results = 'hide'}
SLR_Bayes_fit <- brms::brm(data = CEdata,
                           family = gaussian,
                           formula = LogExpenditure ~ 1 + LogIncome,
                           iter = 5000,
                           warmup = 3000,
                           thin = 1,
                           chains = 1,
                           seed = 720)
```

::: {.panel-tabset}

### Questions
Review the example code above. How many MCMC iterations was run and how many MCMC iterations was saved?
  
  ### Short solution
  There were 5000 MCMC iterations and 3000 for the warmup, meaning the last 2000 iterations were saved. Moreover, since ```thin = 1```, each of the 2000 iterations were saved. 
:::
  
  ## MCMC diagnostics
  
  **Traceplots** show how a model parameter changes along the MCMC chain. A random up-and-down traceplot fully exploring the parameter space indicates little convergence issue.

```{r fig.align = "center"}
bayesplot::mcmc_trace(x = SLR_Bayes_fit,
                      pars = c("b_Intercept", "b_LogIncome", "sigma"))
```

::: {.panel-tabset}

### MCMC convergence issues shown in traceplots

Look for "stickiness" which indicates convergence issues. When it happens, consider running the MCMC chain longer (i.e., increase ```iter```) and using longer burn-in (i.e., increase ```warmup```).

### A figure illustration

![A traceplot showing the parameter ```intercept``` having stickiness across iterations.](figures/traceplot_issue.png)

Reference: [https://stats.stackexchange.com/questions/354981/what-problem-do-these-trace-plots-indicate](https://stats.stackexchange.com/questions/354981/what-problem-do-these-trace-plots-indicate)

:::
  
  **Autocorrelation plots (acf)** show how a model parameter's autocorrelation decreases along the MCMC chain. A fast drop of autocorrelation indicates little convergence issue.

```{r fig.align = "center"}
bayesplot::mcmc_acf(x = SLR_Bayes_fit,
                    pars = c("b_Intercept", "b_LogIncome", "sigma"))
```

::: {.panel-tabset}

### MCMC convergence issues shown in ACF plots

Look for slow decrease of the autocorrelation which indicates convergence issues. When it happens, consider thinning more of the MCMC chain (i.e., increase ```thin```) to get rid of the autocorrelation.

### A figure illustration

![A traceplot showing the parameter ```ts.sim``` having slow decrease of autocorrelation.](figures/ACFplot_issue.png)

Reference: [https://cran.r-project.org/web/packages/BGGM/vignettes/mcmc_diagnostics.html](https://cran.r-project.org/web/packages/BGGM/vignettes/mcmc_diagnostics.html)

:::

# Posterior inference of the fitted simple linear regression (~10 minutes) {#sec-SLR-inference}

We will cover a few inference questions and posterior summaries that can help answer these questions from our fitted simple linear regression on the CE sample.

## A Bayesian credible interval

```{r fig.align = "center"}
post_SLR <- as.data.frame(brms::as_draws_matrix(x = SLR_Bayes_fit))
bayesplot::mcmc_areas(post_SLR, pars = "b_LogIncome", prob = 0.95)
```

```{r}
quantile(post_SLR$b_LogIncome, c(0.025, 0.975))
```

A 95\% **credible interval** for the slope of LogIncome is [0.34, 0.37]. In other words, there is a 95\% posterior probability that the slope of LogIncome is between 0.34 and 0.37. 

```{r}
summary(SLR_Bayes_fit)
```

::: {.panel-tabset}

### Classical inference with the lm package

```{r}
SLR_classical_fit <- lm(formula = LogExpenditure ~ 1 + LogIncome,
                        data = CEdata)

summary(SLR_classical_fit)

confint.default(SLR_classical_fit)
```


A 95\% **confidence interval** of the slope of LogIncome is [0.34, 0.37].

### Questions 
Compare the Bayesian 95\% credible interval vs the classical 95\% confidence interval. Are they similar or are they different? What does each interval represent?

### Short solution
The two intervals are very similar. The classical confidence interval indicates we are 95\% confident that the slope is between 0.34 and 0.37, whereas the Bayesian credible interval represents a 95\% posterior probability that the slope is between 0.34 and 0.37.

:::

## A Bayesian hypothesis test

If we want to evaluate whether the slope of LogIncome is at least 0.35, one can calculate the posterior probability of $Pr$(slope is at least 0.35) as follows. With such a high posterior probability, the results strongly support the claim that the slope of LogIncome is at least 0.35.

```{r}
sum(post_SLR$b_LogIncome >= 0.35)/length(post_SLR$b_LogIncome)
```

## Bayesian prediction and posterior predictive checks

We can perform prediction and posterior predictive checks using the following example R code (refer to Module 2 for some detailed examples).

- **Prediction**: ```predict(SLR_Bayes_fit)```

- **Model checking (posterior predictive check)**: ```pp_check(SLR_Bayes_fit)```

# Highlights in teaching Bayesian simple linear regression {#sec-teaching-highlight}

- This is a multi-parameter model
    - An intercept parameter $\beta_0$
    - A slope parameter $\beta_1$
    - A standard deviation parameter $\sigma$

- The Bayesian inference framework
    - Parameters are **random** variables with distributions that quantify degree of uncertainty
    - Posterior is proportional to likelihood multiplies with prior
    
- Markov chain Monte Carlo (MCMC) estimation
    - Why: The joint posterior of all parameters can be challenging to derive and MCMC estimation helps approximate the joint posterior distribution
    - How: There are many MCMC estimation software to choose from and we focus on the use of the ```brms``` package
    
- MCMC diagnostics
    - Why: We need to evaluate whether the MCMC estimation has done a reasonable job of approximating the joint posterior distribution
    - How: MCMC diagnostics tools (e.g., traceplots and autocorrelation plots) can help us check if there are convergence issues
    
- Inference
    - Credible intervals vs confidence intervals, hypothesis testing
    - Prediction
    - Model checking
    - Comparison to classical inference
    
# Additional regressions {#sec-additiona-regressions}

This section is for later reference.

## Sample script for multiple linear regression

Multiple linear regression is an extension of simple linear regression where now instead of only one predictor variables, we can use multiple predictor variables. The outcome variable is still assumed to be independently normally distributed at its record-indexed mean. Here we include an example where ```LogExpenditure``` is regressed on ```LogIncome``` and binary ```UrbanRural```.

::: {.panel-tabset}

### Code for choosing your own prior

```{r, eval = FALSE}
MLR_Bayes_fit_1 <- brm(data = CEdata,
                       family = gaussian,
                       formula = LogExpenditure ~ 1 + LogIncome + as.factor(UrbanRural),
                       prior = c(prior(normal(0, 10), class = Intercept),
                                 prior(normal(0, 10), class = b),
                                 prior(cauchy(0, 1), class = sigma)),
                       iter = 5000,
                       warmup = 3000,
                       thin = 1,
                       chains = 1,
                       seed = 129)
```

:::

::: {.panel-tabset}


### Code for using default prior

```{r, eval = FALSE}
MLR_Bayes_fit_2 <- brm(data = CEdata,
                       family = gaussian,
                       formula = LogExpenditure ~ 1 + LogIncome + as.factor(UrbanRural),
                       iter = 5000,
                       warmup = 3000,
                       thin = 1,
                       chains = 1,
                       seed = 129)
```

:::

## Logistic regression

Logistic regression is used for binary outcome variables. Check out [Albert and Hu (2019) Chapter 12.4](https://bayesball.github.io/BOOK/bayesian-multiple-regression-and-logistic-models.html#bayesian-logistic-regression) for details of logistic regression (MCMC estimation is done through JAGS).

Here we include sample ```brms``` script for running a logistic regression on:

- **LogExpenditure**: predictor variable

- **UrbanRural**: outcome variable (minus 1)

We also include sample ```glm``` script for running a classical inference with logistic regression.

::: {.panel-tabset}

Default prior choice

### brms script

```{r, eval = FALSE}
logistic_Bayes_fit <- brms::brm(data = CEdata,
                                family = binomial(link = "logit"),
                                formula = (UrbanRural - 1) | trials(n) ~ 1 + LogExpenditure,
                                iter = 5000,
                                warmup = 3000,
                                thin = 1,
                                chains = 1,
                                seed = 257)
```

### Classical inference with the glm package

```{r, eval = FALSE}
logistic_classical_fit <- glm(formula =  (UrbanRural - 1) ~ 1 + LogExpenditure,
                              family = binomial(link = "logit"),
                              data = CEdata)

summary(logistic_classical_fit)

confint.default(logistic_classical_fit)
```

:::

## Multinomial logistic regression

Multinomial logistic regression is an exention of logistic regression where the outcome variable is categorical with more than two levels.

Here we include sample ```brms``` script for running a multinomial logistic regression on:

- **LogIncome**: predictor variable

- **Race**: outcome variable

We also include sample ```nnet``` script for running a classical inference with multinomial logistic regression.

::: {.panel-tabset}

Default prior choice

### brms script

```{r, eval = FALSE}
multi_logistic_Bayes_fit <- brms::brm(data = CEdata,
                                      family = categorical(link = "logit"),
                                      Race ~ 1 + LogIncome,
                                      iter = 5000,
                                      warmup = 3000,
                                      thin = 1,
                                      chains = 1,
                                      seed = 843)
```

### Classical inference with the nnet package

```{r, eval = FALSE}
library(nnet)
multi_logistic_classical_fit <- nnet::multinom(formula =  Race ~ 1 + LogIncome,
                                               data = CEdata)
```
:::

## Poisson regression

Poisson regression is a class of generalized linear models where the outcome variable is count. 

Here we include sample ```brms``` script for running a multinomial logistic regression on:

- **LogIncome** \& **LogExpenditure**: predictor variables

- **KidsCount**: outcome variable

We also include sample ```glm``` script for running a classical inference with Poisson regression.

::: {.panel-tabset}

Default prior choice

### brms script

```{r, eval = FALSE}
Poisson_Bayes_fit <- brms::brm(data = CEdata,
                               family = poisson(link = "log"),
                               formula = KidsCount ~ 1 + LogIncome + LogExpenditure,
                               iter = 5000,
                               warmup = 3000,
                               thin = 1,
                               chains = 1,
                               seed = 853)
```

### Classical inference with the glm package

```{r eval = FALSE}
Poisson_classical_fit <- glm(formula =  KidsCount ~ 1 + LogIncome + LogExpenditure,
                             family = poisson(link = "log"),
                             data = CEdata)
```
:::

# Additional resources {#sec-additional-resources}

- [Simple linear regression with JAGS](https://bayesball.github.io/BOOK/simple-linear-regression.html)

- [Multiple linear regression with JAGS](https://bayesball.github.io/BOOK/bayesian-multiple-regression-and-logistic-models.html#bayesian-multiple-linear-regression)

- [Logistic regression with JAGS](https://bayesball.github.io/BOOK/bayesian-multiple-regression-and-logistic-models.html#bayesian-logistic-regression)

- [Multiple linear regression and logistic regrssion with brms](https://bayesball.github.io/BRMS/multiple-regression-and-logistic-models.html)