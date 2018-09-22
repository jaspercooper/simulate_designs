Gentle introduction to DeclareDesign
================
DVM Bishop
22nd September 2018

For getting set up to use DeclareDesign, see [this
website](https://declaredesign.org/r/declaredesign/).

We’ll be looking at a simplified version of the code that you can find
[here](https://declaredesign.org/library/articles/pretest_posttest.html).
For the sake of clarity, I’ve modified the script in parts.

The script analyses the effectiveness of an experimental design that is
used to assess an intervention effect by comparing the pre- vs post-test
difference score for treated vs untreated cases. For instance, you could
think of this as testing the effectiveness of brain training on a memory
test. We take a sample of people and measure their memory at time 1.
Half of them are then given the brain-training and the others have a
control intervention. Memory is measured again at time 2. If the
training is effective, we expect the improvement in memory from time 1
to time 2 to be higher for the trained than for the untrained people.

Note that using change scores is not the optimal way to analyse
intervention data, but it provides a simple example of how the package
works. Later on we shall see how DeclareDesign makes it easy to compare
different approaches to analysing trial data, in terms of their
efficiency and power.

The user is asked to estimate certain features of the results of the
study, so that data can be simulated with those characteristics. This
can then be used to evaluate how effective the design is in testing the
hypothesis of interest - in this case the treatment effect.

Before we go further, you should include the R packages you need. If you
have not already installed these, you will need to do that first. Go to
Tools|Install Packages on the menu bar and type the names of the three
packages that have ‘library’ commands below. Then run the lines with
those commands.

``` r
library(DeclareDesign)
library(tidyverse)
library(DesignLibrary)
library(beeswarm)
library(knitr)
```

DeclareDesign differentiates four aspects of experimental design that
are not always distinguished. The first of these is the Model of the
world (M). In our example, the Model for the study requires us to
specify a sample size (N), and pre- and post-test scores for those who
do and don’t have the intervention.

It makes life simpler if we think of the test scores as normal random
deviates - i.e., z-scores - where the initial pre-test values have mean
of zero and SD of one.

Regardless of intervention, we expect people’s scores to correlate
between time 1 and time 2: i.e., there will some stability in rank
ordering of people’s memory skills across the two measurement occasions.
You have to provide an estimate of this relationship as a correlation
coefficient.

In addition, an intervention effect will be expected just for the
treated group: you have to provide an estimate of this in standard
deviation units - i.e. as Cohen’s d. In the script this is referred to
as the average treatment effect (`ate`).

In intervention studies, it is very rare for all the participants to
remain in the study up to follow-up. A final parameter to be estimated
is the attrition rate - i.e. the proportion of cases that will be lost
to follow-up.

Let’s start by specifying these sample characteristics for our Model.

``` r
N <- 100 # sample size
ate <- 0.25 # average treatment effect
sd_1 <- 1 # SD for group 1
sd_2 <- 1 # SD for group 2
rho <- 0.5 # correlation between time1 and time2
attrition_rate <- .1 # no outcome data for N*attrition_rate subjects
```

In the next chunk of code, we input these parameters into our model.
When we do this using the `declare_population()` function, we create a
data frame which contains simulated data. Usually, when we run scripts
from DeclareDesign, we would not see the data frame that is gradually
built up as the script proceeds, but to improve understanding of the
functions, we will inspect the simulated data for each chunk, looking
just at the first few rows that are created.

``` r
# M: Model
population <- declare_population(
  N    = N,
  u_t1 = rnorm(N)*sd_1,
  u_t2 = rnorm(N, rho * u_t1, sqrt(1 - rho^2))*sd_2
)

df1 <- population()
kable(head(df1))
```

| ID  |       u\_t1 |       u\_t2 |
| :-- | ----------: | ----------: |
| 001 |   1.1492356 | \-2.2047851 |
| 002 |   0.2321015 |   0.3589651 |
| 003 | \-1.3334377 | \-0.5609404 |
| 004 | \-0.5375723 | \-0.0545800 |
| 005 | \-0.0823834 | \-0.0138572 |
| 006 | \-0.4015819 | \-0.9573552 |

You will see that we have a data frame with four columns. There are 100
rows (or as many as you have specified in `N`), but we just display the
first six. The first column is just a subject identifier. We then have
the two columns `u_t1` and `u_t2`. You can think of these as the
underlying variables from which the observed test scores will be
derived.

We specified that there should be a correlation of `rho` between `u_t1`
and `u_t2`. Let’s see what this looks like:

``` r
mytitle <- paste0('r = ', round(cor(df1$u_t1, df1$u_t2), 3))
plot(df1$u_t1, df1$u_t2, main = mytitle, col = 2, pch = 16)
```

![](DeclareDesignIntro_files/figure-gfm/plotuvars-1.png)<!-- -->

It is unlikely that the correlation will be exactly the one you
specified. That is because we are taking a sample from a population. You
can think of the population as an enormous number of paired observations
of `u_t1` and `u_t2` from which we are just selecting a subset of pairs.
The larger our sample, the closer the observed correlation will be to
the value of `rho` that was specified. Unless you set the random number
generator to give the same answer every time (see [set.seed
command](http://rfunction.com/archives/62)), the specific values of
`u_t1` and `u_t2` will change every time you run the script.

Next, we need to define the potential outcomes of the study – i.e., the
outcome that unit 1 exhibits when assigned to the treatment versus when
they are assigned to the control. We first define the treatment and
control potential outcomes at the pre-test phase, \(Y_t1\). In the chunk
below, you’ll notice that our pretest score potential outcomes,
`Y_t1_Z_0` and `Y_t1_Z_1`, are identical to `u_t1`. Although this may
seem redundant, there are two good reasons to specify the outcomes in
that way. First, specifying the pre-treatment outcomes as identical in
the pre-test period makes explicit our assumption that the treatment
only affects outcomes in the second period. Thus, we’ve explicitly ruled
out anticipation effects.

Second, in general it is helpful to keep separate the original
underlying variables that are generated and the observed variables that
we will enter into the analysis. This provides the flexibility to have
observed variables that are not normally distributed. We could, for
example, generate `Y_t1` as a binary variable through a probit
transformation of `u_t1` by changing the declaration of `Y_t1` in our
population to `draw_binary(N = N,latent = u_t1,link = "probit")`.\[1\]

`Y_t2` is not a simple function of `u_t2`, because it depends on the
intervention. For those who have the intervention, we need to add the
average treatment effect, `ate`, to `u_t2` whenever the treatment
condition, `Z`, is specified as `1`. We have not yet specified which
individuals have treatment, but we create two new columns giving
potential outcomes that tell us what *would* happen if a case were
treated or untreated.

Again, we will inspect the first few rows, and see the new columns
created by the data frame. We’ll also plot the values for the potential
outcomes to observe the treatment effect. Note that the points for the
treated condition are simply increased by an amount equivalent to `ate`.
When we come to the data strategy step below, we will allocate cases as
treated or not treated, but at the model stage, we are just representing
potential outcomes.

``` r
potential_outcomes <- declare_potential_outcomes(
  Y_t1_Z_0 = u_t1,
  Y_t1_Z_1 = u_t1,
  Y_t2_Z_0 = u_t2,
  Y_t2_Z_1 = u_t2 + ate)

df1 <- potential_outcomes(df1)
kable(head(df1))
```

| ID  |       u\_t1 |       u\_t2 | Y\_t1\_Z\_0 | Y\_t1\_Z\_1 | Y\_t2\_Z\_0 | Y\_t2\_Z\_1 |
| :-- | ----------: | ----------: | ----------: | ----------: | ----------: | ----------: |
| 001 |   1.1492356 | \-2.2047851 |   1.1492356 |   1.1492356 | \-2.2047851 | \-1.9547851 |
| 002 |   0.2321015 |   0.3589651 |   0.2321015 |   0.2321015 |   0.3589651 |   0.6089651 |
| 003 | \-1.3334377 | \-0.5609404 | \-1.3334377 | \-1.3334377 | \-0.5609404 | \-0.3109404 |
| 004 | \-0.5375723 | \-0.0545800 | \-0.5375723 | \-0.5375723 | \-0.0545800 |   0.1954200 |
| 005 | \-0.0823834 | \-0.0138572 | \-0.0823834 | \-0.0823834 | \-0.0138572 |   0.2361428 |
| 006 | \-0.4015819 | \-0.9573552 | \-0.4015819 | \-0.4015819 | \-0.9573552 | \-0.7073552 |

We now move to the Inquiry part of the model. Here we just have to
specify the true underlying quantity we want to estimate. In this case,
we specify two. The first is the average treatment effect in the
post-test period (`ATE`), defined as the true mean post-test difference
between the treatment and control potential outcomes of the dependent
variable, `Y`. We also define the true underlying
difference-in-differences (`DiD`) as an estimand, i.e., the true average
difference in pre-post differences between treatment and control. In
fact, because potential outcomes are identical in the pre-test period,
these estimands are equivalent to one another. If we incorporated a
violation of the no-anticipation effects assumption, for example, the
estimands would no longer be equivalent. Specifying our inquiry thus
helps clarify the assumptions on which our inferences are based.

``` r
# I: Inquiry
estimand <- declare_estimand(ATE = mean(Y_t2_Z_1 - Y_t2_Z_0),
                             DiD = mean((Y_t2_Z_1 - Y_t1_Z_1) - (Y_t2_Z_0 - Y_t1_Z_0)))
kable(estimand(df1))
```

| estimand\_label | estimand |
| :-------------- | -------: |
| ATE             |     0.25 |
| DiD             |     0.25 |

The next step is termed the Data Strategy. At this step we determine how
individuals are assigned to groups, and we also take into account any
sample attrition. In effect, we determine which cases will feature in
which groups in the analysis.

By looking at how `df1` changes at each step of the script, we can see
what is being achieved.

First we declare how cases are assigned to treatments (`Z`). The
variable `Z_cond_prob` is the conditional probability that a given case
is in the treatment group, given it was assigned to treatment (and vice
versa for control).\[2\] The default is for 50% to be assigned to each
of two groups. In effect, the script acts to toss a coin that determines
assignment of each case to `Z = 1` or `Z = 0`.

``` r
# D: Data Strategy
assignment <- declare_assignment(prob = .5)

df1 <- assignment(df1)
kable(head(df1))
```

| ID  |       u\_t1 |       u\_t2 | Y\_t1\_Z\_0 | Y\_t1\_Z\_1 | Y\_t2\_Z\_0 | Y\_t2\_Z\_1 | Z | Z\_cond\_prob |
| :-- | ----------: | ----------: | ----------: | ----------: | ----------: | ----------: | -: | ------------: |
| 001 |   1.1492356 | \-2.2047851 |   1.1492356 |   1.1492356 | \-2.2047851 | \-1.9547851 | 0 |           0.5 |
| 002 |   0.2321015 |   0.3589651 |   0.2321015 |   0.2321015 |   0.3589651 |   0.6089651 | 1 |           0.5 |
| 003 | \-1.3334377 | \-0.5609404 | \-1.3334377 | \-1.3334377 | \-0.5609404 | \-0.3109404 | 1 |           0.5 |
| 004 | \-0.5375723 | \-0.0545800 | \-0.5375723 | \-0.5375723 | \-0.0545800 |   0.1954200 | 1 |           0.5 |
| 005 | \-0.0823834 | \-0.0138572 | \-0.0823834 | \-0.0823834 | \-0.0138572 |   0.2361428 | 0 |           0.5 |
| 006 | \-0.4015819 | \-0.9573552 | \-0.4015819 | \-0.4015819 | \-0.9573552 | \-0.7073552 | 0 |           0.5 |

Next we take into account the attrition rate. New columns are added to
the simulation: `R` specifies whether or not the case is lost to
follow-up, according to the attrition rate we specified. Since attrition
is set to only 10%, you probably will only see cases where `R = 1` in
the brief part of the data frame that is shown below, but 10% of cases
will have `R` set to zero. (If you wish you can inspect all values by
typing `df1$R` at the console).

``` r
report <- declare_assignment(prob = (1 - attrition_rate),
                             assignment_variable = R)

df1 <- report(df1)
kable(head(df1))
```

| ID  |       u\_t1 |       u\_t2 | Y\_t1\_Z\_0 | Y\_t1\_Z\_1 | Y\_t2\_Z\_0 | Y\_t2\_Z\_1 | Z | Z\_cond\_prob | R | R\_cond\_prob |
| :-- | ----------: | ----------: | ----------: | ----------: | ----------: | ----------: | -: | ------------: | -: | ------------: |
| 001 |   1.1492356 | \-2.2047851 |   1.1492356 |   1.1492356 | \-2.2047851 | \-1.9547851 | 0 |           0.5 | 1 |           0.9 |
| 002 |   0.2321015 |   0.3589651 |   0.2321015 |   0.2321015 |   0.3589651 |   0.6089651 | 1 |           0.5 | 1 |           0.9 |
| 003 | \-1.3334377 | \-0.5609404 | \-1.3334377 | \-1.3334377 | \-0.5609404 | \-0.3109404 | 1 |           0.5 | 1 |           0.9 |
| 004 | \-0.5375723 | \-0.0545800 | \-0.5375723 | \-0.5375723 | \-0.0545800 |   0.1954200 | 1 |           0.5 | 1 |           0.9 |
| 005 | \-0.0823834 | \-0.0138572 | \-0.0823834 | \-0.0823834 | \-0.0138572 |   0.2361428 | 0 |           0.5 | 1 |           0.9 |
| 006 | \-0.4015819 | \-0.9573552 | \-0.4015819 | \-0.4015819 | \-0.9573552 | \-0.7073552 | 0 |           0.5 | 1 |           0.9 |

Next we need to create the actual values of `Y_t1` and `Y_t2` that apply
with this assignment - for cases where `Z` is 0 we use the values of
`Y_t1_Z_0` and `Y_t2_Z_0`, and for cases where `Z` is 1 we use the value
of `Y_t1_Z_1` and `Y_t2_Z_1`. The declare\_reveal function automatically
reveals the potential outcomes corresponding to a given assignment of
treatment.

``` r
reveal_t1 <- declare_reveal(Y_t1) 
reveal_t2 <- declare_reveal(Y_t2) 

df1 <- reveal_t1(df1)
df1 <- reveal_t2(df1)

kable(head(df1))
```

| ID  |       u\_t1 |       u\_t2 | Y\_t1\_Z\_0 | Y\_t1\_Z\_1 | Y\_t2\_Z\_0 | Y\_t2\_Z\_1 | Z | Z\_cond\_prob | R | R\_cond\_prob |       Y\_t1 |       Y\_t2 |
| :-- | ----------: | ----------: | ----------: | ----------: | ----------: | ----------: | -: | ------------: | -: | ------------: | ----------: | ----------: |
| 001 |   1.1492356 | \-2.2047851 |   1.1492356 |   1.1492356 | \-2.2047851 | \-1.9547851 | 0 |           0.5 | 1 |           0.9 |   1.1492356 | \-2.2047851 |
| 002 |   0.2321015 |   0.3589651 |   0.2321015 |   0.2321015 |   0.3589651 |   0.6089651 | 1 |           0.5 | 1 |           0.9 |   0.2321015 |   0.6089651 |
| 003 | \-1.3334377 | \-0.5609404 | \-1.3334377 | \-1.3334377 | \-0.5609404 | \-0.3109404 | 1 |           0.5 | 1 |           0.9 | \-1.3334377 | \-0.3109404 |
| 004 | \-0.5375723 | \-0.0545800 | \-0.5375723 | \-0.5375723 | \-0.0545800 |   0.1954200 | 1 |           0.5 | 1 |           0.9 | \-0.5375723 |   0.1954200 |
| 005 | \-0.0823834 | \-0.0138572 | \-0.0823834 | \-0.0823834 | \-0.0138572 |   0.2361428 | 0 |           0.5 | 1 |           0.9 | \-0.0823834 | \-0.0138572 |
| 006 | \-0.4015819 | \-0.9573552 | \-0.4015819 | \-0.4015819 | \-0.9573552 | \-0.7073552 | 0 |           0.5 | 1 |           0.9 | \-0.4015819 | \-0.9573552 |

Now we compute the difference score for each individual, which is the
statistic we are using in our analysis. This is computed as the
difference between `Y` at time 1 and `Y` at time 2. The beeswarm plot
below shows the scores for all cases in the two groups. Our specified
attrition rate will determine the proportion of cases lost to follow-up,
shown in the plot as black dots. The mean effect size of treatment is
shown both for the complete sample and the available sample after
excluding those lost to attrition. Note this is unlikely to be exactly
equal to the average treatment effect (`ate`), because the estimate is
based on a sample, whereas `ate` is defined in terms of the
population.

``` r
manipulation <- declare_step(difference = (Y_t2 - Y_t1), handler = fabricate)  
df1 <- manipulation(df1)

kable(head(df1))
```

| ID  |       u\_t1 |       u\_t2 | Y\_t1\_Z\_0 | Y\_t1\_Z\_1 | Y\_t2\_Z\_0 | Y\_t2\_Z\_1 | Z | Z\_cond\_prob | R | R\_cond\_prob |       Y\_t1 |       Y\_t2 |  difference |
| :-- | ----------: | ----------: | ----------: | ----------: | ----------: | ----------: | -: | ------------: | -: | ------------: | ----------: | ----------: | ----------: |
| 001 |   1.1492356 | \-2.2047851 |   1.1492356 |   1.1492356 | \-2.2047851 | \-1.9547851 | 0 |           0.5 | 1 |           0.9 |   1.1492356 | \-2.2047851 | \-3.3540206 |
| 002 |   0.2321015 |   0.3589651 |   0.2321015 |   0.2321015 |   0.3589651 |   0.6089651 | 1 |           0.5 | 1 |           0.9 |   0.2321015 |   0.6089651 |   0.3768636 |
| 003 | \-1.3334377 | \-0.5609404 | \-1.3334377 | \-1.3334377 | \-0.5609404 | \-0.3109404 | 1 |           0.5 | 1 |           0.9 | \-1.3334377 | \-0.3109404 |   1.0224973 |
| 004 | \-0.5375723 | \-0.0545800 | \-0.5375723 | \-0.5375723 | \-0.0545800 |   0.1954200 | 1 |           0.5 | 1 |           0.9 | \-0.5375723 |   0.1954200 |   0.7329923 |
| 005 | \-0.0823834 | \-0.0138572 | \-0.0823834 | \-0.0823834 | \-0.0138572 |   0.2361428 | 0 |           0.5 | 1 |           0.9 | \-0.0823834 | \-0.0138572 |   0.0685263 |
| 006 | \-0.4015819 | \-0.9573552 | \-0.4015819 | \-0.4015819 | \-0.9573552 | \-0.7073552 | 0 |           0.5 | 1 |           0.9 | \-0.4015819 | \-0.9573552 | \-0.5557733 |

``` r
mydiff <- round(mean(df1$difference[df1$Z==1])-mean(df1$difference[df1$Z==0]),3)
mytemp <- filter(df1, R==1)
mean1 <- mean(mytemp$difference[mytemp$Z==1])
mean0 <- mean(mytemp$difference[mytemp$Z==0])
mydiff2 <- round((mean1-mean0),3)
mytitle <- paste0('Average effect without attrition= ', mydiff,'\nAverage effect with attrition=',mydiff2)
beeswarm(df1$difference~df1$Z,xlab='Treatment',pwcol=(df1$R+1),pch=16,main=mytitle,ylab='Memory score difference')

segments(.8, mean0, x1=1.2,y1 = mean0,lty=1,col=1,lwd=2) #add line for mean
segments(1.8, mean1, x1=2.2,y1 = mean1,lty=1,col=1,lwd=2)
text(1.5,-1,'Black dots are\n cases lost to \nfollow-up\n(R = 0)')
```

![](DeclareDesignIntro_files/figure-gfm/beeswarm_plot-1.png)<!-- -->

The final stage is the Answer strategy. This specifies the statistical
procedure used to evaluate the inquiry, using the data simulated from
the model. In the chunk below, we can see the output from the Answer
strategy. You can if you wish view the computed variable `mylm`, which
shows that the output from `estimated_effect(df1)` is identical to the
result obtained by running the code for the robust regression analysis
on `df1`, after excluding those with values of `R = 0` (i.e., the cases
lost to follow-up).

``` r
# A: Answer Strategy
# Use difference score between pretest and posttest
estimated_effect <- declare_estimator(
  difference ~ Z,
  model = lm_robust,
  estimand = c("ATE","DiD"),
  subset = R == 1,
  label = "Change score"
)


kable(estimated_effect(df1)) # print the results of the regression
```

| estimator\_label | term | estimate | std.error | statistic |   p.value |    conf.low | conf.high | df | outcome    | estimand\_label |
| :--------------- | :--- | -------: | --------: | --------: | --------: | ----------: | --------: | -: | :--------- | :-------------- |
| Change score     | Z    | 0.287459 | 0.2249311 |  1.277987 | 0.2046148 | \-0.1595443 | 0.7344624 | 88 | difference | ATE             |
| Change score     | Z    | 0.287459 | 0.2249311 |  1.277987 | 0.2046148 | \-0.1595443 | 0.7344624 | 88 | difference | DiD             |

``` r
mylm <- lm_robust(formula = difference ~ Z, data = subset(df1,R == 1))
```

## Evaluating the design

We’ve worked through the analysis building up simulated data for one
run, but the principal purpose of this exercise is to evaluate the
design by repeatedly running the simulation so that the distribution of
statistics can be obtained.

To do this, you specify the design by creating a list of all the
processes that were involved above, joined by `+`.

The `diagnose_design()` function then computes various diagnostic
statistics for that design. The default number of simulations is 500.
This chunk of code takes time to run, with the amount of time dependent
on the number of simulations.

We can then see the diagnostic statistics for the design. In this script
we’ve assigned the diagnostic information to a variable, `diagnosis1`,
so we can display it differently, but usually you’d just have this
displayed automatically on the console, using the command
`diagnose_design(pretest_posttest_design)`.

``` r
# Design
pretest_posttest_design <- population + potential_outcomes + estimand +
  assignment + reveal_t1 + reveal_t2 + report + manipulation +
  estimated_effect 

diagnosis1 <- diagnose_design(pretest_posttest_design)
kable(reshape_diagnosis(diagnosis1))
```

| Design Label              | Estimand Label | Estimator Label | Term | N Sims | Bias   | RMSE   | Power  | Coverage | Mean Estimate | SD Estimate | Mean Se | Type S Rate | Mean Estimand |
| :------------------------ | :------------- | :-------------- | :--- | :----- | :----- | :----- | :----- | :------- | :------------ | :---------- | :------ | :---------- | :------------ |
| pretest\_posttest\_design | ATE            | Change score    | Z    | 500    | \-0.00 | 0.20   | 0.23   | 0.97     | 0.25          | 0.20        | 0.21    | 0.01        | 0.25          |
|                           |                |                 |      |        | (0.01) | (0.01) | (0.02) | (0.01)   | (0.01)        | (0.01)      | (0.00)  | (0.01)      | (0.00)        |
| pretest\_posttest\_design | DiD            | Change score    | Z    | 500    | \-0.00 | 0.20   | 0.23   | 0.97     | 0.25          | 0.20        | 0.21    | 0.01        | 0.25          |
|                           |                |                 |      |        | (0.01) | (0.01) | (0.02) | (0.01)   | (0.01)        | (0.01)      | (0.00)  | (0.01)      | (0.00)        |

With the design in hand, we can also use it to quickly generate data,
estimates, and estimands. You can try running the code below multiple
times. You’ll get different answers because the experiment is being
simulated from start to finish each
time.

``` r
kable(head(draw_data(pretest_posttest_design)))
```

| ID  |       u\_t1 |       u\_t2 | Y\_t1\_Z\_0 | Y\_t1\_Z\_1 | Y\_t2\_Z\_0 | Y\_t2\_Z\_1 | Z | Z\_cond\_prob |       Y\_t1 |       Y\_t2 | R | R\_cond\_prob |  difference |
| :-- | ----------: | ----------: | ----------: | ----------: | ----------: | ----------: | -: | ------------: | ----------: | ----------: | -: | ------------: | ----------: |
| 001 | \-0.3133212 |   0.1103330 | \-0.3133212 | \-0.3133212 |   0.1103330 |   0.3603330 | 1 |           0.5 | \-0.3133212 |   0.3603330 | 1 |           0.9 |   0.6736542 |
| 002 |   1.3986444 |   1.5162740 |   1.3986444 |   1.3986444 |   1.5162740 |   1.7662740 | 1 |           0.5 |   1.3986444 |   1.7662740 | 1 |           0.9 |   0.3676296 |
| 003 |   0.4303731 |   1.3878644 |   0.4303731 |   0.4303731 |   1.3878644 |   1.6378644 | 0 |           0.5 |   0.4303731 |   1.3878644 | 0 |           0.1 |   0.9574913 |
| 004 |   0.9668249 |   0.4235282 |   0.9668249 |   0.9668249 |   0.4235282 |   0.6735282 | 0 |           0.5 |   0.9668249 |   0.4235282 | 1 |           0.9 | \-0.5432967 |
| 005 | \-1.1597461 | \-1.7222977 | \-1.1597461 | \-1.1597461 | \-1.7222977 | \-1.4722977 | 1 |           0.5 | \-1.1597461 | \-1.4722977 | 1 |           0.9 | \-0.3125515 |
| 006 |   0.4386275 |   0.0512262 |   0.4386275 |   0.4386275 |   0.0512262 |   0.3012262 | 1 |           0.5 |   0.4386275 |   0.3012262 | 1 |           0.9 | \-0.1374013 |

``` r
kable(head(get_estimates(pretest_posttest_design)))
```

| estimator\_label | term |  estimate | std.error | statistic |   p.value |    conf.low | conf.high | df | outcome    | estimand\_label |
| :--------------- | :--- | --------: | --------: | --------: | --------: | ----------: | --------: | -: | :--------- | :-------------- |
| Change score     | Z    | 0.0147479 | 0.2292763 | 0.0643236 | 0.9488585 | \-0.4408906 | 0.4703864 | 88 | difference | ATE             |
| Change score     | Z    | 0.0147479 | 0.2292763 | 0.0643236 | 0.9488585 | \-0.4408906 | 0.4703864 | 88 | difference | DiD             |

``` r
kable(head(get_estimands(pretest_posttest_design)))
```

| estimand\_label | estimand |
| :-------------- | -------: |
| ATE             |     0.25 |
| DiD             |     0.25 |

## Interpreting diagnostic information

To understand the output, it helps to read the [working
paper](https://declaredesign.org/declare.pdf).

I’ve transposed the output and added some explanation in the chunk
below.

``` r
myoutput<-data.frame(t(diagnosis1$diagnosands_df[1,]))
#Add a column for explanation
myoutput$explanation<-'.'
colnames(myoutput)[1]<-'Estimate'
myoutput$explanation[5]<-'Expected difference between estimate and estimand'
myoutput$explanation[7]<-'Root mean squared error'
myoutput$explanation[9]<-'Probability of rejecting null hypothesis of no effect'
myoutput$explanation[11]<-'Probability that estimand falls within conf. interval'
myoutput$explanation[19]<-'Probability that a significant estimate has incorrect sign'
kable(myoutput)
```

|                    | Estimate                  | explanation                                                |
| ------------------ | :------------------------ | :--------------------------------------------------------- |
| design\_label      | pretest\_posttest\_design | .                                                          |
| estimand\_label    | ATE                       | .                                                          |
| estimator\_label   | Change score              | .                                                          |
| term               | Z                         | .                                                          |
| bias               | \-0.0009371753            | Expected difference between estimate and estimand          |
| se(bias)           | 0.009004091               | .                                                          |
| rmse               | 0.2012118                 | Root mean squared error                                    |
| se(rmse)           | 0.007171804               | .                                                          |
| power              | 0.226                     | Probability of rejecting null hypothesis of no effect      |
| se(power)          | 0.01771116                | .                                                          |
| coverage           | 0.97                      | Probability that estimand falls within conf. interval      |
| se(coverage)       | 0.00837647                | .                                                          |
| mean\_estimate     | 0.2490628                 | .                                                          |
| se(mean\_estimate) | 0.009004091               | .                                                          |
| sd\_estimate       | 0.2014112                 | .                                                          |
| se(sd\_estimate)   | 0.007151993               | .                                                          |
| mean\_se           | 0.2104363                 | .                                                          |
| se(mean\_se)       | 0.0006913241              | .                                                          |
| type\_s\_rate      | 0.008849558               | Probability that a significant estimate has incorrect sign |
| se(type\_s\_rate)  | 0.009952575               | .                                                          |
| mean\_estimand     | 0.25                      | .                                                          |
| se(mean\_estimand) | 0                         | .                                                          |
| n\_sims            | 500                       | .                                                          |

## Design library

The goal of this explanation is to give you some understanding of the
processes involved in DeclareDesign; this will make it easier for you to
modify scripts for specific purposes. In practice, however, for the main
types of common design, you can run the same processes with minimal
scripting by using the functions in the DesignLibrary package, which
simplifies the process of using DeclareDesign.

In effect, you just have to specify the key parameters, and the rest is
done for you by the `pretest_posttest_designer()` function. For example,
if you wanted to see what would happen to power if the attrition rate
increased to 20%, you could look at
`diagnose_design(pretest_posttest_designer(attrition_rate = .20))`. Or
you could compare power under three different assumptions about how well
your baseline predicts the endline using
`diagnose_design(expand_design(pretest_posttest_designer,rho =
c(.2,.5,.8)))`.

Note that with the `pretest_posttest_designer()` function, you get
diagnoses for three possible ways of running the analysis. The first is
‘Change score’, which is the same as we have used above - you compare
the pretest-posttest difference for the two treatment conditions. An
alternative approach is to compare the groups on the posttest score,
treating the pretest score as a covariate: this is labelled ‘Condition
on pretest’. The final approach is to just use posttest scores, ignoring
pretest data: this is ‘Posttest only’.

The Design Library is very helpful for understanding strengths and
weaknesses of different
approaches.

``` r
mydes <- pretest_posttest_designer(N = 100, ate = 0.25, sd_1 = 1, sd_2 = 1,
                                   rho = 0.5, attrition_rate = 0.1)
mydeslib <- diagnose_design(mydes)

kable(data.frame(t(mydeslib$diagnosands_df)))
```

|                    | X1           | X2                   | X3            |
| ------------------ | :----------- | :------------------- | :------------ |
| design\_label      | mydes        | mydes                | mydes         |
| estimand\_label    | ATE          | ATE                  | ATE           |
| estimator\_label   | Change score | Condition on pretest | Posttest only |
| term               | Z            | Z                    | Z             |
| bias               | 0.001615089  | \-0.001352337        | \-0.002466738 |
| se(bias)           | 0.007595980  | 0.006661688          | 0.008427578   |
| rmse               | 0.2122084    | 0.1842821            | 0.1988304     |
| se(rmse)           | 0.006303207  | 0.006673300          | 0.007317387   |
| power              | 0.220        | 0.260                | 0.214         |
| se(power)          | 0.01586437   | 0.01671991           | 0.01657741    |
| coverage           | 0.950        | 0.946                | 0.960         |
| se(coverage)       | 0.008839100  | 0.010657932          | 0.009321963   |
| mean\_estimate     | 0.2516151    | 0.2486477            | 0.2475333     |
| se(mean\_estimate) | 0.007595980  | 0.006661688          | 0.008427578   |
| sd\_estimate       | 0.2124148    | 0.1844617            | 0.1990142     |
| se(sd\_estimate)   | 0.006352658  | 0.006711988          | 0.007324738   |
| mean\_se           | 0.2096824    | 0.1828330            | 0.1995765     |
| se(mean\_se)       | 0.0006315406 | 0.0005469638         | 0.0005957276  |
| type\_s\_rate      | 0.000000000  | 0.007692308          | 0.009345794   |
| se(type\_s\_rate)  | 0.000000000  | 0.007191753          | 0.008909758   |
| mean\_estimand     | 0.25         | 0.25                 | 0.25          |
| se(mean\_estimand) | 0            | 0                    | 0             |
| n\_sims            | 500          | 500                  | 500           |

# Session information

Session information is included for technical reasons: it can help
diagnose issues if a script does not run in a given environment.

``` r
sessionInfo()
```

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] bindrcpp_0.2.2           knitr_1.20              
    ##  [3] beeswarm_0.2.3           DesignLibrary_0.1.1.9000
    ##  [5] forcats_0.3.0            stringr_1.3.1           
    ##  [7] dplyr_0.7.6              purrr_0.2.5             
    ##  [9] readr_1.1.1              tidyr_0.8.1             
    ## [11] tibble_1.4.2             ggplot2_3.0.0           
    ## [13] tidyverse_1.2.1          DeclareDesign_0.10.0    
    ## [15] estimatr_0.10.0          fabricatr_0.6.0         
    ## [17] randomizr_0.16.1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.4 haven_1.1.2      lattice_0.20-35  colorspace_1.3-2
    ##  [5] htmltools_0.3.6  yaml_2.2.0       rlang_0.2.2      pillar_1.3.0    
    ##  [9] glue_1.3.0       withr_2.1.2      modelr_0.1.2     readxl_1.1.0    
    ## [13] bindr_0.1.1      plyr_1.8.4       munsell_0.5.0    gtable_0.2.0    
    ## [17] cellranger_1.1.0 rvest_0.3.2      evaluate_0.11    highr_0.7       
    ## [21] broom_0.5.0      Rcpp_0.12.18     scales_1.0.0     backports_1.1.2 
    ## [25] jsonlite_1.5     hms_0.4.2        digest_0.6.16    stringi_1.2.4   
    ## [29] grid_3.5.0       rprojroot_1.3-2  cli_1.0.0        tools_3.5.0     
    ## [33] magrittr_1.5     lazyeval_0.2.1   Formula_1.2-3    crayon_1.3.4    
    ## [37] pkgconfig_2.0.2  xml2_1.2.0       lubridate_1.7.4  assertthat_0.2.0
    ## [41] rmarkdown_1.10   httr_1.3.1       rstudioapi_0.7   R6_2.2.2        
    ## [45] nlme_3.1-137     compiler_3.5.0

1.  Note that, generally speaking, non-linear transformations of `u_t1`
    and `u_t2` might lead to correlations between `Y_t1` and `Y_t2` that
    are not equal to `rho` in expectation.

2.  For example, if the probability of assignment was .8, then units
    assigned to the control would have `Z_cond_prob = .2` and those
    assigned to treatment would have `Z_cond_prob = .8`. This is useful
    for designs where you want to weight observations by their
    assignment probabilities.
