
library(mediation)
library(lme4)
library(dplyr)

set.seed(1)

data("student")

student
student_small = filter(student, SCH_ID < 100)

med.fit <- glmer(attachment ~ catholic + gender + income + pared + (1 + catholic | SCH_ID), family = binomial(link = "logit"), data = student_small)

out.fit <- glmer(fight ~ catholic * attachment + gender + income + pared + (1 + attachment | SCH_ID), family = binomial(link = "logit"), data = student_small)


med.out <- mediate(med.fit, out.fit, treat = "catholic", mediator = "attachment", sims = 500)
summary(med.out)
summary(med.out, output="bygroup")



ranef(med.fit, condVar = TRUE)
ranef(med.fit, condVar = FALSE)

attr(ranef(med.fit, condVar = TRUE)[[1]], "postVar")


