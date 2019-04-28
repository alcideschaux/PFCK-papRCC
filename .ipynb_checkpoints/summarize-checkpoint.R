# UNIVARIATE ANALYSIS
# Summary of numeric data
summarize_num <- function(data, x) {
    
    x <- enquo(x)
    
    num_x <- data %>%
        summarize(
            N = n(),
            Mean = mean(!! x, na.rm = TRUE), 
            SD = sd(!! x, na.rm = TRUE),
            Median = quantile(!! x, probs = 0.5, na.rm = TRUE),
            IQR = IQR(!! x, na.rm = TRUE),
            Min = quantile(!! x, probs = 0, na.rm = TRUE),
            Max = quantile(!! x, probs = 1, na.rm = TRUE),
            Missing = sum(is.na(!! x))
        )
    
    print(num_x)

}

# Summary of categorical (factor) data
summarize_fct <- function(data, x) {
    
    x <- enquo(x)
    
    fct_x <- data %>%
        count(!! x) %>%
        mutate(Freq = n / sum(n) * 100) %>%
        rename(Levels = !! x, N = n) %>%
        mutate(Freq = round(Freq,  digits = 1))

    print(fct_x)

}

# BIVARIATE ANALYSIS
# Summary of numeric data by categorical (factor) data
summarize_nums <- function(data, x, y, pairwise = TRUE, caption = NULL) {
    
    # Creating the dataframe
    x <- enquo(x)
    y <- enquo(y)
    df <- data %>% 
        select(x = !! x, y = !! y) %>%
        mutate(y = factor(y))
    
    # Creating the summary table
    tbl <- df %>% 
        group_by(y) %>% 
        summarize(
            N = n(),
            Mean = mean(x, na.rm = TRUE), 
            SD = sd(x, na.rm = TRUE),
            Median = quantile(x, probs = 0.5, na.rm = TRUE),
            IQR = IQR(x, na.rm = TRUE),
            Min = quantile(x, probs = 0, na.rm = TRUE),
            Max = quantile(x, probs = 1, na.rm = TRUE),
            Missing = sum(is.na(x))
        ) %>%
        rename(
            Levels = y
        )
    
    # Testing for associations
    # Mann-Whitney if y == 2, Kruskal-Wallis if y > 2
    if (nlevels(df$y) > 2) {
        test_bi <- kruskal.test(x ~ y, df)
    }
    
    if (nlevels(df$y) == 2) {
        test_bi <- wilcox.test(x ~ y, df)
    }
    
    # Testing for pairwise associations if y > 2
    if (pairwise == TRUE) {
        test_pairwise <- with(df, pairwise.wilcox.test(x, y, p.adjust.method = "bonferroni"))
    }
    
    # Caption
    if (!is.null(caption)) { 
        paste("##", caption, "##\n") %>% cat()
        cat("\n")
    }
    
    # Results
    print(tbl)
    print(test_bi)
    if(nlevels(df$y) > 2 & pairwise == TRUE) print(test_pairwise)

}

# Summary of a contingency table between 2 factor variables x and y
# useNA = c("no", "ifany", "always")
# prop = c("cols", "rows")
summarize_fcts <- function(data, x, y, useNA = "ifany", prop = c("cols", "rows")) {
    
    x <- enquo(x)
    y <- enquo(y)
    df <- data %>% 
        select(x = !! x, y = !! y)
    
    # Creating the crosstables with frequencies and percentages
    tbl <- table(df$x, df$y, useNA = useNA)
    prop_cols <- round(prop.table(tbl, margin = 2) * 100, digits = 1)
    prop_rows <- round(prop.table(tbl, margin = 1) * 100, digits = 1)
    
    # Testing for associations
    test <- with(df, chisq.test(x, y))
    
    # Printing the tables
    print(tbl)
    
    if("cols" %in% prop) {
        cat("\n Percentage by columns")
        print(prop_cols)
    }
    
    if("rows" %in% prop) {
        cat("\n Percentage by rows")
        print(prop_rows)
    }
    
    print(test)

}

# Creating a function for risk ratios and P values
RR <- function(data, x, y){
    
    x <- enquo(x)
    y <- enquo(y)
    df <- data %>% 
        select(x = !! x, y = !! y)
    
    rr_tbl <- epitools::riskratio(df$x, df$y, verbose = TRUE)
    exact <- Exact::exact.test(table(df$x, df$y), to.plot = FALSE)
    print(rr_tbl$data)
    print(exact)
    print(round(rr_tbl$measure, 2))
}

# A function to estimate odds ratios using binomial logistic regression
OR <- function(data, formula, ...) {
  
  mdl <- glm(as.formula(formula), data = data, family = "binomial", ...)
  mdl_coef <- broom::tidy(mdl) %>% 
    mutate(
      estimate = round(estimate, digits = 2),
      std.error = round(std.error, digits = 2),
      statistic = round(statistic, digits = 2),
      or = round(exp(estimate), digits = 2)
    )
  mdl_confint <- exp(broom::confint_tidy(mdl)) %>% round(digits = 2)
  mdl_tbl <- bind_cols(mdl_coef, mdl_confint)
  
  print(mdl_tbl[-1 , ])

}

# PLOTS
# A function to serial plots of a numeric y varible over a factor x variable by a group variable
plot_serial <- function(data, x, y, group, color = NULL) {
    
    x = enquo(x)
    y = enquo(y)
    group = enquo(group)
    color = enquo(color)
    
    g <- ggplot(data, aes(x = !! x, y = !! y, group = !! group, color = !! color)) +
        geom_line(size = 1, alpha = 0.5) +
        geom_point(size = 3, alpha = 0.5) +
        labs(x = "", y = "")
    
    return(g)
}

# A function to compare paired values in sequential measurements (plots and test)
# x = periods of the time series (v.g., pre and post)
# y = values to be compared
# group = grouping variable (v.g., patient's ID)
compare_ts <- function(data, x, y, group) {
    
    x = enquo(x)
    y = enquo(y)
    group = enquo(group)
    df <- data %>% 
        select(x = !! x, y = !! y, group = !! group)
    
    # Plot
    g <- df %>% 
        ggplot(aes(x = x, y = y, group = group)) +
        geom_point(size = 3, alpha = 0.5) +
        geom_line(size = 1, alpha = 0.5) + 
        labs(x = "", y = "")
    
    # Wilcox test for associations
    df_test <- df %>%
        select(group, x, y) %>% 
        spread(key = x, value = y)
    
    colnames(df_test) <- c("Group", "Pre", "Post")
    
    test <- with(df_test, wilcox.test(Pre, Post, paired = TRUE))
    
    # Printing the results
    print(g)
    print(test)
}

# A function for ploting sequential values of the marker
# When 2 periods are compared, the Wilcoxon test is used
# When > 2 periods are provided, the Kruskal-Wallis test is used
# x = periods of the time series (v.g., pre and post)
# y = values to be compared (they will be scaled)
# group = grouping variable (v.g., patient's ID)
plot_ts <- function(data, x, y, group) {
    
    x = enquo(x)
    y = enquo(y)
    group = enquo(group)
    df <- data %>% 
        select(x = !! x, y = !! y, group = !! group)

    # Wilcox test for associations
    if (nlevels(df$x) == 2) {
        df_test <- df %>%
            select(group, x, y) %>% 
            spread(key = x, value = y)
        test <- wilcox.test(df_test[[2]], df_test[[3]], paired = TRUE)
    }
    
    if (nlevels(df$x) > 2) {
        test <- kruskal.test(y ~ x, data = df)
    }
    
    p_test <- test$p.value %>% 
            formatC(digits = 2)
    p <- paste0("P=", p_test)

    # Plot
    g <- df %>% 
        ggplot(aes(x = x, y = y, group = group)) +
        geom_line(size = 1, alpha = 0.5) +
        geom_point(size = 3, alpha = 0.5) +
        labs(x = "", y = "") +
        geom_label(label = p, x = Inf, y = Inf, vjust = 1.5, hjust = 1.5)
        
    # Returning the plot
    return(g)
}

# A function to plot survival curves using the Kaplan-Meier method (P values included)
# and forest plots of hazard ratios (95% confidence intervals included) estimated
# using Cox's proportional hazards regression modeling
summarize_surv <- function(data, time, event, covariate) {
  
  # Building the dataframe for analysis
  time = enquo(time)
  event = enquo(event)
  covariate = enquo(covariate)
  df <- data %>% 
    select(time = !! time, event = !! event, covariate = !! covariate) %>% 
    mutate(event = as.numeric(event))
  
  # Creating models
  surv_mdl <- survival::survfit(Surv(time, event) ~ covariate, data = df)
  coxph_mdl <- survival::coxph(Surv(time, event) ~ covariate, data = df)
  
  # Ploting models
  surv_plot <- survminer::ggsurvplot(fit = surv_mdl, data = df, pval = TRUE) +
    labs(title = "Kaplan-Meier plot for survival curves")
  coxph_plot <- survminer::ggforest(
    model = coxph_mdl, data = df,
    main = "Forest plot for hazard ratios and 95% confidence intervals"
  )
  # 
  # Printing plots
  print(coxph_plot)
  print(surv_plot)

}
