# $${\color{lightgreen}Regression.analyses.r}$$

This app allows for automating the multivariable and univariable Cox and logistic regression analyses and subsequent visualization of the results.

The arguments should be placed in the following order:

                                       1 - a path to a semicolon-separated CSV file containing the variables to be analyzed,
                                       2 - a variable in the aforementioned CSV file where sample names are stored,
                                       3 - are the variables to be analyzed continuous ('TRUE', 'FALSE'),
                                       4 - a path to a semicolon-separated CSV file containing clinico-pathological and follow-up data,
                                       5 - a variable in the aforementioned CSV file where sample names are stored,
                                       6 - a vector listing any number of dependent clinico-pathological/follow-up variables to be used in the following format: 'discrete_var1;discrete_var2;continuous_time_var1:discrete_censoring_var1;continuous_time_var2:discrete_censoring_var2',
                                       7 - a vector listing any number of independent discrete variables to be used in the following format: 'independent_discrete_var1;independent_discrete_var2', 'NA' if none,
                                       8 - a vector listing any number of independent continuous variables to be used in the following format: 'independent_continuous_var1;independent_continuous_var2', 'NA' if none,
                                       9 - a vector listing one or two grouping variables to be used in the following format: 'grouping_var1;grouping_var2', 'NA' if none,
                                       10 - the number of threads to be used for parallel computations,
                                       11 - a vector listing the types of regression analyses to be performed, 'coxph' for Cox regression, 'lrm' for logistic regression, separated with a semicolon,
                                       12 - a path to the working directory,
                                       13 - a statistical significance level (alpha) to be used,
                                       14 - should the Benjamini-Hochberg correction for multiple comparisons be used ('TRUE', 'FALSE')
                                       15 - should the main independent variables be calibrated so that their values range from 0 to 1 ('TRUE', 'FALSE')
