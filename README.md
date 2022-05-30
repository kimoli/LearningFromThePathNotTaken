# Learning from the path not taken
This repository holds the single-trial learning data collected for Kim, Forrence, & McDougle's manuscript "Motor learning without movement", the R script that processes the data and performs the statistical analyses reported in the manuscript, and an "output" folder that contains csv's with all statistical test results and pdf's with related figure panels.


## Dependencies
The R script requires the following packages:

		- rstatix
		- coin
		- MuMIn
		- lmerTest
		- lme4
		- r2glmm
		- emmeans
		- effsize
		- effectsize
		- magrittr
		- ggplot2
		- ggpubr
		- ggeffects


## Compatibility
The R script in this repo was written using R v. 4.0.3 and the RStudio 2021.09.1+372 "Ghost Orchid" Release for Windows.

## Using the script
If you open the R script, you will see a section near the top entitled "USER TODOS". Please go into this section and change the string value of the variable "homedir" to correspond with this repo's directory on your machine.

You can also set the value of "excludeBasedOnInstrRecall" in the "USER TODOS" section to either exclude or include data from participants who did not completely recall the instructions after the online version of the task. Set this variable to 1 to exclude data from participants who did not recall the instructions, or set this variable to any other number to include data from all participants. By default, data from participants who did not fully recall the instructions were excluded, although we note that the key findings of the manuscript are not affected by these exclusions.

After the "USER TODOS" section, there are a couple of other sections that do some legwork (some custom functions for data organization/running repeated tests, load data, etc.) followed by a section devoted to the data analysis for each experiment. The sections/experiments are labelled in the same way as they are in the manuscript. After the data analysis sections, there are a few sections devoted to saving the statistical output and the generated figures. The figures are formatted and the panels are labelled as they are in the manuscript. Please refer to the manuscript for additional details about the experiments and the contents of the figure panels.

I have tried to name the variables and data frame headers as intuitively as possible, but please feel free to reach out to me if anything is unclear or if you have any other questions.


## Contact Info
Olivia Kim

oakim@princeton.edu

kim.olivia.a@gmail.com