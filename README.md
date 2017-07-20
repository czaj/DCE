# DCE
Models for Discrete Choice Experiments

This is a package of Matlab scripts and functions that allow for the estimation of models for Discrete Choice Experiments (i.e., conditional multinomial logit models). 
This package includes the following models:
- Multinomial (conditional) Logit (MNL)
- Mixed (random parameters) Logit (MXL)
- Generalized Multinomial Logit Model (GMXL)
- Latent Class (LC)
- Latent Class Mixed Logit (LCMXL)
- Multiple Indicators Multiple Causes (MIMIC)
- Hybrid Multinomial Logit (HMNL)
- Hybrid Mixed Logit (HMXL)
- Hybrid Latent Class (HLC)

The models are estimated using the maximum likelihood method and work with the following specifications:
- preference or WTP space
- multiple distribution types (for random parameters)
- non-linear transformations of explanatory variables
- covariates of means, scale, and scale variance (where applicable)
- impose equality restrictions or constraints
- flexible data types (panel structure, non-constant number of choice tasks or alternatives per respondent, and missing data)
- various estimation and numerical optimization algorithms and options
- parallel computing
and more.

The codes are published under a Creative Commons Attribution 4.0 License. This means that you are free to use, share, or modify these codes for any purpose, even commercially. What we ask in return is that you acknowledge the source of the codes or reference one of our papers (see czaj.org/research for details).

We provide these codes for two reasons:
- Evolution - feel free to study, apply, extend, and build upon what we have accomplished.
- Efficiency - we have put considerable effort into making the codes fast and efficient. We hope to receive feedback; if you have any suggestions for improving these codes or simply making them more elegant – please let us know

The codes come with no warranty – we try to make them error free and as researcher friendly as possible; however, certain errors may remain. The demos and documentation are rather scant; if you want to use these codes, be prepared to spend a significant amount of time understanding them.

You may want to begin by adding the DCE and Tools folders to your Matlab pats and checking out the DCE_demo folder or refer to http://czaj.org/research/supplementary-materials for software codes that accompany our papers. 

We gratefully acknowledge the help of (in alphabetical order and in addition to registered GitHub contributors): Danny Campbell, Richard Carson, Marek Giergiczny, William Greene, Arnie Hole, Klaus Moeltner, Nada Wasi, Maciej Wilamowski, and Kenneth Train, whose examples, comments or suggestions we followed when working on these codes.