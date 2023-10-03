# AbundanceConditionalOccupancy


![Satellite image of harvest](images/harvest_image_satellite.png)



This model is for a single species, and will be replicated for 6 different bird species. They are all migratory boreal songbirds and I expect they will react differently to the environmental conditions reflected in the model both due to the species distribution and regional density, their habitat preferences, and their vocal/territorial behaviours.

The data was gathered using ARUs (Autonomous Recording Units) and processed during the first few hours around dawn. All detected individuals of each species of interest is tagged (WTSP, TEWA, REVI, OSFL, YRWA, RCKI).

This single-species integrated abundance-occupancy model described below addresses the different ecological processes that affect the measures of occupancy and abundance.

The two observations we have are: detection of presence/non-detection (related to occupancy through the detection probability) and individual (unmarked) count at each site for 10-15 occasions, or visit. The visits are not menat to be independant replicates: they all occur within a few days and there are several per morning.

The question we are attempting to answer with this model is: *Is the recovery rate of forest songbirds affected when there is a small amount of retention in a harvest?* This question arises from the common wisdom that an ylevel of retention in a harvest is beneficial. However, the scale of the retention patches in this study are much smaller than ever assessed before (~1-5% of the harvest area). So, we are seeing if this theory holds up to small-scale retention, which is prevalent across the province.

**Study design**
To assess recovery, we have a gradient of 22 years of harvests (data was gathered in the last ~3 years), and in all major commercial forests (mixedwood, aspen, upland spruce, pine), across all FMA's (related to geographical timber allotments granted by the provincial government to differnet forestry companies).
- 404+ sites sampled, considered indepent (more than minimum 300m away from each other)
- In equal parts no retention and retention of as small as 1 tree, up to 20,000 m^2 (still fairly small). 
- The visits (1 minute recordings) are transcribed in the online sound processing software Wildtrax and each unique individual per species is identified & tagged.

*Occupancy*

Occupancy is related to whether the species is present or not given the environmental conditions. The assumptions are that, because of the regional spread of the data, one of the largest factor in the species occupancy/distribution is latitude and/or longitude. Additionally, we might consider forest type.
Occupancy in all versions but one of this model is a latent discrete variable, where we assume the observed species occurrence is related to but not equal to the true site occupancy. Occupancy (Royle and Nichols, 2003), is generally defined as the probability that *at least one* individual is present.
 - There is a case to be made that becuase we have 10-15 visits, the conditions that satisty occupancy are easily met. I think that is true for most cases, but because we sampled form end of May to early July, we might have arrived too early at some sites.

 *Abundance*

Abundance is related to how many individuals occupy or use the area around the sampling unit (the ARU), given that the species is present. While we have an observed count of individuals during each visit, the count varies between visits. Conventionally, a multi-visit count is often condensed into a max (maximum ever encountered in one visit, integer), a mean (float), or an intensiy of use (total number of times a unique individual is detected summed for all visits, then divided by survey effort). In this model, I have used the max, but that is easily interchangeable.
The factors presumably influencing this measure of density is the suitability of the habitat in supporting multiple individuals. This is assumed to be related to the size of the retention patch, the age of the harvest, and the proportion of each forest type in the sampling area (or a related measure of green-up/recovery)

*Measuring detection error*

Detection plays a large role in this model at both steps. The probability of detection for a species of interest at the occupancy level of the model is related to Julian date and time-of-day. In this model, this stage of detection rate will be the Probability of detection (p).

Detection of individuals, at the abundance level of the model, is related to two factors: the bird's movement and singing within its territory, and the distance at which we are sampling the bird when it sings and whether its cue is audible given the meteorological and forested conditions. At this stage the variance caused by heterogeneous detection probability between individuals can be estimated through the method described in Rossman et al., (2016) (see bottom of this page).


**First step: Occupancy and detection probability:**

$z_i$ is the true estimated occupancy state.

$$
z_i \sim \text{Bernoulli}(\psi_i)
$$

and the estimated parameter $\psi$ can vary by site in relation to spatial covariates, as in this linear expression:

$$
\text{logit}(\psi_i) = \beta_0 + \beta_1 \times \text{latitude}_i + \beta_2 \times \text{longitude}_i + \text{forest types}_i
$$

where $z \in \{0, 1\}$, is either occupied or not and is estimated, and $\text{logit}(\psi)$ is the prior, and $\text{forest type}$ is the sum of the the coefficients and observations of the differnet forest types.

From the observed samples, this true occupancy state is related by:

$$
y_{ij} \sim \text{Bernoulli}(z_i \times p_{ij})
$$

or, can we write it this way:

$$
y | z \sim \text{Bernoulli}(p(z))
$$

where $y_{ij}$ is the **observed** presence/absence for site $i$ at visit $j$, and $p_{ij}$ is the detection probability at site $i$, visit $j$ given $z_i$. $y_{ij} = 1$ if $N_{ij} > 0$.

$$
\text{logit}(p_{ij}) = \beta_0 + \beta_1 \times \text{time of day} + \beta_2 \times \text{Julian date}
$$


The *detection probability* is defined as the probability of detecting the species given that it is present:

$$
\text{p} = Pr(\text{y_i} = 1 | \text{z_i} = 1)
$$



**Second step: Abundance, as it relates to Occupancy:**

The bird count observations $N_{ij}$, is at the visit-level. For now, it will be summarized as a mean count at the site-level, as $M_i$. The local observed abundance, an observed quantity $N_i$, is conditional on occupancy and is used to estimate the unobserved Poisson rate $\lambda_i$:

$$
[M_i | \psi_i] \sim \text{Poisson}(\lambda_i \times z_i)
$$

$$
\log(\lambda_i) = \alpha_0 + \beta_{\text{size}} \times \text{size} + \beta_{\text{age}} \times \text{age}
$$

where the parameter of the Poisson distribution is dependent on the true occupancy state (when the species is present).


**Implementation in Stan:**

Stan does not support sampling discrete parameters (Stan user guide, Chap. 7). The latent occupancy state, z, is an estimated integer. Although Stan does not directly support this (which is possible in BUGS/JAGS), similar sampling is possible through marginalizing out the latent discrete parameters.

*Marginalisation on the dicrete latent parameters*
1. Occupancy

In estimate a discrete variable (occupancy), the HMC algorithm in Stan fundamentally can't sample from a non-continuous distribution. Stan's work around is to marginalise on the discrete parameter and sample from a continuous distribution that is derived from but independent to the discrete parameter. To marginalise over a discrete parameter, sum over the likelihoods or the joint distribution for all its possible values. We are targetting the marginalisation of the Bernoulli-distributed latent discrete occupancy parameter $z$, given the data $y$, and the continuous parameter $(\psi)$. The joint distribution can be written as:

$$
p(data|\psi, z) = p(data|\psi, z) * p(\psi) * p(z)
$$

where p(data| $\psi$, z) is the likelihood of the data given the parameters, p($\psi$) is the prior on $\psi$ (in our case the log-linked linear expression), and p(z) is the Bernoulli distribution of z.

To marginalise out z, the sum of all possible values of z is:

$$
p(data, \psi) = \sum_{z \in \{0, 1\}} p(data| \psi, z) * p(\psi) * p(z)
$$

This will provide Stan a continuous distribution in terms of the data and the parameter $\psi$ from which to sample.

In the model block, Stan implements the marginalisation on a discrete latent variable with the notation target +=. The target keyword represents the logarithm (because stan operates with log-likelihoods) of the joint density (or posterior density, depending on context) being modeled. The += operation increments the value of the target distribution. 


2. Abundance

Stan has a built-in Zero-Inflated Poisson models, described in Chapter 5 of the user manual. A ZIP is commonly used when there are different processes leading the zero's and the positive counts in an abundance (sounds like occupancy and detection error!). Aka, there are multiple mechanisms behind a non-detection:
- species is present but no cue is given
- species is present, cue is given but not detected
- species is not present

The probabilities of observing 0 and non-0 are described by two different distributions. 

$$ y_n ~
  \begin{cases}
    0       & \quad \text{with probability } \theta \\
    Poisson(y_n | \lambda)  & \quad \text{with probability } 1 - \theta
  \end{cases}
$$

Where $\theta$ is the probability of occupancy
Because Stan does not support sampling conditional on some parameter (with ~), we consider the corresponding likelihoods:

$$ p(y_n | \theta, \lambda) =
  \begin{cases}
    \theta + (1 - \theta) * Poisson(\theta | \lambda)       & \quad \text{if } y_n = 0 \\
    (1 - \theta) * Poisson(y_n | \lambda)  & \quad \text{if} y_n > 0
  \end{cases}
$$

In this mixture model where $\lambda \in \{0, 1\}$ (*is this right?*), each component of the mixture will be estimated with effect data sizes of $\theta$ $N$ and (1 - $\theta$) $N$. 

We need to decide whether the local observed abundance, M_i, is equal to the latent abundance or needs to be used to estimate the true abundance. The observations we have is the local abundance at visit j. The parameters of interest in the model are not at the visit-level, but at the site level (retention size, habitat, etc). The counts are not equal across visits for the same sites; this indicates we are better off treating those as an imperfect observations of the true local latent population size. How do we relate the visit-level abundance to an estimate latent site-level abundance?

**Detection**

Detection probability, $p$,  related to the ocucpancy of a site as $y_{ij}$ ~ Bernoulli($p$), is heterogeneous between visits at the same site. Royle and Nichols (2003) note that heterogeneity in p can be induced by variation in abundance between visits. Let $N_i$ be the abundance at site $i$. The net probability of detection of at least one individual at site $i% is

$$
p_i = 1 ( 1 - \theta)^{N_i}
$$

Rossman et al. (2016), rewrite this detection probability as 

$$
p_i = 1 ( 1 - \theta_{i, j})^{N_{i, j}}
$$

where $\theta_{i, j}$ is the per individual detection probability at site i, visit j. In implementing this version, the main difference is that their $N_{ij}$ is estimated, and mine is observed.

