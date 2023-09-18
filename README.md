# AbundanceConditionalOccupancy

This model is for a single species, and will be replicated for 6 different bird species. They are all migratory boreal songbirds and I expect they will react differnetly to the environmental conditions reflected in the model both due to the species distribution, the regional density, their habitat preferences, and their vocal/territorial behaviours.

The data was gathered using ARUs (Autonomous Recording Units) and processed during the first few hours around dawn. All deteceted individuals of each species of interest is tagged.

This single-species integrated abundance-occupancy model addresses the different ecological processes that affect the measures of occupancy and abundance.

The two observations we have are presence/absence and count for each site per recording, or visit.

Occupancy is related to whether the species is present or not given the environmental conditions. The assumptions are that, because of the regional spread of the data, one of the largest factor in the species occupancy/distribution is latitude and longitude. Additionally, we might consider forest type.

Abundance is related to how many individuals occupy or use the area around the sampling unit (the ARU), given that the species is present. The factors presumably influencing this measure of density is the suitability of the habitat is supporting multiple individuals. This is assumed to be related to the size of the retention patch, the age of the harvest, and the proportion of each forest type in the sampling area.

Detection plays a large role in this model at both steps. The probability of detection for a species is of interest at the occupancy level of the model, because our multi-visit design optimises our detection rate. It will take information such as time of day and Julian date. In this model, this stage of detection rate will be the Probability of detection.

Detection of individuals, at the abundance level of the model, is related to two factors: the bird's movement and singing within its territory, and the distance at which we are sampling the bird when it sings and whether its cue is audible given the meteorological and forested conditions. At this stage the variance caused by heterogeneous detection probability between individuals will be addressed by an adjusted distance sampling function.

This model is applied to a multi-visit framework where there are 414 sites and 15 visits per site.

**First step: Occupancy and detection probability:**

$z_i$ is the true estimated occupancy state.

$$
z_i \sim \text{Bernoulli}(\psi_i)
$$

and the estimated parameter $\psi$ is a linear expression:

$$
\text{logit}(\psi_i) = \beta_0 + \beta_1 \times \text{latitude}_i + \beta_2 \times \text{longitude}_i + \text{forest types}_i
$$

where $z \in \{0, 1\}$, is either occupied or not and is estimated, and $\text{logit}(\psi)$ is the prior, and $\text{forest type}$ is the sum of the the coefficients and observations of the differnet forest types.

From the observed samples, this true occupancy state is related by:

$$
y_{ij} \sim \text{Bernoulli}(z_i \times p_{ij})
$$

where $y_{ij}$ is the **observed** presence/absence for site $i$ at visit $j$, and $p_{ij}$ is the detection probability at site $i$, visit $j$. $y_{ij} = 1$ if $N_{ij} > 0$.

$$
\text{logit}(p_{ij}) = \beta_0 + \beta_1 \times \text{time of day} + \beta_2 \times \text{Julian date}
$$

The bird count observations $N_{ij}$, is at the visit-level. For now, it will be summarized as a mean count at the site-level, as $M_i$. The local observed abundance, an observed quantity $N_i$, is conditional on occupancy and is used to estimate the unobserved Poisson rate $\lambda_i$:

$$
[M_i | \psi_i] \sim \text{Poisson}(\lambda_i \times z_i)
$$

$$
\log(\lambda_i) = \alpha_0 + \beta_{\text{size}} \times \text{size} + \beta_{\text{age}} \times \text{age}
$$

where the parameter of the Poisson distribution is dependent on the true occupancy state (when the species is present).
