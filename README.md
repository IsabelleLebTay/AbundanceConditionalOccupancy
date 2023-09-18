# AbundanceConditionalOccupancy

This integrated abundance-occupancy model addresses the different ecological processes that affect the measures of occupancy and abundance.

The two observations we have are presence/absence and count for each site per recording, or visit.

Occupancy is related to whether the species is present or not given the environmental conditions. The assumptions are that, because of the regional spread of the data, one of the largest factor in the species occupancy/distribution is latitude and longitude. Additionally, we might consider forest type.

Abundance is related to how many individuals occupy or use the area around the sampling unit (the ARU), given that the species is present. The factors presumably influencing this measure of density is the suitability of the habitat is supporting multiple individuals. This is assumed to be related to the size of the retention patch, the age of the harvest, and the proportion of each forest type in the sampling area.

Detection plays a large role in this model at both steps. The probability of detection for a species is of interest at the occupancy level of the model, because our multi-visit design optimises our detection rate. It will take information such as time of day and Julian date. In this model, this stage of detection rate will be the Probability of detection.

Detection of individuals, at the abundance level of the model, is related to two factors: the bird's movement and singing within its territory, and the distance at which we are sampling the bird when it sings and whether its cue is audible given the meteorological and forested conditions. At this stage the variance caused by heterogeneous detection probability between individuals will be addressed by an adjusted distance sampling function.

This model is applied to a multi-visit framework where there are 414 sites and 15 visits per site.

First step: Occupancy and detection probability:
{\textit{z}\textsubscript{i}} is the true estimated occupancy state

\begin{equation}
    z_i \sim Bernoulli(\psi_i)
\end{equation}
\begin{equation}
    logit(\psi_i) = \beta_0 + \beta_1 * latitude_i + \beta_2 * longitude_i + forest \: types_i
\end{equation}

where z $\epsilon$ \{0, 1\}, is either occupied or not and is estimated, and \textit{logit($\psi$)} is the prior.

From the observed samples, this true occupancy state is related by:
\begin{equation}
    y_{ij} \sim Bernoulli(z_i * p_{ij})
\end{equation}

where {\textit{y}\textsubscript{ij}} is the \textbf{observed} presence\textbackslash absence for site i at visit j, and {\textit{p}\textsubscript{ij}} is the detection probability at site i, visit j. {\textit{y}\textsubscript{ij}} is resolved from the observed local abundance, where {\textit{y}\textsubscript{ij}} = 1 if {\textit{N}\textsubscript{ij}} > 0.

\begin{equation}
    logit(p_{ij}) =  \beta_0 + \beta_1 * time \; of \; day + \beta_2 * Julian \; date
\end{equation}

The bird count observations \textit{N\textsubscript{ij}}, is at the visit-level. For now, it will be summarised as a mean count at the site-level, as \textit{M\textsubscript{i}}
The local observed abundance, an observed quantity \textit{N\textsubscript{i}}, is conditional on occupancy and is used to estimate the unobserved Poisson rate \textit{lambda\textsubscript{i}}
\begin{equation}
    [M_i | \psi_i]  \sim Poisson(\lambda_i * z_i)
\end{equation}

\begin{equation}
    log(\lambda_i) = \alpha_0 + \beta_{size} * size + \beta_{age} * age 
\end{equation}

where the parameter of the Poisson distribution is dependant on the true occupancy state (when the species is present). 
