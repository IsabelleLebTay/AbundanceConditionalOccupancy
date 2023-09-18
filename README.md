# AbundanceConditionalOccupancy

This integrated abundance-occupancy model addresses the different ecological processes that affect the measures of occupancy and abundance.

The two observations we have are presence/absence and count for each site per recording, or visit.

Occupancy is related to whether the species is present or not given the environmental conditions. The assumptions are that, because of the regional spread of the data, one of the largest factor in the species occupancy/distribution is latitude and longitude. Additionally, we might consider forest type.

Abundance is related to how many individuals occupy or use the area around the sampling unit (the ARU), given that the species is present. The factors presumably influencing this measure of density is the suitability of the habitat is supporting multiple individuals. This is assumed to be related to the size of the retention patch, the age of the harvest, and the proportion of each forest type in the sampling area.

Detection plays a large role in this model at both steps. The probability of detection for a species is of interest at the occupancy level of the model, because our multi-visit design optimises our detection rate. It will take information such as time of day and Julian date. In this model, this stage of detection rate will be the Probability of detection.

Detection of individuals, at the abundance level of the model, is related to two factors: the bird's movement and singing within its territory, and the distance at which we are sampling the bird when it sings and whether its cue is audible given the meteorological and forested conditions. At this stage the variance caused by heterogeneous detection probability between individuals will be addressed by an adjusted distance sampling function.

This model is applied to a multi-visit framework where there are 414 sites and 15 visits per site.

First step: Occupancy and detection probability:

�
�
z 
i
​
 
​
  is the true estimated occupancy state.

�
�
∼
Bernoulli
(
�
�
)
z 
i
​
 ∼Bernoulli(ψ 
i
​
 )
logit
(
�
�
)
=
�
0
+
�
1
×
latitude
�
+
�
2
×
longitude
�
+
forest types
�
logit(ψ 
i
​
 )=β 
0
​
 +β 
1
​
 ×latitude 
i
​
 +β 
2
​
 ×longitude 
i
​
 +forest types 
i
​
 
where 
�
∈
{
0
,
1
}
z∈{0,1}, is either occupied or not and is estimated, and 
logit
(
�
)
logit(ψ) is the prior.

From the observed samples, this true occupancy state is related by:

�
�
�
∼
Bernoulli
(
�
�
×
�
�
�
)
y 
ij
​
 ∼Bernoulli(z 
i
​
 ×p 
ij
​
 )
where 
�
�
�
y 
ij
​
 
​
  is the observed presence/absence for site 
�
i at visit 
�
j, and 
�
�
�
p 
ij
​
 
​
  is the detection probability at site 
�
i, visit 
�
j. 
�
�
�
y 
ij
​
 
​
  is resolved from the observed local abundance, where 
�
�
�
=
1
y 
ij
​
 
​
 =1 if 
�
�
�
>
0
N 
ij
​
 
​
 >0.

logit
(
�
�
�
)
=
�
0
+
�
1
×
time of day
+
�
2
×
Julian date
logit(p 
ij
​
 )=β 
0
​
 +β 
1
​
 ×time of day+β 
2
​
 ×Julian date
The bird count observations 
�
�
�
N 
ij
​
 
​
 , is at the visit-level. For now, it will be summarized as a mean count at the site-level, as 
�
�
M 
i
​
 
​
 
The local observed abundance, an observed quantity 
�
�
N 
i
​
 
​
 , is conditional on occupancy and is used to estimate the unobserved Poisson rate 
�
�
λ 
i
​
 
​
 :

[
�
�
∣
�
�
]
∼
Poisson
(
�
�
×
�
�
)
[M 
i
​
 ∣ψ 
i
​
 ]∼Poisson(λ 
i
​
 ×z 
i
​
 )
log
⁡
(
�
�
)
=
�
0
+
�
size
×
size
+
�
age
×
age
log(λ 
i
​
 )=α 
0
​
 +β 
size
​
 ×size+β 
age
​
 ×age
where the parameter of the Poisson distribution is dependent on the true occupancy state (when the species is present).
