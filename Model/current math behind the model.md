**Data:**

- Total number of sites: \( I \)
- Number of visits per site: \( J \)
- Bird count for each site \( i \) and visit \( j \): \( N_{i,j} \)
- Various covariates including:
  - Latitude of site \( i \): \( \text{latitude}_i \)
  - Longitude of site \( i \): \( \text{longitude}_i \)
  - ... and others.

**Transformed Data:**

- Observed presence/absence derived from the count \( N \): \( y_{i,j} \)
- Mean bird count for each site \( i \): \( M_i \)

**Parameters:**

Various coefficients for occupancy, detection probability, and abundance.

**Model:**

1. **Occupancy Model:**

The estimated log odds of the latent unobserved occupancy for site \( i \) is modeled as:

$$
\text{logit\_psi}_i \sim \beta_{0\_psi} + \beta_{1\_psi} \times \text{latitude}_i + \beta_{2\_psi} \times \text{longitude}_i + \text{other covariates}
$$

2. **Detection Model:**

The log odds for each site \( i \) and visit \( j \) is modeled as:

$$
\text{logit\_p}_{i,j} \sim \beta_{0\_p} + \beta_{1\_p} \times \text{time\_of\_day}_{i,j} + \beta_{2\_p} \times \text{Julian\_date}_{i,j}
$$

3. **Abundance Model:**

The mean bird count for each site \( i \) given occupancy \( \psi_i \) is modeled as:

$$
M_i | \psi_i \sim \text{Poisson}(\lambda_i \times \psi_i)
$$

And the log lambda for site \( i \) is:

$$
\log(\lambda_i) = \beta_{0\_lambda} + \beta_{\text{size}} \times \text{size}_i + \beta_{\text{age}} \times \text{age}_i + \text{other covariates}
$$
