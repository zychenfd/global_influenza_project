<h1> Persistence Regression Analysis </h1>

<p>
	Analysis below examines the relationship between antigenic drift, air traffic (defined as the converse of the total proportion of flight passenger arrivals), and influenza persistence in the context of pre, during and post pandemic periods (epochs) across eleven world regions.  
</p>

<h1> Model </h1>
<p> The directed acyclic graph (DAG) below indicates the pre-established relationships between relevant variables.

<p align="center">
	<img src="dag_graphviz.png" width="550" height="500" />
</p>

Where Antigenic is antigenic drift, and is the main exposure variable affecting persistence. Other variables represent conditioned variables. To estimate these effects, we built the following Bayesian model:  
</p>

<p align="center"> a<sub>s</sub>, b<sub>s</sub>, c<sub>s</sub> ~ HalfNormal(1)  </p>
<p align="center"> a<sub>l</sub>, b<sub>l</sub>, c<sub>l</sub> ~ Normal(0, 1)  </p>
<p align="center"> a<sub>z</sub>, b<sub>z</sub>, c<sub>z</sub> ~ Normal(0, 1)  </p>
<p align="center"> a = a<sub>l</sub> + a<sub>s</sub>a<sub>z</sub>  </p>
<p align="center"> b = b<sub>l</sub> + b<sub>s</sub>b<sub>z</sub> [region,month] </p>
<p align="center"> c = c<sub>l</sub> + c<sub>s</sub>c<sub>z</sub> [region,month] </p>
<p align="center"> &gamma;  = a + cx </p>
<p align="center"> &epsilon; ~ HalfNormal(1) </p>
<p align="center"> &#373; ~ Normal(&gamma;, &epsilon; ) [observed = w]
<p align="center"> &mu;  = a + bw + cx </p>
<p align="center"> &sigma; ~ HalfNormal(1) </p>
<p align="center"> &#375; ~ Gaussian(&mu;, &sigma;) [observed = y] </p>

Where <i>a</i> is an intercept, <i>b</i> and <i>c</i> are varying parameters of region by month size, and w, x and y are the standarised values of air traffic, antigenic-drift and persistence respectively.


<h1> Results </h1>

<p> The model was sampled using HMC NUTS with 2000 tuning steps, 2000 samples and 0.99 target-accept. The model sampled well with ESS > 1000, R-hats ~ 1 and BFMIs >= 0.75. However posterior predictive checks indicate predictive inefficiency. </p>

<p align="center">
	<img src="rank_plots.png" width="700" height="" />
</p>

The model shows a reasonable fit based on posterior predictive checks.

<p align="center">
	<img src="posterior_predictive.png" width="700" height="250" />
</p>

Estimates of total effects seem to indicate that Africa strongly drives the association between air-traffic, antigenic-drift and persistence. See figures below, shadows in figures below indicate 90% highest density intervals (HDIs).

<p align="center">
	<img src="regression_plots/Africa_pre.png" width="600" height="400" />
</p>
<p align="center">
	<img src="regression_plots/Africa_pan.png" width="600" height="400" />
</p>
<p align="center">
	<img src="regression_plots/Africa_post.png" width="600" height="400" />
</p>

The direct effects of air-traffic and antigenic-drift are summarised in the forest_plots directory. Example plots for antigenic-drift in mostly-tropical regions below:

<p align="center">
	<img src="forest_plots/Africa_anti_drift_forest.png" width="330" height="300" />
	<img src="forest_plots/Southern America_anti_drift_forest.png" width="330" height="300" />
</p>

<p align="center">
	<img src="forest_plots/South-eastern Asia_anti_drift_forest.png" width="330" height="300" />
	<img src="forest_plots/Oceania_anti_drift_forest.png" width="330" height="300" />
</p>

<p align="center">
	<img src="forest_plots/Southern Asia_anti_drift_forest.png" width="330" height="300" />

</p>

Note that during pandemic effects are strong for Africa, while other regions remain around zero, maybe with exception of China.

<p align="center">
	<img src="forest_plots/Europe_anti_drift_forest.png" width="330" height="300" />
	<img src="forest_plots/Western Asia_anti_drift_forest.png" width="330" height="300" />
</p>

<p align="center">
	<img src="forest_plots/Northern America_anti_drift_forest.png" width="330" height="300" />
	<img src="forest_plots/China_anti_drift_forest.png" width="330" height="300" />
</p>

Effects of air-traffic are generally positive, with strong effects on Africa for 2020 and 2021 and on South Asia during 2021.
<p align="center">
	<img src="forest_plots/Africa_air_traffic_forest.png" width="330" height="300" />
	<img src="forest_plots/Southern America_air_traffic_forest.png" width="330" height="300" />
</p>

<p align="center">
	<img src="forest_plots/South-eastern Asia_air_traffic_forest.png" width="330" height="300" />
	<img src="forest_plots/Oceania_air_traffic_forest.png" width="330" height="300" />
</p>

<p align="center">
	<img src="forest_plots/Southern Asia_air_traffic_forest.png" width="330" height="300" />
</p>

Effects remain moderate in other regions, maybe with exception of China, which shows slightly higher effects between 2020 and 2022.

<p align="center">
	<img src="forest_plots/Europe_air_traffic_forest.png" width="330" height="300" />
	<img src="forest_plots/Western Asia_air_traffic_forest.png" width="330" height="300" />
</p>

<p align="center">
	<img src="forest_plots/Northern America_air_traffic_forest.png" width="330" height="300" />
	<img src="forest_plots/China_air_traffic_forest.png" width="330" height="300" />
</p>


<h1> Conclusion </h1>

<p> Present results indicate that Africa shows the strongest pandemic-related effects of antigenic-drift, while the effects of inter-regional air traffic remain relatively constant in most regions, with slight decreases in the pandemic period, Africa and South Asia show the opposite pattern.  </p>

<H1> References </H1>
