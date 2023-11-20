# -*- coding: utf-8 -*-
import jax
import pymc as pm
import numpy as np
import arviz as az
import pandas as pd
import pytensor.tensor as pt
import matplotlib.pyplot as plt
import itertools
prod = itertools.product
from tqdm import tqdm
import graphviz

np.random.seed(27)

#####plotting parameters
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'figure.titlesize': 16})
plt.rcParams['font.family'] = "DeJavu Serif"
plt.rcParams['font.serif'] = "Cambria Math"

## Load data
df = pd.read_csv("persistence_data_antigenic_modified.csv")
 
df['year'] = [y[:4] for y in df.year_month]
df['month'] = [y[5:] for y in df.year_month]
df = df.sort_values("year")
df.reset_index(inplace=True, drop=True)

years = df.year.unique().shape[0]
months = df.month.unique().shape[0]
regions = df.region.unique().shape[0]
climates = df.climate.unique().shape[0]
epochs = df.epoch.unique().shape[0]

pre = np.repeat("pre", len(df[df.year_month.between('2017-01','2020-01')]))
pan = np.repeat("pan", len(df[df.year_month.between('2020-02','2021-07')]))
pos = np.repeat("pos", len(df[df.year_month.between('2021-08','2023-05')]))

df['epoch'] = list(pre) + list(pan) + list(pos)

pers = df.persistence.values
anti = df.antigen_drift.values
air = 1 - df.inter_air_traffic.values #converse of total arrivals proportion


reg_idx = pd.factorize(df.region)[0]
cli_idx = pd.factorize(df.climate)[0]
epo_idx = pd.factorize(df.epoch)[0]
yea_idx  = pd.factorize(df.year)[0]

y = (pers - pers.mean()) / pers.std()
w = (anti - anti.mean()) / anti.std()
x = (air - air.mean()) / air.std()


coords = {'region':df.region.unique(), 'climate':df.climate.unique(), 
          'epoch':df.epoch.unique(), 'year':df.year.unique(), 
          'loc':df.index.values, 'feature':['inter', 'slope']}


with pm.Model(coords=coords) as mod:
    reg_idx = pm.ConstantData("reg_idx", reg_idx, dims="loc")
    cli_idx = pm.ConstantData("cli_idx", cli_idx, dims="loc")
    epo_idx = pm.ConstantData("epo_idx", epo_idx, dims="loc")
    yea_idx = pm.ConstantData("year_idx", yea_idx, dims="loc")
    
    a_s = pm.HalfNormal("a_s", 1)
    a_z = pm.Normal("a_z", mu=0, sigma=1)
    a = pm.Deterministic("a", a_s*a_z)
    
    b_s = pm.HalfNormal("b_s", 1)
    b_l = pm.Normal("b_l", 0, 1)
    b_z = pm.Normal("b_z", mu=0, sigma=1, dims=("region", "year", "epoch"))
    b = pm.Deterministic("b", b_l + b_s*b_z, dims=("region", "year", "epoch"))
    
    c_s = pm.HalfNormal("c_s", 1)
    c_l = pm.Normal("c_l", 0, 1)
    c_z = pm.Normal("c_z", mu=0, sigma=1, dims=("region", "year", "epoch"))
    c = pm.Deterministic("c", c_l + c_s*c_z, dims=("region", "year", "epoch"))
    
    epsilon = pm.HalfNormal("epsilon", 1)
    gamma = a + c[reg_idx,yea_idx,epo_idx]*x
    w_hat = pm.Normal("w", mu=gamma, sigma=epsilon, observed=w)
    
    mu = pm.Deterministic("mu", a + b[reg_idx,yea_idx,epo_idx]*w + c[reg_idx,yea_idx,epo_idx]*x )
    
    sigma = pm.HalfNormal("sigma", 1)
    
    y_hat = pm.Normal("y", mu=mu, sigma=sigma, observed=y)
    
    
with mod:
    ppc = pm.sample_prior_predictive()

    
#az.plot_ppc(ppc, kind="cumulative", group="prior")    
az.plot_ppc(ppc, group="prior")    
plt.title("Prior Predictive")
plt.tight_layout()
plt.savefig("prior_predictive.png", dpi=300)
plt.show()
plt.close()

# idata = az.from_netcdf("idata_causal.nc")

with mod:
    #idata = pm.sampling.jax.sample_numpyro_nuts(1000, tune=1000, chains=4, target_accept=0.95)
    idata = pm.sample(2000, tune=2000, chains=4, cores=12, target_accept=0.99, 
                      nuts_sampler="numpyro", random_seed=27)
    pred = pm.sample_posterior_predictive(idata)



# az.plot_ppc(pred, kind="cumulative", group="posterior") 
az.plot_ppc(pred, group="posterior")    
plt.title("Posterior Predictive")
plt.tight_layout()
plt.savefig("posterior_predictive.png", dpi=300)
plt.show()
plt.close()

az.plot_trace(idata, kind="rank_vlines")
plt.tight_layout()
plt.title("Estimates")
plt.savefig("rank_plots.png", dpi=300)
plt.show()
plt.close()


mu_m = az.extract(idata)['mu'].values.mean(axis=1)
mu_5, mu_95 = az.hdi(az.extract(idata)['mu'].values.T).T

inter = az.extract(idata)['a'].values
slope1 = az.extract(idata)['b'].values.mean(axis=1)
slope2 = az.extract(idata)['c'].values.mean(axis=1)

df['mu_mean'] = mu_m
df['mu_hdi5'] = mu_5
df['mu_hdi95'] = mu_95

df['pers_z'] = y
df['anti_z'] = w
df['air_z'] = x

for e,r in tqdm(prod(range(epochs),range(regions))):
    k = r 
    print(k)
    epo = df.epoch.unique()[e]
    dat = df[(df.region==df.region.unique()[r]) & (df.epoch==epo)]
    reg = df.region.unique()[r]
    if "/" in reg:
        reg = reg.replace("/", "-")
    n1 = reg
    if epo=="pan":
        epo=""
        epo2=""
    if epo=="pre":
        epo2="pre-"
    if epo=="pos":
        epo2="pos-"
    v1 = np.sort([inter.mean(axis=0) + slope1[r,e]*i for i in dat.anti_z])
    v1m = np.sort(v1.mean(axis=1))
    v1_5, v1_95 = np.sort(az.hdi(v1.T, hdi_prob=0.9).T)
    v2 = np.sort([inter.mean(axis=0) + slope2[r,e]*i for i in dat.air_z])
    v2_5, v2_95 = np.sort(az.hdi(v2.T, hdi_prob=0.9).T)
    v2m = np.sort(v2.mean(axis=1))  
    plt.scatter(dat.pers_z, dat.anti_z, marker="^", color='r', alpha=0.2, label="Antigenic Drift")
    plt.plot(np.sort(dat.pers_z), v1m, color="crimson", label="Antigenic Drift")
    plt.fill_between(np.sort(dat.pers_z), v1_5, v1_95, color="crimson", alpha=0.2)
    plt.scatter(dat.pers_z, dat.air_z, marker="o", color='b', alpha=0.2, label="Air Traffic")
    plt.plot(np.sort(dat.pers_z), v2m, color="slateblue", linestyle="--", label="Air Traffic")
    plt.fill_between(np.sort(dat.pers_z), v2_5, v2_95, color="slateblue", alpha=0.2)
    plt.ylabel("Antigenic Drift / Air Traffic (z-scores)")
    plt.xlabel("Persistence (z-scores)")
    plt.grid(alpha=0.2)
    plt.legend(loc="upper left")
    plt.title(n1+" "+epo2+"pandemic")
    plt.tight_layout()
    plt.savefig("./regression_plots/"+str(k)+"-"+n1+"-"+epo+"-pandemic_causal_reg.png", dpi=100)
    plt.close()
    


#idata.to_netcdf("inference_data.nc")

az.plot_energy(idata)
plt.savefig("energy_plot.png", dpi=300)
plt.close()

summ = az.summary(idata, hdi_prob=0.9)
summ.to_csv("inference_summary.csv")

# summ_b = summ[summ.index.str.contains('b')]
# summ_b = summ_b[~summ_b.index.str.contains('b_')]

# summ_c = summ[summ.index.str.contains('c')]
# summ_c = summ_c[~summ_c.index.str.contains('c_')]

#colorblind palette
colors = ['#377eb8','#ff7f00','#4daf4a','#f781bf','#a65628',
          '#984ea3','#999999','#e41a1c','#dede00','#9EDAE5','#C5B0D5']

pos_b = az.extract(idata)['b'].values.mean(axis=(2))
for i in tqdm(range(regions)):
    reg = df.region.unique()[i]
    if "/" in reg:
        reg = reg.replace("/", "-")
    for j in np.flip(range(years)):
        yea = df.year.unique()[j]
        m = pos_b[i,j].mean()
        h5, h95 = az.hdi(pos_b[i,j], hdi_prob=0.9) 
        y = j+1
        plt.scatter([m],[y], s=100, facecolor="w", edgecolor=colors[j], zorder=3)
        plt.hlines(y=y, xmin=h5,xmax=h95, color=colors[j], linewidth=3)
    plt.yticks(ticks=np.arange(years)+1, labels=df.year.unique())
    plt.gca().invert_yaxis()
    plt.xlabel("z-score")
    plt.title(reg+"\nAntigenic Drift Effects")
    plt.grid(color="gainsboro")
    plt.tight_layout()
    plt.savefig("./forest_plots/"+reg+"_anti_drift_forest.png", dpi=300)
    plt.show()
    plt.close()

pos_c = az.extract(idata)['c'].values.mean(axis=(2))
for i in tqdm(range(regions)):
    reg = df.region.unique()[i]
    if "/" in reg:
        reg = reg.replace("/", "-")
    for j in np.flip(range(years)):
        yea = df.year.unique()[j]
        m = pos_c[i,j].mean()
        h5, h95 = az.hdi(pos_c[i,j], hdi_prob=0.9) 
        y = j+1
        plt.scatter([m],[y], s=100, facecolor="w", edgecolor=colors[j], zorder=3)
        plt.hlines(y=y, xmin=h5,xmax=h95, color=colors[j], linewidth=3)
    plt.yticks(ticks=np.arange(years)+1, labels=df.year.unique())
    plt.gca().invert_yaxis()
    plt.xlabel("z-score")
    plt.title(reg+"\nAir Traffic Effects")
    plt.grid(color="gainsboro")
    plt.tight_layout()
    plt.savefig("./forest_plots/"+reg+"_air_traffic_forest.png", dpi=300)
    plt.show()
    plt.close()


######## save forest_plot data ########

#antigenic drift effects
pos_b = az.extract(idata)['b'].values.mean(axis=(2))


summ_anti = {"var":[], "year":[], "region":[],
             "mean":[], "sd":[], 
             "hdi_5%":[], "hdi_95%":[], 
             "hdi_2.5%":[], "hdi_97.5%":[]}

for i in tqdm(range(regions)):
    reg = df.region.unique()[i]
    if "/" in reg:
        reg = reg.replace("/", "-")
    for j in np.flip(range(years)):
        yea = df.year.unique()[j]
        m = pos_b[i,j].mean()
        s = pos_b[i,j].std()
        h5, h95 = az.hdi(pos_b[i,j], hdi_prob=0.9) 
        h25, h975 = az.hdi(pos_b[i,j], hdi_prob=0.95)
        summ_anti['var'].append("antigenic-drift") 
        summ_anti['year'].append(yea) 
        summ_anti['region'].append(reg)
        summ_anti['mean'].append(m)
        summ_anti['sd'].append(s)
        summ_anti['hdi_5%'].append(h5)
        summ_anti['hdi_95%'].append(h95)
        summ_anti['hdi_2.5%'].append(h25)
        summ_anti['hdi_97.5%'].append(h975)
        
summ_anti = pd.DataFrame(summ_anti)
summ_anti.to_csv("antigenic_effect_sumary.csv", index=False)


#air traffic effects
pos_c = az.extract(idata)['c'].values.mean(axis=(2))

summ_at = {"var":[], "year":[], "region":[],
             "mean":[], "sd":[], 
             "hdi_5%":[], "hdi_95%":[], 
             "hdi_2.5%":[], "hdi_97.5%":[]}

for i in tqdm(range(regions)):
    reg = df.region.unique()[i]
    if "/" in reg:
        reg = reg.replace("/", "-")
    for j in np.flip(range(years)):
        yea = df.year.unique()[j]
        m = pos_c[i,j].mean()
        s = pos_c[i,j].std()
        h5, h95 = az.hdi(pos_c[i,j], hdi_prob=0.9) 
        h25, h975 = az.hdi(pos_c[i,j], hdi_prob=0.95)
        summ_at ['var'].append("air-traffic") 
        summ_at ['year'].append(yea) 
        summ_at['region'].append(reg)
        summ_at['mean'].append(m)
        summ_at['sd'].append(s)
        summ_at['hdi_5%'].append(h5)
        summ_at['hdi_95%'].append(h95)
        summ_at['hdi_2.5%'].append(h25)
        summ_at['hdi_97.5%'].append(h975)
        
summ_at = pd.DataFrame(summ_at)
summ_at.to_csv("air-traffic_effect_sumary.csv", index=False)



## Plot DAG

g = graphviz.Digraph('G', filename='dag_graphviz.gv')

g.edge('Antigenic Drift', 'Persistence')
g.edge('Air Traffic', 'Persistence')
g.edge('Epoch', 'Persistence')
g.edge('Year', 'Persistence')
g.edge('Region', 'Persistence')

g.edge('Air Traffic', 'Antigenic Drift')
g.edge('Epoch', 'Antigenic Drift')
g.edge('Year', 'Antigenic Drift')
g.edge('Region', 'Antigenic Drift')

g.edge('Epoch', 'Air Traffic')
g.edge('Year', 'Air Traffic')
g.edge('Region', 'Air Traffic')

g.edge('Year', 'Epoch')
g.edge('Region', 'Epoch')

g.render("dag_graphviz")


#save as png
g = graphviz.Digraph('G', filename='dag_graphviz.gv', format="png")
g.graph_attr['dpi'] = '300'

g.edge('Antigenic Drift', 'Persistence')
g.edge('Air Traffic', 'Persistence')
g.edge('Epoch', 'Persistence')
g.edge('Year', 'Persistence')
g.edge('Region', 'Persistence')

g.edge('Air Traffic', 'Antigenic Drift')
g.edge('Epoch', 'Antigenic Drift')
g.edge('Year', 'Antigenic Drift')
g.edge('Region', 'Antigenic Drift')

g.edge('Epoch', 'Air Traffic')
g.edge('Year', 'Air Traffic')
g.edge('Region', 'Air Traffic')

g.edge('Year', 'Epoch')
g.edge('Region', 'Epoch')

g.render("dag_graphviz")