
# coding: utf-8

# # Purpose:
# 
# __FILL ME IN__

# # Introduction: #
# 
# ## Linkage disequilibrium ##
# 
# Briefly, linkage disequilibrium is the co-occurrence of pairs of alleles within chromosomes in a population at frequencies more elevated or depressed than would be expected if given the overall frequencies of each allele in the population assuming completely random assortment between alleles.  However, alleles are __not__ able to assort in a completely random fashion. In reality, the distance between loci is _inversely proportional_ to the rate of recombination ("mixing") between chromosomal loci. So LD values can not be directly compared between arbitrary groups of allele-pairs.  One must either know the species-specific recombination rate or restrict comparisons to pairs with similar distances and thereby roughly fix the rates of recombination being compared.
# 
# ## _Glossina fuscipes fuscipes_ ##
# 
# 
# We do not know the recombination rate of _G. f. fuscipes_.  Therefore the methods that I have employed below have focused on controlling the effect of variable recombination rates by constructing bins of SNP-pairs that are located on contiguous chromosomal contigs (__not__ considering phasing or diploid/haploid membership) and have similar distances.  For example all SNP-pairs separated by the interval on a supercontig of $[i,j)$ where $|i-j| \sim 100$.

# # Sections of interest
# 
# 
# 
# - [Plot distributions of $r^2$ for various bins](#Plot-distributions-of-$r^2$-for-various-bins)
# 
# 
# - [First thoughts and notes about first thoughts](#First-thoughts-and-notes-about-first-thoughts:)
# 
# 
# - [Generate functions to simulate data from the types of distributions we will be using](#Generate-functions-to-simulate-data-from-the-types-of-distributions-we-will-be-using:)
# 
# 
# - [Build models for learning the distribution parameters](#Build-models-for-learning-the-distribution-parameters)  
#     
#     - [Exponential](#Exponential)
#         - [Results of Exponential excersise](#Results-of-Exponential-excersise:)
#     
#     - [Beta](#Beta)

# # Implementation:

# ## Imports:

# In[1]:

# get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
import seaborn as sns
import ggplot as gp


import numpy as np
import pandas as pd
# import tables as h5

import itertools as it
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy
from scikits import bootstrap as bs
import statsmodels.api as sm
import statsmodels.stats.multitest as smm

import munch

import pymc as mc



# In[2]:

# set figure characteristics

# size
sns.set_context("talk")

# style
sns.set_style("whitegrid")


# ## File paths:

# In[84]:

# define paths to files
ld_table = "/home/gus/remote_mounts/louise/data/genomes/glossina_fuscipes/annotations/SNPs/plink_out/tsetseFINAL_14Oct2014_f2_53.recode.renamed_scaffolds.maf0_05.vcf/ld/r_none_freqs_dprime.ld"
out_pickle="/home/gus/Documents/YalePostDoc/project_stuff/g_f_fucipes_uganda/ddrad58/ld_thresholds/post_MAP_calc.plk"


# ### Define some functions to help us out:

# In[4]:

def get_snps_with_same_distance(df, d, win=100):
    lbound = d - win
    rbound = d + win
    
    if lbound < 0:
        lbound = 0
    
    return df[(df.BP_DELTA >= lbound) & (df.BP_DELTA <= rbound)]


# In[5]:

def get_snps_grouped_by_distances(df, ds, win=100):
    assert isinstance(ds, list)
    assert all([isinstance(d, int) for d in ds])
    
    groups = {}
    
    for d in ds:
        groups[d] = get_snps_with_same_distance(df=df, d=d, win=win)
        
    return groups


# In[6]:

def collect_window_stat(df, upper_lim=1000, win=10,
                        stat_name="window_mean_r2", stat_func=np.mean,
                        return_df=True
                       ):
    
    
    data = {stat_name: [], 
            'd': np.array(range(upper_lim))+1}

    
    for d  in data['d']:
        data[stat_name].append(stat_func(get_snps_with_same_distance(df=df, d=d-1, win=win).R2))

    if return_df:
        return pd.DataFrame(data)
    else:
        return data


# In[7]:

def modified_z_scores(observations):
    
    M_i = (observations - np.median(observations)) / smad(observations)
    
    return M_i

def update_mad_z(df):
    assert isinstance(df, pd.DataFrame)
    
    # add/overwrite mad_z column 
    df['mad_z'] = 'loco'
    
    for d in df.distance_bin.unique():
        d_mask = df.distance_bin == d
        df.loc[d_mask, 'mad_z'] = modified_z_scores(df.R2[d_mask])
        
    # sanity check
    assert 'loco' not in df.mad_z


# In[8]:

def get_snps_in_bin_mask(df, d, win=100):
    
    lbound = d
    rbound = d + win

    return (df.BP_DELTA >= lbound) & (df.BP_DELTA <= rbound)


def update_distance_bin(df, win=100):

    assert isinstance(df, pd.DataFrame)

    # get bin definitions
    longest_d = df.BP_DELTA.max()  # teehee

    ds = xrange(0, longest_d, win)

    # add/overwrite distance_bin column 
    df['distance_bin'] = -1
    df['distance_bin_mean_R2'] = -1
    df['distance_bin_median_R2'] = -1
    
    for d in ds:
        bin_mask = get_snps_in_bin_mask(df=df, d=d, win=win)
        df.loc[bin_mask, 'distance_bin'] = d
        df.loc[bin_mask, 'distance_bin_mean_R2'] = df[bin_mask].R2.mean()
        df.loc[bin_mask, 'distance_bin_median_R2'] = df[bin_mask].R2.median()
    
    # sanity check
    assert df['distance_bin'].min() >= 0
    assert df['distance_bin_mean_R2'].min() >= 0
    assert df['distance_bin_median_R2'].min() >= 0


# ## Load table and do preliminary calculations: 

# In[9]:

ld = pd.read_table(ld_table, sep=" +", engine='python')
ld['R2'] = ld.R**2
ld['BP_DELTA'] = abs(ld.BP_A - ld.BP_B)


# In[10]:

# quick preview of data
ld.head()


# In[11]:

ld.head()


# In[12]:

update_distance_bin(ld, win=50)


# In[13]:

ld.head()


# ## Plot distributions of $r^2$ for various bins

# In[14]:

def plot_bin_dists(df, bin_def="distance_bin <= 500"):
    plt.rcParams['figure.figsize'] = np.array([16, 12]) * 0.65

    p = gp.ggplot(gp.aes(x='R2'), data=df.query(bin_def))
    p = p + gp.geom_histogram(fill='coral') +         gp.facet_wrap("distance_bin") +         gp.theme_seaborn(context='talk') + gp.ggtitle(bin_def)
    
    return p


# # In[15]:

# plot_bin_dists(ld, bin_def="(distance_bin >= 0) and (distance_bin <= 500)")


# # In[16]:

# plot_bin_dists(ld, bin_def="(distance_bin >= 10000) and (distance_bin <= 10500)")


# # In[17]:

# plot_bin_dists(ld, bin_def="(distance_bin >= 100000) and (distance_bin <= 100500)")


# ### A more complete view of the global pattern:

# In[18]:

# ld_thinned = ld.iloc[::50,:]

# sns.jointplot("BP_DELTA",
#               "R2",
#               data = ld_thinned,
#               kind='kde',
#               color="coral",
#               dropna=1);


# # In[19]:

# sns.distplot(ld.R2.dropna(), bins=30, color="coral");
# plt.title("Distribution of $r^2$ irrespective of distance");


# # In[20]:

# sns.distplot(ld.BP_DELTA.dropna(), color="coral");
# plt.title("Distribution of SNP-pair distances");


# # First-thoughts and notes about first-thoughts:
# 
# Where these data approximatly normal-looking one might use a z-score or modified z-score (_median absolute devation_ instead of std-dev) to set a semi-principled threshold for a value being considered "interesting".  But these are __not__ even a little "Normal" looking.
# 
# I then thought about using an Exponential distribution to model them but it can not accommodate the types of data we observe in the shorter distances which have concentrations at the high-end of the $r^2$ spectrum as well.  This may not be a big deal in practice, however?  The fact that the support for the Exponential distribution is $[0,\inf)$ while our data __can only__ occur over $[0,1]$ does bother me a bit though.
# 
# This has led me to the Beta distribution, which in its 2 parameter formulation ($\mathrm{Beta}(\alpha, \beta)$) only supports values over $(0,1)$.  However in a four parameter formulation ($\mathrm{Beta}(\alpha, \beta, a, c)$) accommodates values along $[a,c]$ by scaling linearly $[a,c]$ to within the Beta's normal support interval.
# 
# Below, I will go through what I have done to explore these distributions with our data.

# ----
# # Generate functions to simulate data from the types of distributions we will be using:
# 

# ## Simulate Beta dist/data:

# In[21]:

def sim_beta(a,b):
    return scipy.stats.beta(a,b)


# In[22]:

def plot_beta(a,b):
    plt.rcParams['figure.figsize'] = np.array([16, 12]) * 0.65
    with sns.color_palette("Set2",10):
        if isinstance(a, list):
            a = np.array(a)
        if isinstance(a, int) or isinstance(a, float):
            a = np.array([a])
            
        if isinstance(b, list):
            b = np.array(b)
        if isinstance(b, int) or isinstance(b, float):
            b = np.array([b])
        
        beta = scipy.stats.beta
        
        for i, v in enumerate(a):
            x = np.linspace(beta.ppf(0.01, a=a[i], b=b[i]),
                          beta.ppf(0.99, a=a[i], b=b[i]), 10000)
            plt.plot(x, beta.pdf(x, a=a[i], b=b[i]),
                     lw=2, label=r'$\alpha$ = {a}, $\beta$ = {b}'.format(a=a[i],b=b[i]))
            
        plt.title("Beta PDF(s)".format(a=a,b=b))
        plt.legend(loc='best', frameon=False)
        
        
def plot_beta_cdf(a,b):
    plt.rcParams['figure.figsize'] = np.array([16, 12]) * 0.65
    with sns.color_palette("Set2",10):
        if isinstance(a, list):
            a = np.array(a)
        if isinstance(a, int) or isinstance(a, float):
            a = np.array([a])
            
        if isinstance(b, list):
            b = np.array(b)
        if isinstance(b, int) or isinstance(b, float):
            b = np.array([b])
        
        beta = scipy.stats.beta
        
        for i, v in enumerate(a):
            x = np.linspace(beta.ppf(0.01, a=a[i], b=b[i]),
                          beta.ppf(0.99, a=a[i], b=b[i]), 10000)
            plt.plot(x, beta.cdf(x, a=a[i], b=b[i]),
                     lw=2, label=r'$\alpha$ = {a}, $\beta$ = {b}'.format(a=a[i],b=b[i]))
            
        plt.title("Beta CDF(s)".format(a=a,b=b))
        plt.legend(loc='best', frameon=False)


# In[23]:

# plot_beta([0.5,5,1,2,2],[0.5,1,3,2,5])
# plt.ylim(0,3);


# 
# # Build models for learning the distribution parameters

# ### Scale $r^2$ data to avoid `0`s and `1`s

# In[24]:

ld["R2_scaled_for_B"] = ld.R2.apply(lambda x: ((x-0.5)*0.999) + 0.5)


# In[25]:

ld.head()


# ## Beta model of distance bin's $r^2$ for MCMC learning:
# 
# For `pyMC` we must build model components separately and add them to the model proper. 
# 

# In[26]:

test_bin = ld.query("(distance_bin == 100)").copy() # get a dataframe of the test_bin after the scaling


# In[27]:

alpha_of_beta  = mc.Uniform('alpha_of_beta',0.01,10)
beta_of_beta = mc.Uniform('beta_of_beta',0.01,10)


# In[28]:

r2_distribution_beta = mc.Beta(name='r2_distribution_beta',
                               alpha=alpha_of_beta,
                               beta=beta_of_beta,
                               value=test_bin.R2_scaled_for_B,
                               observed=True,
                               verbose=0
                              )


# In[29]:

alpha_of_beta_original, beta_of_beta_original  = alpha_of_beta.value, beta_of_beta.value
alpha_of_beta_original, beta_of_beta_original


# In[30]:

bin_beta_model = mc.Model([r2_distribution_beta, alpha_of_beta, beta_of_beta])


# In[31]:

bin_beta_MAP = mc.MAP(bin_beta_model)


# In[32]:

bin_beta_MAP.fit()


# In[33]:

alpha_of_beta_after_map, beta_of_beta_after_map  = alpha_of_beta.value, beta_of_beta.value
alpha_of_beta_after_map, beta_of_beta_after_map


# In[34]:

a_list = []
b_list = []

alpha_of_beta.value, beta_of_beta.value = alpha_of_beta.random(),beta_of_beta.random()
alpha_of_beta_random, beta_of_beta_random = alpha_of_beta.value, beta_of_beta.value
a_list.append(alpha_of_beta_random)
b_list.append(beta_of_beta_random)
print alpha_of_beta_random, beta_of_beta_random
bin_beta_MCMC = mc.MCMC(bin_beta_model)


# In[35]:

bin_beta_MCMC.sample(iter=40000, burn=20000, thin=1)

alpha_of_beta_after_mcmc, beta_of_beta_after_mcmc  = alpha_of_beta.value, beta_of_beta.value

a_list.append(alpha_of_beta_after_mcmc)
b_list.append(beta_of_beta_after_mcmc)

alpha_of_beta_after_mcmc, beta_of_beta_after_mcmc


# In[36]:

# mc.Matplot.plot(bin_beta_MCMC)


# In[37]:

# plot_beta(a_list, b_list)
# plt.ylim(0,50);
# plt.text(.5,42, "*The green function is the random starting params.");


# # In[38]:

# plot_beta(alpha_of_beta_after_mcmc.item(), beta_of_beta_after_mcmc.item())
# plt.text(0.21, 275,'Learned params alone on an auto-scaled canvas rather than scaled to allow view of both curves.')


# # In[39]:

# plot_beta_cdf(alpha_of_beta_after_mcmc.item(), beta_of_beta_after_mcmc.item())
# plt.title("Cumulative distribution function of learned Beta params");


# ## Record the probabilities of obtaining each R2 value (or smaller) in the test bin:
# 
# __[QUESTION:]__ I am not sure whether multiple testing comes in to play here in the filtering step or not so I calculate both uncorrected and corrected values.  Which should I use for the filter?

# In[40]:

# Record the probabilities of obtaining each R2 value (or smaller) in the bin after 
test_bin['cdf'] = scipy.stats.beta.cdf(test_bin.R2_scaled_for_B, alpha_of_beta_after_mcmc, beta_of_beta_after_mcmc)
test_bin['one_minus_cdf'] = 1 - test_bin['cdf']
test_bin['one_minus_cdf_BH'] = smm.multipletests(test_bin.one_minus_cdf, method='fdr_bh')[1]


# ### Filter for $r^2$ values with corrected p-values $ \le 0.05$

# In[41]:

len(test_bin.query("one_minus_cdf_BH <= 0.05"))


# In[42]:

# bin_beta_model_dag = mc.graph.dag(model=bin_beta_model, format='jpg', prog='dot',
#                                          path=None, name="bin_beta_model", consts=True, legend=True, 
#                                          collapse_deterministics=False, collapse_potentials=False, label_edges=True)


# __Network graph of Beta Model__
# 
# 
# ![](bin_beta_model.jpg)
# 

# # Begin function development

# In[43]:

# function to set up and yield models

def yield_models(dataframe):
    
    verbose=True

    bin_ids = dataframe.distance_bin.unique()
    for bin_id in bin_ids:
        
        assert isinstance(bin_id, int)
        
        # get our distance binned r^2 data in a nice dataframe
        data = dataframe.query("(distance_bin == {bin_id})".format(bin_id=bin_id))

        # generate names for stocastics
        alpha_name = "{bin_id}_alpha".format(bin_id=bin_id)
        beta_name = "{bin_id}_beta".format(bin_id=bin_id)
        r2_dist_name = "{bin_id}_r2_distribution_beta".format(bin_id=bin_id)

        # set priors for parameters
        alpha_of_beta  = mc.Uniform(alpha_name,0.01,10)
        beta_of_beta = mc.Uniform(beta_name,0.01,10)

        # set the data
        r2_distribution_beta = mc.Beta(name=r2_dist_name,
                                       alpha=alpha_of_beta,
                                       beta=beta_of_beta,
                                       value=data.R2_scaled_for_B.dropna(),
                                       observed=True,
                                       verbose=0
                                      )

        # create and yield the model object tagged with its bin_id
        model = mc.Model([r2_distribution_beta, alpha_of_beta, beta_of_beta])
        model.bin_id_tag = bin_id
        model.MCMC_run = None  # allow us to know how many times we had to do the experiments rather than MAP
        if verbose:
            print bin_id
        yield model


# In[44]:

def record_parameters_and_probabilities(df, model):
    
    # set up empty slots for our upcoming values
    # but only we havent already
    try:
        df.alpha_param[0]
    except AttributeError:
        df['alpha_param'] = np.nan
        df['beta_param'] = np.nan
        df['cdf'] = np.nan
        df['one_minus_cdf'] = np.nan
        df['one_minus_cdf_BH'] = np.nan
        df['MAP_succeeded'] = False
    
    # set up a mask for this distance_bin
    bin_mask = df.distance_bin == model.bin_id_tag
    
    # perform Max A Posteriori 
    model_runner = mc.MAP(model)
    
    
    # Parameter names for access through the model
    a_name = "{bin_id}_alpha".format(bin_id=model.bin_id_tag)
    b_name = "{bin_id}_beta".format(bin_id=model.bin_id_tag)
    
    try:
        # if we can take the short cut, great!
        model_runner.fit()
        model.MCMC_run = False
    except RuntimeError as exc:
        if "Posterior probability optimization converged to value with zero probability." in exc.message:
            # Otherwise lets run the MCMC experiment 
            model_runner = mc.MCMC(model)
            model.MCMC_run = model_runner
            
            
            # re-initialize the alpha and beta starting values randomly
            # to avoid non-allowed values that seem to slip in from the MAP.fit()
            current_alpha,current_beta = model.get_node(a_name),model.get_node(b_name)
            current_alpha.random()
            current_beta.random()
            
            # do the learning
            model_runner.sample(iter=40000, burn=20000, thin=1)
            
            
    
    # get and store parameters
    alpha, beta = model.get_node(a_name),model.get_node(b_name)
    
    df[bin_mask].alpha_param = alpha
    df[bin_mask].beta_param = beta
    
    # get and store probabilities
    # Record the probabilities of obtaining each R2 value (or smaller) in the bin afer scaling
    df.loc[bin_mask, 'cdf'] = scipy.stats.beta.cdf(df[bin_mask].R2_scaled_for_B,
                                               alpha.value.item(),
                                               beta.value.item())
    df.loc[bin_mask, 'one_minus_cdf'] = 1 - df[bin_mask]['cdf']
    df.loc[bin_mask, 'one_minus_cdf_BH'] = smm.multipletests(df[bin_mask].one_minus_cdf,
                                                         method='fdr_bh')[1]


# # Do the calculations

# In[45]:

# pdb


# In[46]:

models = {}
for model in yield_models(ld):
    models[model.bin_id_tag] = model  # Save models in case we need to plot MCMC convergence plots
    record_parameters_and_probabilities(ld, model)
ld.head()


# In[ ]:

ld.to_pickle(out_pickle)


# # In[47]:

# # a_23800_alpha = models[23800].get_node("23800_alpha").value


# # # In[48]:

# # models[23800].get_node("23800_beta").value


# # In[49]:

# # plot_beta(4.109947507870923,2.3511229452031075)


# # # Examine results that needed MCMC

# # In[58]:

# mcmc_runs = [m for m in models.values() if m.MCMC_run]


# # In[61]:

# mc.Matplot.plot(mcmc_runs[2].MCMC_run)


# # In[62]:

# plot_bin_dists(ld, bin_def='distance_bin == 9150')


# # In[74]:

# # mod_x = models[9150]

# # sns.distplot(mod_x.MCMC_run.get_node('9150_r2_distribution_beta').value)


# # In[69]:

# mod_x_b = mod_x.MCMC_run.get_node('9150_r2_distribution_beta')


# # In[73]:

# mod_x.MCMC_run.get_node('9150_r2_distribution_beta').value


# # In[76]:

# r9150 = ld.query('distance_bin == 9150')


# # In[83]:

# r9150.R2_scaled_for_B.dropna()


# In[ ]:



