
# #### Initial library imports:

# In[1]:

from pylab import *
from gfunc.dev import devel as dev
from gfunc import maths as m
from gfunc.scripts.gfunc_build_n_way_one2one import reset_random_edges
from gfunc.xpermutations import xuniqueCombinations

# from IPython.display import display, SVG, Image
# from IPython.display import display_pretty, display_html, display_jpeg, display_png, display_json, display_latex, display_svg


# Out[1]:

#     bestChoose is 'choose' from 'rSeq'.
# 

#     /home/gus/virtualenvs/py275/lib/python2.7/site-packages/pytz/__init__.py:31: UserWarning: Module argparse was already imported from /home/gus/Dropbox/repos/git/khmer/python/argparse.pyc, but /home/gus/virtualenvs/py275/lib/python2.7/site-packages is being added to sys.path
#       from pkg_resources import resource_stream
# 

# In[2]:

a = 1
rcParams['figure.figsize'] = 8*a, 6*a

rcParams['font.size'] = 14


# # PTCI Review:
# 
# **PTCI ~ Phylogentic Transcriptional Conservation Index**
# 
# - a measure of the similarity (putative conservation) of the RNA expression (abundance) profiles between orthologous genes between species.  
# 
# 
# - **The original score combined:**
# 
#     - the Pearson correlation value ($r$) calulated between orthologous genes
# 
#     - the corresponding p-value ($p$) of that correlation
#     
#     - and a scaled weight derived from the evolutionary distance between the orthologs' species ($w(d)$)
#     
#     - the $w(d)$ value is a scaled version of the million years since a last common ancestor between the species
#     
# 

# 
# ## $PTCI = (r_{xprn} + (\frac{r_{tfbs}}{2})) \cdot w(d)$

# ## the $p-value$ of the correlation did not add much to score so it has been removed:
# 
# 
# ## $r_{expn} \cdot w(d)$

# In[3]:

#Image(url='https://raw.github.com/xguse/uci-thesis-latex/my_text/figures/figs/mosqPhyloTree.png')


# Out[3]:

#     <IPython.core.display.Image at 0x4037410>

# # Incorperating similarity between orthologs' 5' flanking regions

# This will be represented by the Pearson correlation between TFBS scores between the orthologous promoters.
# 
# The new score must then combine:
# 
# - ## $r_{xprn}$
# - ## $r_{tfbs}$
# - ## $w(d)$

# # Lets explore the posible scoring landscapes with all posible combinations of xprn and tfbs correlations using addition or multiplication operations
# 
# ## $r_{xprn} \cdot r_{tfbs}$
# ## $r_{xprn} + r_{tfbs}$
# 
# [score evaluations](http://nbviewer.ipython.org/urls/raw.github.com/xguse/ipy_notebooks/master/NewPTCI.ipynb)

# ### Setting the PTCI parameters:

# In[36]:

#set the parameters of ptci to calculate
master_kind = 'rsrd'
master_w_min = 1.0
master_w_max = 1.1

null_reps = 60
ptci_threshold = 1

save_figs = True
save_dir = '/home/gus/Dropbox/repos/git/uci-thesis-latex/figures/figs/ecr_team_ptci_20130918_orthodb7_20131124'


# In[5]:

def save_figure(save_as):
    if save_figs:
        out_path = '%s/%s' % (save_dir.rstrip('/'),save_as)
        print 'saving figure as: %s' % (out_path)
        savefig(out_path)


# # 1: Analyze 1-to-1 ortholog correlations (pairwise only):

# #### Run external script to load expressin/divergence/etc data into the gFunc graph:

# In[6]:

# run gfunc_build_n_way_one2one script and capture output in variable `trap`
get_ipython().magic(u'run -i ../../../Dropbox/repos/git/gfunc/src/gfunc/scripts/gfunc_build_n_way_one2one.py ../../../Dropbox/common/gfunc/Aa_Ag_Cq_sort_pci.new_data.ecr_team_orthodb7.yaml')

graphHandler,graphBuilder,n_way_ortho_table,ortho_parser = trap

# extract edge data from the graphs
edges = graphHandler.edge_dict.values() 


# Out[6]:

#     main() completed.
# 

# #### Function to calculate and store the orthologous expression profile correlation values in the graph edges:

# In[7]:

def get_edge_r_and_p_vals(edges,quiet=True):
    """
    set and get r and p vals from list of edges
    """
    # collect all the results using edge_correlation()
    edge_r_and_p_values = [dev.edge_correlation(edge) for edge in edges]
    
    if not quiet:
        print "r_vals before cleaning: %s" % (len(edge_r_and_p_values))

    # get rid of any results that equal None
    edge_r_and_p_values = [x for x in edge_r_and_p_values if not dev.is_none_or_nan(x)]
    
    if not quiet:
        print "Returning %s r_vals." % (len(edge_r_and_p_values))
        
    return edge_r_and_p_values


# #### Set the expression correlation values and also store them in external variable for easy access after cleaning:

# In[8]:

r_and_p_values = get_edge_r_and_p_vals(edges,quiet=False)


# Out[8]:

#     r_vals before cleaning: 15303
#     Returning 13536 r_vals.
# 

# In[9]:

r_values = [r_and_p_values[i][0] for i in range(len(r_and_p_values))]


# #### Histogram of r-values:

# In[10]:

hist(r_values,bins=50,histtype='stepfilled',cumulative=False, color='b')
xlabel('correlation values')
ylabel('number of edges in each bin')
title('r values for pairwise edge comparisons')


# Out[10]:

#     <matplotlib.text.Text at 0x17c39990>

# image file:

# #### Use the correlation values to calculate and store the REAL DATA z-score statistics for later:

# In[11]:

z_stats = m.get_z_score_stats(r_values)


# In[12]:

print "mean:\t%s\nmedian:\t%s\nstdv:\t%s" % (z_stats[0],z_stats[1],z_stats[2])


# Out[12]:

#     mean:	0.250232161388
#     median:	0.363913598987
#     stdv:	0.579497665578
# 

# #### Function to use z-score stats to calculate and store z-score converted r-values in the gFunc graph:

# In[13]:

def set_z_vals(graphHandler,z_stats,use_center='median'):
    z_stats = {'mean':z_stats[0],'median':z_stats[1],'stdv':z_stats[2]}
    
    center = z_stats[use_center]
    stdv   = z_stats['stdv']
    
    def z_val(r_val,center,stdv):
        return  (r_val - center) / stdv
    
    edges = graphHandler.edge_dict.values() 
    for edge in edges:
        try:
            edge.data.z_val = z_val(edge.data.r_val,center,stdv)
            
        except (TypeError,AttributeError) as exc:
            if 'TypeError' in str(exc):
                edge.data.z_val = None
            elif 'AttributeError' in str(exc):
                dev.edge_correlation(edge)
                if edge.data.r_val == None:
                    edge.data.z_val = None
                else:
                    edge.data.z_val = z_val(edge.data.r_val,center,stdv)
            


# In[14]:

set_z_vals(graphHandler,z_stats,use_center='median')


# #### Function to calculate and store the 1-to-1 pairwise PTCI values in the graph edges

# In[15]:

def get_pairwise_ptci_vals(edges,kind='rpd',quiet=True,w_min=1.0,w_max=1.1):
    if not quiet:
        print "kind: %s" % (kind)
    pairwise_ptci_vals = [dev.get_ptci(edge,kind,w_min,w_max) for edge in edges]
    if not quiet:
        print "ptci_vals before cleaning: %s" % (len(pairwise_ptci_vals))
    # remove any None values
    pairwise_ptci_vals = [ptci for ptci in pairwise_ptci_vals if not dev.is_none_or_nan(ptci)]
    if not quiet:
        print "Returning %s ptci_vals." % (len(pairwise_ptci_vals))
        
    return pairwise_ptci_vals


# #### Set and collect the REAL DATA 1-to-1 PTCIs

# In[16]:

pairwise_ptci_vals = get_pairwise_ptci_vals(edges,kind=master_kind,quiet=True,w_min=master_w_min,w_max=master_w_max)


# #### Function to calculate and store RANDOMIZED 1-to-1 pairwise PTCI values in the graph edges to generate many NULL distributions

# In[17]:

def get_null_pairwise_ptci_distributions(graphHandler,graphBuilder,n_way_ortho_table,ortho_parser,reps=50,kind='rpd',quiet=True,w_min=1.0,w_max=1.1):
    """
    """
    null_paired_ptci_distributions = []

    for rep in range(reps):
        # scramble edges for this rep and set new r&p vals
        reset_random_edges(graphHandler,graphBuilder,n_way_ortho_table,ortho_parser)
        graphHandler.measure_relations()
        
        # do prep
        null_edges = graphHandler.edge_dict.values()
        null_r_and_p_values = get_edge_r_and_p_vals(null_edges,quiet)
        null_r_values = [null_r_and_p_values[i][0] for i in range(len(null_r_and_p_values))]
        null_z_stats = m.get_z_score_stats(null_r_values)
        set_z_vals(graphHandler,null_z_stats,use_center='median')
        
        # calculate null ptci vals
        null_pairwise_ptci_vals = get_pairwise_ptci_vals(null_edges,kind,quiet,w_min,w_max)
        
        # collect null ptci distribution
        null_paired_ptci_distributions.append(null_pairwise_ptci_vals)
        
    
    return null_paired_ptci_distributions


# #### Set and collect the NULL DATA 1-to-1 PTCIs

# In[18]:

null_pairwise_ptci_distributions = get_null_pairwise_ptci_distributions(graphHandler,
                                                                        graphBuilder,
                                                                        n_way_ortho_table,
                                                                        ortho_parser,
                                                                        reps=null_reps,
                                                                        kind=master_kind,
                                                                        quiet=True,
                                                                        w_min=master_w_min,
                                                                        w_max=master_w_max)


# #### Histogram overlays contrasting the REAL and NULL distributions of 1-to-1 orthologous PTCI values

# In[19]:

# Show what the actual data looks like for comparison
real_hist_data_p = hist(pairwise_ptci_vals,bins=50,histtype='stepfilled',cumulative=False, color='c',alpha=.7, label='Real Data')
real_data_bins_p = real_hist_data_p[1]

# Graph null distributions as grey slightly transparent histograms
null_label = 'Null Data'

null_hists_data_p = []
for null_dist in null_pairwise_ptci_distributions:
    nhd = hist(null_dist,bins=real_data_bins_p,histtype='step',cumulative=False, color='k',alpha=.1,label=null_label)
    null_label = None
    null_hists_data_p.append(nhd)

null_bin_members_counts_p = dev.get_null_bin_members_counts_dataframe(null_hists_data_p)
bincenters_p = 0.5*(real_data_bins_p[1:]+real_data_bins_p[:-1])
    

plot(bincenters_p,null_bin_members_counts_p.median(),'-r', label='Null Median (%s reps)' % (len(null_pairwise_ptci_distributions)))
#plot(bincenters_p,null_bin_members_counts_p.mean(),'-b', label='Null Mean')


xlabel(r'PTCI  $w(d)$ = [%s, %s]' % (master_w_min, master_w_max))
ylabel('Pairwise comparisons per bin')
#title('null distributions for pairwise edge correlations (%s reps)' % (len(null_pairwise_ptci_distributions)))
legend(loc=2)

save_figure('pairwise_ptci_hist.pdf')


# Out[19]:

#     saving figure as: /home/gus/Dropbox/repos/git/uci-thesis-latex/figures/figs/ecr_team_ptci_20130918_orthodb7/pairwise_ptci_hist.pdf
# 

# image file:

# In[20]:

# Show what the actual data looks like for comparison
real_hist_data_cum_p = hist(pairwise_ptci_vals,bins=50,histtype='stepfilled',cumulative=-1, color='c',alpha=.7, label='Real Data')
real_data_bins_cum_p = real_hist_data_cum_p[1]

# Graph null distributions as grey slightly transparent histograms
null_label = 'Null Data'

null_hists_data_cum_p = []
for null_dist in null_pairwise_ptci_distributions:
    nhd = hist(null_dist,bins=real_data_bins_cum_p,histtype='step',cumulative=-1, color='k',alpha=.1,label=null_label)
    null_label = None
    null_hists_data_cum_p.append(nhd)

null_bin_members_counts_cum_p = dev.get_null_bin_members_counts_dataframe(null_hists_data_cum_p)
bincenters_cum_p = 0.5*(real_data_bins_cum_p[1:]+real_data_bins_cum_p[:-1])
    

plot(bincenters_cum_p,null_bin_members_counts_cum_p.median(),'-r',
     label='Null Median (%s reps)' % (len(null_pairwise_ptci_distributions)))
#plot(bincenters_cum_p,null_bin_members_counts_cum_p.mean(),'-b', label='Null Mean')


xlabel(r'PTCI  $w(d)$ = [%s, %s]' % (master_w_min, master_w_max))
ylabel('Pairwise comparisons per bin')

legend(loc=0)

save_figure('pairwise_ptci_cum_hist.pdf')


# Out[20]:

#     saving figure as: /home/gus/Dropbox/repos/git/uci-thesis-latex/figures/figs/ecr_team_ptci_20130918_orthodb7/pairwise_ptci_cum_hist.pdf
# 

# image file:

# ## Plot FDR vs PTCI value

# In[21]:

null_divided_cum_p = null_bin_members_counts_cum_p.median()/real_hist_data_cum_p[0]

plot(real_data_bins_cum_p[:-1],null_divided_cum_p,'c')
axvline(ls='--', c='k',x=ptci_threshold)
axhline(ls='--',c='g',y=.25, label='FDR = 0.25')
axhline(ls='--',c='m',y=.15, label='FDR = 0.15')
axhline(ls='--',c='r',y=.05, label='FDR = 0.05')
xlabel(r'PTCI  $w(d)$ = [%s, %s]' % (master_w_min, master_w_max))
ylabel('FDR (%s reps)' % (len(null_pairwise_ptci_distributions)))
title('False discovery rate vs pairwise PTCI thrsholds')
legend()


# Out[21]:

#     <matplotlib.legend.Legend at 0x22338790>

# image file:

# In[22]:

null_divided_cum_p = null_bin_members_counts_cum_p.median()/real_hist_data_cum_p[0]

plot(real_data_bins_cum_p[:-1],null_divided_cum_p,'c')
axvline(ls='--', c='k',x=ptci_threshold)
axhline(ls='--',c='g',y=.25, label='FDR = 0.25')
axhline(ls='--',c='m',y=.15, label='FDR = 0.15')
axhline(ls='--',c='r',y=.05, label='FDR = 0.05')
xlabel(r'PTCI  $w(d)$ = [%s, %s]' % (master_w_min, master_w_max))
ylabel('FDR (%s reps)' % (len(null_pairwise_ptci_distributions)))
title('False discovery rate vs pairwise PTCI thrsholds')
legend()

save_figure('pairwise_ptci_fdr.pdf')


# Out[22]:

#     saving figure as: /home/gus/Dropbox/repos/git/uci-thesis-latex/figures/figs/ecr_team_ptci_20130918_orthodb7/pairwise_ptci_fdr.pdf
# 

# image file:

# # 2: Analyze composite $N$-way ortholog correlations (1-to-1 orthologs in all $N$ species):

# ### We need to regenrate the graph database after the randomization step above

# In[23]:

# run gfunc_build_n_way_one2one script and capture output in variable `trap`
get_ipython().magic(u'run -i ../../../Dropbox/repos/git/gfunc/src/gfunc/scripts/gfunc_build_n_way_one2one.py ../../../Dropbox/common/gfunc/Aa_Ag_Cq_sort_pci.new_data.ecr_team_orthodb7.yaml')

graphHandler,graphBuilder,n_way_ortho_table,ortho_parser = trap

edges = graphHandler.edge_dict.values() 
r_and_p_values = get_edge_r_and_p_vals(edges,quiet=False)
set_z_vals(graphHandler, z_stats, use_center='median')


# Out[23]:

#     main() completed.
#     r_vals before cleaning: 15303
#     Returning 13536 r_vals.
# 

# #### Function to calculate and return the mean PTCI values for $N$-way ortholog subgraphs if and only if all edges successfully generated non-None value PTCI results.

# In[24]:

def calc_mean_ptcis(graphHandler,n_way_ortho_table,kind='rpd',quiet=True,w_min=1.0,w_max=1.1):
    """
    returns list of n-way averaged PTCI values for N-way ortholog subgraphs if and only if 
    all edges successfully generated non-None value PTCI results.
    """
    
    # dictionary of nodes/edges indexed by gene names
    node_dict = graphHandler.node_dict
    edge_dict = graphHandler.edge_dict
    
    graph = graphHandler.graph 
    
    mean_ptcis = []
    
    # calculate all pairwise combinations of indexes
    # so that each ortho-edge of n-way orthos are obtained
    
    index_combos = [ x for x in xuniqueCombinations(range(len(n_way_ortho_table.columns)),2)]
    
    
    for node_list in n_way_ortho_table.itertuples():
        
        node_list = node_list[1:]
        
        ortho_edges = []
        for i in index_combos:
            key = tuple(sorted([node_list[i[0]],node_list[i[1]]]))
            
            try:
                ortho_edges.append(edge_dict[key])
            except KeyError:
                break
                
        ptcis = [dev.get_ptci(edge,kind,w_min,w_max) for edge in ortho_edges]
        
        try:
            mean_ptci = np.mean(ptcis)
            mean_ptcis.append(mean_ptci)
        except TypeError:
            pass

    return mean_ptcis


# #### Calulate the REAL DATA mean PTCI values

# In[25]:

n_way_mean_ptcis = calc_mean_ptcis(graphHandler,n_way_ortho_table,kind=master_kind,quiet=True,w_min=master_w_min,w_max=master_w_max)


# #### Function to calculate RANDOMIZED mean PTCI values to generate many NULL distributions

# In[26]:

def get_null_mean_ptci_distributions(graphHandler,graphBuilder,n_way_ortho_table,ortho_parser,reps=50,kind='rpd',quiet=True,w_min=1.0,w_max=1.1):
    """
    """
    null_mean_ptci_distributions = []

    for rep in range(reps):
        # scramble edges for this rep and set new r&p vals
        reset_random_edges(graphHandler,graphBuilder,n_way_ortho_table,ortho_parser)
        graphHandler.measure_relations()
        
        
        # do prep
        null_edges = graphHandler.edge_dict.values()
        null_r_and_p_values = get_edge_r_and_p_vals(null_edges,quiet)
        null_r_values = [null_r_and_p_values[i][0] for i in range(len(null_r_and_p_values))]
        null_z_stats = m.get_z_score_stats(null_r_values)
        set_z_vals(graphHandler,null_z_stats,use_center='median')
        
        # calculate null ptci vals
        null_mean_ptci_vals = calc_mean_ptcis(graphHandler,n_way_ortho_table,kind,quiet,w_min,w_max)
        
        # collect null ptci distribution
        null_mean_ptci_distributions.append(null_mean_ptci_vals)
        
    
    return null_mean_ptci_distributions


# #### Calculate RANDOMIZED mean PTCI values

# In[27]:

null_mean_ptci_distributions = get_null_mean_ptci_distributions(graphHandler,
                                                                graphBuilder,
                                                                n_way_ortho_table,
                                                                ortho_parser,
                                                                reps=null_reps,
                                                                kind=master_kind,
                                                                quiet=True,
                                                                w_min=master_w_min,
                                                                w_max=master_w_max)


# #### Histogram overlays contrasting the REAL and NULL distributions of $N$-way orthologous PTCI values

# In[28]:

# Show what the actual data looks like for comparison
real_hist_data = hist(n_way_mean_ptcis,bins=50,histtype='stepfilled',cumulative=False, color='c',alpha=.7, label='Real Data')
real_data_bins = real_hist_data[1]

# Graph null distributions as grey slightly transparent histograms
null_label = 'Null Data'
null_hists_data = []
for null_dist in null_mean_ptci_distributions:
    nhd = hist(null_dist,bins=real_data_bins,histtype='step',cumulative=False, color='k',alpha=.1,label=null_label)
    null_label = None
    null_hists_data.append(nhd)

null_bin_members_counts = dev.get_null_bin_members_counts_dataframe(null_hists_data)
bincenters = 0.5*(real_data_bins[1:]+real_data_bins[:-1])
    

plot(bincenters,null_bin_members_counts.median(),'-r',
     label='Null Median (%s reps)' % (len(null_mean_ptci_distributions)))
#plot(bincenters,null_bin_members_counts.mean(),'-b', label='Null Mean')

#axvline(ls='--',x=.65)

xlabel(r'PTCI  $w(d)$ = [%s, %s]' % (master_w_min, master_w_max))
ylabel('Mean comparisons per bin')
#title('null distributions for mean edge correlations (%s reps)' % (len(null_mean_ptci_distributions)))
legend()

save_figure('mean_ptci_hist.pdf')


# Out[28]:

#     saving figure as: /home/gus/Dropbox/repos/git/uci-thesis-latex/figures/figs/ecr_team_ptci_20130918_orthodb7/mean_ptci_hist.pdf
# 

# image file:

# In[29]:

# Show what the actual data looks like for comparison
real_hist_data_cum = hist(n_way_mean_ptcis,bins=50,histtype='stepfilled',cumulative=-1, color='c',alpha=.7, label='Real Data')
real_data_bins_cum = real_hist_data_cum[1]

# Graph null distributions as grey slightly transparent histograms
null_label = 'Null Data'
null_hists_data_cum = []
for null_dist in null_mean_ptci_distributions:
    nhd = hist(null_dist,bins=real_data_bins_cum,histtype='step',cumulative=-1, color='k',alpha=.1,label=null_label)
    null_label = None
    null_hists_data_cum.append(nhd)

null_bin_members_counts_cum = dev.get_null_bin_members_counts_dataframe(null_hists_data_cum)
bincenters_cum = 0.5*(real_data_bins_cum[1:]+real_data_bins_cum[:-1])
    

plot(bincenters_cum,null_bin_members_counts_cum.median(),'-r',
     label='Null Median (%s reps)' % (len(null_mean_ptci_distributions)))
#plot(bincenters_cum,null_bin_members_counts_cum.mean(),'-b', label='Null Mean')

xlabel(r'PTCI  $w(d)$ = [%s, %s]' % (master_w_min, master_w_max))
ylabel('Mean comparisons per bin')
#title('null distributions for mean edge correlations (%s reps)' % (len(null_mean_ptci_distributions)))
legend(loc=0)

save_figure('mean_ptci_cum_hist.pdf')


# Out[29]:

#     saving figure as: /home/gus/Dropbox/repos/git/uci-thesis-latex/figures/figs/ecr_team_ptci_20130918_orthodb7/mean_ptci_cum_hist.pdf
# 

# image file:

# ## Plot FDR vs mean(PTCI)

# In[30]:

null_subtracted_cum = real_hist_data_cum[0] - null_bin_members_counts_cum.median()
null_divided_cum = null_bin_members_counts_cum.median()/real_hist_data_cum[0]

plot(real_data_bins_cum[:-1],null_divided_cum,'c')
axvline(ls='--', c='k',x=ptci_threshold)
axhline(ls='--',c='g',y=.25, label='FDR = 0.25')
axhline(ls='--',c='m',y=.15, label='FDR = 0.15')
axhline(ls='--',c='r',y=.05, label='FDR = 0.05')
xlabel(r'PTCI  $w(d)$ = [%s, %s]' % (master_w_min, master_w_max))
ylabel('FDR')
title('Cum FDR for mean N-way edge correlations (%s reps)' % (len(null_mean_ptci_distributions)))
legend()


# Out[30]:

#     <matplotlib.legend.Legend at 0x13f63e10>

# image file:

# In[35]:

null_subtracted_cum = real_hist_data_cum[0] - null_bin_members_counts_cum.median()
null_divided_cum = null_bin_members_counts_cum.median()/real_hist_data_cum[0]

plot(real_data_bins_cum[:-1],null_divided_cum,'c')
axvline(ls='--', c='k',x=ptci_threshold)
axhline(ls='--',c='g',y=.25, label='FDR = 0.25')
axhline(ls='--',c='m',y=.15, label='FDR = 0.15')
axhline(ls='--',c='r',y=.05, label='FDR = 0.05')
xlabel(r'PTCI  $w(d)$ = [%s, %s]' % (master_w_min, master_w_max))
ylabel('FDR (%s reps)' % (len(null_mean_ptci_distributions)))
title('False discovery rate vs mean PTCI thrsholds')
legend()

save_figure('mean_ptci_fdr.pdf')


# Out[35]:

#     saving figure as: /home/gus/Dropbox/repos/git/uci-thesis-latex/figures/figs/ecr_team_ptci_20130918_orthodb7/mean_ptci_fdr.pdf
# 

# image file:

# In[32]:

#dev.calc_ptci_rsrd??


# In[33]:

#get_ipython().system(u"fembot -t 'ECR team is done' -r3")


# Out[33]:

#     
#     -: (wav)
#     
#       Encoding: Signed PCM    
#       Channels: 1 @ 16-bit   
#     Samplerate: 22050Hz      
#     Replaygain: off         
#       Duration: unknown      
#     
#     play WARN echo: echo: warning >>> gain-out can cause saturation of output <<<
#     



