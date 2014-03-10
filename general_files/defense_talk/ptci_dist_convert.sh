DEFTALK=/home/gus/Dropbox/common/ipy_notebooks/general_files/defense_talk/dist_hist
INSECT=/home/gus/Dropbox/repos/git/uci-thesis-latex/figures/figs/jaspar_insect_ptci_20130918_orthodb7
ECR=/home/gus/Dropbox/repos/git/uci-thesis-latex/figures/figs/ecr_team_ptci_20130918_orthodb7

ins="ins_"
ecr="ecr_"

pdf2svg $INSECT/mean_ptci_hist.pdf $DEFTALK/${ins}mean_ptci_hist.svg 
pdf2svg $INSECT/mean_ptci_cum_hist.pdf $DEFTALK/${ins}mean_ptci_cum_hist.svg
pdf2svg $INSECT/mean_ptci_fdr.pdf $DEFTALK/${ins}mean_ptci_fdr.svg 

pdf2svg $INSECT/pairwise_ptci_hist.pdf $DEFTALK/${ins}pairwise_ptci_hist.svg 
pdf2svg $INSECT/pairwise_ptci_cum_hist.pdf $DEFTALK/${ins}pairwise_ptci_cum_hist.svg
pdf2svg $INSECT/pairwise_ptci_fdr.pdf $DEFTALK/${ins}pairwise_ptci_fdr.svg 

pdf2svg $ECR/mean_ptci_hist.pdf $DEFTALK/${ecr}mean_ptci_hist.svg 
pdf2svg $ECR/mean_ptci_cum_hist.pdf $DEFTALK/${ecr}mean_ptci_cum_hist.svg
pdf2svg $ECR/mean_ptci_fdr.pdf $DEFTALK/${ecr}mean_ptci_fdr.svg 

pdf2svg $ECR/pairwise_ptci_hist.pdf $DEFTALK/${ecr}pairwise_ptci_hist.svg 
pdf2svg $ECR/pairwise_ptci_cum_hist.pdf $DEFTALK/${ecr}pairwise_ptci_cum_hist.svg
pdf2svg $ECR/pairwise_ptci_fdr.pdf $DEFTALK/${ecr}pairwise_ptci_fdr.svg 