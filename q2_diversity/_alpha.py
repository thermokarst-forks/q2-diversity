# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import biom
import scipy
import skbio
import qiime
import pandas as pd
import seaborn as sns
import skbio.diversity
from statsmodels.sandbox.stats.multicomp import multipletests


# We should consider moving these functions to scikit-bio. They're part of
# the private API here for now.
def phylogenetic_metrics():
    return {'faith_pd'}


def non_phylogenetic_metrics():
    return {'ace', 'chao1', 'chao1_ci', 'berger_parker_d', 'brillouin_d',
            'dominance', 'doubles', 'enspie', 'esty_ci', 'fisher_alpha',
            'goods_coverage', 'heip_e', 'kempton_taylor_q', 'margalef',
            'mcintosh_d', 'mcintosh_e', 'menhinick', 'michaelis_menten_fit',
            'observed_otus', 'osd', 'pielou_e', 'robbins', 'shannon',
            'simpson', 'simpson_e', 'singles', 'strong', 'gini_index',
            'lladser_pe', 'lladser_ci'}


def alpha_phylogenetic(table: biom.Table, phylogeny: skbio.TreeNode,
                       metric: str) -> pd.Series:
    if metric not in phylogenetic_metrics():
        raise ValueError("Unknown phylogenetic metric: %s" % metric)

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')
    feature_ids = table.ids(axis='observation')

    result = skbio.diversity.alpha_diversity(
                metric=metric,
                counts=counts,
                ids=sample_ids,
                otu_ids=feature_ids,
                tree=phylogeny)
    result.name = metric
    return result


def alpha(table: biom.Table, metric: str) -> pd.Series:
    if metric not in non_phylogenetic_metrics():
        raise ValueError("Unknown metric: %s" % metric)

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')

    result = skbio.diversity.alpha_diversity(metric=metric, counts=counts,
                                             ids=sample_ids)
    result.name = metric
    return result


def alpha_compare(output_dir: str, alpha_diversity: pd.Series,
                  metadata: qiime.MetadataCategory) -> None:
    metadata = metadata.to_series()
    data = pd.concat([alpha_diversity, metadata], axis=1)
    names = []
    groups = []
    for name, group in data.groupby(metadata.name):
        names.append(name)
        groups.append(list(group[alpha_diversity.name]))
    ax = sns.boxplot(data=groups)
    ax.set_xticklabels(names)
    ax.set_ylabel(alpha_diversity.name)
    ax.get_figure().savefig(os.path.join(output_dir, 'boxplots.png'))
    ax.get_figure().savefig(os.path.join(output_dir, 'boxplots.pdf'))

    # perform Kruskal-Wallis across all groups
    kw_H_all, kw_p_all = scipy.stats.mstats.kruskalwallis(*groups)

    # perform pairwise Kruskal-Wallis across all pairs of groups and
    # correct for multiple comparisons
    kw_H_pairwise = []
    for i in range(len(names)):
        for j in range(i):
            H, p = scipy.stats.mstats.kruskalwallis(groups[i], groups[j])
            kw_H_pairwise.append([names[j], names[i], H, p])
    kw_H_pairwise = pd.DataFrame(
        kw_H_pairwise, columns=['Group 1', 'Group 2', 'H', 'p-value'])
    kw_H_pairwise.set_index(['Group 1', 'Group 2'], inplace=True)
    kw_H_pairwise['q-value'] = multipletests(
                            kw_H_pairwise['p-value'], method='fdr_bh')[1]
    kw_H_pairwise.sort_index(inplace=True)
    kw_H_pairwise.to_csv(os.path.join(output_dir,
                                      'kruskal-wallis-pairwise.csv'))

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        fh.write('<h1>Alpha diversity boxplots</h1>\n')
        fh.write('<img src="./boxplots.png" width="700" height="500"><p>\n')
        fh.write('Image source (<a href="./boxplots.png">png</a> | '
                 '<a href="./boxplots.pdf">pdf</a>)<p>\n')
        fh.write('<h1>Kruskal-Wallis (all groups)</h1>\n')
        fh.write('H: %f<br>\n' % kw_H_all)
        fh.write('p-value: %e<p>\n' % kw_p_all)

        fh.write('<h1>Kruskal-Wallis (pairwise)</h1>\n')
        fh.write('Kruskal-Wallis test applied to all pairs of groups. q-'
                 'values are multiple comparisons adjusted p-values, corrected'
                 ' with Benjamini/Hochberg False Discovery Rate correction.<p>'
                 '\n')
        fh.write(kw_H_pairwise.to_html())
        fh.write('<a href="./kruskal-wallis-pairwise.csv">csv</a>\n')
        fh.write('</body></html>')
