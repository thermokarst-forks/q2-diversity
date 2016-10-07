# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import biom
import skbio
import skbio.diversity


# We should consider moving these functions to scikit-bio. They're part of
# the private API here for now.
def phylogenetic_metrics():
    return {'unweighted_unifrac', 'weighted_unifrac'}


def non_phylogenetic_metrics():
    return {'cityblock', 'euclidean', 'seuclidean', 'sqeuclidean', 'cosine',
            'correlation', 'hamming', 'jaccard', 'chebyshev', 'canberra',
            'braycurtis', 'mahalanobis', 'yule', 'matching', 'dice',
            'kulsinski', 'rogerstanimoto', 'russellrao', 'sokalmichener',
            'sokalsneath', 'wminkowski'}


def beta_phylogenetic(table: biom.Table, phylogeny: skbio.TreeNode,
                      metric: str)-> skbio.DistanceMatrix:
    if metric not in phylogenetic_metrics():
        raise ValueError("Unknown phylogenetic metric: %s" % metric)

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')
    feature_ids = table.ids(axis='observation')

    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids,
        otu_ids=feature_ids,
        tree=phylogeny
    )


def beta(table: biom.Table, metric: str)-> skbio.DistanceMatrix:
    if metric not in non_phylogenetic_metrics():
        raise ValueError("Unknown metric: %s" % metric)

    counts = table.matrix_data.toarray().astype(int).T
    sample_ids = table.ids(axis='sample')

    return skbio.diversity.beta_diversity(
        metric=metric,
        counts=counts,
        ids=sample_ids
    )


_beta_group_significance_fns = {'permanova': skbio.stats.distance.permanova,
                                'anosim': skbio.stats.distance.anosim}
