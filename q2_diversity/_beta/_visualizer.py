# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import collections
import urllib.parse
import pkg_resources
import itertools
import functools

import qiime2
import skbio
import biom
import skbio.diversity
from scipy import spatial, cluster
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import q2templates
from q2_feature_table import rarefy
from ._method import beta, beta_phylogenetic, phylogenetic_metrics
from statsmodels.sandbox.stats.multicomp import multipletests


TEMPLATES = pkg_resources.resource_filename('q2_diversity', '_beta')


def bioenv(output_dir: str, distance_matrix: skbio.DistanceMatrix,
           metadata: qiime2.Metadata) -> None:
    # convert metadata to numeric values where applicable, drop the non-numeric
    # values, and then drop samples that contain NaNs
    df = metadata.to_dataframe()
    df = df.apply(lambda x: pd.to_numeric(x, errors='ignore'))

    # filter categorical columns
    pre_filtered_cols = set(df.columns)
    df = df.select_dtypes([np.number]).dropna()
    filtered_categorical_cols = pre_filtered_cols - set(df.columns)

    # filter 0 variance numerical columns
    pre_filtered_cols = set(df.columns)
    df = df.loc[:, df.var() != 0]
    filtered_zero_variance_cols = pre_filtered_cols - set(df.columns)

    # filter the distance matrix to exclude samples that were dropped from
    # the metadata, and keep track of how many samples survived the filtering
    # so that information can be presented to the user.
    initial_dm_length = distance_matrix.shape[0]
    distance_matrix = distance_matrix.filter(df.index, strict=False)
    filtered_dm_length = distance_matrix.shape[0]

    result = skbio.stats.distance.bioenv(distance_matrix, df)
    result = result.to_html(classes='table table-striped table-hover').replace(
        'border="1"', 'border="0"')

    index = os.path.join(TEMPLATES, 'bioenv_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'initial_dm_length': initial_dm_length,
        'filtered_dm_length': filtered_dm_length,
        'filtered_categorical_cols': ', '.join(filtered_categorical_cols),
        'filtered_zero_variance_cols': ', '.join(filtered_zero_variance_cols),
        'result': result})


_beta_group_significance_fns = {'permanova': skbio.stats.distance.permanova,
                                'anosim': skbio.stats.distance.anosim}


def _get_distance_boxplot_data(distance_matrix, group_id, groupings):
    x_ticklabels = []
    all_group_distances = []

    # extract the within group distances
    within_group_distances = []
    group = groupings[group_id]
    for i, sid1 in enumerate(group):
        for sid2 in group[:i]:
            within_group_distances.append(distance_matrix[sid1, sid2])
    x_ticklabels.append('%s (n=%d)' %
                        (group_id, len(within_group_distances)))
    all_group_distances.append(within_group_distances)

    # extract between group distances for group to each other group
    for other_group_id, other_group in groupings.items():
        between_group_distances = []
        if group_id == other_group_id:
            continue
        for sid1 in group:
            for sid2 in other_group:
                between_group_distances.append(distance_matrix[sid1, sid2])
        x_ticklabels.append('%s (n=%d)' %
                            (other_group_id, len(between_group_distances)))
        all_group_distances.append(between_group_distances)
    return all_group_distances, x_ticklabels


def _get_pairwise_group_significance_stats(
        distance_matrix, group1_id, group2_id, groupings, metadata,
        beta_group_significance_fn, permutations):
    group1_group2_samples = groupings[group1_id] + groupings[group2_id]
    metadata = metadata[group1_group2_samples]
    distance_matrix = distance_matrix.filter(group1_group2_samples)
    return beta_group_significance_fn(distance_matrix, metadata,
                                      permutations=permutations)


def beta_group_significance(output_dir: str,
                            distance_matrix: skbio.DistanceMatrix,
                            metadata: qiime2.MetadataCategory,
                            method: str='permanova',
                            pairwise: bool=False,
                            permutations: int=999) -> None:
    try:
        beta_group_significance_fn = _beta_group_significance_fns[method]
    except KeyError:
        raise ValueError('Unknown group significance method %s. The available '
                         'options are %s.' %
                         (method,
                          ', '.join(_beta_group_significance_fns)))

    # Cast metadata to numeric (if applicable), which gives better sorting
    # in boxplots. Then filter any samples that are not in the distance matrix,
    # and drop samples with have no data for this metadata
    # category, including those with empty strings as values.
    metadata = pd.to_numeric(metadata.to_series(), errors='ignore')
    metadata = metadata.loc[list(distance_matrix.ids)]
    metadata = metadata.replace(r'', np.nan).dropna()

    # filter the distance matrix to exclude samples that were dropped from
    # the metadata, and keep track of how many samples survived the filtering
    # so that information can be presented to the user.
    initial_dm_length = distance_matrix.shape[0]
    distance_matrix = distance_matrix.filter(metadata.index)
    filtered_dm_length = distance_matrix.shape[0]

    # Run the significance test
    result = beta_group_significance_fn(distance_matrix, metadata,
                                        permutations=permutations)

    # Generate distance boxplots
    sns.set_style("white")
    # Identify the groups, then compute the within group distances and the
    # between group distances, and generate one boxplot per group.
    # groups will be an OrderedDict mapping group id to the sample ids in that
    # group. The order is used both on the x-axis, and in the layout of the
    # boxplots in the visualization.
    groupings = collections.OrderedDict(
        [(id, list(series.index))
         for id, series in sorted(metadata.groupby(metadata))])

    for group_id in groupings:
        group_distances, x_ticklabels = \
            _get_distance_boxplot_data(distance_matrix, group_id, groupings)

        ax = sns.boxplot(data=group_distances, flierprops={
            'marker': 'o', 'markeredgecolor': 'black', 'markeredgewidth': 0.5,
            'alpha': 0.5})
        ax.set_xticklabels(x_ticklabels, rotation=90)
        ax.set_xlabel('Group')
        ax.set_ylabel('Distance')
        ax.set_title('Distances to %s' % group_id)
        # change the color of the boxes to white
        for box in ax.artists:
            box.set_facecolor('white')
        sns.despine()
        plt.tight_layout()
        fig = ax.get_figure()
        fig.savefig(os.path.join(output_dir, '%s-boxplots.png' %
                                 urllib.parse.quote_plus(str(group_id))))
        fig.savefig(os.path.join(output_dir, '%s-boxplots.pdf' %
                                 urllib.parse.quote_plus(str(group_id))))
        fig.clear()

    result_html = result.to_frame().to_html(classes=("table table-striped "
                                                     "table-hover"))
    result_html = result_html.replace('border="1"', 'border="0"')

    if pairwise:
        pairwise_results = []
        for group1_id, group2_id in itertools.combinations(groupings, 2):
            pairwise_result = \
                _get_pairwise_group_significance_stats(
                    distance_matrix=distance_matrix,
                    group1_id=group1_id,
                    group2_id=group2_id,
                    groupings=groupings,
                    metadata=metadata,
                    beta_group_significance_fn=beta_group_significance_fn,
                    permutations=permutations)
            pairwise_results.append([group1_id,
                                     group2_id,
                                     pairwise_result['sample size'],
                                     permutations,
                                     pairwise_result['test statistic'],
                                     pairwise_result['p-value']])
        columns = ['Group 1', 'Group 2', 'Sample size', 'Permutations',
                   result['test statistic name'], 'p-value']
        pairwise_results = pd.DataFrame(pairwise_results, columns=columns)
        pairwise_results.set_index(['Group 1', 'Group 2'], inplace=True)
        pairwise_results['q-value'] = multipletests(
            pairwise_results['p-value'], method='fdr_bh')[1]
        pairwise_results.sort_index(inplace=True)
        pairwise_path = os.path.join(
            output_dir, '%s-pairwise.csv' % method)
        pairwise_results.to_csv(pairwise_path)

        pairwise_results_html = pairwise_results.to_html(
            classes=("table table-striped table-hover"))
        pairwise_results_html = pairwise_results_html.replace(
            'border="1"', 'border="0"')
    else:
        pairwise_results_html = None

    index = os.path.join(
        TEMPLATES, 'beta_group_significance_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'initial_dm_length': initial_dm_length,
        'filtered_dm_length': filtered_dm_length,
        'method': method,
        'groupings': groupings,
        'result': result_html,
        'pairwise_results': pairwise_results_html
    })


def beta_rarefaction(output_dir: str, table: biom.Table, sampling_depth: int,
                     metric: str, num_iterations: int,
                     phylogeny: skbio.TreeNode = None,
                     color_scheme: str='BrBG',
                     method: str='spearman') -> None:
        rare_trees = []
        dms = []
        test_statistics = {'spearman': 'rho', 'pearson': 'r'}

        if metric in phylogenetic_metrics():
            if phylogeny is None:
                raise ValueError("A phylogenetic metric (%s) was requested, "
                                 "but a phylogenetic tree was not provided. "
                                 "Phylogeny must be provided when using a "
                                 "phylogenetic diversity metric." % metric)
            beta_metric = functools.partial(beta_phylogenetic,
                                            phylogeny=phylogeny)
        else:
            beta_metric = beta

        for i in range(num_iterations):
            rt = rarefy(table, sampling_depth)
            dm = beta_metric(table=rt, metric=metric)
            dms.append(dm)
            condensed_dm = dm.condensed_form()
            upgma_tree = cluster.hierarchy.average(condensed_dm)
            skbio_tree = skbio.tree.TreeNode.from_linkage_matrix(upgma_tree,
                                                                 dm.ids)
            leaves = _get_leaves(skbio_tree)
            rare_trees.append(leaves)

        dm = beta_metric(table=rt, metric=metric)
        condensed_dm = dm.condensed_form()
        upgma_tree = cluster.hierarchy.average(condensed_dm)
        master_tree = skbio.tree.TreeNode.from_linkage_matrix(upgma_tree,
                                                              dm.ids)
        master_leaves = _get_leaves(master_tree, master=True)

        for leaves in master_leaves.keys():
            nodes_w_leaves = 0
            for tree in rare_trees:
                if leaves in tree:
                    nodes_w_leaves += 1
            master_leaves[leaves] = nodes_w_leaves / num_iterations

        for node in master_tree.non_tips():
            leaves = tuple([tip.name for tip in node.tips()])
            if leaves in master_leaves.keys():
                node.name = str(master_leaves[leaves])

        master_tree.write(os.path.join(output_dir, 'master-tree.nwk'),
                          'newick')

        similarity_mtx = np.ones(shape=(num_iterations, num_iterations))

        for i in range(num_iterations):
            for j in range(i):
                r, p, n = skbio.stats.distance.mantel(
                    dms[i], dms[j], method=method, permutations=0, strict=True)
                similarity_mtx[i, j] = similarity_mtx[j, i] = r

        plt.figure()
        sns.heatmap(similarity_mtx, cmap=color_scheme,
                    vmin=-1.0, vmax=1.0, annot=False,
                    cbar_kws={'ticks': [1, 0.5, 0, -0.5, -1]}).set(
                        xlabel=test_statistics[method],
                        ylabel=test_statistics[method])
        frame = plt.gca()
        frame.axes.get_xaxis().set_ticks([])
        frame.axes.get_yaxis().set_ticks([])
        plt.savefig(os.path.join(output_dir, 'heatmap.svg'))
        similarity_mtx_fp = \
            os.path.join(output_dir, 'rarefaction-iteration-similarities.tsv')
        np.savetxt(similarity_mtx_fp, similarity_mtx, fmt='%1.3f',
                   delimiter='\t')

        heatmap_template = os.path.join(
            TEMPLATES, 'beta_rarefaction_assets', 'heatmap.html')
        master_tree_template = os.path.join(
            TEMPLATES, 'beta_rarefaction_assets', 'master-tree.html')
        index = os.path.join(
            TEMPLATES, 'beta_rarefaction_assets', 'index.html')
        templates = [index, heatmap_template, master_tree_template]

        q2templates.render(templates, output_dir, context={
            'tabs': [{'url': 'heatmap.html',
                      'title': 'Heatmap'},
                     {'url': 'master-tree.html',
                      'title': 'Master Tree'}]})


def _get_leaves(tree, master=False):
    nodes = dict() if master else list()
    for node in tree.non_tips():
        leaves = tuple(sorted([leaf.name for leaf in node.tips()]))

        if master:
            nodes[leaves] = 0
        else:
            nodes.append(leaves)

    return nodes


def _metadata_distance(metadata: pd.Series)-> skbio.DistanceMatrix:
    # This code is derived from @jairideout's scikit-bio cookbook recipe,
    # "Exploring Microbial Community Diversity"
    # https://github.com/biocore/scikit-bio-cookbook
    distances = spatial.distance.pdist(
        metadata.values[:, np.newaxis], metric='euclidean')
    return skbio.DistanceMatrix(distances, ids=metadata.index)


def beta_correlation(output_dir: str,
                     distance_matrix: skbio.DistanceMatrix,
                     metadata: qiime2.MetadataCategory,
                     method: str='spearman',
                     permutations: int=999) -> None:
    test_statistics = {'spearman': 'rho', 'pearson': 'r'}
    alt_hypothesis = 'two-sided'
    try:
        metadata = pd.to_numeric(metadata.to_series(), errors='raise')
    except ValueError as e:
        raise ValueError('Only numeric data can be used with the Mantel test. '
                         'Non-numeric data was encountered in the sample '
                         'metadata. Orignal error message follows:\n%s' %
                         str(e))

    initial_metadata_length = len(metadata)
    metadata = metadata.loc[list(distance_matrix.ids)]
    metadata = metadata.replace(r'', np.nan).dropna()
    filtered_metadata_length = len(metadata)

    ids_with_missing_metadata = set(distance_matrix.ids) - set(metadata.index)
    if len(ids_with_missing_metadata) > 0:
        raise ValueError('All samples in distance matrix must be present '
                         'and contain data in the sample metadata. The '
                         'following samples were present in the distance '
                         'matrix, but were missing from the sample metadata '
                         'or had no data: %s' %
                         ', '.join(ids_with_missing_metadata))

    metadata_distances = _metadata_distance(metadata)
    r, p, n = skbio.stats.distance.mantel(
        distance_matrix, metadata_distances, method=method,
        permutations=permutations, alternative=alt_hypothesis, strict=True)

    result = pd.Series([method.title(), n, permutations, alt_hypothesis,
                        metadata.name, r, p],
                       index=['Method', 'Sample size', 'Permutations',
                              'Alternative hypothesis', 'Metadata category',
                              '%s %s' % (method.title(),
                                         test_statistics[method]),
                              'p-value'],
                       name='Mantel test results')
    result_html = result.to_frame().to_html(classes=("table table-striped "
                                                     "table-hover"))
    result_html = result_html.replace('border="1"', 'border="0"')

    scatter_data = []
    for id1, id2 in itertools.combinations(distance_matrix.ids, 2):
        scatter_data.append((distance_matrix[id1, id2],
                             metadata_distances[id1, id2]))
    x = 'Input distance'
    y = 'Euclidean distance of\n%s' % metadata.name
    plt.figure()
    scatter_data = pd.DataFrame(scatter_data, columns=[x, y])
    sns.regplot(x=x, y=y, data=scatter_data, fit_reg=False)
    plt.savefig(os.path.join(output_dir, 'beta-correlation-scatter.png'))
    plt.savefig(os.path.join(output_dir, 'beta-correlation-scatter.pdf'))

    index = os.path.join(
        TEMPLATES, 'beta_correlation_assets', 'index.html')
    q2templates.render(index, output_dir, context={
        'initial_metadata_length': initial_metadata_length,
        'filtered_metadata_length': filtered_metadata_length,
        'result': result_html
    })
