# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import io
import biom
import skbio
import numpy as np
import pandas as pd
import pandas.util.testing as pdt

from q2_diversity import alpha, alpha_phylogenetic


class AlphaTests(unittest.TestCase):

    def test_alpha(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        actual = alpha(table=t, metric='observed_otus')
        # expected computed by hand
        expected = pd.Series({'S1': 1, 'S2': 2, 'S3': 2}, name='observed_otus')
        pdt.assert_series_equal(actual, expected)

    def test_alpha_phylo_metric(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            alpha(table=t, metric='faith_pd')

    def test_alpha_unknown_metric(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            alpha(table=t, metric='not-a-metric')

    def test_alpha_phylogenetic(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        actual = alpha_phylogenetic(table=t, phylogeny=tree, metric='faith_pd')
        # expected computed with skbio.diversity.alpha_diversity
        expected = pd.Series({'S1': 0.75, 'S2': 1.0, 'S3': 1.0},
                             name='faith_pd')
        pdt.assert_series_equal(actual, expected)

    def test_alpha_phylogenetic_non_phylo_metric(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        with self.assertRaises(ValueError):
            alpha_phylogenetic(table=t, phylogeny=tree, metric='observed_otus')

    def test_alpha_phylogenetic_unknown_metric(self):
        t = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                       ['O1', 'O2'],
                       ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        with self.assertRaises(ValueError):
            alpha_phylogenetic(table=t, phylogeny=tree, metric='not-a-metric')
