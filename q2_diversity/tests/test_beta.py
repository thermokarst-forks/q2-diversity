# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import io
import os
import tempfile
import glob
import collections

import skbio
import numpy as np
import numpy.testing as npt
from biom.table import Table
import pandas as pd
import qiime

from q2_diversity import (beta, beta_phylogenetic, bioenv,
                          beta_group_significance)
from q2_diversity._beta._visualizer import _get_distance_boxplot_data


class BetaDiversityTests(unittest.TestCase):

    def test_beta(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        actual = beta(table=t, metric='braycurtis')
        # expected computed with scipy.spatial.distance.braycurtis
        expected = skbio.DistanceMatrix([[0.0000000, 0.3333333, 0.6666667],
                                         [0.3333333, 0.0000000, 0.4285714],
                                         [0.6666667, 0.4285714, 0.0000000]],
                                        ids=['S1', 'S2', 'S3'])

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_beta_phylo_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            beta(table=t, metric='unweighted_unifrac')

    def test_beta_unknown_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            beta(table=t, metric='not-a-metric')

    def test_beta_phylogenetic(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        actual = beta_phylogenetic(
            table=t, phylogeny=tree, metric='unweighted_unifrac')
        # expected computed with skbio.diversity.beta_diversity
        expected = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                         [0.25, 0.00, 0.00],
                                         [0.25, 0.00, 0.00]],
                                        ids=['S1', 'S2', 'S3'])

        self.assertEqual(actual.ids, expected.ids)
        for id1 in actual.ids:
            for id2 in actual.ids:
                npt.assert_almost_equal(actual[id1, id2], expected[id1, id2])

    def test_beta_phylogenetic_non_phylo_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        with self.assertRaises(ValueError):
            beta_phylogenetic(table=t, phylogeny=tree, metric='braycurtis')

    def test_beta_phylogenetic_unknown_metric(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        tree = skbio.TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        with self.assertRaises(ValueError):
            beta_phylogenetic(table=t, phylogeny=tree, metric='not-a-metric')


class BioenvTests(unittest.TestCase):

    def test_bioenv(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.Metadata(
            pd.DataFrame([['1.0', 'a'], ['2.0', 'b'], ['3.0', 'c']],
                         index=['sample1', 'sample2', 'sample3'],
                         columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:
            bioenv(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('metadata1' in open(index_fp).read())
            self.assertFalse('metadata2' in open(index_fp).read())
            self.assertFalse('Warning' in open(index_fp).read())

    def test_bioenv_exclude_missing_data(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.Metadata(
            pd.DataFrame([['1.0', '2.0'], ['2.0', ''], ['3.0', '42.0']],
                         index=['sample1', 'sample2', 'sample3'],
                         columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:
            bioenv(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('metadata1' in open(index_fp).read())
            self.assertTrue('metadata2' in open(index_fp).read())
            self.assertTrue('Warning' in open(index_fp).read())
            self.assertTrue('contained 3 samples' in open(index_fp).read())
            self.assertTrue('only 2 sample' in open(index_fp).read())

    def test_bioenv_extra_metadata(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.Metadata(
            pd.DataFrame([['1.0', 'a'], ['2.0', 'b'], ['3.0', 'c'],
                          ['4.0', 'd']],
                         index=['sample1', 'sample2', 'sample3', 'sample4'],
                         columns=['metadata1', 'metadata2']))
        with tempfile.TemporaryDirectory() as output_dir:
            bioenv(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            self.assertTrue('metadata1' in open(index_fp).read())
            self.assertFalse('metadata2' in open(index_fp).read())
            self.assertFalse('Warning' in open(index_fp).read())


class BetaGroupSignificanceTests(unittest.TestCase):

    def test_permanova(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.MetadataCategory(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=['sample1', 'sample2', 'sample3']))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            # all expected boxplots are generated
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.png')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.png')))
            # no extra boxplots are generated
            self.assertEqual(len(glob.glob('%s/*-boxplots.pdf' % output_dir)),
                             2)
            self.assertEqual(len(glob.glob('%s/*-boxplots.png' % output_dir)),
                             2)
            self.assertTrue('PERMANOVA results' in open(index_fp).read())
            self.assertFalse('Warning' in open(index_fp).read())

    def test_anosim(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.MetadataCategory(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=['sample1', 'sample2', 'sample3']))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md, method='anosim',
                                    permutations=42)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))
            # all expected boxplots are generated
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'a-boxplots.png')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.pdf')))
            self.assertTrue(os.path.exists(
                            os.path.join(output_dir, 'b-boxplots.png')))
            # no extra boxplots are generated
            self.assertEqual(len(glob.glob('%s/*-boxplots.pdf' % output_dir)),
                             2)
            self.assertEqual(len(glob.glob('%s/*-boxplots.png' % output_dir)),
                             2)
            self.assertTrue('ANOSIM results' in open(index_fp).read())
            self.assertTrue('<td>42</td>' in open(index_fp).read())
            self.assertFalse('Warning' in open(index_fp).read())

    def test_alt_permutations(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.MetadataCategory(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=['sample1', 'sample2', 'sample3']))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md, permutations=42)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue('<td>42</td>' in open(index_fp).read())

    def test_invalid_method(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.MetadataCategory(
            pd.Series(['a', 'b', 'b'], name='a or b',
                      index=['sample1', 'sample2', 'sample3']))

        with self.assertRaises(ValueError):
            with tempfile.TemporaryDirectory() as output_dir:
                beta_group_significance(output_dir, dm, md, method='bad!')

    def test_filtered_samples_numeric_metadata(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25, 0.66],
                                   [0.25, 0.00, 0.00, 0.66],
                                   [0.25, 0.00, 0.00, 0.66],
                                   [0.66, 0.66, 0.66, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3',
                                       'sample4'])
        md = qiime.MetadataCategory(
            pd.Series(['1.0', '2.0', '2.0', ''], name='a or b',
                      index=['sample1', 'sample2', 'sample3', 'sample4']))
        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue('Warning' in open(index_fp).read())

    def test_filtered_samples_str_metadata(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25, 0.66],
                                   [0.25, 0.00, 0.00, 0.66],
                                   [0.25, 0.00, 0.00, 0.66],
                                   [0.66, 0.66, 0.66, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3',
                                       'sample4'])
        md = qiime.MetadataCategory(
            pd.Series(['a', 'b', 'b', ''], name='a or b',
                      index=['sample1', 'sample2', 'sample3', 'sample4']))
        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue('Warning' in open(index_fp).read())

    def test_extra_metadata(self):
        dm = skbio.DistanceMatrix([[0.00, 0.25, 0.25],
                                   [0.25, 0.00, 0.00],
                                   [0.25, 0.00, 0.00]],
                                  ids=['sample1', 'sample2', 'sample3'])
        md = qiime.MetadataCategory(
            pd.Series(['a', 'b', 'b', 'c'], name='a or b',
                      index=['sample1', 'sample2', 'sample3', 'sample4']))

        with tempfile.TemporaryDirectory() as output_dir:
            beta_group_significance(output_dir, dm, md, permutations=42)
            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue('<td>2</td>' in open(index_fp).read())

    def test_get_distance_boxplot_data_two_groups(self):
        dm = skbio.DistanceMatrix([[0.00, 0.12, 0.13, 0.14, 0.15],
                                   [0.12, 0.00, 0.22, 0.23, 0.24],
                                   [0.13, 0.22, 0.00, 0.31, 0.32],
                                   [0.14, 0.23, 0.31, 0.00, 0.44],
                                   [0.15, 0.24, 0.32, 0.44, 0.00]],
                                  ids=['s1', 's2', 's3', 's4', 's5'])

        groupings = collections.OrderedDict(
            [('g1', ['s1', 's2']), ('g2', ['s3', 's4', 's5'])])
        obs = _get_distance_boxplot_data(dm, 'g1', groupings)
        exp_data = [[0.12], [0.13, 0.14, 0.15, 0.22, 0.23, 0.24]]
        exp_labels = ['g1 (n=1)', 'g2 (n=6)']
        self.assertEqual(obs[0], exp_data)
        self.assertEqual(obs[1], exp_labels)

    def test_get_distance_boxplot_data_within_always_first(self):
        dm = skbio.DistanceMatrix([[0.00, 0.12, 0.13, 0.14, 0.15],
                                   [0.12, 0.00, 0.22, 0.23, 0.24],
                                   [0.13, 0.22, 0.00, 0.31, 0.32],
                                   [0.14, 0.23, 0.31, 0.00, 0.44],
                                   [0.15, 0.24, 0.32, 0.44, 0.00]],
                                  ids=['s1', 's2', 's3', 's4', 's5'])

        groupings = collections.OrderedDict(
            [('g2', ['s3', 's4', 's5']), ('g1', ['s1', 's2'])])
        obs = _get_distance_boxplot_data(dm, 'g1', groupings)
        exp_data = [[0.12], [0.13, 0.14, 0.15, 0.22, 0.23, 0.24]]
        exp_labels = ['g1 (n=1)', 'g2 (n=6)']
        self.assertEqual(obs[0], exp_data)
        self.assertEqual(obs[1], exp_labels)

    def test_get_distance_boxplot_data_three_groups(self):
        dm = skbio.DistanceMatrix([[0.00, 0.12, 0.13, 0.14, 0.15],
                                   [0.12, 0.00, 0.22, 0.23, 0.24],
                                   [0.13, 0.22, 0.00, 0.31, 0.32],
                                   [0.14, 0.23, 0.31, 0.00, 0.44],
                                   [0.15, 0.24, 0.32, 0.44, 0.00]],
                                  ids=['s1', 's2', 's3', 's4', 's5'])

        groupings = collections.OrderedDict(
            [('g1', ['s1', 's2']), ('g2', ['s3', 's5']), ('g3', ['s4'])])
        obs = _get_distance_boxplot_data(dm, 'g1', groupings)
        exp_data = [[0.12], [0.13, 0.15, 0.22, 0.24], [0.14, 0.23]]
        exp_labels = ['g1 (n=1)', 'g2 (n=4)', 'g3 (n=2)']
        self.assertEqual(obs[0], exp_data)
        self.assertEqual(obs[1], exp_labels)

    def test_get_distance_boxplot_data_between_order_retained(self):
        dm = skbio.DistanceMatrix([[0.00, 0.12, 0.13, 0.14, 0.15],
                                   [0.12, 0.00, 0.22, 0.23, 0.24],
                                   [0.13, 0.22, 0.00, 0.31, 0.32],
                                   [0.14, 0.23, 0.31, 0.00, 0.44],
                                   [0.15, 0.24, 0.32, 0.44, 0.00]],
                                  ids=['s1', 's2', 's3', 's4', 's5'])

        groupings = collections.OrderedDict(
            [('g1', ['s1', 's2']), ('g3', ['s4']), ('g2', ['s3', 's5'])])
        obs = _get_distance_boxplot_data(dm, 'g1', groupings)
        exp_data = [[0.12], [0.14, 0.23], [0.13, 0.15, 0.22, 0.24]]
        exp_labels = ['g1 (n=1)', 'g3 (n=2)', 'g2 (n=4)']
        self.assertEqual(obs[0], exp_data)
        self.assertEqual(obs[1], exp_labels)

if __name__ == "__main__":
    unittest.main()
