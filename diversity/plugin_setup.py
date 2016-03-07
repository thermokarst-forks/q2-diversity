# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.plugin import Plugin, Str

import diversity
from diversity.artifact_types import DistanceMatrix, Phylogeny
from feature_table.artifact_types import FeatureTable, Frequency

plugin = Plugin(
    name='diversity',
    version=diversity.__version__,
    website='https://github.com/qiime2-plugins/diversity',
    package='diversity'
)

# TODO create decorator for promoting functions to workflows. This info would
# be moved to the decorator calls.
plugin.register_function(
    function=diversity.beta_diversity,
    # TODO require a uniform sampling effort FeatureTable when predicates exist
    inputs={'metric': Str,
            'feature_table': FeatureTable[Frequency],
            # TODO this is optional; how do we handle that here?
            'phylogeny': Phylogeny},
    outputs=[('distance_matrix', DistanceMatrix)],
    name='Beta diversity',
    doc="Let's compute some pairwise distances!"
)

plugin.register_workflow('workflows/feature_table_to_pcoa.md')
