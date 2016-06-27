# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.plugin import Plugin, Str

import q2_diversity
from q2_types import FeatureTable, Frequency, DistanceMatrix, Phylogeny

plugin = Plugin(
    name='diversity',
    version=q2_diversity.__version__,
    website='https://github.com/qiime2-plugins/q2-diversity',
    package='q2_diversity'
)

# TODO create decorator for promoting functions to methods. This info would
# be moved to the decorator calls.
plugin.methods.register_function(
    function=q2_diversity.beta_diversity,
    # TODO require a uniform sampling effort FeatureTable when predicates exist
    inputs={'feature_table': FeatureTable[Frequency],
            # TODO this is optional; how do we handle that here?
            'phylogeny': Phylogeny},
    parameters={'metric': Str},
    outputs=[('distance_matrix', DistanceMatrix)],
    name='Beta diversity',
    description="Let's compute some pairwise distances!"
)

plugin.methods.register_markdown('markdown/feature_table_to_pcoa.md')
