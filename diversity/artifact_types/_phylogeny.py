# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from qiime.plugin import Type
import skbio


class Phylogeny(Type, variant_of=Type.Artifact):
    def load(self, data_reader):
        fh = data_reader.get_file('phylogeny.nwk')
        return skbio.TreeNode.read(fh)

    def save(self, data, data_writer):
        fh = data_writer.create_file('phylogeny.nwk')
        data.write(fh)
