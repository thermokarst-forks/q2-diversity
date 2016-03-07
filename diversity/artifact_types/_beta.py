# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from qiime.plugin import Type
import pandas as pd
import skbio


class DistanceMatrix(Type, variant_of=(Type.Artifact, Type.Metadata)):
    def load(self, data_reader):
        fh = data_reader.get_file('distance-matrix.tsv')
        return skbio.DistanceMatrix.read(fh)

    def save(self, data, data_writer):
        fh = data_writer.create_file('distance-matrix.tsv')
        data.write(fh)

    def get_columns(self, data):
        columns = data.ids
        # TODO implement controlled vocabulary for column types
        return pd.Series([None] * len(columns), index=columns)

    def get_values(self, data, column):
        return pd.Series(data[column], index=data.ids)
