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


class PCoAResults(Type, variant_of=(Type.Artifact, Type.Metadata)):
    def load(self, data_reader):
        fh = data_reader.get_file('pcoa-results.txt')
        return skbio.OrdinationResults.read(fh)

    def save(self, data, data_writer):
        fh = data_writer.create_file('pcoa-results.txt')
        data.write(fh)

    def get_columns(self, data):
        columns = ['PC%d' % i for i in range(1, len(data.samples.columns) + 1)]
        # TODO implement controlled vocabulary for column types
        return pd.Series([None] * len(columns), index=columns)

    def get_values(self, data, column):
        column_idx = int(column.split('PC')[1]) - 1
        return pd.Series(data.samples.columns.iloc[column_idx],
                         index=data.samples.index)
