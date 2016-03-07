# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="diversity",
    # TODO stop duplicating version string
    version="0.0.0-dev",
    packages=find_packages(),
    install_requires=['feature_table', 'scikit-bio >= 0.4.2, < 0.5.0',
                      'qiime >= 2.0.0'],
    package_data={'diversity': ['workflows/*md']},
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Core diversity analyses.",
    license="BSD",
    url="http://www.qiime.org",
    entry_points={
        'qiime.plugin': ['diversity=diversity.plugin_setup:plugin']
    }
)
