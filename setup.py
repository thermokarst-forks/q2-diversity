# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="q2-diversity",
    # TODO stop duplicating version string
    version="0.0.0-dev",
    packages=find_packages(),
    install_requires=['q2-feature-table', 'scikit-bio', 'qiime >= 2.0.0',
                      'q2-types'],
    package_data={'q2_diversity': ['markdown/*md']},
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Core diversity analyses.",
    license="BSD",
    url="http://www.qiime.org",
    entry_points={
        'qiime.plugin': ['q2-diversity=q2_diversity.plugin_setup:plugin']
    }
)
