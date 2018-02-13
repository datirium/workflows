#! /usr/bin/env python3
"""
****************************************************************************

 Copyright (C) 2018 Datirium. LLC.
 All rights reserved.
 Contact: Datirium, LLC (datirium@datirium.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 ****************************************************************************"""


from setuptools import setup, find_packages
from os import symlink
from subprocess import check_output
from time import strftime, gmtime


def get_git_tag():
    return check_output(['git', 'describe', '--contains']).strip()


def get_git_timestamp():
    gitinfo = check_output(
        ['git', 'log', '--first-parent', '--max-count=1',
         '--format=format:%ct', '.']).strip()
    return strftime('%Y%m%d%H%M%S', gmtime(int(gitinfo)))


try:
    symlink('.', './biowardrobe_cwl_workflows', target_is_directory=True)
except:
    pass


def get_version():
    """
    Tries to get package version with following order:
    1. from git_version file - when installing from pip, this is the only source to get version
    2. from tag
    3. from commit timestamp
    :return: package version
    """
    version = '1.0.0'                                       # set default version
    try:
        version = get_git_tag()                             # try to get version info from the closest tag
    except Exception:
        try:
            version = "1.0.{}".format(get_git_timestamp())  # try to get version info from commit date
        except Exception:
            pass
    return version


setup(
    name='biowardrobe-cwl-workflows',
    description="The wrapped BioWardrobe's CWL files for python packaging",
    version=get_version(),
    url='https://github.com/datirium/workflows',
    # download_url='https://github.com/datirium/workflows/archive/v1.0.2.zip',
    author='Datirium, LLC',
    author_email='porter@datirium.com',
    license='Apache-2.0',
    packages=find_packages(
        exclude=[
            "biowardrobe_cwl_workflows.biowardrobe_cwl_workflows",
            "biowardrobe_cwl_workflows.biowardrobe_cwl_workflows.*"]),
    install_requires=[
        'cwltool',
        'jsonmerge',
        'ruamel.yaml < 0.15',
        'apache-airflow >= 1.9.0, < 2'
    ],
    zip_safe=True,
    include_package_data=True,
    package_data={
        'biowardrobe_cwl_workflows': ['workflows/*.cwl',
                                      'expressiontools/*.cwl',
                                      'metadata/*.cwl',
                                      'tools/*.cwl',
                                      'tools/metadata/*.yaml',
                                      'tools/metadata/*.yml']
    }
)
