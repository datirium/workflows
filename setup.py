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
import os


try:
    os.symlink('.', './biowardrobe_cwl_workflows', target_is_directory=True)
except:
    pass


setup(
    name='biowardrobe-cwl-workflows',
    description='Python package to extend BioWardrobe functionality with CWL Airflow',
    version='1.0.0',
    url='https://github.com/datirium/biowardrobe-airflow-analysis',
    download_url='https://github.com/datirium/biowardrobe-airflow-analysis',
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
    zip_safe=False,
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
