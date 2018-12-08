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


import os
from sys import platform
from shutil import copyfile
from sqlparse import split

from airflow.utils.db import merge_conn
from airflow import models
from airflow.settings import DAGS_FOLDER, AIRFLOW_HOME
from airflow.bin.cli import api_client
from airflow import configuration as conf

from warnings import filterwarnings
from MySQLdb import Warning

from biowardrobe_cwl_workflows import available

filterwarnings('ignore', category=Warning)


# system_folder = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "system"))


def generate_biowardrobe_workflow():

    _template = u"""#!/usr/bin/env python3
from airflow import DAG
from biowardrobe_cwl_workflows import workflow
dag = workflow("{}")
"""
    all_workflows = available()
    for workflow in all_workflows:
        if not workflow:
            continue

        _filename = os.path.abspath(os.path.join(
            DAGS_FOLDER,
            os.path.basename(os.path.splitext(workflow)[0]) + '.py')
        )
        print(_filename)
        with open(_filename, 'w') as generated_workflow_stream:
            generated_workflow_stream.write(_template.format(workflow))

    try:
        api_client.get_pool(name='basic_analysis')
    except Exception as e:
        api_client.create_pool(name='basic_analysis',
                               slots=1,
                               description="pool to run basic analysis")

    if not conf.has_option('cwl', 'tmp_folder'):
        if not os.path.exists(conf.AIRFLOW_CONFIG+'.orig'):
            copyfile(conf.AIRFLOW_CONFIG, conf.AIRFLOW_CONFIG+'.orig')
        with open(conf.AIRFLOW_CONFIG, 'w') as fp:
            # for s in ['mesos', 'kerberos', 'celery', 'smtp', 'email', 'dask', 'ldap']:
            #     conf.conf.remove_section(s)

            conf.conf.add_section('cwl')
            conf.set('cwl', 'tmp_folder', os.path.join(AIRFLOW_HOME, 'tmp'))

            conf.set('core', 'logging_level', 'WARNING')
            conf.set('core', 'load_examples', 'False')
            conf.set('webserver', 'dag_default_view', 'graph')
            conf.set('webserver', 'dag_orientation', 'TB')
            conf.set('webserver', 'web_server_worker_timeout', '120')
            conf.set('scheduler', 'job_heartbeat_sec', '20')
            conf.set('scheduler', 'scheduler_heartbeat_sec', '20')
            conf.set('scheduler', 'min_file_process_interval', '30')
            conf.conf.write(fp)

    # startup_scripts = ['com.datirium.airflow-scheduler.plist', 'com.datirium.airflow-webserver.plist']
    # if platform == "darwin":
    #     _sys_dir = os.path.expanduser('~/Library/LaunchAgents')
    #     for scripts in startup_scripts:
    #         with open(os.path.join(system_folder, 'macosx', scripts), 'r') as s:
    #             data = s.read()
    #             # OS X
    #         dst = os.path.join(_sys_dir, scripts)
    #
    #         if os.path.exists(dst):
    #             with open(dst + '.new', 'w') as w:
    #                 w.write(data.format(AIRFLOW_HOME=AIRFLOW_HOME))
    #         else:
    #             with open(dst, 'w') as w:
    #                 w.write(data.format(AIRFLOW_HOME=AIRFLOW_HOME))

    # if platform == "linux" or platform == "linux2":
    # linux
    # elif platform == "win32":
    # Windows...

    # TODO: tmp, dags do not exist ???

# generate_biowardrobe_workflow()