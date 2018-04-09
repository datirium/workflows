#!/usr/bin/env python3
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
import logging
from functools import lru_cache
from json import dumps

_logger = logging.getLogger(__name__)


@lru_cache(maxsize=256)
def available(workflow=None):
    workflows_folder = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+"/workflows")

    all_workflows = {}
    for root, dirs, files in os.walk(workflows_folder):
        all_workflows.update(
            {filename: os.path.join(root, filename) for filename in files if os.path.splitext(filename)[1] == '.cwl'}
        )
    _logger.debug("all_workflows: {0}".format(dumps(all_workflows, indent=4)))

    if workflow and workflow not in all_workflows:
        raise Exception("Can't find workflow %s" % workflow)

    return all_workflows[workflow] if workflow else all_workflows
