#!/usr/bin/env python3

import os
import logging
from functools import lru_cache

_logger = logging.getLogger(__name__)


@lru_cache(maxsize=256)
def available_workflows(workflow=None):
    workflows_folder = os.path.abspath(os.path.dirname(os.path.abspath(__file__))+"/workflows")

    all_workflows = {}
    for root, dirs, files in os.walk(workflows_folder):
        all_workflows.update(
            {filename: os.path.join(root, filename)
             for filename in files
             if os.path.splitext(filename)[1] == '.cwl'
             and (filename not in all_workflows
                  or os.path.getctime(os.path.join(root, filename)) <
                  os.path.getctime(all_workflows[filename]))
             }
        )
    _logger.debug("all_workflows: {0}".format(all_workflows))

    if workflow and workflow not in all_workflows:
        raise Exception("Can't find workflow %s" % workflow)

    return all_workflows[workflow] if workflow else all_workflows
