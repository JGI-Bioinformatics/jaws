#!/usr/bin/env python

import os
import sys
import json
import jinja2
import pprint

def get_config_json_data(config_json):
    with open(config_json) as json_file:
        jsondata = json.load(json_file)
    return jsondata

if len(sys.argv) != 3:
    print("Usage: jinja2_translate.py deploy.json template.json")
    sys.exit()

config_file = sys.argv[1]
template_file = sys.argv[2]
config_data = get_config_json_data(config_file)

with open(template_file) as fh:
    tp = jinja2.Template(fh.read())

print(tp.render(config_data))
