# -*- coding: utf-8 -*-

"""Top-level package for ROTLA."""

__author__ = """Christopher Andrew Lavender, Adam Burkholder"""
__email__ = 'c.andrew.lavender@gmail.com, adam.burkholder@nih.gov'
__version__ = '0.1.0'

import os
import ConfigParser

PATHS = dict()
config = ConfigParser.ConfigParser()
config.read(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '../config', 'paths.cfg'))
for key, value in config.items('paths'):
    PATHS[key] = value
