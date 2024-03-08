#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:47:41 2023

@author: bjorn
"""

from jinja2 import Environment, FileSystemLoader
from jinja2 import nodes
from jinja2.ext import Extension
from jinja2.utils import markupsafe # Markup

class IncludeRawExtension(Extension):
    tags = {"include_sequence_file"}

    def parse(self, parser):
        lineno = parser.stream.expect("name:include_sequence_file").lineno
        template = parser.parse_expression()
        result = self.call_method("_render", [template], lineno=lineno)
        return nodes.Output([result], lineno=lineno)

    def _render(self, filename):
        return markupsafe.Markup(self.environment.loader.get_source(self.environment, filename)[0])

# Use the extension when setting up Jinja
environment = Environment(loader=FileSystemLoader("."),
                          extensions=[IncludeRawExtension])

template = environment.get_template("jinjatemplate.txt")

print(template.render())
