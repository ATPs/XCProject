#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

from nose.tools import *
from grep import utils
import re,os

def test_convert_patterns():
    # test that I can convert a list of patterns to expressions
    results = utils.convert_patterns(['.*.py'])
    # assert that they are equal to my expectations
    assert_equal(results,[re.compile(".*.py")])
    
def test_troll_directories():
    # given a directory, return all its contents
    results = utils.troll_directories('.')
    # assert that we have the same contents
    assert_true(os.path.join(".","NOTES") in results)

def test_apply_patterns():
    # get a list of directories
    files = utils.troll_directories('tests')
    patterns = utils.convert_patterns(['test_'])
    print(files,patterns)
    # apply a simple pattern on them
    utils.apply_patterns(files, patterns)
    # assert that we get the right results