#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

tests for markerUtils

"""

# required so that we can find the central model code from here!

import sys
sys.path.append("../") 

import unittest
import numpy as np
from solver.physics.markerUtils import *

class TestApplyMarkerContrib(unittest.TestCase):

    def setUp(self):
        # Create a 3x3 field for testing
        self.field = np.zeros((3, 3))

    def test_default_width(self):
        """Test with default width (-1), all nodes should be updated."""
        applyMarkerContrib(self.field, 1.0, 0.5, 0.5, 0, 0, 1.0)
        expected = np.array([
            [0.25, 0.25, 0.0],
            [0.25, 0.25, 0.0],
            [0.0,  0.0,  0.0]
        ])
        np.testing.assert_array_almost_equal(self.field, expected)
        np.testing.assert_almost_equal(np.sum(self.field), 1.0)

    def test_default_width_non_zero_ij(self):
            """Test with default width, different position, all nodes should be updated."""
            self.field = np.zeros((3, 3))
            applyMarkerContrib(self.field, 1.0, 0.5, 0.5, 1, 1, 1.0)
            expected = np.array([
                [0.0, 0.0, 0.0],
                [0.0, 0.25, 0.25],
                [0.0,  0.25,  0.25]
            ])
            np.testing.assert_array_almost_equal(self.field, expected)
            np.testing.assert_almost_equal(np.sum(self.field), 1.0)

    ############ not working ########################
    # this function doesn't actually handle this case!
    def test_default_width_edge_ij(self):
                """Test with default width, marker is off grid, 
                   only 2 nodes should be updated."""
                self.field = np.zeros((3, 3))
                applyMarkerContrib(self.field, 1.0, 0.5, 0.5, 2, 1, 1.0)
                expected = np.array([
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.25],
                    [0.0,  0.0,  0.25]
                ])
                np.testing.assert_array_almost_equal(self.field, expected)

    def test_custom_width_all_inside(self):
        """Test with custom width, all nodes inside width."""
        self.field = np.zeros((3, 3))
        applyMarkerContrib(self.field, 1.0, 0.2, 0.2, 0, 0, 1.0, width=0.9)
        expected = np.array([
            [0.64, 0.16, 0.0],
            [0.16,  0.04, 0.0],
            [0.0,  0.0, 0.0]
        ])
        np.testing.assert_array_almost_equal(self.field, expected)
        np.testing.assert_almost_equal(np.sum(self.field), 1.0)

    def test_custom_width_1_inside(self):
        """Test with custom width, 1 node in range."""
        self.field = np.zeros((3, 3))
        applyMarkerContrib(self.field, 1.0, 0.6, 0.2, 0, 0, 1.0, width=0.5)
        expected = np.array([
            [0.0, 0.48, 0.0],
            [0.0, 0.0,  0.0],
            [0.0, 0.0,  0.0]
        ])
        np.testing.assert_array_almost_equal(self.field, expected)

    def test_custom_width_all_outside(self):
        """Test with custom width, all nodes outside width."""
        self.field = np.zeros((3, 3))
        applyMarkerContrib(self.field, 1.0, 0.6, 0.6, 0, 0, 1.0, width=0.2)
        expected = np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0,  0.0]
        ])
        np.testing.assert_array_almost_equal(self.field, expected)

    def test_zero_distances(self):
        """Test with zero distances (dxm=0, dym=0)."""
        self.field = np.zeros((3, 3))
        applyMarkerContrib(self.field, 1.0, 0.0, 0.0, 0, 0, 1.0)
        expected = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0]
        ])
        np.testing.assert_array_almost_equal(self.field, expected)
        np.testing.assert_almost_equal(np.sum(self.field), 1.0)

    def test_zero_marker_value(self):
        """Test with zero marker value, field should remain unchanged."""
        self.field = np.ones((3, 3))
        applyMarkerContrib(self.field, 0.0, 0.5, 0.5, 0, 0, 1.0)
        expected = np.ones((3, 3))
        np.testing.assert_array_almost_equal(self.field, expected)

    def test_zero_mwt(self):
        """Test with zero marker weight, field should remain unchanged."""
        self.field = np.ones((3, 3))
        applyMarkerContrib(self.field, 1.0, 0.5, 0.5, 0, 0, 0.0)
        expected = np.ones((3, 3))
        np.testing.assert_array_almost_equal(self.field, expected)


# run the tests!
if __name__ == '__main__':
    unittest.main()