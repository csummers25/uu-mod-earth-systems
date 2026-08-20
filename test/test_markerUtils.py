#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

tests for markerUtils functions

currently tested:
 - applyMarkerContrib
 - applyGribContrib
 - findNearestNode
 - 

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

###################################################################
class TestApplyGridContrib(unittest.TestCase):

    def test_basic_case(self):
        """ Simple 2x2 grid, all values 1 """
        field = np.ones((2, 2))
        xn, yn = 0, 0
        dxm, dym = 0.5, 0.5
        result = applyGridContrib(field, xn, yn, dxm, dym)
        self.assertAlmostEqual(result, 1.0)

    def test_zero_fractional_distance(self):
        """ dxm and dym are 0, should return top-left node value """
        field = np.array([[1, 2], [3, 4]])
        xn, yn = 0, 0
        dxm, dym = 0.0, 0.0
        result = applyGridContrib(field, xn, yn, dxm, dym)
        self.assertAlmostEqual(result, field[yn, xn])

    def test_one_fractional_distance(self):
        """ dxm and dym are 1, should return bottom-right node value """
        field = np.array([[1, 2], [3, 4]])
        xn, yn = 0, 0
        dxm, dym = 1.0, 1.0
        result = applyGridContrib(field, xn, yn, dxm, dym)
        self.assertAlmostEqual(result, field[yn+1, xn+1])

    def test_edge_case_top_right(self):
        """ Test at top-right corner of a 2x2 grid """
        field = np.array([[1, 2], [3, 4]])
        xn, yn = 0, 0
        dxm, dym = 1.0, 0.0
        result = applyGridContrib(field, xn, yn, dxm, dym)
        self.assertAlmostEqual(result, field[yn, xn+1])

    def test_edge_case_bottom_left(self):
        """ Test at bottom-left corner of a 2x2 grid """
        field = np.array([[1, 2], [3, 4]])
        xn, yn = 0, 0
        dxm, dym = 0.0, 1.0
        result = applyGridContrib(field, xn, yn, dxm, dym)
        self.assertAlmostEqual(result, field[yn+1, xn])

    def test_basic_interpolation(self):
        """ Test middle of nodes """
        field = np.array([[0, 10], [20, 30]])
        xn, yn = 0, 0
        dxm, dym = 0.5, 0.5
        result = applyGridContrib(field, xn, yn, dxm, dym)
        expected = 0.25*0 + 0.25*10 + 0.25*20 + 0.25*30
        self.assertAlmostEqual(result, expected)

    def test_off_centre(self):
        """ Test with off centre position """
        field = np.array([[1, 3], [5, 7]])
        xn, yn = 0, 0
        dxm, dym = 0.25, 0.75
        result = applyGridContrib(field, xn, yn, dxm, dym)
        expected = (0.75*0.25*1) + (0.25*0.25*3) + (0.75*0.75*5) + (0.25*0.75*7)
        self.assertAlmostEqual(result, expected)

    def test_larger_grid(self):
        """ Test with a larger grid, not just 2x2 """
        field = np.arange(12).reshape(3, 4)
        xn, yn = 1, 1
        dxm, dym = 0.5, 0.5
        result = applyGridContrib(field, xn, yn, dxm, dym)
        expected = (0.25*field[1,1]) + (0.25*field[2,1]) + (0.25*field[1,2]) + (0.25*field[2,2])
        self.assertAlmostEqual(result, expected)

    def test_out_of_bounds(self):
        """ Test with indices that would cause out of bounds access """
        field = np.ones((2, 2))
        xn, yn = 1, 1
        dxm, dym = 0.5, 0.5
        with self.assertRaises(IndexError):
            applyGridContrib(field, xn, yn, dxm, dym)

###################################################################
class TestFindNearestNode(unittest.TestCase):

    def test_marker_at_node(self):
        gridx = np.array([0, 1, 2, 3])
        gridy = np.array([0, 1, 2, 3])
        xnum, ynum = len(gridx), len(gridy) 
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 1, 1)
        self.assertEqual((xn, yn), (1, 1))

    def test_marker_between_nodes_x(self):
        gridx = np.array([0, 1, 2, 3])
        gridy = np.array([0, 1, 2, 3])
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 1.5, 1)
        self.assertEqual((xn, yn), (1, 1))

    def test_marker_between_nodes_y(self):
        gridx = np.array([0, 1, 2, 3])
        gridy = np.array([0, 1, 2, 3])
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 1, 1.5)
        self.assertEqual((xn, yn), (1, 1))

    def test_marker_between_nodes_xy(self):
        gridx = np.array([0, 1, 2, 3])
        gridy = np.array([0, 1, 2, 3])
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 1.5, 1.5)
        self.assertEqual((xn, yn), (1, 1))

    def test_marker_left_of_all_nodes(self):
        gridx = np.array([1, 2, 3, 4])
        gridy = np.array([1, 2, 3, 4])
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 0.5, 2)
        self.assertEqual((xn, yn), (0, 1))

    def test_marker_right_of_all_nodes(self):
        gridx = np.array([1, 2, 3, 4])
        gridy = np.array([1, 2, 3, 4])
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 4.5, 2)
        self.assertEqual((xn, yn), (2, 1))

    def test_marker_below_all_nodes(self):
        gridx = np.array([1, 2, 3, 4])
        gridy = np.array([1, 2, 3, 4])
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 2, 0.5)
        self.assertEqual((xn, yn), (1, 0))

    def test_marker_above_all_nodes(self):
        gridx = np.array([1, 2, 3, 4])
        gridy = np.array([1, 2, 3, 4])
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 2, 4.5)
        self.assertEqual((xn, yn), (1, 2))

    def test_non_uniform_grid_x(self):
        gridx = np.array([0, 0.5, 2, 3])
        gridy = np.array([0, 1, 2, 3])
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 1, 1)
        self.assertEqual((xn, yn), (1, 1))

    def test_non_uniform_grid_y(self):
        gridx = np.array([0, 1, 2, 3])
        gridy = np.array([0, 0.5, 2, 3])
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 1, 1)
        self.assertEqual((xn, yn), (1, 1))

    def test_large_grid(self):
        gridx = np.linspace(0, 10, 11)
        gridy = np.linspace(0, 10, 11)
        xnum, ynum = len(gridx), len(gridy)
        xn, yn = findNearestNode(gridx, gridy, xnum, ynum, 5.3, 3.7)
        self.assertEqual((xn, yn), (5, 3))


# run the tests!
if __name__ == '__main__':
    unittest.main()