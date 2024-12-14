import unittest
import numpy as np
from Toolbox.Airfoil import Airfoil
from Toolbox.Member import Member
from Data.Panels import Laminates, Sandwiches

class TestAirfoil(unittest.TestCase):
    def setUp(self):
        """
        Create a basic Airfoil object for testing.
        This uses mock data for simplicity.
        """
        sparlocations = [80, 250]  # Example spar positions
        chordlength = 300  # Example chord length
        thickness = 15  # Example thickness
        sparmembers = [Member(Laminates['CFQIezcomposites_spreadtow']) for _ in sparlocations]
        topmembers = [Member(Laminates['CFQIezcomposites_spreadtow']) for _ in range(len(sparlocations) + 1)]
        botmembers = [Member(Laminates['CFQIezcomposites_spreadtow']) for _ in range(len(sparlocations) + 1)]

        reinforcement_top = Laminates['Reinforcement']
        reinforcement_bottom = Laminates['Reinforcement']
        trstart, trend = 60, 100
        brstart, brend = 60, 100

        self.airfoil = Airfoil(
            airfoilname='NACA2410',
            thickness=thickness,
            chordlength=chordlength,
            sparlocations=sparlocations,
            topmembers=topmembers,
            botmembers=botmembers,
            sparmembers=sparmembers,
            trpanel=reinforcement_top,
            trstart=trstart,
            trend=trend,
            brpanel=reinforcement_bottom,
            brstart=brstart,
            brend=brend
        )

    def test_initialization(self):
        """Test initialization of the Airfoil object."""
        self.assertEqual(self.airfoil.airfoilname, 'NACA2410')
        self.assertEqual(self.airfoil.thickness, 15)
        self.assertEqual(len(self.airfoil.sparlocations), 2)

        # Check length of the top- and botmember lists, they should be one longer than sparlocations list
        self.assertEqual(len(self.airfoil.topmembers), 3)  # len(sparlocations) + 1
        self.assertEqual(len(self.airfoil.botmembers), 3)  # len(sparlocations) + 1

    def test_top_botmemberlist_order(self):
        """Test the order of the top and bot member lists from leading to trailing edge."""

        # Extract the x-coordinates of the start points for topmembers
        top_start_x_coords = [member.startcoord[0] for member in self.airfoil.topmembers]
        # Ensure the x-coordinates are in ascending order
        print('top_start_x_coords',top_start_x_coords)
        self.assertTrue(
            all(x1 <= x2 for x1, x2 in zip(top_start_x_coords, top_start_x_coords[1:])),
            f"Top members are not ordered correctly: {top_start_x_coords}"
        )

        # Extract the x-coordinates of the start points for botmembers
        bot_start_x_coords = [member.startcoord[0] for member in self.airfoil.botmembers]
        # Ensure the x-coordinates are in ascending order
        self.assertTrue(
            all(x1 <= x2 for x1, x2 in zip(bot_start_x_coords, bot_start_x_coords[1:])),
            f"Bottom members are not ordered correctly: {bot_start_x_coords}"
        )

    def test_top_height_at(self):
        """Test interpolation of top height at a given chordwise position."""
        x = 100
        height = self.airfoil.Top_height_at(x)
        self.assertTrue(height > 0, "Top height should be positive.")
        self.assertIsInstance(height, float, "Top height should be a float.")

    def test_bot_height_at(self):
        """Test interpolation of bottom height at a given chordwise position."""
        x = 150
        height = self.airfoil.Bot_height_at(x)
        self.assertTrue(height <= 0, "Bottom height should be non-positive.")
        self.assertIsInstance(height, float, "Bottom height should be a float.")

    def test_neutral_points(self):
        """Test calculation of neutral points."""
        xbar, ybar = self.airfoil.Neutralpoints()
        self.assertIsInstance(xbar, float, "Neutral x-coordinate should be a float.")
        self.assertIsInstance(ybar, float, "Neutral y-coordinate should be a float.")

    def test_neutral_points_within_bounds(self):
        """Test that the neutral points (xbar, ybar) are within reasonable bounds, namely within the airfoil section"""
        # Call the method to calculate neutral points
        xbar, ybar = self.airfoil.Neutralpoints()

        max_y = max(self.airfoil.topcoordinates, key=lambda coord: coord[1])[1]
        min_y = min(self.airfoil.botcoordinates, key=lambda coord: coord[1])[1]

        ybar_min = min_y
        ybar_max = max_y
        # Define acceptable bounds for xbar and ybar
        xbar_min, xbar_max = 0, self.airfoil.chordlength  # xbar should be within the chord length

        # Assert that xbar and ybar fall within the bounds
        self.assertTrue(xbar_min <= xbar <= xbar_max, f"xbar {xbar} is out of bounds ({xbar_min}, {xbar_max}).")
        self.assertTrue(ybar_min <= ybar <= ybar_max, f"ybar {ybar} is out of bounds ({ybar_min}, {ybar_max}).")

    def test_EI_calculation(self):
        """Test calculation of bending stiffness (EI)."""
        # Ensure neutral points are calculated first
        self.airfoil.Neutralpoints()

        # Now calculate bending stiffness
        EIxx, EIyy, EIxy = self.airfoil.CalculateEI()
        self.assertTrue(EIxx > 0, "EIxx should be positive.")
        self.assertTrue(EIyy > 0, "EIyy should be positive.")
        self.assertTrue(EIxy >= 0, "EIxy should be non-negative.")

    def test_curvatures(self):
        """Test curvature calculations."""
        self.airfoil.Neutralpoints()

        # Now calculate bending stiffness
        self.airfoil.CalculateEI()
        kx, ky = self.airfoil.curvatures(1000, 500)  # Example moment values
        self.assertIsInstance(kx, float, "kx should be a float.")
        self.assertIsInstance(ky, float, "ky should be a float.")

    def test_section_shear_flow(self):
        """Test shear flow calculations."""

        """Test curvature calculations."""
        self.airfoil.Neutralpoints()

        # Now calculate bending stiffness
        self.airfoil.CalculateEI()

        self.airfoil.Mx = 1000
        self.airfoil.My = 200
        self.airfoil.SectionShearFlows(100, 50)  # Example shear force values
        self.assertTrue(hasattr(self.airfoil, 'sectionlist'), "Airfoil should have a sectionlist attribute after shear flow calculation.")
        self.assertGreater(len(self.airfoil.sectionlist), 0, "Section list should not be empty.")

if __name__ == "__main__":
    unittest.main()
