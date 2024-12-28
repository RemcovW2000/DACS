import unittest
import numpy as np
from Toolbox.airfoil import Airfoil
from Toolbox.member import Member
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

    def test_segment_I(self):
        t = 1
        P1_v = [0, 0]
        P2_v = [0, 1]

        P1_h = [0, 0]
        P2_h = [1, 0]

        P1_d = [0, 0]
        P2_d = [1, 1]

        P1_vneg = [0, 0]
        P2_vneg = [0, -1]

        P1_hneg = [0, 0]
        P2_hneg = [-1, 0]

        P1_dneg = [0, 0]
        P2_dneg = [-1, -1]

        Iv = self.airfoil.segment_I(t, P1_v, P2_v)
        Ih = self.airfoil.segment_I(t, P1_h, P2_h)
        Id = self.airfoil.segment_I(t, P1_d, P2_d)

        Ivneg = self.airfoil.segment_I(t, P1_vneg, P2_vneg)
        Ihneg = self.airfoil.segment_I(t, P1_hneg, P2_hneg)
        Idneg = self.airfoil.segment_I(t, P1_dneg, P2_dneg)

        self.assertEqual(Iv[2], -Ivneg[2])
        self.assertEqual(Ih[2], -Ihneg[2])
        self.assertEqual(Id[2], -Idneg[2])

        self.assertTrue(len(Iv) == 3)
        self.assertTrue(len(Ih) == 3)
        self.assertTrue(len(Id) == 3)

        self.assertTrue(Iv[0] > 0 and Iv[1] > 0)
        self.assertTrue(Ih[0] > 0 and Ih[1] > 0)
        self.assertTrue(Id[0] > 0 and Id[1] > 0)

    def test_top_botmemberlist_order(self):
        """Test the order of the top and bot member lists from leading to trailing edge."""

        # Extract the x-coordinates of the start points for topmembers
        top_start_x_coords = [member.startcoord[0] for member in self.airfoil.topmembers]
        # Ensure the x-coordinates are in ascending order
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
        height = self.airfoil.top_height_at(x)
        self.assertTrue(height > 0, "Top height should be positive.")
        self.assertIsInstance(height, float, "Top height should be a float.")

    def test_bot_height_at(self):
        """Test interpolation of bottom height at a given chordwise position."""
        x = 150
        height = self.airfoil.bot_height_at(x)
        self.assertTrue(height <= 0, "Bottom height should be non-positive.")
        self.assertIsInstance(height, float, "Bottom height should be a float.")

    def test_neutral_points(self):
        """Test calculation of neutral points."""
        xbar, ybar = self.airfoil.neutral_points()
        self.assertIsInstance(xbar, float, "Neutral x-coordinate should be a float.")
        self.assertIsInstance(ybar, float, "Neutral y-coordinate should be a float.")

    def test_neutral_points_within_bounds(self):
        """Test that the neutral points (xbar, ybar) are within reasonable bounds, namely within the airfoil section"""
        # Call the method to calculate neutral points
        xbar, ybar = self.airfoil.neutral_points()

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
        self.airfoil.neutral_points()

        # Now calculate bending stiffness
        EIxx, EIyy, EIxy = self.airfoil.calculate_EI()
        self.assertTrue(EIxx > 0, "EIxx should be positive.")
        self.assertTrue(EIyy > 0, "EIyy should be positive.")
        self.assertTrue(EIxy >= 0, "EIxy should be non-negative.")

    def test_curvatures(self):
        """Test curvature calculations."""
        self.airfoil.neutral_points()

        # Now calculate bending stiffness
        self.airfoil.calculate_EI()
        kx, ky = self.airfoil.curvatures(1000, 500)  # Example moment values
        self.assertIsInstance(kx, float, "kx should be a float.")
        self.assertIsInstance(ky, float, "ky should be a float.")

    def test_section_shear_flow(self):
        """Test shear flow calculations"""
        self.airfoil.neutral_points()

        # Now calculate bending stiffness
        self.airfoil.calculate_EI()

        self.airfoil.Mx = 1000
        self.airfoil.My = 200
        self.airfoil.section_shear_flows(100, 50)  # Example shear force values
        self.assertTrue(hasattr(self.airfoil, 'sectionlist'), "Airfoil should have a sectionlist attribute after shear flow calculation.")
        self.assertTrue(len(self.airfoil.sectionlist) == len(self.airfoil.sparlocations) + 1, "Section list should be size of sparlocations + 1.")

    def test_nr_of_booms(self):
        self.airfoil.neutral_points()

        # Now calculate bending stiffness
        self.airfoil.calculate_EI()

        self.airfoil.Mx = 1000
        self.airfoil.My = 200
        self.airfoil.section_shear_flows(100, 50)

    def test_no_shared_subobjects(self):
        """
        Test to ensure no objects share sub-objects (deep copy verification).
        """

        self.airfoil.Mx = 1000
        self.airfoil.My = 200
        self.Sx, self.Sy = 100, 50

        # TODO: fix these dimensions so they are in the global (wing) FOR?
        self.E0, self.N0 = 0, 0 # relative to the local airfoil FOR

        # Functions for the normal stresses/deformations:
        self.airfoil.neutral_points() # calculates neutral point
        self.airfoil.calculate_EI() # then calculates EI around neutral point

        # Functions for shear stress solving
        self.airfoil.section_shear_flows(self.Sx, self.Sy)
        self.airfoil.section_shear_correction()
        self.airfoil.shear_superposition()
        # we want to check the total shear force exerted by the section:
        self.airfoil.shear_force_analysis()
        # Collect all sub-objects into a list for comparison
        subobjects = []

        # Traverse through the hierarchy
        for member in self.airfoil.topmembers + self.airfoil.botmembers + self.airfoil.sparmembers:
            subobjects.append(member)
            subobjects.extend(member.booms)  # Add booms to the list if present
            subobjects.extend(member.segments)

        # Ensure no two sub-objects have the same ID
        ids = [id(obj) for obj in subobjects]
        self.assertEqual(len(ids), len(set(ids)), "Sub-objects are shared between parent objects!")

    def test_section_order(self):
        self.airfoil.Mx = 1000
        self.airfoil.My = 200
        self.Sx, self.Sy = 100, 50

        self.E0, self.N0 = 0, 0  # relative to the local airfoil FOR

        # Functions for the normal stresses/deformations:
        self.airfoil.neutral_points()  # calculates neutral point
        self.airfoil.calculate_EI()  # then calculates EI around neutral point

        # Functions for shear stress solving
        self.airfoil.section_shear_flows(self.Sx, self.Sy)
        self.airfoil.section_shear_correction()
        self.airfoil.shear_superposition()
        # we want to check the total shear force exerted by the section:
        self.airfoil.shear_force_analysis()

        xstarts = []
        for i, section in enumerate(self.airfoil.sectionlist):
            x_start_firstmember = self.airfoil.sectionlist[i].members[0].startcoord[0]
            xstarts.append(x_start_firstmember)

        # order of start points of first members should increase:
        is_ascending = all(xstarts[i] <= xstarts[i + 1] for i in range(len(xstarts) - 1))
        self.assertTrue(is_ascending)


if __name__ == "__main__":
    unittest.main()
