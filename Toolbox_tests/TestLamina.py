import unittest
import numpy as np
from Toolbox.lamina import Lamina

class TestLamina(unittest.TestCase):
    def setUp(self):
        # Setup a test Lamina object with arbitrary properties
        self.elastic_properties = [135000, 11200, 5000, 0.3]  # Example: E1, E2, G12, v12
        self.failure_properties = [230000, 0.2, 1.1, 2550, 1470, 69, 300, 100]  # Example values
        self.theta = 45  # Fiber orientation in degrees
        self.t = 0.125  # Thickness of the lamina
        self.test_lamina = Lamina(t=self.t, theta=self.theta,
                                  elasticproperties=self.elastic_properties,
                                  failureproperties=self.failure_properties)

    def test_calculate_Q_matrix(self):
        """Test the calculation of the stiffness matrix Q."""
        self.test_lamina.calculate_QS()
        Q = self.test_lamina.Q
        self.assertEqual(Q.shape, (3, 3))  # Check if Q is a 3x3 matrix
        self.assertTrue(np.all(Q > 0), "All entries in Q should be positive.")

    def test_stress_analysis(self):
        """Test stress analysis with a given strain state."""
        # Assign a strain state (example values)
        self.test_lamina.Epsilon = np.array([0.001, 0.0005, 0.0001])
        stresses = self.test_lamina.stress_analysis()
        self.assertEqual(stresses.shape, (3, ))  # Check if stress is a 3, array
        self.assertTrue(np.all(stresses >= 0), "All stresses should be non-negative in this test case.")

    def test_failure_analysis(self):
        """Test the failure analysis under a stress state."""
        sigma = np.array([200, 50, 30])  # Example stress state
        failure, IFFfactor, FFfactor = self.test_lamina.failure_analysis(sigma)
        self.assertTrue(failure in [0, 1, 2], "Failure state should be 0, 1, or 2.")
        self.assertGreaterEqual(IFFfactor, 0, "IFF factor should be non-negative.")
        self.assertGreaterEqual(FFfactor, 0, "FF factor should be non-negative.")

if __name__ == "__main__":
    unittest.main()
