import os
import glob
import numpy as np
import subprocess
import unittest
import shutil

class TestFilterAntenna(unittest.TestCase):
    
    def setUp(self):
        """Set up a temporary workspace and generate dummy antenna data."""
        self.test_dir = "test_filter_workspace"
        os.makedirs(self.test_dir, exist_ok=True)
        
        # Get absolute path to the script before changing directories
        self.script_path = os.path.abspath("filter_antenna.py")
        
        self.orig_dir = os.getcwd()
        os.chdir(self.test_dir)

        # Generate dummy data mimicking {prefix}_{antno}_{pol}_f.npy
        
        # Antenna 0: Valid X, Valid Y (Should pass X, Y, and Dual tests)
        np.save("test_0_X_f.npy", np.ones(10))
        np.save("test_0_Y_f.npy", np.ones(10))

        # Antenna 1: Valid X, Missing Y (Should pass X test, fail Y and Dual tests)
        np.save("test_1_X_f.npy", np.ones(10))

        # Antenna 2: Valid X, Zero-size Y (Should pass X test, fail Y and Dual tests)
        np.save("test_2_X_f.npy", np.ones(10))
        np.save("test_2_Y_f.npy", np.array([]))

        # Antenna 3: Missing X, Valid Y (Should fail X and Dual tests, pass Y test)
        np.save("test_3_Y_f.npy", np.ones(10))

    def tearDown(self):
        """Clean up the temporary workspace after tests."""
        os.chdir(self.orig_dir)
        shutil.rmtree(self.test_dir)

    def run_filter(self, pols):
        """Helper to run the filter_antenna script."""
        cmd = ["python3", self.script_path, "--pols"] + pols
        subprocess.run(cmd, check=True, capture_output=True)

    def get_filtered_ants(self):
        """Helper to check which antennas successfully passed and generated _filtered.npy files."""
        filtered_files = glob.glob("*_filtered.npy")
        # Extract antenna numbers from the filenames
        ants = set([int(f.split('_')[1]) for f in filtered_files])
        return ants

    def test_dual_pol(self):
        """Test with --pols X Y"""
        self.run_filter(["X", "Y"])
        ants = self.get_filtered_ants()
        # Only Antenna 0 has both valid X and Y files
        self.assertEqual(ants, {0}, f"Expected only Ant 0 to pass Dual Pol, got {ants}")

    def test_single_pol_x(self):
        """Test with --pols X"""
        self.run_filter(["X"])
        ants = self.get_filtered_ants()
        # Antennas 0, 1, and 2 have valid X files. (Ant 3 is missing X entirely).
        self.assertEqual(ants, {0, 1, 2}, f"Expected Ants 0, 1, 2 to pass Single Pol X, got {ants}")

    def test_single_pol_y(self):
        """Test with --pols Y"""
        self.run_filter(["Y"])
        ants = self.get_filtered_ants()
        # Antennas 0 and 3 have valid Y files. (Ant 1 is missing Y, Ant 2 has zero-size Y).
        self.assertEqual(ants, {0, 3}, f"Expected Ants 0, 3 to pass Single Pol Y, got {ants}")

if __name__ == '__main__':
    # Run the tests
    print("Running filter_antenna.py tests...")
    unittest.main(verbosity=2)
