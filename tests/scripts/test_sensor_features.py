import unittest
import hashlib
import pandas as pd
import utils
import yaml
import os

class TestSensorFeatures(unittest.TestCase):

    # Hack to run code in positional order (not 100% full proof)
    unittest.TestLoader.sortTestMethodsUsing = lambda self, a, b: (a < b) - (a > b)
    
    @classmethod
    def setUpClass(cls):
        # Runs once to Setup env
        global configs 
        with open(r'tests/settings/testing_config.yaml') as file:
            configs = yaml.full_load(file)


    def test_sensors_files_exist(self):
        # Loop through the file_list dictionary and check if the files exist. 

        file_lists = utils.generate_sensor_file_lists(configs)
        for each in file_lists:
            for out_file, _ in file_lists[each]:
                self.assertEqual(os.path.exists(out_file), 1)


    def test_sensors_features_calculations(self):
        calc_files = utils.generate_sensor_file_lists(configs)
        for each in calc_files:
            for act_result, exp_result in calc_files[each]:
                df_act = pd.read_csv(act_result)
                df_exp = pd.read_csv(exp_result)
                if df_act.empty:
                    self.assertTrue(df_exp.empty)
                else:
                    pd.testing.assert_frame_equal(df_exp, df_act, obj=df_exp)


if __name__ == '__main__':

    unittest.main()