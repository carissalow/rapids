from snakemake.io import expand
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
            #for out_file, _ in file_lists[each]:
            self.assertEqual(os.path.exists( each[0]), 1)


    def test_sensors_features_calculations(self):

                       
        sensor_file_list = utils.generate_sensor_file_lists(configs)        
        for each in sensor_file_list:
            for act_result, exp_result in sensor_file_list:
                df_act = pd.read_csv(act_result)
                df_exp = pd.read_csv(exp_result)
                if df_act.empty:
                    self.assertTrue(df_exp.empty)
                else:
                    # The order in which the columns varies from time to time so
                    # the columns are sorted before doing the comparision 
                    df_exp = df_exp.reindex(sorted(df_exp.columns), axis=1)
                    df_act = df_act.reindex(sorted(df_act.columns), axis=1)
                    pd.testing.assert_frame_equal(df_exp, df_act, obj=df_exp)


if __name__ == '__main__':

    unittest.main()