import unittest
from blind_delta import search
import os
import pandas as pd


class TestBoincWorker(unittest.TestCase):
	def test_worker_vs_mp_version(self):
		# Run original version
		search(2000, 2, [3, 3], 0, 1, 100000, -1010, -2020, 1000000000, 3)
		mp_output_filename = 'BlindDelta[3, 3]_0_1.csv'

		# Run worker, as it would in BOINC
		test_config_file = os.path.join('tests_data', 'test_job.json')
		worker_output_filename = os.path.join('tests_data', 'test_output.csv')
		os.system(f'execute_from_json.py {test_config_file} {worker_output_filename}')

		mp_output = pd.read_csv(mp_output_filename)
		worker_output = pd.read_csv(worker_output_filename)

		# There is a mismatch in the way we strore the coefficents, in one we use () and the other []
		worker_output['Coefficients'] = worker_output['Coefficients'].apply(lambda x: x.strip('()[]'))
		mp_output['Coefficients'] = mp_output['Coefficients'].apply(lambda x: x.strip('()[]'))

		# In more complicated search spaces, the two versions might not be sortted. 
		# Here it is not a problem.
		assert all(worker_output == mp_output)


if __name__ == '__main__':
    unittest.main()
