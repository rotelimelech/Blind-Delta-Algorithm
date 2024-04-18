from blind_delta import calc_individual
import json
import sys

DEFAULT_OUTPUT_FILENAME = 'blind_delta_output.csv'


def main():
	if len(sys.argv) not in (2, 3) or sys.argv[1] in ('-h', '--help', '-?', '/?'):
		print('Usage:')
		print('execute_from_json.py coefficent_ranges.json [output_filename.csv]')

		exit()


if __name__ == '__main__':
	main()