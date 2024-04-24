"""
This file is the build configuration for cx_Freeze to generate an executable binary across platforms
"""
from cx_Freeze import setup, Executable

build_options = {'packages': [],
                 'excludes': ['tkinter'],
                 'build_exe': 'build/executable' # sets name of folder under build folder in which executables end up
                 }

base = 'console'

executables = [
    Executable('execute_from_json.py', base='console', target_name='BlindDeltaAlgo')
]

setup(name='BlindDeltaAlgo',
      version='0.1',
      description='Blind Delta Algorithm',
      options={'build_exe': build_options},
      executables=executables)
