from cx_Freeze import setup, Executable

# Dependencies are automatically detected, might need fine-tuning
build_options = {'packages': [], 'excludes': ['tkinter']}

base = 'console'

# The script being bundled is just randomly selected for the initial GitHub Action POC
executables = [
    Executable('scripts/boinc/blind_delta.py', base='console', target_name='BlindDeltaAlgo')
]

setup(name='BlindDeltaAlgo',
      version='0.1',
      description='Blind Delta Algorithm',
      options={'build_exe': build_options},
      executables=executables)
