import os
import subprocess

folder_path = 'scripts/analysis'

for filename in os.listdir(folder_path):
    if filename.endswith('.py'):
        full_path = os.path.join(folder_path, filename)
        print(f'Running {filename}...')
        subprocess.run(['python', full_path])
