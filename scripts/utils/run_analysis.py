import os
import subprocess
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("analysis.log"),
        logging.StreamHandler()
    ]
)

# Analysis script file path
folder_path = 'scripts/analysis'

# Create results directory if it doesn't exist
os.makedirs('analysis', exist_ok=True)

# Define the order in which scripts should run for dependency reasons
# This ensures that dependencies are generated before they're needed
ordered_scripts = [
    'mutation_spectrum_analysis.py',
    'genomic_context_analysis.py',     # Process genomic context first
    'statistical_pattern_analysis.py',  # Statistical analysis next
    'mutational_signature_analysis.py', # Generate mutation signatures 
    'regional_enrichment_analysis.py',  # Regional analysis
    'variation.py',                    # Process variation data
    'scaffold_distribution_analysis.py', # Scaffold analysis
    'population_spectrum_analysis.py',  # Population analysis
    'TC_Visualization.py'             # Generate visualizations     # Mutation spectrum analysis
]

# Run scripts in specified order
for script in ordered_scripts:
    full_path = os.path.join(folder_path, script)
    if os.path.exists(full_path):
        try:
            logging.info(f'Running {script}...')
            subprocess.run(['python3', full_path], check=True)
            logging.info(f'Successfully completed {script}')
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running {script}: {e}")
    else:
        logging.warning(f"Script {script} not found in {folder_path}")

# Run any additional scripts that weren't in the ordered list
for filename in os.listdir(folder_path):
    if filename.endswith('.py') and filename not in ordered_scripts:
        full_path = os.path.join(folder_path, filename)
        try:
            logging.info(f'Running additional script {filename}...')
            subprocess.run(['python3', full_path], check=True)
            logging.info(f'Successfully completed {filename}')
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running {filename}: {e}")

logging.info("All analysis scripts completed")
