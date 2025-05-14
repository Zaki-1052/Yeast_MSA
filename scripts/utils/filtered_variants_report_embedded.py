import base64
import os
import re
from pathlib import Path

def embed_images(html_file_path):
    # Read the HTML file
    with open(html_file_path, 'r') as file:
        html_content = file.read()
    
    # Get directory of the HTML file
    html_dir = os.path.dirname(html_file_path)
    
    # Find all image references
    img_pattern = re.compile(r'<img src="([^"]+)"')
    matches = img_pattern.findall(html_content)
    
    # Replace each image reference with base64 data
    for img_path in matches:
        try:
            # Get absolute path to the image
            absolute_img_path = os.path.join(html_dir, img_path)
            
            # Get the image MIME type
            extension = Path(img_path).suffix.lower()
            mime_type = {
                '.png': 'image/png',
                '.jpg': 'image/jpeg',
                '.jpeg': 'image/jpeg',
                '.gif': 'image/gif',
                '.svg': 'image/svg+xml'
            }.get(extension, 'image/png')
            
            # Read and encode the image
            with open(absolute_img_path, 'rb') as img_file:
                img_data = base64.b64encode(img_file.read()).decode('utf-8')
            
            # Create data URL
            data_url = f'data:{mime_type};base64,{img_data}'
            
            # Replace in HTML
            html_content = html_content.replace(f'src="{img_path}"', f'src="{data_url}"')
            print(f"Embedded: {img_path}")
            
        except Exception as e:
            print(f"Error embedding {img_path}: {e}")
    
    # Write the modified HTML to a new file
    output_path = html_file_path.replace('.html', '_embedded.html')
    with open(output_path, 'w') as file:
        file.write(html_content)
    
    return output_path

if __name__ == "__main__":
    html_file = "/Users/zakiralibhai/Documents/GitHub/Yeast_MSA/results/filtered_scaffold_variants/visualizations/filtered_variants_report.html"
    output_file = embed_images(html_file)
    print(f"Created embedded HTML file: {output_file}")