#!/usr/bin/env python3
"""
find_genes_of_interest.py - Find genes of interest in W303 yeast strain annotation
Using multiple search strategies to identify S. cerevisiae genes in the W303 annotation
"""

import os
import sys
import re
from collections import defaultdict
from Bio import SeqIO

# Define known information about our genes of interest (from SGD database)
GENES_INFO = {
    "YHR190W": {
        "name": "ERG9",
        "aliases": ["ERG9", "SRT1"],
        "description": "Squalene synthase; catalyzes the conversion of two farnesyl pyrophosphate molecules to squalene",
        "keywords": ["squalene", "synthase", "farnesyl", "pyrophosphate", "ergosterol"]
    },
    "YGR175C": {
        "name": "ERG1",
        "aliases": ["ERG1"],
        "description": "Squalene epoxidase; catalyzes the epoxidation of squalene to 2,3-oxidosqualene",
        "keywords": ["squalene", "epoxidase", "oxidosqualene", "ergosterol"]
    },
    "YHR072W": {
        "name": "ERG7",
        "aliases": ["ERG7"],
        "description": "Lanosterol synthase; catalyzes cyclization of 2,3-oxidosqualene to lanosterol",
        "keywords": ["lanosterol", "synthase", "oxidosqualene", "ergosterol"]
    },
    "YHR007C": {
        "name": "ERG11",
        "aliases": ["ERG11", "CYP51"],
        "description": "Lanosterol 14-alpha-demethylase; catalyzes the C-14 demethylation of lanosterol",
        "keywords": ["lanosterol", "demethylase", "ergosterol", "CYP51"]
    },
    "YNL280C": {
        "name": "ERG24",
        "aliases": ["ERG24"],
        "description": "C-14 sterol reductase; catalyzes the reduction of the C-14 double bond of 4,4-dimethylcholesta-8,14,24-triene-3-beta-ol",
        "keywords": ["sterol", "reductase", "ergosterol"]
    },
    "YGR060W": {
        "name": "ERG25",
        "aliases": ["ERG25"],
        "description": "C-4 methyl sterol oxidase; catalyzes the removal of C-4 methyl groups from the sterol ring",
        "keywords": ["sterol", "oxidase", "ergosterol", "methyl"]
    },
    "YML008C": {
        "name": "ERG6",
        "aliases": ["ERG6", "ISE1", "LIS1", "SED6"],
        "description": "Delta(24)-sterol C-methyltransferase; converts zymosterol to fecosterol",
        "keywords": ["sterol", "methyltransferase", "zymosterol", "fecosterol", "ergosterol"]
    },
    "YMR202W": {
        "name": "ERG2",
        "aliases": ["ERG2"],
        "description": "C-8 sterol isomerase; catalyzes the isomerization of the delta-8 double bond to the delta-7 position",
        "keywords": ["sterol", "isomerase", "ergosterol"]
    },
    "YLR056W": {
        "name": "ERG3",
        "aliases": ["ERG3", "SYR1"],
        "description": "C-5 sterol desaturase; catalyzes the introduction of a C-5(6) double bond into episterol",
        "keywords": ["sterol", "desaturase", "episterol", "ergosterol"]
    },
    "YMR015C": {
        "name": "ERG5",
        "aliases": ["ERG5"],
        "description": "C-22 sterol desaturase; catalyzes the formation of the C-22(23) double bond in the sterol side chain",
        "keywords": ["sterol", "desaturase", "ergosterol"]
    },
    "YGL012W": {
        "name": "ERG4",
        "aliases": ["ERG4"],
        "description": "C-24(28) sterol reductase; reduces the C-24(28) double bond during ergosterol biosynthesis",
        "keywords": ["sterol", "reductase", "ergosterol"]
    }
}

def search_genes_of_interest(genbank_dir, output_file):
    """
    Search for genes of interest using multiple strategies
    
    Args:
        genbank_dir (str): Directory containing GenBank files
        output_file (str): Output file to write results
    
    Returns:
        dict: Dictionary of matches found
    """
    print(f"Searching for genes of interest in {genbank_dir}...")
    
    # Initialize results container
    results = defaultdict(list)
    exact_matches = set()
    all_annotations = []
    
    # Process each GenBank file
    files_processed = 0
    for filename in sorted(os.listdir(genbank_dir)):
        if not filename.endswith(('.genbank', '.gb', '.gbk')):
            continue
            
        filepath = os.path.join(genbank_dir, filename)
        files_processed += 1
        
        if files_processed % 50 == 0:
            print(f"Processed {files_processed} files...")
        
        try:
            # Parse the GenBank file
            record = SeqIO.read(filepath, "genbank")
            
            # Extract scaffold ID
            scaffold_id = record.id
            if scaffold_id == "unknown" and "_scaffold_" in filename:
                scaffold_num = filename.split("_scaffold_")[1].split(".")[0]
                scaffold_id = f"scaffold_{scaffold_num}"
                
            # Extract JRIU ID from source feature
            jriu_id = None
            for feature in record.features:
                if feature.type == "source":
                    for note in feature.qualifiers.get('note', []):
                        if note.startswith('JRIU'):
                            jriu_id = note
                            break
            
            # Process each feature to look for genes of interest
            for feature in record.features:
                # Skip non-gene features
                if feature.type not in ["gene", "CDS"]:
                    continue
                
                # Extract basic feature information
                gene_id = None
                if 'gene' in feature.qualifiers:
                    gene_id = feature.qualifiers['gene'][0]
                elif 'locus_tag' in feature.qualifiers:
                    gene_id = feature.qualifiers['locus_tag'][0]
                
                if not gene_id:
                    continue
                
                location = feature.location
                start = int(location.start)
                end = int(location.end)
                strand = '+' if location.strand == 1 else '-'
                
                # Collect all annotation text to search
                annotation_text = []
                
                # Extract qualifiers
                product = feature.qualifiers.get('product', [''])[0]
                note = ' '.join(feature.qualifiers.get('note', ['']))
                inference = ' '.join(feature.qualifiers.get('inference', ['']))
                
                # Add all text to the annotation
                annotation_text = [product, note, inference]
                annotation_text = ' '.join([t for t in annotation_text if t]).lower()
                
                # Store full annotation for later analysis
                all_annotations.append({
                    'gene_id': gene_id,
                    'scaffold_id': scaffold_id,
                    'jriu_id': jriu_id,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'product': product,
                    'note': note,
                    'inference': inference,
                    'annotation_text': annotation_text
                })
                
                # Strategy 1: Direct gene ID match
                for sc_gene_id, info in GENES_INFO.items():
                    if sc_gene_id in inference or sc_gene_id in note:
                        if f"{gene_id}_{sc_gene_id}" not in exact_matches:
                            results[sc_gene_id].append({
                                'match_type': 'direct_id',
                                'gene_id': gene_id,
                                'scaffold_id': scaffold_id,
                                'jriu_id': jriu_id,
                                'start': start,
                                'end': end,
                                'strand': strand,
                                'product': product,
                                'note': note,
                                'confidence': 'high'
                            })
                            exact_matches.add(f"{gene_id}_{sc_gene_id}")
                
                # Strategy 2: Gene name/alias match
                for sc_gene_id, info in GENES_INFO.items():
                    for alias in info['aliases']:
                        alias_pattern = r'\b' + re.escape(alias) + r'\b'
                        if re.search(alias_pattern, annotation_text, re.IGNORECASE) and f"{gene_id}_{sc_gene_id}" not in exact_matches:
                            results[sc_gene_id].append({
                                'match_type': 'alias',
                                'gene_id': gene_id,
                                'scaffold_id': scaffold_id,
                                'jriu_id': jriu_id,
                                'start': start,
                                'end': end,
                                'strand': strand,
                                'product': product,
                                'note': note,
                                'confidence': 'high',
                                'matched_alias': alias
                            })
                            exact_matches.add(f"{gene_id}_{sc_gene_id}")
        
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
    
    # Strategy 3: Keyword/function matching for remaining unmatched genes
    unmatched_genes = [g for g in GENES_INFO.keys() if g not in results]
    if unmatched_genes:
        print(f"\nSearching for {len(unmatched_genes)} unmatched genes using keyword matching...")
        
        for sc_gene_id in unmatched_genes:
            info = GENES_INFO[sc_gene_id]
            keyword_matches = []
            
            # Count keyword matches in each annotation
            for annotation in all_annotations:
                annotation_text = annotation['annotation_text']
                match_score = 0
                matched_keywords = []
                
                # Check each keyword
                for keyword in info['keywords']:
                    if keyword.lower() in annotation_text:
                        match_score += 1
                        matched_keywords.append(keyword)
                
                if match_score >= 2:  # Require at least 2 keyword matches for confidence
                    keyword_matches.append({
                        'gene_id': annotation['gene_id'],
                        'scaffold_id': annotation['scaffold_id'],
                        'jriu_id': annotation['jriu_id'],
                        'start': annotation['start'],
                        'end': annotation['end'],
                        'strand': annotation['strand'],
                        'product': annotation['product'],
                        'note': annotation['note'],
                        'match_score': match_score,
                        'matched_keywords': matched_keywords
                    })
            
            # Sort by match score and add top matches
            if keyword_matches:
                keyword_matches.sort(key=lambda x: x['match_score'], reverse=True)
                for match in keyword_matches[:3]:  # Take top 3 matches
                    confidence = 'medium' if match['match_score'] >= 3 else 'low'
                    results[sc_gene_id].append({
                        'match_type': 'keyword',
                        'gene_id': match['gene_id'],
                        'scaffold_id': match['scaffold_id'],
                        'jriu_id': match['jriu_id'],
                        'start': match['start'],
                        'end': match['end'],
                        'strand': match['strand'],
                        'product': match['product'],
                        'note': match['note'],
                        'confidence': confidence,
                        'matched_keywords': ', '.join(match['matched_keywords']),
                        'match_score': match['match_score']
                    })
    
    # Write results to file
    with open(output_file, 'w') as f:
        # Write header
        f.write("sc_gene_id\tsc_gene_name\tmatch_type\tconfidence\tw303_gene_id\tscaffold_id\tjriu_id\tstart\tend\tstrand\tproduct\tmatched_terms\n")
        
        # Write data rows
        for sc_gene_id, matches in sorted(results.items()):
            sc_gene_name = GENES_INFO[sc_gene_id]['name']
            
            for match in matches:
                match_type = match['match_type']
                gene_id = match['gene_id']
                scaffold_id = match['scaffold_id']
                jriu_id = match['jriu_id'] or 'unknown'
                start = match['start']
                end = match['end']
                strand = match['strand']
                product = match.get('product', '').replace('\t', ' ').replace('\n', ' ')
                confidence = match.get('confidence', 'unknown')
                
                # Get matched terms based on match type
                if match_type == 'direct_id':
                    matched_terms = sc_gene_id
                elif match_type == 'alias':
                    matched_terms = match.get('matched_alias', '')
                elif match_type == 'keyword':
                    matched_terms = match.get('matched_keywords', '')
                else:
                    matched_terms = ''
                
                f.write(f"{sc_gene_id}\t{sc_gene_name}\t{match_type}\t{confidence}\t{gene_id}\t{scaffold_id}\t{jriu_id}\t{start}\t{end}\t{strand}\t{product}\t{matched_terms}\n")
    
    return results

def print_results_summary(results):
    """Print a summary of the search results"""
    print("\n=== Genes of Interest Search Results ===\n")
    
    found_genes = 0
    total_matches = 0
    
    for sc_gene_id, matches in sorted(results.items()):
        if matches:
            found_genes += 1
            total_matches += len(matches)
            
            print(f"{sc_gene_id} ({GENES_INFO[sc_gene_id]['name']}): {len(matches)} potential match(es)")
            
            for i, match in enumerate(matches):
                confidence = match.get('confidence', 'unknown')
                match_type = match['match_type']
                gene_id = match['gene_id']
                product = match.get('product', '')[:50] + ('...' if len(match.get('product', '')) > 50 else '')
                
                if match_type == 'direct_id':
                    match_detail = f"Direct ID match"
                elif match_type == 'alias':
                    match_detail = f"Alias match: {match.get('matched_alias', '')}"
                elif match_type == 'keyword':
                    match_detail = f"Keyword match (score: {match.get('match_score', 0)}): {match.get('matched_keywords', '')}"
                else:
                    match_detail = ""
                
                print(f"  {i+1}. {gene_id} - {product} [{confidence} confidence] - {match_detail}")
            
            print("")
    
    print(f"Summary: Found {found_genes}/{len(GENES_INFO)} genes of interest with {total_matches} total matches")
    
    if found_genes < len(GENES_INFO):
        missing_genes = [f"{g} ({GENES_INFO[g]['name']})" for g in GENES_INFO if g not in results or not results[g]]
        print(f"Missing genes: {', '.join(missing_genes)}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python find_genes_of_interest.py <genbank_directory> [output_file]")
        sys.exit(1)
    
    genbank_dir = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "genes_of_interest_search.tsv"
    
    if not os.path.isdir(genbank_dir):
        print(f"Error: {genbank_dir} is not a directory")
        sys.exit(1)
    
    results = search_genes_of_interest(genbank_dir, output_file)
    print_results_summary(results)
    print(f"\nComplete results written to: {output_file}")

if __name__ == "__main__":
    main()