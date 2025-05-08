def load_gene_mapping():
    """
    Load gene mapping data from reference files.
    
    This function loads gene data from the reference directory, including:
    1. Gene coordinates and information from gene_mapping_full.tsv
    2. Genes of interest (ergosterol pathway) from genes_of_interest_mapping.tsv
    
    Returns:
        bool: True if data was loaded successfully, False otherwise
    """
    global GENE_DATA, SCAFFOLD_GENES, GENES_OF_INTEREST
    
    # Clear existing data
    GENE_DATA.clear()
    SCAFFOLD_GENES.clear()
    GENES_OF_INTEREST.clear()
    
    # Define possible file paths for gene mapping data, prioritizing gene_mapping_full.tsv
    gene_mapping_paths = [
        "reference/gene_mapping_full.tsv",  # New comprehensive mapping
        "reference/gene_mapping.tsv",       # Fallback to original mapping if full doesn't exist
    ]
    
    # Define path for ergosterol pathway genes (genes of interest)
    genes_of_interest_path = "reference/genes_of_interest_mapping.tsv"
    
    # Load gene mapping data
    gene_mapping_file = None
    for path in gene_mapping_paths:
        if os.path.exists(path):
            gene_mapping_file = path
            break
    
    if gene_mapping_file:
        try:
            # Load gene mapping data
            gene_df = pd.read_csv(gene_mapping_file, sep='\t')
            print(f"Loaded {len(gene_df)} genes from {gene_mapping_file}")
            logging.info(f"Loaded {len(gene_df)} genes from {gene_mapping_file}")
            
            # Process each gene to build GENE_DATA and SCAFFOLD_GENES
            for _, row in gene_df.iterrows():
                gene_id = row['w303_gene_id']
                scaffold = row['w303_scaffold']
                
                # Store gene data
                gene_data = {
                    'gene_id': gene_id,
                    'locus_tag': row['locus_tag'] if 'locus_tag' in row else None,
                    'sc_gene_id': row['sc_gene_id'] if 'sc_gene_id' in row else None,
                    'erg_name': row['erg_name'] if 'erg_name' in row else None,
                    'scaffold': scaffold,
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'strand': row['strand'] if 'strand' in row else None,
                    'product': row['product'] if 'product' in row else None
                }
                
                # Add chromosome_id if it's in the data
                if 'chromosome_id' in row and not pd.isna(row['chromosome_id']):
                    gene_data['chromosome_id'] = row['chromosome_id']
                
                GENE_DATA[gene_id] = gene_data
                
                # Map scaffold to genes
                if scaffold not in SCAFFOLD_GENES:
                    SCAFFOLD_GENES[scaffold] = []
                SCAFFOLD_GENES[scaffold].append(gene_id)
                
                # Also map chromosome_id to genes if available
                if 'chromosome_id' in gene_data and gene_data['chromosome_id']:
                    chromosome_id = gene_data['chromosome_id']
                    if chromosome_id not in SCAFFOLD_GENES:
                        SCAFFOLD_GENES[chromosome_id] = []
                    SCAFFOLD_GENES[chromosome_id].append(gene_id)
            
            # Now, load the genes of interest (ergosterol pathway genes) from a separate file
            if os.path.exists(genes_of_interest_path):
                try:
                    # Load genes of interest data
                    erg_df = pd.read_csv(genes_of_interest_path, sep='\t')
                    print(f"Loaded {len(erg_df)} ergosterol pathway genes from {genes_of_interest_path}")
                    logging.info(f"Loaded {len(erg_df)} ergosterol pathway genes from {genes_of_interest_path}")
                    
                    # Add these genes to GENES_OF_INTEREST
                    for _, row in erg_df.iterrows():
                        if 'w303_gene_id' in row and not pd.isna(row['w303_gene_id']):
                            GENES_OF_INTEREST.add(row['w303_gene_id'])
                except Exception as e:
                    print(f"Error loading ergosterol genes: {e}")
                    logging.error(f"Error loading ergosterol genes: {e}")
                    # If this fails, try to identify them from the main gene data
                    for gene_id, gene_data in GENE_DATA.items():
                        if gene_data['erg_name'] and (gene_data['erg_name'].startswith('ERG') or 'ERG' in gene_data['erg_name']):
                            GENES_OF_INTEREST.add(gene_id)
            else:
                # If the genes_of_interest file doesn't exist, try to identify them from the main gene data
                print(f"No specific ergosterol genes file found. Identifying based on gene names...")
                logging.warning(f"No specific ergosterol genes file found. Identifying based on gene names...")
                
                for gene_id, gene_data in GENE_DATA.items():
                    if gene_data['erg_name'] and (gene_data['erg_name'].startswith('ERG') or 'ERG' in gene_data['erg_name']):
                        GENES_OF_INTEREST.add(gene_id)
            
            # Print statistics about loaded data
            genes_with_chroms = sum(1 for gene_data in GENE_DATA.values() if 'chromosome_id' in gene_data)
            print(f"Loaded {len(GENE_DATA)} genes, {genes_with_chroms} have chromosome IDs")
            print(f"Mapped genes to {len(SCAFFOLD_GENES)} scaffolds/chromosomes")
            print(f"Identified {len(GENES_OF_INTEREST)} ergosterol pathway genes")
            
            # Print sample of genes of interest
            print(f"Sample ergosterol genes: {list(GENES_OF_INTEREST)[:5]}")
            
            # Print sample of chromosome IDs and scaffolds
            chrom_ids = set(gene_data.get('chromosome_id') for gene_data in GENE_DATA.values() 
                        if 'chromosome_id' in gene_data)
            scaffold_ids = set(gene_data['scaffold'] for gene_data in GENE_DATA.values())
            
            print(f"Sample chromosome IDs: {list(chrom_ids)[:5]}")
            print(f"Sample scaffold IDs: {list(scaffold_ids)[:5]}")
            print(f"Sample SCAFFOLD_GENES keys: {list(SCAFFOLD_GENES.keys())[:10]}")
            
            logging.info(f"Loaded {len(GENE_DATA)} genes, {genes_with_chroms} have chromosome IDs")
            logging.info(f"Mapped genes to {len(SCAFFOLD_GENES)} scaffolds/chromosomes")
            logging.info(f"Identified {len(GENES_OF_INTEREST)} ergosterol pathway genes")
            
            return True
        except Exception as e:
            print(f"Error loading gene mapping data: {e}")
            logging.error(f"Error loading gene mapping data: {e}")
            return False
    else:
        print("No gene mapping file found. Gene mapping will not be available.")
        logging.warning("No gene mapping file found. Gene mapping will not be available.")
        return False