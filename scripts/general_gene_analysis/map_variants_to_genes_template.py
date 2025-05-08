def map_variants_to_genes(variant_df):
    """
    Map variants to genes based on their genomic coordinates.
    
    Args:
        variant_df (pandas.DataFrame): DataFrame containing variant information
    
    Returns:
        pandas.DataFrame: The original DataFrame with additional gene-related columns
    """
    if not GENE_DATA or not SCAFFOLD_GENES:
        print("Gene mapping data not loaded. Cannot map variants to genes.")
        logging.warning("Gene mapping data not loaded. Cannot map variants to genes.")
        return variant_df
    
    # Create a copy of the input dataframe
    result_df = variant_df.copy()
    
    # Initialize gene-related columns
    result_df['in_gene'] = False
    result_df['gene_id'] = None
    result_df['gene_name'] = None
    result_df['gene_type'] = None
    result_df['gene_product'] = None
    
    # Define fields we expect in variant data
    expected_chrom_field = 'CHROM'
    expected_pos_field = 'POS'
    
    # Check and extract chromosome and position if not already present
    if expected_chrom_field not in result_df.columns or expected_pos_field not in result_df.columns:
        # Try to extract from Variant_ID if available (common format: scaffold_position_ref_alt)
        if 'Variant_ID' in result_df.columns:
            try:
                result_df[['CHROM', 'POS', 'REF', 'ALT']] = result_df['Variant_ID'].str.split('_', expand=True)
                result_df['POS'] = pd.to_numeric(result_df['POS'], errors='coerce')
            except Exception as e:
                print(f"Error extracting position from Variant_ID: {e}")
                logging.error(f"Error extracting position from Variant_ID: {e}")
                return variant_df
    
    # Count matches and misses for debugging
    matches = 0
    misses = 0
    chrom_misses = {}
    
    # Build a chromosome ID to scaffold mapping for more flexible lookups
    chrom_to_scaffold = {}
    scaffold_to_chrom = {}
    
    # First pass: build the chromosome/scaffold mapping from gene data
    for gene_id, gene_data in GENE_DATA.items():
        scaffold = gene_data['scaffold']
        if 'chromosome_id' in gene_data and gene_data['chromosome_id']:
            chrom_id = gene_data['chromosome_id']
            chrom_to_scaffold[chrom_id] = scaffold
            scaffold_to_chrom[scaffold] = chrom_id
    
    # Print some diagnostics about the mapping
    print(f"Built chromosome-to-scaffold mapping with {len(chrom_to_scaffold)} entries")
    print(f"Sample chromosome IDs: {list(chrom_to_scaffold.keys())[:5]}")
    print(f"Sample scaffold IDs: {list(scaffold_to_chrom.keys())[:5]}")
    
    # Map each variant to genes
    for idx, row in result_df.iterrows():
        # Skip rows with missing chromosome or position
        if expected_chrom_field not in row or expected_pos_field not in row:
            continue
        if pd.isna(row[expected_chrom_field]) or pd.isna(row[expected_pos_field]):
            continue
        
        chrom = str(row[expected_chrom_field])
        position = int(row[expected_pos_field])
        gene_found = False
        
        # Try direct match first
        if chrom in SCAFFOLD_GENES:
            for gene_id in SCAFFOLD_GENES[chrom]:
                gene_data = GENE_DATA[gene_id]
                
                # Check if position falls within gene coordinates
                if gene_data['start'] <= position <= gene_data['end']:
                    result_df.at[idx, 'in_gene'] = True
                    result_df.at[idx, 'gene_id'] = gene_id
                    
                    # Add gene name if available
                    if gene_data['erg_name']:
                        result_df.at[idx, 'gene_name'] = gene_data['erg_name']
                    elif gene_data['sc_gene_id']:
                        result_df.at[idx, 'gene_name'] = gene_data['sc_gene_id']
                    
                    # Set gene type based on presence in genes of interest
                    if gene_id in GENES_OF_INTEREST:
                        result_df.at[idx, 'gene_type'] = 'ergosterol'
                    else:
                        result_df.at[idx, 'gene_type'] = 'other'
                    
                    # Add gene product description if available
                    if gene_data['product']:
                        result_df.at[idx, 'gene_product'] = gene_data['product']
                    
                    gene_found = True
                    matches += 1
                    break
        
        # If not found by direct match, try using the mapping
        if not gene_found:
            # Try to map chromosome ID to scaffold
            if chrom in chrom_to_scaffold:
                scaffold = chrom_to_scaffold[chrom]
                if scaffold in SCAFFOLD_GENES:
                    for gene_id in SCAFFOLD_GENES[scaffold]:
                        gene_data = GENE_DATA[gene_id]
                        
                        # Check if position falls within gene coordinates
                        if gene_data['start'] <= position <= gene_data['end']:
                            result_df.at[idx, 'in_gene'] = True
                            result_df.at[idx, 'gene_id'] = gene_id
                            
                            # Add gene name if available
                            if gene_data['erg_name']:
                                result_df.at[idx, 'gene_name'] = gene_data['erg_name']
                            elif gene_data['sc_gene_id']:
                                result_df.at[idx, 'gene_name'] = gene_data['sc_gene_id']
                            
                            # Set gene type based on presence in genes of interest
                            if gene_id in GENES_OF_INTEREST:
                                result_df.at[idx, 'gene_type'] = 'ergosterol'
                            else:
                                result_df.at[idx, 'gene_type'] = 'other'
                            
                            # Add gene product description if available
                            if gene_data['product']:
                                result_df.at[idx, 'gene_product'] = gene_data['product']
                            
                            gene_found = True
                            matches += 1
                            break
            
            # Try to map scaffold to chromosome ID
            elif chrom in scaffold_to_chrom:
                chrom_id = scaffold_to_chrom[chrom]
                if chrom_id in SCAFFOLD_GENES:
                    for gene_id in SCAFFOLD_GENES[chrom_id]:
                        gene_data = GENE_DATA[gene_id]
                        
                        # Check if position falls within gene coordinates
                        if gene_data['start'] <= position <= gene_data['end']:
                            result_df.at[idx, 'in_gene'] = True
                            result_df.at[idx, 'gene_id'] = gene_id
                            
                            # Add gene name if available
                            if gene_data['erg_name']:
                                result_df.at[idx, 'gene_name'] = gene_data['erg_name']
                            elif gene_data['sc_gene_id']:
                                result_df.at[idx, 'gene_name'] = gene_data['sc_gene_id']
                            
                            # Set gene type based on presence in genes of interest
                            if gene_id in GENES_OF_INTEREST:
                                result_df.at[idx, 'gene_type'] = 'ergosterol'
                            else:
                                result_df.at[idx, 'gene_type'] = 'other'
                            
                            # Add gene product description if available
                            if gene_data['product']:
                                result_df.at[idx, 'gene_product'] = gene_data['product']
                            
                            gene_found = True
                            matches += 1
                            break
        
        # Track chromosome misses for diagnostics if no match found
        if not gene_found:
            if chrom not in chrom_misses:
                chrom_misses[chrom] = 0
            chrom_misses[chrom] += 1
            misses += 1
    
    # Log the mapping results
    in_gene_count = sum(result_df['in_gene'])
    ergosterol_count = sum(result_df['gene_type'] == 'ergosterol')
    
    print(f"Mapped {in_gene_count} out of {len(result_df)} variants to genes ({in_gene_count/len(result_df)*100:.1f}%)")
    print(f"Found {ergosterol_count} variants in ergosterol pathway genes")
    print(f"Chromosome matches: {matches}, misses: {misses}")
    
    # Print top missed chromosomes for diagnostics
    if misses > 0:
        print("Top 5 chromosomes with missing mappings:")
        for chrom, count in sorted(chrom_misses.items(), key=lambda x: x[1], reverse=True)[:5]:
            print(f"  {chrom}: {count} misses")
    
    logging.info(f"Mapped {in_gene_count} out of {len(result_df)} variants to genes ({in_gene_count/len(result_df)*100:.1f}%)")
    logging.info(f"Found {ergosterol_count} variants in ergosterol pathway genes")
    
    return result_df