#!/usr/bin/env python3
"""
Build SQLite database for transcripts and analyze downstream AUG codons.
Caches transcript sequences and structures for future use.
"""

import pandas as pd
import pysam
import sqlite3
import re
from Bio.Seq import Seq
import argparse
import logging
from typing import List, Tuple, Dict, Optional, Set
import gzip
import os
import json
from datetime import datetime

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    return logging.getLogger(__name__)

class TranscriptDatabase:
    def __init__(self, db_file: str):
        self.logger = setup_logging()
        self.db_file = db_file
        self.conn = sqlite3.connect(db_file)
        self.conn.row_factory = sqlite3.Row
        self._create_tables()
    
    def _create_tables(self):
        """Create tables for transcript data and sequence cache"""
        cursor = self.conn.cursor()
        
        # Transcript structure table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS transcripts (
                transcript_id TEXT PRIMARY KEY,
                gene_symbol TEXT,
                chromosome TEXT,
                strand TEXT,
                cds_regions TEXT,  -- JSON array of [start, end] pairs
                exon_regions TEXT,  -- JSON array of [start, end] pairs
                cds_length INTEGER,
                cds_start_position INTEGER,  -- Position where CDS starts in full transcript
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        # Sequence cache table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS sequence_cache (
                transcript_id TEXT PRIMARY KEY,
                cds_sequence TEXT,
                sequence_length INTEGER,
                cds_offset INTEGER DEFAULT 0,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (transcript_id) REFERENCES transcripts (transcript_id)
            )
        """)
        
        # Update database table schema for new columns
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS aug_analysis (
                variant_id TEXT,
                gene TEXT,
                transcript_id TEXT,
                chromosome TEXT,
                genomic_position INTEGER,
                ptc_codon_position INTEGER,
                ptc_annotation TEXT,
                downstream_aug_codon_position INTEGER,
                downstream_aug_frame TEXT,
                distance_codons INTEGER,
                distance_nucleotides INTEGER,
                distance_start_to_aug_codons INTEGER,
                distance_start_to_aug_nt INTEGER,
                nearest_plus1_frame_distance_nt INTEGER,
                nearest_plus2_frame_distance_nt INTEGER,
                has_plus1_frame_aug BOOLEAN,
                has_plus2_frame_aug BOOLEAN,
                distance_from_start_codons INTEGER,
                distance_from_start_nt INTEGER,
                distance_to_end_codons INTEGER,
                distance_to_end_nt INTEGER,
                transcript_length_codons INTEGER,
                transcript_length_nt INTEGER,
                ptc_position_percent REAL,
                original_kozak_sequence TEXT,
                original_kozak_score INTEGER,
                original_kozak_strength TEXT,
                original_kozak_features TEXT,
                downstream_kozak_sequence TEXT,
                downstream_kozak_score INTEGER,
                downstream_kozak_strength TEXT,
                downstream_kozak_features TEXT,
                analysis_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (variant_id, transcript_id, downstream_aug_codon_position)
            )
        """)
        
        # Create indexes for faster queries
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_transcripts_gene ON transcripts(gene_symbol)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_aug_analysis_gene ON aug_analysis(gene)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_aug_analysis_distance ON aug_analysis(distance_codons)")
        
        self.conn.commit()
        self.logger.info("Database tables created/verified")
    
    def add_transcript(self, transcript_id: str, gene_symbol: str, chromosome: str, 
                      strand: str, cds_regions: List[Tuple[int, int]], 
                      exon_regions: List[Tuple[int, int]] = None):
        """Add transcript structure to database"""
        cursor = self.conn.cursor()
        
        cds_regions_json = json.dumps(cds_regions)
        exon_regions_json = json.dumps(exon_regions) if exon_regions else json.dumps([])
        cds_length = sum(end - start + 1 for start, end in cds_regions)
        
        # Calculate where CDS starts in the full transcript
        cds_start_position = 0
        if exon_regions and cds_regions:
            # Find cumulative length up to first CDS region
            first_cds_start = cds_regions[0][0]
            for exon_start, exon_end in exon_regions:
                if exon_end < first_cds_start:
                    # This exon is entirely before CDS (5'UTR)
                    cds_start_position += exon_end - exon_start + 1
                elif exon_start < first_cds_start <= exon_end:
                    # CDS starts within this exon
                    cds_start_position += first_cds_start - exon_start
                    break
        
        cursor.execute("""
            INSERT OR REPLACE INTO transcripts 
            (transcript_id, gene_symbol, chromosome, strand, cds_regions, exon_regions, cds_length, cds_start_position)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (transcript_id, gene_symbol, chromosome, strand, cds_regions_json, exon_regions_json, cds_length, cds_start_position))
        
        self.conn.commit()
    
    def add_sequence(self, transcript_id: str, cds_sequence: str, cds_offset: int = 0):
        """Add CDS sequence to cache with optional 5'UTR offset"""
        cursor = self.conn.cursor()
        
        # Store offset as metadata in a new column
        cursor.execute("""
            INSERT OR REPLACE INTO sequence_cache 
            (transcript_id, cds_sequence, sequence_length, cds_offset)
            VALUES (?, ?, ?, ?)
        """, (transcript_id, cds_sequence, len(cds_sequence), cds_offset))
        
        self.conn.commit()
    
    def get_transcript(self, transcript_id: str) -> Optional[Dict]:
        """Get transcript structure from database"""
        cursor = self.conn.cursor()
        
        # Try exact match first, then without version
        base_id = transcript_id.split('.')[0]
        
        for tid in [transcript_id, base_id]:
            cursor.execute("SELECT * FROM transcripts WHERE transcript_id = ?", (tid,))
            row = cursor.fetchone()
            
            if row:
                return {
                    'transcript_id': row['transcript_id'],
                    'gene_symbol': row['gene_symbol'],
                    'chromosome': row['chromosome'],
                    'strand': row['strand'],
                    'cds_regions': json.loads(row['cds_regions']),
                    'exon_regions': json.loads(row['exon_regions']) if row['exon_regions'] else [],
                    'cds_length': row['cds_length'],
                    'cds_start_position': row['cds_start_position'] if 'cds_start_position' in row.keys() else 0
                }
        
        return None
    
    def get_sequence(self, transcript_id: str) -> Optional[Tuple[str, int]]:
        """Get CDS sequence from cache with offset
        
        Returns:
            Tuple of (sequence, cds_offset) or None
        """
        cursor = self.conn.cursor()
        
        base_id = transcript_id.split('.')[0]
        
        for tid in [transcript_id, base_id]:
            cursor.execute("SELECT cds_sequence, cds_offset FROM sequence_cache WHERE transcript_id = ?", (tid,))
            row = cursor.fetchone()
            
            if row:
                cds_offset = row['cds_offset'] if 'cds_offset' in row.keys() else 0
                return (row['cds_sequence'], cds_offset)
        
        return None
    
    def save_aug_results(self, results_df: pd.DataFrame):
        """Save AUG analysis results to database"""
        if len(results_df) == 0:
            self.logger.warning("No results to save - empty DataFrame")
            return
            
        # Clean up the DataFrame to avoid SQL issues
        results_df = results_df.copy()
        
        # Replace any problematic characters in column names
        results_df.columns = [col.replace(' ', '_').replace('-', '_') for col in results_df.columns]
        
        # Ensure all text columns are properly formatted
        for col in results_df.select_dtypes(include=['object']).columns:
            results_df[col] = results_df[col].astype(str)
        
        try:
            results_df.to_sql('aug_analysis', self.conn, if_exists='replace', index=False)
            self.conn.commit()
            self.logger.info(f"Saved {len(results_df)} analysis results to database")
        except Exception as e:
            self.logger.error(f"Failed to save results to database: {e}")
            # Try to save as CSV backup
            backup_file = "aug_analysis_backup.csv"
            results_df.to_csv(backup_file, index=False)
            self.logger.info(f"Saved backup results to {backup_file}")
    
    def get_database_stats(self):
        """Get statistics about the database"""
        cursor = self.conn.cursor()
        
        cursor.execute("SELECT COUNT(*) FROM transcripts")
        transcript_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM sequence_cache")
        sequence_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM aug_analysis")
        analysis_count = cursor.fetchone()[0]
        
        return {
            'transcripts': transcript_count,
            'cached_sequences': sequence_count,
            'analysis_results': analysis_count
        }
    
    def close(self):
        """Close database connection"""
        self.conn.close()

class GTFParser:
    def __init__(self, gtf_file: str):
        self.logger = setup_logging()
        self.gtf_file = gtf_file
    
    def extract_transcripts_from_variants(self, variant_file: str) -> Set[str]:
        """Extract unique transcript IDs from variant file"""
        self.logger.info(f"Extracting transcript IDs from {variant_file}")
        
        df = pd.read_csv(variant_file, sep=None, engine='python')  # Auto-detect separator
        transcript_ids = set()
        
        # Check if we have txnames column (direct transcript IDs)
        if 'txnames' in df.columns:
            self.logger.info("Using 'txnames' column for transcript IDs")
            for transcript_id in df['txnames'].dropna():
                transcript_ids.add(str(transcript_id))
                # Also add base ID without version
                base_id = str(transcript_id).split('.')[0]
                transcript_ids.add(base_id)
        
        # Fallback: try to extract from annotation columns (V3 or V8)
        elif 'V3' in df.columns:
            self.logger.info("Extracting transcript IDs from 'V3' annotation column")
            for annotation in df['V3'].dropna():
                parts = str(annotation).split(':')
                if len(parts) >= 2:
                    transcript_id = parts[1]
                    transcript_ids.add(transcript_id)
                    # Also add base ID without version
                    base_id = transcript_id.split('.')[0]
                    transcript_ids.add(base_id)
        elif 'V8' in df.columns:
            self.logger.info("Extracting transcript IDs from 'V8' annotation column")
            for annotation in df['V8'].dropna():
                parts = str(annotation).split(':')
                if len(parts) >= 2:
                    transcript_id = parts[1]
                    transcript_ids.add(transcript_id)
                    # Also add base ID without version
                    base_id = transcript_id.split('.')[0]
                    transcript_ids.add(base_id)
        else:
            self.logger.error("Could not find transcript ID column. Expected 'txnames', 'V3', or 'V8' column.")
            return set()
        
        self.logger.info(f"Found {len(transcript_ids)} unique transcript IDs")
        return transcript_ids
    
    def parse_gtf_for_transcripts(self, transcript_ids: Set[str], db: TranscriptDatabase):
        """Parse GTF file and extract data for specific transcripts"""
        self.logger.info(f"Parsing GTF file for {len(transcript_ids)} transcripts")
        
        open_func = gzip.open if self.gtf_file.endswith('.gz') else open
        
        transcript_data = {}
        processed_count = 0
        
        with open_func(self.gtf_file, 'rt') as f:
            for line_num, line in enumerate(f):
                if line.startswith('#'):
                    continue
                
                if line_num % 100000 == 0:
                    self.logger.info(f"Processed {line_num} GTF lines, found {processed_count} target transcripts")
                
                parts = line.strip().split('\t')
                # Extract both exon and CDS features
                if len(parts) < 9 or parts[2] not in ['CDS', 'exon']:
                    continue
                
                # Parse attributes
                attributes = {}
                for attr in parts[8].split(';'):
                    if ' "' in attr:
                        match = re.match(r'(\w+)\s+"([^"]+)"', attr.strip())
                        if match:
                            key, value = match.groups()
                            attributes[key] = value.strip('"')
                
                if 'transcript_id' not in attributes:
                    continue
                
                transcript_id = attributes['transcript_id']
                base_id = transcript_id.split('.')[0]
                
                # Check if this is one of our target transcripts
                if transcript_id in transcript_ids or base_id in transcript_ids:
                    if transcript_id not in transcript_data:
                        transcript_data[transcript_id] = {
                            'gene_symbol': attributes.get('gene_name', attributes.get('gene_id', 'unknown')),
                            'chromosome': parts[0],
                            'strand': parts[6],
                            'cds_regions': [],
                            'exon_regions': []
                        }
                        processed_count += 1
                    
                    start = int(parts[3])
                    end = int(parts[4])
                    feature_type = parts[2]
                    
                    if feature_type == 'CDS':
                        transcript_data[transcript_id]['cds_regions'].append((start, end))
                    elif feature_type == 'exon':
                        transcript_data[transcript_id]['exon_regions'].append((start, end))
        
        # Sort CDS regions and add to database
        self.logger.info(f"Adding {len(transcript_data)} transcripts to database")
        for transcript_id, data in transcript_data.items():
            data['cds_regions'] = sorted(data['cds_regions'])
            data['exon_regions'] = sorted(data['exon_regions'])
            db.add_transcript(
                transcript_id,
                data['gene_symbol'],
                data['chromosome'],
                data['strand'],
                data['cds_regions'],
                data['exon_regions']
            )
        
        self.logger.info(f"Successfully processed {len(transcript_data)} transcripts")
        return transcript_data

class TranscriptAnalyzer:
    def __init__(self, fasta_file: str, db: TranscriptDatabase):
        self.logger = setup_logging()
        self.fasta = pysam.FastaFile(fasta_file)
        self.db = db
    
    def parse_variant_id(self, variant_id: str) -> Dict:
        """Parse variant ID like: chr19_58353186_G_A"""
        parts = variant_id.split('_')
        if len(parts) >= 4:
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[2]
            alt = parts[3]
            return {
                'chromosome': chrom,
                'position': pos,
                'ref_allele': ref,
                'alt_allele': alt
            }
        return {}
    
    def extract_protein_position_from_annotation(self, annotation: str) -> int:
        """Extract protein position from annotation like: A1BG:ENST00000263100.7:exon3:c.C82T:p.Q28X"""
        # Remove trailing comma if present
        annotation = annotation.rstrip(',')
        
        # Look for the protein change part (p.Q28X)
        match = re.search(r'p\.[A-Z*](\d+)[X*]', annotation)
        if match:
            return int(match.group(1))
        return None
    
    def extract_cds_sequence(self, transcript_id: str) -> Optional[Tuple[str, int]]:
        """Extract CDS sequence using database info
        
        Returns:
            Tuple of (sequence, cds_offset) where cds_offset is the position where CDS starts
        """
        # Check cache first
        cached = self.db.get_sequence(transcript_id)
        if cached:
            return cached
        
        # Get transcript structure from database
        transcript_info = self.db.get_transcript(transcript_id)
        if not transcript_info:
            return None
        
        chromosome = transcript_info['chromosome']
        strand = transcript_info['strand']
        exon_regions = transcript_info.get('exon_regions', [])
        cds_regions = transcript_info['cds_regions']
        cds_start_pos = transcript_info.get('cds_start_position', 0)
        
        # If we have exon regions, extract full transcript then slice to CDS with context
        if exon_regions:
            # Extract full transcript sequence
            full_transcript = ""
            for start, end in exon_regions:
                region_seq = self.fasta.fetch(chromosome, start - 1, end)
                full_transcript += region_seq
            
            # Reverse complement if on negative strand
            if strand == '-':
                seq_obj = Seq(full_transcript)
                full_transcript = str(seq_obj.reverse_complement())
            
            full_transcript = full_transcript.upper()
            
            # Calculate CDS end position
            cds_length = sum(end - start + 1 for start, end in cds_regions)
            cds_end_pos = cds_start_pos + cds_length
            
            # Extract with upstream context (at least 15nt for Kozak)
            context_start = max(0, cds_start_pos - 15)
            cds_with_context = full_transcript[context_start:cds_end_pos]
            
            # Calculate offset (how much 5'UTR we included)
            offset = cds_start_pos - context_start
            
            # Store the sequence with offset
            self.db.add_sequence(transcript_id, cds_with_context, offset)
            
            return (cds_with_context, offset)
        else:
            # Fallback to CDS-only extraction (original method)
            cds_sequence = ""
            for start, end in cds_regions:
                region_seq = self.fasta.fetch(chromosome, start - 1, end)
                cds_sequence += region_seq
            
            if strand == '-':
                seq_obj = Seq(cds_sequence)
                cds_sequence = str(seq_obj.reverse_complement())
            
            cds_sequence = cds_sequence.upper()
            self.db.add_sequence(transcript_id, cds_sequence, 0)
            return (cds_sequence, 0)
    
    def extract_kozak_sequence_at_position(self, sequence: str, aug_position: int) -> Dict:
        """Extract Kozak sequence around a specific AUG position"""
        # Convert to 0-based nucleotide position
        aug_nt_pos = (aug_position - 1) * 3
        
        # Check if we have enough sequence (need at least ATG + some context)
        if aug_nt_pos < 0 or aug_nt_pos + 3 > len(sequence):
            return {
                'kozak_sequence': '',
                'kozak_score': 0,
                'kozak_strength': 'insufficient_sequence',
                'kozak_features': 'none'
            }
        
        # Kozak sequence: -9 to +4 around ATG (but work with what we have)
        kozak_start = max(0, aug_nt_pos - 9)
        kozak_end = min(len(sequence), aug_nt_pos + 7)  # ATG + 4 nucleotides
        kozak_sequence = sequence[kozak_start:kozak_end]
        
        # Need at least the ATG itself
        if len(kozak_sequence) < 3:
            return {
                'kozak_sequence': kozak_sequence,
                'kozak_score': 0,
                'kozak_strength': 'insufficient_sequence',
                'kozak_features': 'none'
            }
        
        # Calculate Kozak strength
        kozak_score = 0
        kozak_features = []
        
        if len(kozak_sequence) >= 10:
            # Position -3 relative to ATG start
            pos_minus3_idx = (aug_nt_pos - 3) - kozak_start
            if 0 <= pos_minus3_idx < len(kozak_sequence):
                pos_minus3 = kozak_sequence[pos_minus3_idx]
                if pos_minus3 == 'A':
                    kozak_score += 3
                    kozak_features.append('-3A')
                elif pos_minus3 == 'G':
                    kozak_score += 2
                    kozak_features.append('-3G')
            
            # Position +4 after ATG
            pos_plus4_idx = (aug_nt_pos + 4) - kozak_start
            if 0 <= pos_plus4_idx < len(kozak_sequence):
                pos_plus4 = kozak_sequence[pos_plus4_idx]
                if pos_plus4 == 'G':
                    kozak_score += 3
                    kozak_features.append('+4G')
            
            # Position -6
            pos_minus6_idx = (aug_nt_pos - 6) - kozak_start
            if 0 <= pos_minus6_idx < len(kozak_sequence):
                pos_minus6 = kozak_sequence[pos_minus6_idx]
                if pos_minus6 == 'A':
                    kozak_score += 1
                    kozak_features.append('-6A')
            
            # Position +1 (first position after ATG)
            pos_plus1_idx = (aug_nt_pos + 3) - kozak_start
            if 0 <= pos_plus1_idx < len(kozak_sequence):
                pos_plus1 = kozak_sequence[pos_plus1_idx]
                if pos_plus1 == 'C':
                    kozak_score += 1
                    kozak_features.append('+1C')
        
        # Classify strength
        if kozak_score >= 5:
            strength = 'strong'
        elif kozak_score >= 3:
            strength = 'moderate' 
        elif kozak_score >= 1:
            strength = 'weak'
        else:
            strength = 'very_weak'
        
        return {
            'kozak_sequence': kozak_sequence,
            'kozak_score': kozak_score,
            'kozak_strength': strength,
            'kozak_features': ','.join(kozak_features) if kozak_features else 'none'
        }
    
    def extract_original_start_kozak(self, sequence: str, cds_offset: int = 0) -> Dict:
        """Extract Kozak sequence around the original start codon
        
        Args:
            sequence: The transcript sequence (may include 5'UTR context)
            cds_offset: Number of nucleotides before the CDS start (5'UTR length)
        """
        # If we have 5'UTR context, the ATG is at position cds_offset
        # Otherwise, it's at position 0
        if cds_offset > 0:
            # Extract Kozak: -9 to +4 around ATG at position cds_offset
            kozak_start = max(0, cds_offset - 9)
            kozak_end = min(len(sequence), cds_offset + 7)  # ATG (3) + 4
            kozak_sequence = sequence[kozak_start:kozak_end]
            
            if len(kozak_sequence) < 3:
                return {
                    'kozak_sequence': '',
                    'kozak_score': 0,
                    'kozak_strength': 'insufficient_sequence',
                    'kozak_features': 'none'
                }
            
            # Calculate Kozak strength
            kozak_score = 0
            kozak_features = []
            
            # Position -3 relative to ATG (should be A or G)
            pos_minus3_idx = cds_offset - 3 - kozak_start
            if 0 <= pos_minus3_idx < len(kozak_sequence):
                pos_minus3 = kozak_sequence[pos_minus3_idx]
                if pos_minus3 == 'A':
                    kozak_score += 3
                    kozak_features.append('-3A')
                elif pos_minus3 == 'G':
                    kozak_score += 2
                    kozak_features.append('-3G')
            
            # Position +4 after ATG (should be G)
            pos_plus4_idx = cds_offset + 4 - kozak_start
            if 0 <= pos_plus4_idx < len(kozak_sequence):
                pos_plus4 = kozak_sequence[pos_plus4_idx]
                if pos_plus4 == 'G':
                    kozak_score += 3
                    kozak_features.append('+4G')
            
            # Position -6
            pos_minus6_idx = cds_offset - 6 - kozak_start
            if 0 <= pos_minus6_idx < len(kozak_sequence):
                pos_minus6 = kozak_sequence[pos_minus6_idx]
                if pos_minus6 == 'A':
                    kozak_score += 1
                    kozak_features.append('-6A')
            
            # Position +1 (first position after ATG)
            pos_plus1_idx = cds_offset + 3 - kozak_start
            if 0 <= pos_plus1_idx < len(kozak_sequence):
                pos_plus1 = kozak_sequence[pos_plus1_idx]
                if pos_plus1 == 'C':
                    kozak_score += 1
                    kozak_features.append('+1C')
            
            # Classify strength
            if kozak_score >= 5:
                strength = 'strong'
            elif kozak_score >= 3:
                strength = 'moderate' 
            elif kozak_score >= 1:
                strength = 'weak'
            else:
                strength = 'very_weak'
            
            return {
                'kozak_sequence': kozak_sequence,
                'kozak_score': kozak_score,
                'kozak_strength': strength,
                'kozak_features': ','.join(kozak_features) if kozak_features else 'none'
            }
        else:
            # Fallback to position 1 method (no 5'UTR context)
            return self.extract_kozak_sequence_at_position(sequence, 1)
    
    def calculate_distances_to_start(self, ptc_codon_position: int, sequence_length: int) -> Dict:
        """Calculate distances from PTC to start and end of transcript"""
        # Distance from start codon (codon 1) to PTC
        distance_from_start_codons = ptc_codon_position - 1  # PTC position minus start position
        distance_from_start_nt = distance_from_start_codons * 3
        
        # Distance from PTC to end of transcript
        total_codons = sequence_length // 3
        distance_to_end_codons = total_codons - ptc_codon_position
        distance_to_end_nt = distance_to_end_codons * 3
        
        # Position as percentage of transcript
        ptc_position_percent = (ptc_codon_position / total_codons) * 100 if total_codons > 0 else 0
        
        return {
            'distance_from_start_codons': distance_from_start_codons,
            'distance_from_start_nt': distance_from_start_nt,
            'distance_to_end_codons': distance_to_end_codons,
            'distance_to_end_nt': distance_to_end_nt,
            'transcript_length_codons': total_codons,
            'transcript_length_nt': sequence_length,
            'ptc_position_percent': round(ptc_position_percent, 1)
        }
    
    def calculate_start_to_downstream_aug_distance(self, downstream_aug_position: int) -> Dict:
        """Calculate distance from original start codon to downstream AUG"""
        # Distance from start codon (position 1) to downstream AUG
        distance_start_to_aug_codons = downstream_aug_position - 1  # AUG position minus start position
        distance_start_to_aug_nt = distance_start_to_aug_codons * 3
        
        return {
            'distance_start_to_aug_codons': distance_start_to_aug_codons,
            'distance_start_to_aug_nt': distance_start_to_aug_nt
        }
    
    def find_all_downstream_augs(self, sequence: str, ptc_codon_position: int) -> Dict:
        """Find downstream AUG codons in all reading frames"""
        # Convert to 0-based indexing
        ptc_nucleotide_pos = (ptc_codon_position - 1) * 3
        
        # Start searching from the nucleotide after the PTC
        search_start = ptc_nucleotide_pos + 3
        
        # Find AUGs in all three reading frames
        augs_by_frame = {0: [], 1: [], 2: []}  # 0 = in-frame, 1 = +1 frame, 2 = +2 frame
        
        # Search through the rest of the sequence
        for i in range(search_start, len(sequence) - 2):
            if sequence[i:i+3] == 'ATG':
                # Determine reading frame relative to original start codon (position 0)
                # Original start is at nucleotide position 0, so frame = i % 3
                frame = i % 3
                
                # Calculate positions and distances
                aug_nt_position = i + 1  # Convert to 1-based
                
                if frame == 0:  # In-frame with original start codon
                    aug_codon_position = (i // 3) + 1
                    distance_codons = aug_codon_position - ptc_codon_position
                    distance_nucleotides = distance_codons * 3
                    augs_by_frame[0].append((aug_codon_position, distance_codons, distance_nucleotides, aug_nt_position))
                else:
                    # Out-of-frame - calculate distance in nucleotides from PTC
                    distance_nucleotides = i - ptc_nucleotide_pos
                    augs_by_frame[frame].append((None, None, distance_nucleotides, aug_nt_position))
        
        # Find the nearest AUG in each frame
        nearest_augs = {}
        for frame in [0, 1, 2]:
            if augs_by_frame[frame]:
                # Get the first (nearest) AUG in this frame
                nearest_augs[frame] = augs_by_frame[frame][0]
            else:
                nearest_augs[frame] = None
        
        return {
            'in_frame_augs': augs_by_frame[0],
            'plus1_frame_augs': augs_by_frame[1], 
            'plus2_frame_augs': augs_by_frame[2],
            'nearest_in_frame': nearest_augs[0],
            'nearest_plus1_frame': nearest_augs[1],
            'nearest_plus2_frame': nearest_augs[2],
            'has_in_frame_aug': len(augs_by_frame[0]) > 0,
            'has_any_downstream_aug': any(len(augs_by_frame[f]) > 0 for f in [0, 1, 2])
        }
    
    def process_variants(self, input_file: str, output_file: str):
        """Process variants and find downstream AUGs"""
        self.logger.info(f"Loading variants from {input_file}")
        
        # Read the data - auto-detect separator
        with open(input_file, 'r') as f:
            first_line = f.readline()
            if ',' in first_line and first_line.count(',') > first_line.count('\t'):
                separator = ','
            else:
                separator = '\t'
        
        df = pd.read_csv(input_file, sep=separator)
        self.logger.info(f"Loaded {len(df)} variants")
        
        # Check for required columns
        required_cols = ['variantID', 'txnames', 'gene']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            self.logger.error(f"Missing required columns: {missing_cols}")
            self.logger.info(f"Available columns: {list(df.columns)}")
            return pd.DataFrame()
        
        results = []
        cache_hits = 0
        db_hits = 0
        failures = 0
        
        for idx, row in df.iterrows():
            if idx % 500 == 0:
                self.logger.info(f"Processing variant {idx + 1}/{len(df)} (Cache: {cache_hits}, DB: {db_hits}, Failed: {failures})")
            
            try:
                # Get basic variant info from direct columns
                variant_id = row['variantID']
                gene = row['gene']
                transcript_id = row['txnames']
                
                # Parse variant ID for genomic coordinates
                variant_info = self.parse_variant_id(variant_id)
                if not variant_info:
                    if idx < 5:
                        self.logger.warning(f"Could not parse variant ID: {variant_id}")
                    failures += 1
                    continue
                
                chrom = variant_info['chromosome']
                pos = variant_info['position']
                
                # Get protein position from matched_annotation column
                annotation = None
                if 'matched_annotation' in row.index and pd.notna(row['matched_annotation']):
                    annotation = str(row['matched_annotation']).strip()
                elif 'V3' in row.index and pd.notna(row['V3']):
                    annotation = str(row['V3']).strip()
                elif 'V8' in row.index and pd.notna(row['V8']):
                    annotation = str(row['V8']).strip()
                else:
                    if idx < 5:
                        self.logger.warning(f"No annotation found for variant {idx}")
                    failures += 1
                    continue
                
                if idx < 5:
                    self.logger.info(f"Processing: {variant_id}, transcript: {transcript_id}, annotation: {annotation}")
                
                # Extract PTC position from annotation
                ptc_position = self.extract_protein_position_from_annotation(annotation)
                if not ptc_position:
                    if idx < 5:
                        self.logger.warning(f"Could not extract protein position from: {annotation}")
                    failures += 1
                    continue
                
                # Check if sequence is cached
                if self.db.get_sequence(transcript_id):
                    cache_hits += 1
                else:
                    db_hits += 1
                
                # Get CDS sequence with offset
                seq_result = self.extract_cds_sequence(transcript_id)
                if not seq_result:
                    failures += 1
                    continue
                
                sequence, cds_offset = seq_result
                
                # Extract original start codon Kozak sequence (with proper offset)
                original_kozak_info = self.extract_original_start_kozak(sequence, cds_offset)
                
                # For downstream AUG analysis, we need to work with the CDS portion
                # Adjust positions to account for the 5'UTR offset
                if cds_offset > 0:
                    # Extract just the CDS for downstream analysis
                    cds_only = sequence[cds_offset:]
                else:
                    cds_only = sequence
                
                # Calculate distances to start and end (using CDS length)
                distance_info = self.calculate_distances_to_start(ptc_position, len(cds_only))
                
                # Find all downstream AUGs (in-frame and out-of-frame) - use CDS only
                all_augs = self.find_all_downstream_augs(cds_only, ptc_position)
                
                # Process ONLY the nearest in-frame AUG (if any)
                if all_augs['in_frame_augs']:
                    # Get only the first (nearest) in-frame AUG
                    aug_codon_pos, dist_codons, dist_nt, aug_nt_pos = all_augs['in_frame_augs'][0]
                    
                    # Calculate distance from start to this downstream AUG
                    start_to_aug_info = self.calculate_start_to_downstream_aug_distance(aug_codon_pos)
                    
                    # Extract Kozak sequence for this downstream AUG (use CDS only)
                    downstream_kozak_info = self.extract_kozak_sequence_at_position(cds_only, aug_codon_pos)
                    
                    # Get info about nearest AUGs in other frames
                    nearest_plus1 = all_augs['nearest_plus1_frame']
                    nearest_plus2 = all_augs['nearest_plus2_frame']
                    
                    results.append({
                        'variant_id': variant_id,
                        'gene': gene,
                        'transcript_id': transcript_id,
                        'chromosome': chrom,
                        'genomic_position': pos,
                        'ptc_codon_position': ptc_position,
                        'ptc_annotation': annotation,
                        
                        # Nearest in-frame downstream AUG
                        'nearest_inframe_aug_codon_position': aug_codon_pos,
                        'nearest_inframe_aug_distance_codons': dist_codons,
                        'nearest_inframe_aug_distance_nt': dist_nt,
                        'distance_start_to_inframe_aug_codons': start_to_aug_info['distance_start_to_aug_codons'],
                        'distance_start_to_inframe_aug_nt': start_to_aug_info['distance_start_to_aug_nt'],
                        'nearest_inframe_kozak_sequence': downstream_kozak_info['kozak_sequence'],
                        'nearest_inframe_kozak_score': downstream_kozak_info['kozak_score'],
                        'nearest_inframe_kozak_strength': downstream_kozak_info['kozak_strength'],
                        'nearest_inframe_kozak_features': downstream_kozak_info['kozak_features'],
                        
                        # Nearest out-of-frame AUG distances
                        'nearest_plus1_frame_distance_nt': nearest_plus1[2] if nearest_plus1 else None,
                        'nearest_plus2_frame_distance_nt': nearest_plus2[2] if nearest_plus2 else None,
                        'has_plus1_frame_aug': nearest_plus1 is not None,
                        'has_plus2_frame_aug': nearest_plus2 is not None,
                        
                        # PTC context in transcript
                        'distance_from_start_codons': distance_info['distance_from_start_codons'],
                        'distance_from_start_nt': distance_info['distance_from_start_nt'],
                        'distance_to_end_codons': distance_info['distance_to_end_codons'],
                        'distance_to_end_nt': distance_info['distance_to_end_nt'],
                        'transcript_length_codons': distance_info['transcript_length_codons'],
                        'transcript_length_nt': distance_info['transcript_length_nt'],
                        'ptc_position_percent': distance_info['ptc_position_percent'],
                        
                        # Original start codon Kozak
                        'original_kozak_sequence': original_kozak_info['kozak_sequence'],
                        'original_kozak_score': original_kozak_info['kozak_score'],
                        'original_kozak_strength': original_kozak_info['kozak_strength'],
                        'original_kozak_features': original_kozak_info['kozak_features']
                    })
                
                # If no in-frame AUGs, record the nearest out-of-frame AUG
                elif all_augs['has_any_downstream_aug']:
                    # Find the nearest AUG in any out-of-frame
                    nearest_aug = None
                    nearest_frame = None
                    min_distance = float('inf')
                    
                    for frame_name, frame_data in [('plus1_frame', all_augs['nearest_plus1_frame']), 
                                                  ('plus2_frame', all_augs['nearest_plus2_frame'])]:
                        if frame_data and frame_data[2] < min_distance:
                            nearest_aug = frame_data
                            nearest_frame = frame_name
                            min_distance = frame_data[2]
                    
                    if nearest_aug:
                        # Get info about all frame AUGs
                        nearest_plus1 = all_augs['nearest_plus1_frame']
                        nearest_plus2 = all_augs['nearest_plus2_frame']
                        
                        results.append({
                            'variant_id': variant_id,
                            'gene': gene,
                            'transcript_id': transcript_id,
                            'chromosome': chrom,
                            'genomic_position': pos,
                            'ptc_codon_position': ptc_position,
                            'ptc_annotation': annotation,
                            
                            # No in-frame AUG available
                            'nearest_inframe_aug_codon_position': None,
                            'nearest_inframe_aug_distance_codons': None,
                            'nearest_inframe_aug_distance_nt': None,
                            'distance_start_to_inframe_aug_codons': None,
                            'distance_start_to_inframe_aug_nt': None,
                            'nearest_inframe_kozak_sequence': None,
                            'nearest_inframe_kozak_score': None,
                            'nearest_inframe_kozak_strength': None,
                            'nearest_inframe_kozak_features': None,
                            
                            # Nearest out-of-frame AUG info
                            'nearest_outframe_aug_frame': nearest_frame,
                            'nearest_outframe_aug_distance_nt': nearest_aug[2],
                            'nearest_plus1_frame_distance_nt': nearest_plus1[2] if nearest_plus1 else None,
                            'nearest_plus2_frame_distance_nt': nearest_plus2[2] if nearest_plus2 else None,
                            'has_plus1_frame_aug': nearest_plus1 is not None,
                            'has_plus2_frame_aug': nearest_plus2 is not None,
                            
                            # PTC context in transcript
                            'distance_from_start_codons': distance_info['distance_from_start_codons'],
                            'distance_from_start_nt': distance_info['distance_from_start_nt'],
                            'distance_to_end_codons': distance_info['distance_to_end_codons'],
                            'distance_to_end_nt': distance_info['distance_to_end_nt'],
                            'transcript_length_codons': distance_info['transcript_length_codons'],
                            'transcript_length_nt': distance_info['transcript_length_nt'],
                            'ptc_position_percent': distance_info['ptc_position_percent'],
                            
                            # Original start codon Kozak
                            'original_kozak_sequence': original_kozak_info['kozak_sequence'],
                            'original_kozak_score': original_kozak_info['kozak_score'],
                            'original_kozak_strength': original_kozak_info['kozak_strength'],
                            'original_kozak_features': original_kozak_info['kozak_features']
                        })
                
                else:
                    # No downstream AUGs in any frame
                    results.append({
                        'variant_id': variant_id,
                        'gene': gene,
                        'transcript_id': transcript_id,
                        'chromosome': chrom,
                        'genomic_position': pos,
                        'ptc_codon_position': ptc_position,
                        'ptc_annotation': annotation,
                        
                        # No downstream AUGs
                        'nearest_inframe_aug_codon_position': None,
                        'nearest_inframe_aug_distance_codons': None,
                        'nearest_inframe_aug_distance_nt': None,
                        'distance_start_to_inframe_aug_codons': None,
                        'distance_start_to_inframe_aug_nt': None,
                        'nearest_inframe_kozak_sequence': None,
                        'nearest_inframe_kozak_score': None,
                        'nearest_inframe_kozak_strength': None,
                        'nearest_inframe_kozak_features': None,
                        
                        # No out-of-frame AUGs
                        'nearest_outframe_aug_frame': None,
                        'nearest_outframe_aug_distance_nt': None,
                        'nearest_plus1_frame_distance_nt': None,
                        'nearest_plus2_frame_distance_nt': None,
                        'has_plus1_frame_aug': False,
                        'has_plus2_frame_aug': False,
                        
                        # PTC context in transcript
                        'distance_from_start_codons': distance_info['distance_from_start_codons'],
                        'distance_from_start_nt': distance_info['distance_from_start_nt'],
                        'distance_to_end_codons': distance_info['distance_to_end_codons'],
                        'distance_to_end_nt': distance_info['distance_to_end_nt'],
                        'transcript_length_codons': distance_info['transcript_length_codons'],
                        'transcript_length_nt': distance_info['transcript_length_nt'],
                        'ptc_position_percent': distance_info['ptc_position_percent'],
                        
                        # Original start codon Kozak
                        'original_kozak_sequence': original_kozak_info['kozak_sequence'],
                        'original_kozak_score': original_kozak_info['kozak_score'],
                        'original_kozak_strength': original_kozak_info['kozak_strength'],
                        'original_kozak_features': original_kozak_info['kozak_features']
                    })
                
            except Exception as e:
                self.logger.error(f"Error processing variant {idx}: {e}")
                failures += 1
                continue
        
        self.logger.info(f"Processing complete: Cache hits: {cache_hits}, DB lookups: {db_hits}, Failures: {failures}")
        
        # Save results
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_file, sep='\t', index=False)
        
        # Also save to database for future reference
        self.db.save_aug_results(results_df)
        
        self.logger.info(f"Results saved to {output_file} and database")
        return results_df
    
    def analyze_results(self, results_df: pd.DataFrame):
        """Analyze results and show statistics"""
        print("\n" + "="*60)
        print("DOWNSTREAM AUG ANALYSIS RESULTS")
        print("="*60)
        
        total_variants = len(results_df)
        variants_with_inframe_augs = results_df['nearest_inframe_aug_codon_position'].notna().sum()
        variants_with_outframe_only = (
            (results_df['nearest_inframe_aug_codon_position'].isna()) & 
            (results_df.get('nearest_outframe_aug_distance_nt', pd.Series([None]*len(results_df))).notna())
        ).sum()
        variants_without_augs = (
            (results_df['nearest_inframe_aug_codon_position'].isna()) & 
            (results_df.get('nearest_outframe_aug_distance_nt', pd.Series([None]*len(results_df))).isna())
        ).sum()
        
        print(f"\nSummary:")
        print(f"  Total variants analyzed: {total_variants}")
        print(f"  Variants with in-frame downstream AUGs: {variants_with_inframe_augs} ({variants_with_inframe_augs/total_variants*100:.1f}%)")
        print(f"  Variants with only out-of-frame AUGs: {variants_with_outframe_only} ({variants_with_outframe_only/total_variants*100:.1f}%)")
        print(f"  Variants without any downstream AUGs: {variants_without_augs} ({variants_without_augs/total_variants*100:.1f}%)")
        
        # Database statistics
        db_stats = self.db.get_database_stats()
        print(f"\nDatabase Cache:")
        print(f"  Transcripts in database: {db_stats['transcripts']}")
        print(f"  Sequences cached: {db_stats['cached_sequences']}")
        print(f"  Analysis results stored: {db_stats['analysis_results']}")
        
        # Additional analysis for reading frames and position context
        if len(results_df) > 0:
            print(f"\nReading Frame Analysis:")
            inframe_count = results_df['nearest_inframe_aug_codon_position'].notna().sum()
            outframe_only = (
                (results_df['nearest_inframe_aug_codon_position'].isna()) & 
                (results_df.get('nearest_outframe_aug_distance_nt', pd.Series([None]*len(results_df))).notna())
            ).sum()
            none_count = (
                (results_df['nearest_inframe_aug_codon_position'].isna()) & 
                (results_df.get('nearest_outframe_aug_distance_nt', pd.Series([None]*len(results_df))).isna())
            ).sum()
            
            print(f"  Has in-frame AUG: {inframe_count} ({inframe_count/len(results_df)*100:.1f}%)")
            print(f"  Only out-of-frame AUGs: {outframe_only} ({outframe_only/len(results_df)*100:.1f}%)")
            print(f"  No downstream AUGs: {none_count} ({none_count/len(results_df)*100:.1f}%)")
            
            print(f"\nTranscript Position Analysis:")
            if 'ptc_position_percent' in results_df.columns:
                positions = results_df['ptc_position_percent'].dropna()
                print(f"  PTC position in transcript (%):")
                print(f"    Mean: {positions.mean():.1f}%")
                print(f"    Median: {positions.median():.1f}%")
                print(f"    Range: {positions.min():.1f}% - {positions.max():.1f}%")
            
            if 'distance_from_start_nt' in results_df.columns:
                start_distances = results_df['distance_from_start_nt'].dropna()
                print(f"  Distance from start codon:")
                print(f"    Mean: {start_distances.mean():.1f} nt")
                print(f"    Median: {start_distances.median():.1f} nt")
                
                # Start-proximal analysis (key finding from paper)
                within_200nt = (start_distances <= 200).sum()
                within_100nt = (start_distances <= 100).sum()
                print(f"    Within 200nt of start: {within_200nt} ({within_200nt/len(start_distances)*100:.1f}%)")
                print(f"    Within 100nt of start: {within_100nt} ({within_100nt/len(start_distances)*100:.1f}%)")
            
            # Kozak strength analysis
            if 'original_kozak_strength' in results_df.columns:
                kozak_counts = results_df['original_kozak_strength'].value_counts()
                print(f"\nOriginal Kozak Sequence Strength:")
                for strength, count in kozak_counts.items():
                    percentage = (count / len(results_df)) * 100
                    print(f"  {strength}: {count} ({percentage:.1f}%)")
            
            # Out-of-frame AUG analysis
            if 'has_plus1_frame_aug' in results_df.columns:
                plus1_count = results_df['has_plus1_frame_aug'].sum()
                plus2_count = results_df['has_plus2_frame_aug'].sum()
                print(f"\nOut-of-Frame AUG Availability:")
                print(f"  Has +1 frame AUG: {plus1_count} ({plus1_count/len(results_df)*100:.1f}%)")
                print(f"  Has +2 frame AUG: {plus2_count} ({plus2_count/len(results_df)*100:.1f}%)")
        
        # Distance analysis for variants with in-frame downstream AUGs
        aug_data = results_df[results_df['nearest_inframe_aug_codon_position'].notna()]
        if len(aug_data) > 0:
            print(f"\nDownstream In-Frame AUG Distance Statistics:")
            distances = aug_data['nearest_inframe_aug_distance_codons']
            print(f"  Variants with in-frame AUGs: {len(aug_data)}")
            print(f"  Distance range: {distances.min()} - {distances.max()} codons")
            print(f"  Mean distance: {distances.mean():.1f} codons ({distances.mean()*3:.1f} nucleotides)")
            print(f"  Median distance: {distances.median():.1f} codons ({distances.median()*3:.1f} nucleotides)")
            
            # Distance distribution
            print(f"\nDistance Distribution:")
            ranges = [(1, 10), (11, 50), (51, 100), (101, 200), (201, float('inf'))]
            for start, end in ranges:
                if end == float('inf'):
                    count = len(aug_data[distances >= start])
                    print(f"  {start}+ codons: {count} ({count/len(aug_data)*100:.1f}%)")
                else:
                    count = len(aug_data[(distances >= start) & (distances <= end)])
                    print(f"  {start}-{end} codons: {count} ({count/len(aug_data)*100:.1f}%)")
            
            # Kozak strength of downstream AUGs
            if 'nearest_inframe_kozak_strength' in aug_data.columns:
                downstream_kozak_counts = aug_data['nearest_inframe_kozak_strength'].value_counts()
                print(f"\nNearest In-Frame AUG Kozak Strength:")
                for strength, count in downstream_kozak_counts.items():
                    percentage = (count / len(aug_data)) * 100
                    print(f"  {strength}: {count} ({percentage:.1f}%)")

def main():
    parser = argparse.ArgumentParser(description="Build transcript database and analyze downstream AUGs")
    parser.add_argument("--input", required=True, help="Input CSV file with variant data")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA file")
    parser.add_argument("--gtf", required=True, help="Gene annotation GTF file")
    parser.add_argument("--db", required=True, help="SQLite database file (will be created/updated)")
    parser.add_argument("--output", required=True, help="Output TSV file for results")
    parser.add_argument("--build-db", action="store_true", help="Build/update database from GTF")
    parser.add_argument("--analyze", action="store_true", help="Perform summary analysis")
    
    args = parser.parse_args()
    
    # Initialize database
    db = TranscriptDatabase(args.db)
    
    # Build database if requested
    if args.build_db:
        gtf_parser = GTFParser(args.gtf)
        transcript_ids = gtf_parser.extract_transcripts_from_variants(args.input)
        gtf_parser.parse_gtf_for_transcripts(transcript_ids, db)
    
    # Run analysis
    analyzer = TranscriptAnalyzer(args.fasta, db)
    results = analyzer.process_variants(args.input, args.output)
    
    if args.analyze:
        analyzer.analyze_results(results)
    
    db.close()

if __name__ == "__main__":
    main()