#!/usr/bin/env python3
"""
Downstream AUG analysis with Kozak strength calculation
Uses txnames to select the correct annotation from V3
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
import json

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
        
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS transcripts (
                transcript_id TEXT PRIMARY KEY,
                gene_symbol TEXT,
                chromosome TEXT,
                strand TEXT,
                cds_regions TEXT,
                exon_regions TEXT,
                cds_length INTEGER,
                cds_start_position INTEGER,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
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
        
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS aug_analysis (
                variant_id TEXT,
                gene TEXT,
                transcript_id TEXT,
                chromosome TEXT,
                genomic_position INTEGER,
                ptc_codon_position INTEGER,
                ptc_annotation TEXT,
                nearest_inframe_aug_codon_position INTEGER,
                nearest_inframe_aug_distance_codons INTEGER,
                nearest_inframe_aug_distance_nt INTEGER,
                distance_start_to_inframe_aug_codons INTEGER,
                distance_start_to_inframe_aug_nt INTEGER,
                nearest_inframe_kozak_sequence TEXT,
                nearest_inframe_kozak_score INTEGER,
                nearest_inframe_kozak_strength TEXT,
                nearest_inframe_kozak_features TEXT,
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
                analysis_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        """)
        
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_transcripts_gene ON transcripts(gene_symbol)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_aug_analysis_gene ON aug_analysis(gene)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_aug_analysis_variant ON aug_analysis(variant_id)")
        
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
        
        cds_start_position = 0
        if exon_regions and cds_regions:
            first_cds_start = cds_regions[0][0]
            for exon_start, exon_end in exon_regions:
                if exon_end < first_cds_start:
                    cds_start_position += exon_end - exon_start + 1
                elif exon_start < first_cds_start <= exon_end:
                    cds_start_position += first_cds_start - exon_start
                    break
        
        cursor.execute("""
            INSERT OR REPLACE INTO transcripts 
            (transcript_id, gene_symbol, chromosome, strand, cds_regions, exon_regions, cds_length, cds_start_position)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (transcript_id, gene_symbol, chromosome, strand, cds_regions_json, exon_regions_json, cds_length, cds_start_position))
        
        self.conn.commit()
    
    def add_sequence(self, transcript_id: str, cds_sequence: str, cds_offset: int = 0):
        """Add CDS sequence to cache"""
        cursor = self.conn.cursor()
        cursor.execute("""
            INSERT OR REPLACE INTO sequence_cache 
            (transcript_id, cds_sequence, sequence_length, cds_offset)
            VALUES (?, ?, ?, ?)
        """, (transcript_id, cds_sequence, len(cds_sequence), cds_offset))
        self.conn.commit()
    
    def get_transcript(self, transcript_id: str) -> Optional[Dict]:
        """Get transcript structure from database"""
        cursor = self.conn.cursor()
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
        """Get CDS sequence from cache"""
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
        
        results_df = results_df.copy()
        results_df.columns = [col.replace(' ', '_').replace('-', '_') for col in results_df.columns]
        
        for col in results_df.select_dtypes(include=['object']).columns:
            results_df[col] = results_df[col].astype(str)
        
        try:
            results_df.to_sql('aug_analysis', self.conn, if_exists='replace', index=False)
            self.conn.commit()
            self.logger.info(f"Saved {len(results_df)} analysis results to database")
        except Exception as e:
            self.logger.error(f"Failed to save results to database: {e}")
            backup_file = "aug_analysis_backup.csv"
            results_df.to_csv(backup_file, index=False)
            self.logger.info(f"Saved backup results to {backup_file}")
    
    def close(self):
        """Close database connection"""
        self.conn.close()

class GTFParser:
    def __init__(self, gtf_file: str):
        self.logger = setup_logging()
        self.gtf_file = gtf_file
    
    def extract_transcripts_from_variants(self, variant_file: str) -> Set[str]:
        """Extract unique transcript IDs from txnames column"""
        self.logger.info(f"Extracting transcript IDs from {variant_file}")
        
        df = pd.read_csv(variant_file, sep=None, engine='python')
        transcript_ids = set()
        
        if 'txnames' in df.columns:
            self.logger.info("Using 'txnames' column for transcript IDs")
            for transcript_id in df['txnames'].dropna():
                transcript_ids.add(str(transcript_id))
                base_id = str(transcript_id).split('.')[0]
                transcript_ids.add(base_id)
        else:
            self.logger.error("Could not find 'txnames' column")
            return set()
        
        self.logger.info(f"Found {len(transcript_ids)} unique transcript IDs from {len(df):,} variants")
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
                    self.logger.info(f"Processed {line_num:,} GTF lines, found {processed_count} target transcripts")
                
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] not in ['CDS', 'exon']:
                    continue
                
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
                    
                    start, end = int(parts[3]), int(parts[4])
                    feature_type = parts[2]
                    
                    if feature_type == 'CDS':
                        transcript_data[transcript_id]['cds_regions'].append((start, end))
                    elif feature_type == 'exon':
                        transcript_data[transcript_id]['exon_regions'].append((start, end))
        
        self.logger.info(f"Adding {len(transcript_data)} transcripts to database")
        for transcript_id, data in transcript_data.items():
            data['cds_regions'] = sorted(data['cds_regions'])
            data['exon_regions'] = sorted(data['exon_regions'])
            db.add_transcript(
                transcript_id, data['gene_symbol'], data['chromosome'],
                data['strand'], data['cds_regions'], data['exon_regions']
            )
        
        self.logger.info(f"Successfully processed {len(transcript_data)} transcripts")

class TranscriptAnalyzer:
    def __init__(self, fasta_file: str, db: TranscriptDatabase):
        self.logger = setup_logging()
        self.fasta = pysam.FastaFile(fasta_file)
        self.db = db
    
    def parse_variant_id(self, variant_id: str) -> Dict:
        """Parse variant ID"""
        parts = variant_id.split('_')
        if len(parts) >= 4:
            return {
                'chromosome': parts[0],
                'position': int(parts[1]),
                'ref_allele': parts[2],
                'alt_allele': parts[3]
            }
        return {}
    
    def find_matching_annotation(self, txname: str, v3_annotations: str) -> Optional[str]:
        """Find the annotation in V3 that matches the transcript in txnames
        
        Args:
            txname: Transcript ID from txnames (e.g., ENST00000370128.8)
            v3_annotations: Comma-separated annotations from V3
        
        Returns:
            The matching annotation or None
        """
        if pd.isna(v3_annotations):
            return None
        
        # Split V3 by comma to get individual annotations
        annotations = [a.strip() for a in str(v3_annotations).split(',')]
        
        # Extract base transcript ID (without version)
        txname_base = txname.split('.')[0]
        
        # Look for matching transcript in each annotation
        for annotation in annotations:
            # Annotation format: GENE:TRANSCRIPT:exon:c.change:p.change
            if ':' in annotation:
                parts = annotation.split(':')
                if len(parts) >= 2:
                    annot_transcript = parts[1]
                    annot_transcript_base = annot_transcript.split('.')[0]
                    
                    # Match either exact or base transcript ID
                    if annot_transcript == txname or annot_transcript_base == txname_base:
                        return annotation
        
        return None
    
    def extract_protein_position_from_annotation(self, annotation: str) -> Optional[int]:
        """Extract protein position from annotation like: GENE:ENST:exon:c.change:p.Q28X"""
        if not annotation:
            return None
        
        annotation = annotation.rstrip(',')
        match = re.search(r'p\.[A-Z*](\d+)[X*]', annotation)
        return int(match.group(1)) if match else None
    
    def extract_cds_sequence(self, transcript_id: str) -> Optional[Tuple[str, int]]:
        """Extract CDS sequence using database info"""
        cached = self.db.get_sequence(transcript_id)
        if cached:
            return cached
        
        transcript_info = self.db.get_transcript(transcript_id)
        if not transcript_info:
            return None
        
        chromosome = transcript_info['chromosome']
        strand = transcript_info['strand']
        exon_regions = transcript_info.get('exon_regions', [])
        cds_regions = transcript_info['cds_regions']
        cds_start_pos = transcript_info.get('cds_start_position', 0)
        
        if exon_regions:
            full_transcript = ""
            for start, end in exon_regions:
                region_seq = self.fasta.fetch(chromosome, start - 1, end)
                full_transcript += region_seq
            
            if strand == '-':
                seq_obj = Seq(full_transcript)
                full_transcript = str(seq_obj.reverse_complement())
            
            full_transcript = full_transcript.upper()
            cds_length = sum(end - start + 1 for start, end in cds_regions)
            cds_end_pos = cds_start_pos + cds_length
            context_start = max(0, cds_start_pos - 15)
            cds_with_context = full_transcript[context_start:cds_end_pos]
            offset = cds_start_pos - context_start
            
            self.db.add_sequence(transcript_id, cds_with_context, offset)
            return (cds_with_context, offset)
        else:
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
        aug_nt_pos = (aug_position - 1) * 3
        
        if aug_nt_pos < 0 or aug_nt_pos + 3 > len(sequence):
            return {
                'kozak_sequence': '',
                'kozak_score': 0,
                'kozak_strength': 'insufficient_sequence',
                'kozak_features': 'none'
            }
        
        kozak_start = max(0, aug_nt_pos - 9)
        kozak_end = min(len(sequence), aug_nt_pos + 7)
        kozak_sequence = sequence[kozak_start:kozak_end]
        
        if len(kozak_sequence) < 3:
            return {
                'kozak_sequence': kozak_sequence,
                'kozak_score': 0,
                'kozak_strength': 'insufficient_sequence',
                'kozak_features': 'none'
            }
        
        kozak_score = 0
        kozak_features = []
        
        if len(kozak_sequence) >= 10:
            pos_minus3_idx = (aug_nt_pos - 3) - kozak_start
            if 0 <= pos_minus3_idx < len(kozak_sequence):
                pos_minus3 = kozak_sequence[pos_minus3_idx]
                if pos_minus3 == 'A':
                    kozak_score += 3
                    kozak_features.append('-3A')
                elif pos_minus3 == 'G':
                    kozak_score += 2
                    kozak_features.append('-3G')
            
            pos_plus4_idx = (aug_nt_pos + 4) - kozak_start
            if 0 <= pos_plus4_idx < len(kozak_sequence):
                pos_plus4 = kozak_sequence[pos_plus4_idx]
                if pos_plus4 == 'G':
                    kozak_score += 3
                    kozak_features.append('+4G')
            
            pos_minus6_idx = (aug_nt_pos - 6) - kozak_start
            if 0 <= pos_minus6_idx < len(kozak_sequence):
                pos_minus6 = kozak_sequence[pos_minus6_idx]
                if pos_minus6 == 'A':
                    kozak_score += 1
                    kozak_features.append('-6A')
            
            pos_plus1_idx = (aug_nt_pos + 3) - kozak_start
            if 0 <= pos_plus1_idx < len(kozak_sequence):
                pos_plus1 = kozak_sequence[pos_plus1_idx]
                if pos_plus1 == 'C':
                    kozak_score += 1
                    kozak_features.append('+1C')
        
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
        """Extract Kozak sequence around the original start codon"""
        if cds_offset > 0:
            kozak_start = max(0, cds_offset - 9)
            kozak_end = min(len(sequence), cds_offset + 7)
            kozak_sequence = sequence[kozak_start:kozak_end]
            
            if len(kozak_sequence) < 3:
                return {
                    'kozak_sequence': '',
                    'kozak_score': 0,
                    'kozak_strength': 'insufficient_sequence',
                    'kozak_features': 'none'
                }
            
            kozak_score = 0
            kozak_features = []
            
            pos_minus3_idx = cds_offset - 3 - kozak_start
            if 0 <= pos_minus3_idx < len(kozak_sequence):
                pos_minus3 = kozak_sequence[pos_minus3_idx]
                if pos_minus3 == 'A':
                    kozak_score += 3
                    kozak_features.append('-3A')
                elif pos_minus3 == 'G':
                    kozak_score += 2
                    kozak_features.append('-3G')
            
            pos_plus4_idx = cds_offset + 4 - kozak_start
            if 0 <= pos_plus4_idx < len(kozak_sequence):
                pos_plus4 = kozak_sequence[pos_plus4_idx]
                if pos_plus4 == 'G':
                    kozak_score += 3
                    kozak_features.append('+4G')
            
            pos_minus6_idx = cds_offset - 6 - kozak_start
            if 0 <= pos_minus6_idx < len(kozak_sequence):
                pos_minus6 = kozak_sequence[pos_minus6_idx]
                if pos_minus6 == 'A':
                    kozak_score += 1
                    kozak_features.append('-6A')
            
            pos_plus1_idx = cds_offset + 3 - kozak_start
            if 0 <= pos_plus1_idx < len(kozak_sequence):
                pos_plus1 = kozak_sequence[pos_plus1_idx]
                if pos_plus1 == 'C':
                    kozak_score += 1
                    kozak_features.append('+1C')
            
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
            return self.extract_kozak_sequence_at_position(sequence, 1)
    
    def calculate_distances_to_start(self, ptc_codon_position: int, sequence_length: int) -> Dict:
        """Calculate distances from PTC to start and end of transcript"""
        distance_from_start_codons = ptc_codon_position - 1
        distance_from_start_nt = distance_from_start_codons * 3
        
        total_codons = sequence_length // 3
        distance_to_end_codons = total_codons - ptc_codon_position
        distance_to_end_nt = distance_to_end_codons * 3
        
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
        distance_start_to_aug_codons = downstream_aug_position - 1
        distance_start_to_aug_nt = distance_start_to_aug_codons * 3
        
        return {
            'distance_start_to_aug_codons': distance_start_to_aug_codons,
            'distance_start_to_aug_nt': distance_start_to_aug_nt
        }
    
    def find_all_downstream_augs(self, sequence: str, ptc_codon_position: int) -> Dict:
        """Find downstream AUG codons in all reading frames"""
        ptc_nucleotide_pos = (ptc_codon_position - 1) * 3
        search_start = ptc_nucleotide_pos + 3
        
        augs_by_frame = {0: [], 1: [], 2: []}
        
        for i in range(search_start, len(sequence) - 2):
            if sequence[i:i+3] == 'ATG':
                frame = i % 3
                
                if frame == 0:
                    aug_codon_position = (i // 3) + 1
                    distance_codons = aug_codon_position - ptc_codon_position
                    distance_nucleotides = distance_codons * 3
                    augs_by_frame[0].append((aug_codon_position, distance_codons, distance_nucleotides, i))
                else:
                    distance_nucleotides = i - ptc_nucleotide_pos
                    augs_by_frame[frame].append((None, None, distance_nucleotides, i))
        
        nearest_augs = {}
        for frame in [0, 1, 2]:
            nearest_augs[frame] = augs_by_frame[frame][0] if augs_by_frame[frame] else None
        
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
        """Process variants and find downstream AUGs with Kozak analysis"""
        self.logger.info(f"Loading variants from {input_file}")
        
        with open(input_file, 'r') as f:
            first_line = f.readline()
            separator = ',' if (',' in first_line and first_line.count(',') > first_line.count('\t')) else '\t'
        
        df = pd.read_csv(input_file, sep=separator)
        self.logger.info(f"Loaded {len(df)} variants")
        
        required_cols = ['variantID', 'txnames', 'V3']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            self.logger.error(f"Missing required columns: {missing_cols}")
            self.logger.info(f"Available columns: {list(df.columns)}")
            return pd.DataFrame()
        
        results = []
        cache_hits = 0
        db_hits = 0
        failures = 0
        no_matching_annotation = 0
        
        for idx, row in df.iterrows():
            if idx % 500 == 0:
                self.logger.info(f"Processing variant {idx + 1}/{len(df)} (Cache: {cache_hits}, DB: {db_hits}, Failed: {failures}, No match: {no_matching_annotation})")
            
            try:
                variant_id = row['variantID']
                gene = (row.get('gene') or row.get('GENE_ID') or row.get('Gene') or 
                       row.get('gene_id') or row.get('gene_name') or 'Unknown')
                transcript_id = str(row['txnames']).strip()
                v3_annotations = row['V3']
                
                # Find the matching annotation in V3 for this transcript
                matched_annotation = self.find_matching_annotation(transcript_id, v3_annotations)
                
                if not matched_annotation:
                    if idx < 5:
                        self.logger.warning(f"No matching annotation found for transcript {transcript_id} in V3: {v3_annotations}")
                    no_matching_annotation += 1
                    failures += 1
                    continue
                
                if idx < 5:
                    self.logger.info(f"Processing: {variant_id}")
                    self.logger.info(f"  Transcript: {transcript_id}")
                    self.logger.info(f"  Matched annotation: {matched_annotation}")
                
                variant_info = self.parse_variant_id(variant_id)
                if not variant_info:
                    if idx < 5:
                        self.logger.warning(f"Could not parse variant ID: {variant_id}")
                    failures += 1
                    continue
                
                chrom = variant_info['chromosome']
                pos = variant_info['position']
                
                # Extract PTC position from the matched annotation
                ptc_position = self.extract_protein_position_from_annotation(matched_annotation)
                if not ptc_position:
                    if idx < 5:
                        self.logger.warning(f"Could not extract protein position from: {matched_annotation}")
                    failures += 1
                    continue
                
                if self.db.get_sequence(transcript_id):
                    cache_hits += 1
                else:
                    db_hits += 1
                
                seq_result = self.extract_cds_sequence(transcript_id)
                if not seq_result:
                    failures += 1
                    continue
                
                sequence, cds_offset = seq_result
                
                # Extract original start codon Kozak
                original_kozak_info = self.extract_original_start_kozak(sequence, cds_offset)
                
                # Work with CDS portion for downstream analysis
                cds_only = sequence[cds_offset:] if cds_offset > 0 else sequence
                
                # Calculate distances to start and end
                distance_info = self.calculate_distances_to_start(ptc_position, len(cds_only))
                
                # Find all downstream AUGs
                all_augs = self.find_all_downstream_augs(cds_only, ptc_position)
                
                # Process in-frame AUG
                if all_augs['in_frame_augs']:
                    aug_codon_pos, dist_codons, dist_nt, aug_nt_pos = all_augs['in_frame_augs'][0]
                    
                    start_to_aug_info = self.calculate_start_to_downstream_aug_distance(aug_codon_pos)
                    downstream_kozak_info = self.extract_kozak_sequence_at_position(cds_only, aug_codon_pos)
                    
                    nearest_plus1 = all_augs['nearest_plus1_frame']
                    nearest_plus2 = all_augs['nearest_plus2_frame']
                    
                    results.append({
                        'variant_id': variant_id,
                        'gene': gene,
                        'transcript_id': transcript_id,
                        'chromosome': chrom,
                        'genomic_position': pos,
                        'ptc_codon_position': ptc_position,
                        'ptc_annotation': matched_annotation,
                        
                        'nearest_inframe_aug_codon_position': aug_codon_pos,
                        'nearest_inframe_aug_distance_codons': dist_codons,
                        'nearest_inframe_aug_distance_nt': dist_nt,
                        'distance_start_to_inframe_aug_codons': start_to_aug_info['distance_start_to_aug_codons'],
                        'distance_start_to_inframe_aug_nt': start_to_aug_info['distance_start_to_aug_nt'],
                        'nearest_inframe_kozak_sequence': downstream_kozak_info['kozak_sequence'],
                        'nearest_inframe_kozak_score': downstream_kozak_info['kozak_score'],
                        'nearest_inframe_kozak_strength': downstream_kozak_info['kozak_strength'],
                        'nearest_inframe_kozak_features': downstream_kozak_info['kozak_features'],
                        
                        'nearest_plus1_frame_distance_nt': nearest_plus1[2] if nearest_plus1 else None,
                        'nearest_plus2_frame_distance_nt': nearest_plus2[2] if nearest_plus2 else None,
                        'has_plus1_frame_aug': nearest_plus1 is not None,
                        'has_plus2_frame_aug': nearest_plus2 is not None,
                        
                        'distance_from_start_codons': distance_info['distance_from_start_codons'],
                        'distance_from_start_nt': distance_info['distance_from_start_nt'],
                        'distance_to_end_codons': distance_info['distance_to_end_codons'],
                        'distance_to_end_nt': distance_info['distance_to_end_nt'],
                        'transcript_length_codons': distance_info['transcript_length_codons'],
                        'transcript_length_nt': distance_info['transcript_length_nt'],
                        'ptc_position_percent': distance_info['ptc_position_percent'],
                        
                        'original_kozak_sequence': original_kozak_info['kozak_sequence'],
                        'original_kozak_score': original_kozak_info['kozak_score'],
                        'original_kozak_strength': original_kozak_info['kozak_strength'],
                        'original_kozak_features': original_kozak_info['kozak_features']
                    })
                
                else:
                    # No in-frame AUG
                    nearest_plus1 = all_augs['nearest_plus1_frame']
                    nearest_plus2 = all_augs['nearest_plus2_frame']
                    
                    results.append({
                        'variant_id': variant_id,
                        'gene': gene,
                        'transcript_id': transcript_id,
                        'chromosome': chrom,
                        'genomic_position': pos,
                        'ptc_codon_position': ptc_position,
                        'ptc_annotation': matched_annotation,
                        
                        'nearest_inframe_aug_codon_position': None,
                        'nearest_inframe_aug_distance_codons': None,
                        'nearest_inframe_aug_distance_nt': None,
                        'distance_start_to_inframe_aug_codons': None,
                        'distance_start_to_inframe_aug_nt': None,
                        'nearest_inframe_kozak_sequence': None,
                        'nearest_inframe_kozak_score': None,
                        'nearest_inframe_kozak_strength': None,
                        'nearest_inframe_kozak_features': None,
                        
                        'nearest_plus1_frame_distance_nt': nearest_plus1[2] if nearest_plus1 else None,
                        'nearest_plus2_frame_distance_nt': nearest_plus2[2] if nearest_plus2 else None,
                        'has_plus1_frame_aug': nearest_plus1 is not None,
                        'has_plus2_frame_aug': nearest_plus2 is not None,
                        
                        'distance_from_start_codons': distance_info['distance_from_start_codons'],
                        'distance_from_start_nt': distance_info['distance_from_start_nt'],
                        'distance_to_end_codons': distance_info['distance_to_end_codons'],
                        'distance_to_end_nt': distance_info['distance_to_end_nt'],
                        'transcript_length_codons': distance_info['transcript_length_codons'],
                        'transcript_length_nt': distance_info['transcript_length_nt'],
                        'ptc_position_percent': distance_info['ptc_position_percent'],
                        
                        'original_kozak_sequence': original_kozak_info['kozak_sequence'],
                        'original_kozak_score': original_kozak_info['kozak_score'],
                        'original_kozak_strength': original_kozak_info['kozak_strength'],
                        'original_kozak_features': original_kozak_info['kozak_features']
                    })
                
            except Exception as e:
                self.logger.error(f"Error processing variant {idx}: {e}")
                failures += 1
                continue
        
        self.logger.info(f"Processing complete:")
        self.logger.info(f"  Successful: {len(results)}")
        self.logger.info(f"  Cache hits: {cache_hits}")
        self.logger.info(f"  DB lookups: {db_hits}")
        self.logger.info(f"  No matching annotation: {no_matching_annotation}")
        self.logger.info(f"  Total failures: {failures}")
        
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_file, sep='\t', index=False)
        self.db.save_aug_results(results_df)
        
        self.logger.info(f"Results saved to {output_file} and database")
        return results_df

def main():
    parser = argparse.ArgumentParser(description="Downstream AUG analysis with Kozak strength (txnames-matched)")
    parser.add_argument("--input", required=True, help="Input CSV file with variant data")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA file")
    parser.add_argument("--gtf", required=True, help="Gene annotation GTF file")
    parser.add_argument("--db", required=True, help="SQLite database file")
    parser.add_argument("--output", required=True, help="Output TSV file for results")
    parser.add_argument("--build-db", action="store_true", help="Build/update database from GTF")
    
    args = parser.parse_args()
    
    db = TranscriptDatabase(args.db)
    
    if args.build_db:
        gtf_parser = GTFParser(args.gtf)
        transcript_ids = gtf_parser.extract_transcripts_from_variants(args.input)
        gtf_parser.parse_gtf_for_transcripts(transcript_ids, db)
    
    analyzer = TranscriptAnalyzer(args.fasta, db)
    results = analyzer.process_variants(args.input, args.output)
    
    db.close()

if __name__ == "__main__":
    main()