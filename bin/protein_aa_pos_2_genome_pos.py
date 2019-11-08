#! /bin/env python3
# ------------------
# * CDS feature is used for CDS region for each exon
# * Stop coden is not included in CDS
# * frame is the coden position of first nt
# ------------------
####################

import os
import sys
import numpy as np
import pandas as pd
import argparse
from Bio import SeqIO
sys.path.append('/gpfs/share/home/1501111485/Code/pycas')
import pycas
from pycas.utils.gff import Gff, GffFile

####################

def amino_acid(x, information):
    infocol = ['name', 'triletter', 'uniletter', 'class', 'polarity', 'charge', 'codon']
    assert information in infocol, 'Wrong information: {0}. options should be in [{1}]'.format(
        information,
        ', '.join(infocol)
    )
    info = [
        ['Alanine', 'Ala', 'A', 'aliphatic', 'nonpolar', 'neutral', ['GCU', 'GCC', 'GCA', 'GCG']],
        ['Arginine', 'Arg', 'R', 'basic', 'basic polar', 'positive', ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG']],
        ['Asparagine', 'Asn', 'N', 'amide', 'polar', 'neutral', ['AAU', 'AAC']],
        ['Aspartic acid', 'Asp', 'D', 'acid', 'acidic polar', 'negative', ['GAU', 'GAC']],
        ['Cysteine', 'Cys', 'C', 'sulfur-containing', 'nonpolar', 'neutral', ['UGU', 'UGC']],
        ['Glutamic acid', 'Glu', 'E', 'acid', 'acidic polar', 'negative', ['GAA', 'GAG']],
        ['Glutamine', 'Gln', 'Q', 'amide', 'polar', 'neutral', ['CAA', 'CAG']],
        ['Glycine', 'Gly', 'G', 'aliphatic', 'nonpolar', 'neutral', ['GGU', 'GGC', 'GGA', 'GGG']],
        ['Histidine', 'His', 'H', 'basic aromatic', 'basic polar', 'positive(10%) neutral(90%)', ['CAU', 'CAC']],
        ['Isoleucine', 'Ile', 'I', 'aliphatic', 'nonpolar', 'neutral', ['AUU', 'AUC', 'AUA']],
        ['Leucine', 'Leu', 'L', 'aliphatic', 'nonpolar', 'neutral', ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG']],
        ['Lysine', 'Lys', 'K', 'basic', 'basic polar', 'positive', ['AAA', 'AAG']],
        ['Methionine', 'Met', 'M', 'sulfur-containing', 'nonpolar', 'neutral', ['AUG']],
        ['Phenylalanine', 'Phe', 'F', 'aromatic', 'nonpolar', 'neutral', ['UUU', 'UUC']],
        ['Proline', 'Pro', 'P', 'cyclic', 'nonpolar', 'neutral', ['CCU', 'CCC', 'CCA', 'CCG']],
        ['Serine', 'Ser', 'S', 'hydroxyl-containing', 'polar', 'neutral', ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC']],
        ['Threonine', 'Thr', 'T', 'hydroxyl-containing', 'polar', 'neutral', ['ACU', 'ACC', 'ACA', 'ACG']],
        ['Tryptophan', 'Trp', 'W', 'aromatic', 'nonpolar', 'neutral', ['UGG']],
        ['Tyrosine', 'Tyr', 'Y', 'aromatic', 'polar', 'neutral', ['UAU', 'UAC']],
        ['Valine', 'Val', 'V', 'aliphatic', 'nonpolar', 'neutral', ['GUU', 'GUC', 'GUA', 'GUG']],
        ['STOP', 'STP', '*', '', '', '', ['UAA', 'UGA', 'UAG']]
    ]
    infodict = dict(
        zip(infocol, [ [x[i] for x in info ] for i in range(len(infocol))])
    )
    for i in range(len(infocol)):
        if i == 6:
            for codons in infodict['codon']:
                if x in codons:
                    return infodict[information][infodict['codon'].index(codons)]
        else:
            if x in infodict[infocol[i]]:
                return infodict[information][infodict[infocol[i]].index(x)]
    raise ValueError('Wrong input x: {}'.format(x))


def complement(seq):
    # calculate the complement sequence string
    # input a string
    # return a string
    iupac = {
        "A": "T", "G": "C", "C": "G", "T": "A", "Y": "R", "R": "Y",
        "W": "W", "S": "S", "K": "M", "M": "K", "D": "H", "V": "B",
        "H": "D", "B": "V", "N": "N", "X": "X", "-": "-",
        "a": "t", "g": "c", "c": "g", "t": "a", "y": "r", "r": "y",
        "w": "w", "s": "s", "k": "m", "m": "k", "d": "h", "v": "b",
        "h": "d", "b": "v", "n": "n", "x": "x", "-": "-",
    }
    return ''.join(iupac[i] for i in seq)


def reverse_complement(seq):
    # calculate the reverse complement sequence string
    # input a string
    # return a string
    return complement(seq)[::-1]

####################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='protein_aa_pos_2_genome_pos')
    parser.add_argument(
        '-f', '--fasta', required=True,
        help='reference genome fasta file.'
    )
    parser.add_argument(
        '-g', '--gtf', required=True,
        help='gtf file for the reference genome.'
    )
    parser.add_argument(
        '-p', '--protein-pos',
        help='input tab separated file for the protein AA position, [Gene Name, Amino Acid (one character), Protein Position].'
    )
    parser.add_argument(
        '-o', '--output', help='output file for the AA pos on genome without splicing.'
    )
    parser.add_argument(
        '-s', '--splice-output', help='output file for the AA pos on genome with splicing.'
    )
    args = vars(parser.parse_args())

    ####################
    # __console__ = sys.stdout
    # if 'output' in args.keys():
    #     sys.stdout = open(args['output'], 'w')

    vseqname = [''.join(['chr', str(i)]) for i in list(range(1, 23)) + ['X', 'Y']]

    protein_pos = pd.read_table(args['protein_pos'], header=None)
    protein_pos.columns = ['uniprot_1', 'uniprot_2', 'gene_id', 'aa', 'aa_prot_pos']
    protein_pos = protein_pos.loc[np.logical_not(protein_pos.duplicated())]
    protein_pos.loc[:,'aa_cds_pos_3'] = protein_pos['aa_prot_pos'] * 3
    protein_pos.loc[:,'aa_cds_pos_2'] = protein_pos['aa_cds_pos_3'] - 1
    protein_pos.loc[:,'aa_cds_pos_1'] = protein_pos['aa_cds_pos_3'] - 2

    # ------------------
    # genome sequence
    fa = dict()
    for record in SeqIO.parse(args["fasta"], "fasta"):
        if record.id in vseqname:
            fa[record.id] = record
    # ------------------

    # ------------------
    # GFF

    gfffile = GffFile(
        args['gtf'], header=None
    )

    gff = pd.concat([gfffile.gff.gff, gfffile.gffattr], axis=1)
    gff.loc[:,'gene_id'] = gff.loc[:,'gene_id'].str.replace('"', '')
    gff.loc[:,'gene_name'] = gff.loc[:,'gene_name'].str.replace('"', '')
    gff.loc[:,'transcript_id'] = gff.loc[:,'transcript_id'].str.replace('"', '')
    gff = gff.loc[gff['seqname'].isin(vseqname).values]
    # ------------------

    # CDS information:
    #     seqname: name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
    #     source: name of the program that generated this feature, or the data source (database or project name)
    #     feature: feature type name, e.g. Gene, Variation, Similarity
    #     start: Start position of the feature, with sequence numbering starting at 1.
    #     end: End position of the feature, with sequence numbering starting at 1.
    #     score: A floating point value.
    #     strand: defined as + (forward) or - (reverse).
    #     frame: One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    #     attribute: A semicolon-separated list of tag-value pairs, providing additional information about each feature.
    #     gene_id: gene ID
    #     gene_name: gene name
    #     p_id: protein ID
    #     transcript_id: transcript ID
    #     tss_id: transcript starting site ID
    #     cds_no: CDS number in the transcript
    #     cds_len: CDS length
    #     cds_cumlen: CDS cumulated length, all length including the CDS and CDSs before
    #     cds_ntstart: the start length of nucleotide in the CDS, cds_cumlen - cds_len
    # Amino acid information:
    #     uniprot_1: uniprot accession number
    #     uniprot_2: uniprot ID
    #     aa: uni-letter of amino acid
    #     aa_prot_pos: amino acid position in the protein sequence
    #     aa_cds_pos_1: amino acid 1st nucleotide of codon position in the whole CDS, aa_prot_pos * 3 - 2
    #     aa_cds_pos_2: amino acid 2nd nucleotide of codon position in the whole CDS, aa_prot_pos * 3 - 1
    #     aa_cds_pos_3: amino acid 3rd nucleotide of codon position in the whole CDS, aa_prot_pos * 3
    #     aa_exon_pos_1: amino acid 1st nucleotide of codon position in the current CDS exon, aa_cds_pos_1 - cds_ntstart
    #     aa_exon_pos_2: amino acid 2nd nucleotide of codon position in the current CDS exon, aa_cds_pos_2 - cds_ntstart
    #     aa_exon_pos_3: amino acid 3rd nucleotide of codon position in the current CDS exon, aa_cds_pos_3 - cds_ntstart
    #     aa_genome_pos_start: amino acid codon genome posiiton start
    #     aa_genome_pos_end: amino acid codon genome posiiton end
    #     aa_genome_seq: amino acid codon genome sequence
    #     aa_genome_pos_1: amino acid 1st nucleotide of codon genome posiiton
    #     aa_genome_pos_2: amino acid 2nd nucleotide of codon genome posiiton
    #     aa_genome_pos_3: amino acid 3rd nucleotide of codon genome posiiton
    #     aa_genome_seq_1: amino acid 1st nucleotide of codon genome seq
    #     aa_genome_seq_2: amino acid 2nd nucleotide of codon genome seq
    #     aa_genome_seq_3: amino acid 3rd nucleotide of codon genome seq

    resultcol = [
        'seqname', 'source', 'feature', 'start', 'end',
        'score', 'strand', 'frame', 'attribute',
        'gene_id', 'gene_name', 'p_id', 'transcript_id',
        'tss_id', 'cds_no', 'cds_len', 'cds_cumlen', 'cds_ntstart',
        'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos',
        'aa_cds_pos_1', 'aa_cds_pos_2', 'aa_cds_pos_3',
        'aa_exon_pos_1', 'aa_exon_pos_2', 'aa_exon_pos_3',
        'aa_genome_pos_start', 'aa_genome_pos_end', 'aa_genome_seq',
        'aa_genome_pos_1', 'aa_genome_pos_2', 'aa_genome_pos_3',
        'aa_genome_seq_1', 'aa_genome_seq_2', 'aa_genome_seq_3'
    ]

    result = dict()
    result['nosplice'] = list()
    result['splice'] = list()

    for x in gff['gene_id'].str[0].unique():
        gffx = gff.loc[gff['gene_id'].str[0] == x]
        # ------------------
        # CDS
        # watson strand +
        cdsw = gffx.loc[
            np.logical_and(gffx.feature == 'CDS', gffx.strand=='+')
        ].copy()
        cdsw.loc[:,'cds_len'] = cdsw['end'] - cdsw['start'] + 1
        cdsw.sort_values(
            ['transcript_id', 'start'], ascending=True, inplace=True
        )
        cdsw.reset_index(drop=True, inplace=True)
        # CDS number in transcripts
        cdsw.loc[:, 'cds_no'] = cdsw.groupby(
            'transcript_id'
        ).rank(method='first', ascending=True)['start']
        # calculate the cumulate length of nucleotides in the CDS
        cdsw.loc[:,'cds_cumlen'] = cdsw.groupby('transcript_id')['cds_len'].cumsum()
        # the nucleotides start position in the CDS
        cdsw.loc[:,'cds_ntstart'] = cdsw['cds_cumlen'] - cdsw['cds_len'] + 1
        # merge the CDS DF with protein aa position DF
        cdswp = pd.merge(cdsw, protein_pos, on=['gene_id'], how='left')
        # keep aa and exon pair with codon position inside the exon
        cdswp = cdswp.loc[
            np.logical_or(
                np.logical_and(         # condition for the end of codon in exon
                    cdswp['aa_cds_pos_1'] >= cdswp['cds_ntstart'],
                    cdswp['aa_cds_pos_1'] <= cdswp['cds_cumlen']
                ),
                np.logical_and(         # condition for the start of codon in exon
                    cdswp['aa_cds_pos_3'] >= cdswp['cds_ntstart'],
                    cdswp['aa_cds_pos_3'] <= cdswp['cds_cumlen']
                )
            )
        ]
        # position of aa relative to exon start
        cdswp.loc[:, 'aa_exon_pos_1'] = cdswp['aa_cds_pos_1'] - cdswp['cds_ntstart'] + 1
        cdswp.loc[:, 'aa_exon_pos_2'] = cdswp['aa_cds_pos_2'] - cdswp['cds_ntstart'] + 1
        cdswp.loc[:, 'aa_exon_pos_3'] = cdswp['aa_cds_pos_3'] - cdswp['cds_ntstart'] + 1
        # calculate the position of last nt in codon of aa
        cdswp.loc[:, 'aa_genome_pos_end'] = (
            cdswp['start'] + cdswp['aa_exon_pos_3'] - 1
        ).astype(int)
        # calculate the position of first nt in codon of aa
        cdswp.loc[:, 'aa_genome_pos_start'] = (cdswp.loc[:, 'aa_genome_pos_end'] - 2).astype(int)
        # get the genome sequence of codon
        # fa dict Bio.Seq sequence start from 0
        cdswp.loc[:,'aa_genome_seq'] = cdswp.apply(
            lambda x: str(fa[x['seqname']].seq[x['aa_genome_pos_start'] - 1: x['aa_genome_pos_end']]).upper(),
            axis=1
        )
        cdswp.loc[:, 'aa_genome_pos_1'] = cdswp['aa_genome_pos_start']
        cdswp.loc[:, 'aa_genome_pos_2'] = cdswp['aa_genome_pos_start'] + 1
        cdswp.loc[:, 'aa_genome_pos_3'] = cdswp['aa_genome_pos_start'] + 2
        cdswp.loc[:, 'aa_genome_seq_1'] = cdswp.apply(
            lambda x: str(fa[x['seqname']].seq[x['aa_genome_pos_1'] - 1]).upper(),
            axis=1
        )
        cdswp.loc[:, 'aa_genome_seq_2'] = cdswp.apply(
            lambda x: str(fa[x['seqname']].seq[x['aa_genome_pos_2'] - 1]).upper(),
            axis=1
        )
        cdswp.loc[:, 'aa_genome_seq_3'] = cdswp.apply(
            lambda x: str(fa[x['seqname']].seq[x['aa_genome_pos_3'] - 1]).upper(),
            axis=1
        )
        # keep exon and aa pair with right codon
        # no splicing aa
        result['nosplice'].append(
            cdswp.loc[
                cdswp.apply(
                    lambda x: (
                        x['aa_exon_pos_3'] >= 3 and (
                            x['aa_genome_seq'] in
                            [p.replace('U', 'T') for p in amino_acid(x['aa'], 'codon')]
                        )
                    ),
                    axis=1
                )
            ].copy()
        )
        # aa with splicing
        # 1|23 split
        cdswp_spl_a_1 = cdswp.loc[
            cdswp.apply(
                lambda x: (
                    x['aa_cds_pos_1'] == x['cds_cumlen'] and (
                        x['aa_genome_seq_1'] in [
                            p.replace('U', 'T')[0]
                            for p in amino_acid(x['aa'], 'codon')
                        ]
                    )
                ),
                axis=1
            ),
            [
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos',
                'aa_cds_pos_1', 'aa_genome_pos_1', 'aa_genome_seq_1'
            ]
        ].copy()
        cdswp_spl_a_23 = cdswp.loc[
            cdswp.apply(
                lambda x:(
                    x['frame'] == '2' and x['aa_exon_pos_3'] == 2 and
                    (
                        x['aa_genome_seq'][1:] in [
                            p.replace('U', 'T')[1:]
                            for p in amino_acid(x['aa'], 'codon')
                        ]
                    )
                ),
                axis=1
            ),
            [
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos',
                'aa_cds_pos_2', 'aa_genome_pos_2', 'aa_genome_seq_2',
                'aa_cds_pos_3', 'aa_genome_pos_3', 'aa_genome_seq_3'
            ]
        ].copy()
        cdswp_spl_a = pd.merge(
            cdswp_spl_a_23, cdswp_spl_a_1, on=[
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos'
            ], how='inner'
        )
        # 12|3 split
        cdswp_spl_b_12 = cdswp.loc[
            cdswp.apply(
                lambda x: (
                    x['aa_cds_pos_2'] == x['cds_cumlen'] and (
                        ''.join([x['aa_genome_seq_1'],x['aa_genome_seq_2']]) in [
                            p.replace('U', 'T')[0:2]
                            for p in amino_acid(x['aa'], 'codon')
                        ]
                    )
                ),
                axis=1
            ),
            [
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos',
                'aa_cds_pos_1', 'aa_genome_pos_1', 'aa_genome_seq_1',
                'aa_cds_pos_2', 'aa_genome_pos_2', 'aa_genome_seq_2'
            ]
        ].copy()
        cdswp_spl_b_3 = cdswp.loc[
            cdswp.apply(
                lambda x: (
                    x['frame'] == '1' and x['aa_exon_pos_3'] == 1 and
                    (
                        x['aa_genome_seq'][2] in [
                            p.replace('U', 'T')[2]
                            for p in amino_acid(x['aa'], 'codon')
                        ]
                    )
                ),
                axis=1
            ),
            [
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos',
                'aa_cds_pos_3', 'aa_genome_pos_3', 'aa_genome_seq_3'
            ]
        ].copy()
        cdswp_spl_b = pd.merge(
            cdswp_spl_b_3, cdswp_spl_b_12, on=[
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos'
            ], how='inner'
        )
        result['splice'].append(
            pd.concat([cdswp_spl_a, cdswp_spl_b], axis=0).reset_index(drop=True)
        )
        del cdswp

        ####################
        # crick strand -
        cdsc = gffx.loc[
            np.logical_and(gffx.feature == 'CDS', gffx.strand=='-')
        ].copy()
        cdsc.loc[:,'cds_len'] = cdsc['end'] - cdsc['start'] + 1
        cdsc.sort_values(['transcript_id', 'end'], ascending=False, inplace=True)
        cdsc.reset_index(drop=True, inplace=True)
        # CDS number in transcripts
        cdsc.loc[:, 'cds_no'] = cdsc.groupby(
            'transcript_id'
        ).rank(method='first', ascending=False)['end']
        # calculate the cumulate length of nucleotides in the CDS
        cdsc.loc[:,'cds_cumlen'] = cdsc.groupby('transcript_id')['cds_len'].cumsum()
        # the nucleotides start position in the CDS
        cdsc.loc[:,'cds_ntstart'] = cdsc['cds_cumlen'] - cdsc['cds_len'] + 1
        # merge the CDS DF with protein aa position DF
        cdscp = pd.merge(cdsc, protein_pos, on=['gene_id'], how='left')
        # keep aa and exon pair with codon position inside the exon
        cdscp = cdscp.loc[
            np.logical_or(
                np.logical_and(         # condition for the end of codon in exon
                    cdscp['aa_cds_pos_1'] >= cdscp['cds_ntstart'],
                    cdscp['aa_cds_pos_1'] <= cdscp['cds_cumlen']
                ),
                np.logical_and(         # condition for the start of codon in exon
                    cdscp['aa_cds_pos_3'] >= cdscp['cds_ntstart'],
                    cdscp['aa_cds_pos_3'] <= cdscp['cds_cumlen']
                )
            )
        ]
        # position of aa relative to exon start
        cdscp.loc[:, 'aa_exon_pos_1'] = cdscp['aa_cds_pos_1'] - cdscp['cds_ntstart'] + 1
        cdscp.loc[:, 'aa_exon_pos_2'] = cdscp['aa_cds_pos_2'] - cdscp['cds_ntstart'] + 1
        cdscp.loc[:, 'aa_exon_pos_3'] = cdscp['aa_cds_pos_3'] - cdscp['cds_ntstart'] + 1
        # calculate the position of last nt in codon of aa
        cdscp.loc[:, 'aa_genome_pos_start'] = (
            cdscp['end'] - cdscp['aa_exon_pos_3'] + 1
        ).astype(int)
        # calculate the position of first nt in codon of aa
        cdscp.loc[:, 'aa_genome_pos_end'] = (cdscp.loc[:, 'aa_genome_pos_start'] + 2).astype(int)
        # get the genome sequence of codon
        # fa dict Bio.Seq sequence start from 0
        cdscp.loc[:,'aa_genome_seq'] = cdscp.apply(
            lambda x: str(fa[x['seqname']].seq[x['aa_genome_pos_start'] - 1: x['aa_genome_pos_end']]).upper(),
            axis=1
        )
        cdscp.loc[:, 'aa_genome_pos_1'] = cdscp['aa_genome_pos_end']
        cdscp.loc[:, 'aa_genome_pos_2'] = cdscp['aa_genome_pos_end'] - 1
        cdscp.loc[:, 'aa_genome_pos_3'] = cdscp['aa_genome_pos_end'] - 2
        cdscp.loc[:, 'aa_genome_seq_1'] = cdscp.apply(
            lambda x: str(fa[x['seqname']].seq[x['aa_genome_pos_1'] - 1]).upper(),
            axis=1
        )
        cdscp.loc[:, 'aa_genome_seq_2'] = cdscp.apply(
            lambda x: str(fa[x['seqname']].seq[x['aa_genome_pos_2'] - 1]).upper(),
            axis=1
        )
        cdscp.loc[:, 'aa_genome_seq_3'] = cdscp.apply(
            lambda x: str(fa[x['seqname']].seq[x['aa_genome_pos_3'] - 1]).upper(),
            axis=1
        )
        # keep exon and aa pair with right codon
        # no splicing aa
        result['nosplice'].append(
            cdscp.loc[
                cdscp.apply(
                    lambda x: (
                        x['aa_exon_pos_3'] >= 3 and (
                            x['aa_genome_seq'] in
                            [
                                reverse_complement(p.replace('U', 'T'))
                                for p in amino_acid(x['aa'], 'codon')
                            ]
                        )
                    ),
                    axis=1
                )
            ].copy()
        )
        # aa with splicing
        # 1|23 split
        cdscp_spl_a_1 = cdscp.loc[
            cdscp.apply(
                lambda x: (
                    x['aa_cds_pos_1'] == x['cds_cumlen'] and (
                        x['aa_genome_seq_1'] in [
                            reverse_complement(p.replace('U', 'T')[0])
                            for p in amino_acid(x['aa'], 'codon')
                        ]
                    )
                ),
                axis=1
            ),
            [
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos',
                'aa_cds_pos_1', 'aa_genome_pos_1', 'aa_genome_seq_1'
            ]
        ].copy()
        cdscp_spl_a_23 = cdscp.loc[
            cdscp.apply(
                lambda x:(
                    x['frame'] == '2' and x['aa_exon_pos_3'] == 2 and
                    (
                        x['aa_genome_seq'][0:2] in [
                            reverse_complement(p.replace('U', 'T')[1:])
                            for p in amino_acid(x['aa'], 'codon')
                        ]
                    )
                ),
                axis=1
            ),
            [
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos',
                'aa_cds_pos_2', 'aa_genome_pos_2', 'aa_genome_seq_2',
                'aa_cds_pos_3', 'aa_genome_pos_3', 'aa_genome_seq_3'
            ]
        ].copy()
        cdscp_spl_a = pd.merge(
            cdscp_spl_a_23, cdscp_spl_a_1, on=[
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos'
            ], how='inner'
        )
        # 12|3 split
        cdscp_spl_b_12 = cdscp.loc[
            cdscp.apply(
                lambda x: (
                    x['aa_cds_pos_2'] == x['cds_cumlen'] and (
                        ''.join([x['aa_genome_seq_2'],x['aa_genome_seq_1']]) in [
                            reverse_complement(p.replace('U', 'T')[0:2])
                            for p in amino_acid(x['aa'], 'codon')
                        ]
                    )
                ),
                axis=1
            ),
            [
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos',
                'aa_cds_pos_1', 'aa_genome_pos_1', 'aa_genome_seq_1',
                'aa_cds_pos_2', 'aa_genome_pos_2', 'aa_genome_seq_2'
            ]
        ].copy()
        cdscp_spl_b_3 = cdscp.loc[
            cdscp.apply(
                lambda x: (
                    x['frame'] == '1' and x['aa_exon_pos_3'] == 1 and
                    (
                        x['aa_genome_seq'][0] in [
                            reverse_complement(p.replace('U', 'T')[2])
                            for p in amino_acid(x['aa'], 'codon')
                        ]
                    )
                ),
                axis=1
            ),
            [
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos',
                'aa_cds_pos_3', 'aa_genome_pos_3', 'aa_genome_seq_3'
            ]
        ].copy()
        cdscp_spl_b = pd.merge(
            cdscp_spl_b_3, cdscp_spl_b_12, on=[
                'seqname', 'strand', 'gene_id', 'gene_name', 'p_id', 'transcript_id',
                'uniprot_1', 'uniprot_2', 'aa', 'aa_prot_pos'
            ], how='inner'
        )
        result['splice'].append(
            pd.concat([cdscp_spl_a, cdscp_spl_b], axis=0).reset_index(drop=True)
        )
        del cdscp

    nosplice = pd.concat(result['nosplice'], axis=0)
    nosplice[resultcol].to_csv(args['output'], index=False, quoting=0, sep='\t')
    splice = pd.concat(result['splice'], axis=0)
    splice.to_csv(args['splice_output'], index=False, quoting=0, sep='\t')
    # ------------------
    # sys.stdout = __console__
################################################################################
