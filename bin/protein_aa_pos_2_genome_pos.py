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

####################
# gff

# -*-coding:utf-8-*-

# GFF format
# 1. seqname: 'chr1', 'chr2', ..., 'chr22', 'chrX', 'chrY'
# 2. source:
# 3. feature: 'gene', 'exon', 'intron', 'CDS', 'region', ...
# 4. start: number of start site
# 5. end: number of end site
# 6. score: score number
# 7. strand: '+', '-', '.'
# 8. frame: The frame of the coding sequence. '0', '1', '2', '.'
# 9. attribute:


# ------------------
# Functions
# ------------------

def regionseq(seqname, seq, start, end):
    # get the sequences by start and end for the input.
    region = pd.concat(
        [pd.Series(start), pd.Series(end)],
        axis=1
    )
    region.columns = ['start', 'end']
    region['seq'] = region.apply(
        lambda x: ''.join(list(seq[x['start']:x['end']])),
        axis=1
    )
    region['seqname'] = [seqname] * region.shape[0]
    return region[['seqname', 'start', 'end', 'seq']]

def gff_overlap(gff1, gff2):
    # gff1 and gff2 should be pd.DataFrame
    def apinbr(ap, bs, be):
        # a point in b region
        # ap, bs, be should be pd.Series
        if len(bs) != len(be):
            raise ValueError('bs should have same length with be')
        ap = pd.Series(ap)
        bs = pd.Series(bs)
        be = pd.Series(be)
        dim = (len(ap), len(bs))
        one = np.array([np.array(x) for x in ap.map(lambda x: x >= bs)])
        one.shape = dim
        two = np.array([np.array(x) for x in ap.map(lambda x: x <= be)])
        two.shape = dim
        return np.where(np.logical_and(one, two))

    def chr_overlap(gff1, gff2):
        # get the index and columns of two gff file for same chromosome
        gff1colnames = gff1.columns + '.1'
        gff2colnames = gff2.columns + '.2'
        gff1index = gff1.index      # gff1 index
        gff2index = gff2.index      # gff2 index
        # get the pairs of gff1 overlap gff2
        # gff1 left overlap gff2
        one_left = apinbr(gff1['start'], gff2['start'], gff2['end'])
        # gff1 right overlap gff2
        one_right = apinbr(gff1['end'], gff2['start'], gff2['end'])
        # gff2 left overlap gff1
        two_left = apinbr(gff2['start'], gff1['start'], gff1['end'])
        # gff2 right overlap gff1
        two_right = apinbr(gff2['end'], gff1['start'], gff1['end'])
        # merge two kinds of overlaps
        allpairs = (
            np.concatenate(
                (one_left[0], one_right[0], two_left[1], two_right[1])
            ),
            np.concatenate(
                (one_left[1], one_right[1], two_left[0], two_right[0])
            )
        )
        pairindex = [
            '{0}.{1}'.format(
                gff1index[allpairs[0][i]],
                gff2index[allpairs[1][i]]
            )
            for i in range(len(allpairs[0]))
        ]
        notduplicated = np.where(
            np.logical_not(pd.Series(pairindex).duplicated())
        )
        allpair_rd = (allpairs[0][notduplicated], allpairs[1][notduplicated])
        # get the gff1 in overlaps
        gff1side = gff1.loc[gff1index[allpair_rd[0]]].reset_index(drop=True)
        # get the gff1 in overlaps
        gff2side = gff2.loc[gff2index[allpair_rd[1]]].reset_index(drop=True)
        # merge two gff data
        gffpaired = pd.concat([gff1side, gff2side], axis=1)
        gffpaired.columns = gff1colnames.append(gff2colnames)
        return gffpaired

    resultlist = list()
    for seqname in gff1['seqname'].unique():
        resultlist.append(
            chr_overlap(
                gff1.loc[gff1['seqname'] == seqname],
                gff2.loc[gff2['seqname'] == seqname]
            )
        )
    result = pd.concat(resultlist, axis=0)
    return result


def gff_not_overlap(gff1, gff2):
    overlap = gff_overlap(gff1, gff2)
    gff1str = gff1['seqname'] + ':' + gff1['start'].map(str) + '-' + gff1['end'].map(str)
    gff2str = gff2['seqname'] + ':' + gff2['start'].map(str) + '-' + gff2['end'].map(str)
    overlap1str = overlap['seqname.1'] + ':' + overlap['start.1'].map(str) + '-' + overlap['end.1'].map(str)
    overlap2str = overlap['seqname.2'] + ':' + overlap['start.2'].map(str) + '-' + overlap['end.2'].map(str)
    gff1notoverlap = gff1.loc[np.logical_not(gff1str.isin(overlap1str))]
    gff2notoverlap = gff2.loc[np.logical_not(gff2str.isin(overlap2str))]
    return (gff1notoverlap, gff2notoverlap)


def gff_inmerge(gff1, gff2, method='intersect'):
    # get the merged common region in two gffs
    # method: intersect, union
    overlap = gff_overlap(gff1, gff2)
    colname2del = list(overlap.columns)
    overlap['seqname'] = overlap['seqname.1']
    overlap['source'] = overlap['source.1'] + '.' + overlap['source.2']
    overlap['feature'] = overlap['feature.1'] + '.' + overlap['feature.2']
    overlap['frame'] = overlap['frame.1'] + '.' + overlap['frame.2']
    if method == 'intersect':
        overlap['start'] = overlap[['start.1', 'start.2']].max(axis=1)
        overlap['end'] = overlap[['end.1', 'end.2']].min(axis=1)
    elif method == 'union':
        overlap['start'] = overlap[['start.1', 'start.2']].min(axis=1)
        overlap['end'] = overlap[['end.1', 'end.2']].max(axis=1)
    else:
        raise ValueError('Method should be intersect or union.')
    overlap['score'] = overlap['score.1']
    overlap['strand'] = overlap[['strand.1', 'strand.2']].apply(
        lambda x: x['strand.1'] if x['strand.1'] == x['strand.2'] else '.',
        axis=1
    )
    overlap['attribute'] = overlap['attribute.1'] + '.' + overlap['attribute.2']
    # remove old columns
    for x in colname2del:
        del overlap[x]
    # return results
    return overlap.reset_index(drop=True)


def gff_collapse(gff):
    # merge overlaped region in one gff
    selfmerge = gff_inmerge(gff, gff, 'union')
    # remove small end
    selfmerge2 = selfmerge.groupby(
        [
            'seqname', 'source', 'feature', 'frame',
            'start', 'score', 'strand', 'attribute'
        ]
    ).apply(
        lambda x: x['end'].max()
    )
    selfmerge2.name = 'end'
    selfmerge2 = selfmerge2.reset_index()[selfmerge.columns]
    # remove large start
    selfmerge3 = selfmerge2.groupby(
        [
            'seqname', 'source', 'feature', 'frame',
            'end', 'score', 'strand', 'attribute'
        ]
    ).apply(
        lambda x: x['start'].min()
    )
    selfmerge3.name = 'start'
    selfmerge3 = selfmerge3.reset_index()[selfmerge.columns]
    return selfmerge3


def gff_intersect(gff1, gff2):
    # get the intersect of gff
    result = gff_inmerge(gff1, gff2, 'intersect')
    return result


def gff_union(gff1, gff2):
    # get the union of gff
    result = gff_inmerge(gff1, gff2, 'union')
    result = gff_collapse(result)
    return result


def gff_merge(gff1, gff2, method='union'):
    # get merged gff
    overlap = gff_inmerge(gff1, gff2, method=method)
    notoverlap = pd.concat(
        gff_not_overlap(gff1, gff2)
    ).reset_index(drop=True)
    result = pd.concat([overlap, notoverlap]).reset_index(drop=True)
    return result


def gff_subtract(gff1, gff2):
    # subtract the intersect region from gff1
    pass


# ------------------
# Variables
# ------------------



# ------------------
# Errors
# ------------------


class GFFError(ValueError):
    pass


# ------------------
# Classes
# ------------------



class Gff:
    def __init__(self,
                 seqname,
                 source,
                 feature,
                 start,
                 end,
                 score,
                 strand,
                 frame,
                 attribute,
                 log_level='info'):
        # make gff DataFrame
        self.gff = pd.DataFrame(
            {
                'seqname': seqname,
                'source': source,
                'feature': feature,
                'start': start,
                'end': end,
                'score': score,
                'strand': strand,
                'frame': frame,
                'attribute': attribute
            }
        )
        self.gff = self.gff[
            ['seqname', 'source', 'feature', 'start', 'end',
             'score', 'strand', 'frame', 'attribute']
        ]
        self.loc = self.gff.loc
        self.index = self.gff.index
        self.columns = self.gff.columns
        self.shape = self.gff.shape

    def __getitem__(self, key):
        return self.gff[key]

    def __len__(self):
        return len(self.gff)

    def __repr__(self):
        return self.gff.__repr__()

    def head(self, n=5):
        return self.loc[0:n,]

    def tail(self, n=5):
        rows=self.gff.shape[0]
        return self.loc[list(range(rows - n, rows)),]

    def intersect(self, gff):
        outtype = 'DF'
        if isinstance(gff, Gff):
            inputgff = gff.gff
            outtype = 'Gff'
        elif isinstance(gff, pd.DataFrame):
            inputgff = gff
        else:
            raise TypeError('Only accept Gff or DataFrame class.')
        result = gff_intersect(self.gff, gff)
        if outtype == 'Gff':
            result = Gff(
                result['seqname'],
                result['source'],
                result['feature'],
                result['start'],
                result['end'],
                result['score'],
                result['strand'],
                result['frame'],
                result['attribute']
            )
        return result

    def union(self, gff):
        outtype = 'DF'
        if isinstance(gff, Gff):
            inputgff = gff.gff
            outtype = 'Gff'
        elif isinstance(gff, pd.DataFrame):
            inputgff = gff
        else:
            raise TypeError('Only accept Gff or DataFrame class.')
        result = gff_union(self.gff, gff)
        if outtype == 'Gff':
            result = Gff(
                result['seqname'],
                result['source'],
                result['feature'],
                result['start'],
                result['end'],
                result['score'],
                result['strand'],
                result['frame'],
                result['attribute']
            )
        return result

    def merge(self, gff, method='union'):
        outtype = 'DF'
        if isinstance(gff, Gff):
            inputgff = gff.gff
            outtype = 'Gff'
        elif isinstance(gff, pd.DataFrame):
            inputgff = gff
        else:
            raise TypeError('Only accept Gff or DataFrame class.')
        result = gff_merge(self.gff, inputgff, method)
        if outtype == 'Gff':
            result = Gff(
                result['seqname'],
                result['source'],
                result['feature'],
                result['start'],
                result['end'],
                result['score'],
                result['strand'],
                result['frame'],
                result['attribute']
            )
        return result

    def overlap(self, gff):
        outtype = 'DF'
        if isinstance(gff, Gff):
            inputgff = gff.gff
            outtype = 'Gff'
        elif isinstance(gff, pd.DataFrame):
            inputgff = gff
        else:
            raise TypeError('Only accept Gff or DataFrame class.')
        result = gff_overlap(self.gff, inputgff)
        if outtype == 'Gff':
            result = Gff(
                result['seqname'],
                result['source'],
                result['feature'],
                result['start'],
                result['end'],
                result['score'],
                result['strand'],
                result['frame'],
                result['attribute']
            )
        return result

    def not_overlap(self, gff):
        outtype = 'DF'
        if isinstance(gff, Gff):
            inputgff = gff.gff
            outtype = 'Gff'
        elif isinstance(gff, pd.DataFrame):
            inputgff = gff
        else:
            raise TypeError('Only accept Gff or DataFrame class.')
        result = gff_not_overlap(self.gff, inputgff)
        if outtype == 'Gff':
            result = Gff(
                result['seqname'],
                result['source'],
                result['feature'],
                result['start'],
                result['end'],
                result['score'],
                result['strand'],
                result['frame'],
                result['attribute']
            )
        return result

    def get_seq(self, seqdict):
        # return sequences from seq dict, with seq name as key and
        # sequence string as value
        if not isinstance(seqdict, dict):
            raise TypeError(
                'seqdict should be a dict contains sequence string, with seqname as dict keys.'
            )
        result = dict()
        for x in self.gff['seqname'].unique():
            result[x] = regionseq(
                x, seqdict[x],
                self.gff['start'][self.gff['seqname'] == x].tolist(),
                self.gff['end'][self.gff['seqname'] == x].tolist()
            )
        return pd.concat(result.values(), axis=0)


class GffFile:
    '''
    Manage gff or gtf file
    '''
    def __init__(self,
                 gfffile,
                 header=None,
                 log_level='info'):

        self.__colnames = ['seqname', 'source', 'feature',
                           'start', 'end', 'score',
                           'strand', 'frame', 'attribute']
        self.gfffile = gfffile
        self.__format = self.gfffile.split('.')[-1]
        parsegff = self.__parsegff(
            self.gfffile,
            self.__format,
            header
        )
        self.gff = Gff(
            parsegff['seqname'],
            parsegff['source'],
            parsegff['feature'],
            parsegff['start'],
            parsegff['end'],
            parsegff['score'],
            parsegff['strand'],
            parsegff['frame'],
            parsegff['attribute']
        )
        attrcol = [x for x in parsegff.columns if x not in self.__colnames]
        self.gffattr = parsegff[attrcol]

    def __parsegff(self, filename, fileformat, header):
        '''
        Parse gff3 or gtf file into pandas DataFrame.
        Parameters
        ----------
        filename: the file path of the gff3 or gtf file.
        fileformat: the format of file, should be "gtf" or "gff".
        columns: the column names of the file.
        Returns
        -------
        DataFrame stored parsed information from gff or gtf
        Examples
        --------
        >>> import pandas as pd
        >>> mygff = parsegff(
                './hsa.gtf',
                'gtf',
                ['seqname', 'source', 'feature', 'start', 'end',
                 'score', 'strand', 'frame', 'attribute']
            )
        '''
        if fileformat == 'gtf':
            # make different separation according to fileformat.
            sep1, sep2 = ' ', '; '
        elif fileformat == 'gff3' or fileformat == 'gff':
            sep1, sep2 = '=', ';'
        data = pd.read_table(
            filename,
            sep='\t',
            header=header,
            index_col=False,
            comment='#'
        )                           # read gtf/gff data
        data.columns = self.__colnames
        attr_split = data[self.__colnames[-1]].apply(
            lambda row: [
                [y.replace(';"', '') for y in x.split(sep1)]
                for x in row.split(sep2) if len(x) > 2
            ]
        )
        # split the attributes in attribute columns
        attr_dict = attr_split.apply(dict)
        # make splited attributes into dicts
        attr_columns = attr_dict.apply(
            lambda row: list(row.keys())
        ).tolist()                  # get attr columns names
        attr_names = list(
            set([
                attr_columns[i][j] for i in range(len(attr_columns))
                for j in range(len(attr_columns[i]))
            ])
        )
        # get attr columns names
        attr = pd.DataFrame(
            dict(
                zip(
                    attr_names,
                    [
                        pd.Series(
                            [x[attr_name] if attr_name in x else np.NaN
                             for x in attr_dict]
                        ) for attr_name in attr_names
                    ]
                )
            )
        )                           # make attr columns
        data = data.join(
            attr
        )
        # link attr columns with annotation.
        return data


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
