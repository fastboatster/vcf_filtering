#! /usr/bin/env python3
import vcf
import sys


def get_alt(ann_field, var_type, var_subtype, ref, alts):
    """
    Return ALT value for a record in ANN list
    :param ann_field: one subrecord in ANN list in the INFO field in the vcf record
    :param var_type: snp or indel
    :param var_subtype: record subtype, i.e., single snp or indel or multiple indels (insertion/dup or deletions)
    :param ref: reference allele
    :param alts: list of ALTs
    :return: ALT value for the record in ANN list
    """
    # case when record is SNV:
    if var_type == 'snp':
        alt = ann_field[0]  # get Allele value
    # single indel or single ins/dup:
    if var_type == 'indel':
        if var_subtype == 'ins':
            alt = ref + ann_field[0]  # concat ref and allele value in the ANN record
        if var_subtype == 'del':
            # get the first value in alt list
            alt = alts[0].sequence
    # multiple deletions or insertions in the same record:
        if var_subtype == 'unknown':
            # check if record is ins/dup or del
            # get hgvs.c and check
            hgvs_c = ann_field[9]
            if 'ins' in hgvs_c or 'dup' in hgvs_c:
                alt = ref + ann_field[0]
            # this only handles the case when there's only one kind of deletion in the record
            # but all the cases when there's 2 or more don't affect protein structure in this particular vcf file
            if 'del' in hgvs_c:
                # get the shortest string in the list of ALTs
                alt = min(alts, key=len)
    return alt


def process_record(record):
    """
    Prints out chromosome, position, ref. allele, alt. allele, gene name, transcript id and resulting amino acid change
    for vcf file record
    :param record: pyvcf3 record object
    :return: None
    """
    # get chr, pos and ref:
    chrom = record.CHROM
    pos = int(record.POS)
    ref = record.REF
    alts = list(record.ALT)
    # get ANN list (list of annotations in the INFO field) for the record
    ann = record.INFO['ANN']
    # go over the items in ANN list
    for field in ann:
        # tokenize the field, turn it into a list:
        fields = field.split('|')
        # get biotype
        tr_biotyp = fields[7]
        # and hgvs.p values
        hgvs_p = fields[10]
        # only work with "protein_coding" fields which affect protein seq:
        if tr_biotyp == 'protein_coding' and hgvs_p:
            gene_name = fields[3]
            transcript_id = fields[6]
            # get ALT based on the variant type and subrecord type:
            alt = get_alt(fields, record.var_type, record.var_subtype, ref, alts)
            # now print transcript description:
            print(chrom, pos, ref, alt, gene_name, transcript_id, hgvs_p)


def main():
    vcf_file = sys.argv[1]
    # read vcf file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        # process each record:
        process_record(record)
    pass


if __name__ == "__main__":
    main()

