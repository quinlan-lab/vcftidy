"""
This implements:

1. allele splitting: ref=A, alt=AT,AAT,ACT becomes 3 separate variants.
   Adjusts PL (or GL) and removes AF and DP fields as needed.

2. left alignment as described in 'Unified Representation of Genetic Variants',
   Bioinformatics Bioinformatics (2015): btv112. From Tan et al.

3. simplify variants: pos=123, ref=AAT, alt=AT becomes pos=124, ref=AT, alt=T

4. When splitting alleles, it pulls out the correct subset of annotations from
   SnpEff (ANN= or EFF= tags) and VEP.


It adjusts AD for GATK and RO,AO for freebayes. It should also parse the header and
keep fields with Number=R or Number=A or Number=G, adjust those, and discard
others from the format field.
"""

from collections import OrderedDict
import sys

INFO_EXCLUDE = ('AF',)
GT_EXCLUDE = ('DP',)

def snpeff_ann(ann, alt, alts):
    # in these, the annos are separated by "," and the fields within each anno
    # are the alleles. They are just the allele itself
    ann = ann.split(",")
    return ",".join(x for x in ann if x.startswith("%s|" % alt))

def snpeff_eff(eff, alt, alts):
    """
    >>> eff = 'CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE28E|196|MICAL3|protein_coding|CODING|ENST00000498573|1|1|WARNING_TRANSCRIPT_NO_START_CODON),CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE4E|188|MICAL3|protein_coding|CODING|ENST00000578984|1|1|WARNING_TRANSCRIPT_INCOMPLETE),CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE998E|2002|MICAL3|protein_coding|CODING|ENST00000441493|21|1),CODON_INSERTION(MODERATE||gag/gaGGAg|E28EE|196|MICAL3|protein_coding|CODING|ENST00000498573|1|2|WARNING_TRANSCRIPT_NO_START_CODON),CODON_INSERTION(MODERATE||gag/gaGGAg|E4EE|188|MICAL3|protein_coding|CODING|ENST00000578984|1|2|WARNING_TRANSCRIPT_INCOMPLETE),CODON_INSERTION(MODERATE||gag/gaGGAg|E998EE|2002|MICAL3|protein_coding|CODING|ENST00000441493|21|2),DOWNSTREAM(MODIFIER||2412||966|MICAL3|protein_coding|CODING|ENST00000400561||2),DOWNSTREAM(MODIFIER||2415||966|MICAL3|protein_coding|CODING|ENST00000400561||1),DOWNSTREAM(MODIFIER||2422||966|MICAL3|protein_coding|CODING|ENST00000444520||2)'
    >>> snpeff_eff(eff, 'TTT', ['CTTCT', 'TTT'])
    'CODON_INSERTION(MODERATE||gag/gaGGAg|E28EE|196|MICAL3|protein_coding|CODING|ENST00000498573|1|2|WARNING_TRANSCRIPT_NO_START_CODON),CODON_INSERTION(MODERATE||gag/gaGGAg|E4EE|188|MICAL3|protein_coding|CODING|ENST00000578984|1|2|WARNING_TRANSCRIPT_INCOMPLETE),CODON_INSERTION(MODERATE||gag/gaGGAg|E998EE|2002|MICAL3|protein_coding|CODING|ENST00000441493|21|2),DOWNSTREAM(MODIFIER||2412||966|MICAL3|protein_coding|CODING|ENST00000400561||2),DOWNSTREAM(MODIFIER||2422||966|MICAL3|protein_coding|CODING|ENST00000444520||2)'

    >>> snpeff_eff(eff, 'CTTCT', ['CTTCT', 'TTT'])
    'CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE28E|196|MICAL3|protein_coding|CODING|ENST00000498573|1|1|WARNING_TRANSCRIPT_NO_START_CODON),CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE4E|188|MICAL3|protein_coding|CODING|ENST00000578984|1|1|WARNING_TRANSCRIPT_INCOMPLETE),CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE998E|2002|MICAL3|protein_coding|CODING|ENST00000441493|21|1),DOWNSTREAM(MODIFIER||2415||966|MICAL3|protein_coding|CODING|ENST00000400561||1)'

    """
    # http://snpeff.sourceforge.net/SnpEff_manual.version_4_0.html#input
    eff = eff.split(",")

    # this version of snpeff uses the 1-based index of the alt.
    # it stores that index in the 10th position of the EFF field.
    idx = str(alts.index(alt) + 1)
    subset = [e for e in eff if e.split("|")[10].rstrip(")") == idx]
    
    return ",".join(subset)


def vep_csq(csq, alt, alts):
    """
    >>> csq = 'intron_variant|||ENSG00000047056|WDR37|ENST00000263150||||-/494|protein_coding|1,downstream_gene_variant|||ENSG00000047056|WDR37|ENST00000436154||||-/209|protein_coding|1,intron_variant|||ENSG00000047056|WDR37|ENST00000358220||||-/494|protein_coding|1,stop_lost|Tga/Cga|*/R|ENSG00000047056|WDR37|ENST00000381329|9/9|||250/249|protein_coding|1'
    >>> vep_csq(csq, 'T', 'T') == csq
    True
    >>> vep_csq(csq, 'T', ['T', 'A']) == csq
    True

    >>> vep_csq(csq, 'A', ['T', 'A']) == ''
    True
    >>> extra = csq.split(",")[-1]

    # change its source
    >>> extra = extra[:-1] + "2"
    >>> vep_csq(csq + "," + extra, 'A', ['T', 'A']) == extra
    True
    >>> vep_csq(csq + "," + extra, 'T', ['T', 'A']) == csq
    True
    """

    csq = csq.split(',')
    idx = str(alts.index(alt) + 1)
    subset = [c for c in csq if c.split('|')[-1] == idx]
    return ",".join(subset)

# these have the signature fn(info_dict, alt, alts).
# e.g. fn(info, 'A', ['A', 'C', 'T'])
EFFECT_PARSERS = {'ANN': snpeff_ann,
                  'EFF': snpeff_eff,
                  'CSQ': vep_csq}

def fix_info(info, alt, alts):
    """
    this just dispatches to the correct function from EFFECT_PARSERS
    and sets the field
    """
    info = info.copy()
    assert isinstance(info, dict)
    try:
        # check if e.g. EFF/ANN/SQR are present in INFO and use the
        # correct func.
        field, fixer_fn = next((k, EFFECT_PARSERS[k]) for k in EFFECT_PARSERS if k in info)
    except StopIteration:
        return info
    info[field] = fixer_fn(info[field], alt, alts)
    return info

def infostr(infodict, exclude=INFO_EXCLUDE):
    if isinstance(infodict, basestring):
        return infodict
    for k in infodict:
        if isinstance(infodict[k], list):
            infodict[k] = ",".join(map(str, infodict[k]))
    s = ";".join(k + (('=' + str(v)) if v is not None else '') for
                    k, v in infodict.iteritems() if not k in exclude)
    return s

import copy

def fix_genos(genos, alt_idx, exclude=GT_EXCLUDE):
    """
    This example is taken from the Vt wiki for documenting the vt decompose
    tools which performs a similar task.

    >>> genos = [{'GT': ['1', '2'], 'PL': '281,5,9,58,0,115,338,46,116,809', 'DP': '81'},
    ...          {'GT': ['0', '0'], 'PL': '0,30,323,31,365,483,38,291,325,567', 'DP': '86'}]

    >>> fix_genos(genos, 3)
    [{'GT': ['.', '.'], 'PL': '281,338,809'}, {'GT': ['0', '0'], 'PL': '0,38,567'}]

    >>> fix_genos(genos, 2)
    [{'GT': ['.', '1'], 'PL': '281,58,115'}, {'GT': ['0', '0'], 'PL': '0,31,483'}]

    >>> fix_genos(genos, 1)
    [{'GT': ['1', '.'], 'PL': '281,5,9'}, {'GT': ['0', '0'], 'PL': '0,30,323'}]

    >>> genos = [{'GT': ['0', '1'], 'PL': '281,5,9,58,0,115,338,46,116,809', 'DP': '81', 'AO': ['22', '33']},
    ...          {'GT': ['0', '2'], 'PL': '0,30,323,31,365,483,38,291,325,567', 'DP': '86', 'AO': ['122', '133']}]

    >>> fix_genos(genos, 1)
    [{'GT': ['0', '1'], 'AO': '22', 'PL': '281,5,9'}, {'GT': ['0', '.'], 'AO': '122', 'PL': '0,30,323'}]
    
    >>> fix_genos(genos, 2)
    [{'GT': ['0', '.'], 'AO': '33', 'PL': '281,58,115'}, {'GT': ['0', '1'], 'AO': '133', 'PL': '0,31,483'}]


    >>> genos = [{'GT': ['0', '0'], 'PL': '281,5,9,58,0,115,338,46,116,809', 'GL': '281,5,9,58,0,115,338,46,116,809', 'DP': '81', 'AO': '22'}]
    >>> fix_genos(genos, 1)
    [{'GT': ['0', '0'], 'GL': '281,5,9', 'AO': '22', 'PL': '281,5,9'}]

    >>> fix_genos(genos, 2)
    [{'GT': ['0', '0'], 'GL': '281,58,115', 'AO': '-1', 'PL': '281,58,115'}]

    """

    assert int(alt_idx) > 0
    # this copy needed since we are modifying in place.
    genos = [copy.deepcopy(x) for x in genos]

    alt_idx = str(alt_idx)
    for n, geno in enumerate(genos):
        # this copy *also* needed since we are modifying in place.
        geno['GT'] = geno['GT'][:]
        for ex in exclude:
            geno.pop(ex, None)

        # pull out only the appropriate likelihoods.
        # according to VCF spec, these will be in the order:
        # 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
        # ordering given by: Order(j/k) = (k*(k+1)/2)+j.

        a, b = 0, int(alt_idx)
        for key in ('PL', 'GL'):
            if not key in geno: continue

            vals = geno[key].split(",")
            # so do this for (a, a), (a, b), (b, b)
            likelihood = []
            for j, k in ((a, a), (a, b), (b, b)):
                order = (k * (k + 1) / 2) + j
                #print >>sys.stderr, a, b, (j, k), 'order:' + str(order), vals#, 'value:' + str(vals[order])
                try:
                    likelihood.append(vals[order])
                except IndexError:
                    # needed for older freebayes.
                    likelihood.append('-1')
            geno[key] = ",".join(likelihood)

        for key in ('AO', 'PAO', 'TYPE', 'QA', 'PQA', 'SAF'):
            if key in geno:
                #assert "," in geno[key] or isinstance(geno[key], list), (key, geno)
                if not isinstance(geno[key], list):
                    vals = geno[key].split(",")
                else:
                    vals = map(str, geno[key])
                try:
                    geno[key] = vals[int(alt_idx) - 1]
                except IndexError:
                    # required for older versions of freebayes which didn't
                    # correctly support variants with number of alts > 2.
                    geno[key] = '-1'

        # these are split by comma and include the ref as the first.
        for key in ('AD',):
            if key in geno:
                vals = geno[key]
                if isinstance(vals, basestring):
                    vals = vals.split(",")
                # don't subtract 1 from vals index since the ref is 0 index.
                geno[key] = ",".join((vals[0], vals[int(alt_idx)]))

        # set the current allele to 1, keep 0 at 0 and set all others to '.'
        for i, allele in enumerate(geno['GT']):

            if allele == alt_idx:
                geno['GT'][i] = '1'
            elif allele != '0':
                geno['GT'][i] = '.'

        # adjust the depth to be alt-specific.
        if 'RO' in geno and 'AO' in geno:
            if '-1' in (geno['RO'], geno['AO']):
                geno['DP'] = '-1'
            else:
                geno['DP'] = str(int(geno['AO']) + int(geno['RO']))

        if 'AD' in geno:
            vals = map(int, geno['AD'].split(","))
            geno['DP'] = str(sum(vals))

    return genos

def fmt_genos(genos):
    for g in genos:
        if 'GT' in g:
            sep = g.pop('__sep', '/')
            g['GT'] = sep.join(g['GT'])
    vals = [":".join(str(v) for v in f.values()) for f in genos]
    return vals
    

def varsplit(ref, alts, info, frmt, gts):
    """
    Split a variant with multiple alternate alleles into separate records to aid
    in annotation.

    >>> n = varsplit('TA', ['TAA', 'TAAA', 'T'], 'AF=0.2,0.3,0.1;KKKK',
    ...              'GT:DP:PL', ['1/2:81:281,5,9,58,0,115,338,46,116,809',
    ...              '0/0:86:0,30,323,31,365,483,38,291,325,567'])
    >>> next(n)
    ('TA', 'TAA', 'KKKK;ORIG_ALLELES=TA/TAA,TAAA,T;ORIG_ALLELE_i=1', 'GT:PL', ['1/.:281,5,9', '0/0:0,30,323'])

    >>> next(n)
    ('TA', 'TAAA', 'KKKK;ORIG_ALLELES=TA/TAA,TAAA,T;ORIG_ALLELE_i=2', 'GT:PL', ['./1:281,58,115', '0/0:0,31,483'])

    >>> next(n)
    ('TA', 'T', 'KKKK;ORIG_ALLELES=TA/TAA,TAAA,T;ORIG_ALLELE_i=3', 'GT:PL', ['./.:281,338,809', '0/0:0,38,567'])


    # test new snp eff ANN= field
    >>> n = varsplit('G', ['A', 'C', 'T'],
    ...              'ANN=A|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.|||||||', 
    ...              'GT:DP:PL', ['1/2:81:281,5,9,58,0,115,338,46,116,809'])

    >>> next(n)
    ('G', 'A', 'ANN=A|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||,A|.|.|.|.|.|.|.|.|.||||||;ORIG_ALLELES=G/A,C,T;ORIG_ALLELE_i=1', 'GT:PL', ['1/.:281,5,9'])

    >>> next(n)
    ('G', 'C', 'ANN=C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||,C|.|.|.|.|.|.|.|.|.||||||;ORIG_ALLELES=G/A,C,T;ORIG_ALLELE_i=2', 'GT:PL', ['./1:281,58,115'])

    >>> next(n)
    ('G', 'T', 'ANN=T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.||||||,T|.|.|.|.|.|.|.|.|.|||||||;ORIG_ALLELES=G/A,C,T;ORIG_ALLELE_i=3', 'GT:PL', ['./.:281,338,809'])

    >>> next(n)
    Traceback (most recent call last):
    ...
    StopIteration

    # test older snpeff EFF= field
    >>> n = varsplit('TTCC', ['T', 'TTCCTCC'], 'EFF=CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE28E|196|MICAL3|protein_coding|CODING|ENST00000498573|1|1|WARNING_TRANSCRIPT_NO_START_CODON),CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE4E|188|MICAL3|protein_coding|CODING|ENST00000578984|1|1|WARNING_TRANSCRIPT_INCOMPLETE),CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE998E|2002|MICAL3|protein_coding|CODING|ENST00000441493|21|1),CODON_INSERTION(MODERATE||gag/gaGGAg|E28EE|196|MICAL3|protein_coding|CODING|ENST00000498573|1|2|WARNING_TRANSCRIPT_NO_START_CODON),CODON_INSERTION(MODERATE||gag/gaGGAg|E4EE|188|MICAL3|protein_coding|CODING|ENST00000578984|1|2|WARNING_TRANSCRIPT_INCOMPLETE),CODON_INSERTION(MODERATE||gag/gaGGAg|E998EE|2002|MICAL3|protein_coding|CODING|ENST00000441493|21|2),DOWNSTREAM(MODIFIER||2412||966|MICAL3|protein_coding|CODING|ENST00000400561||2),DOWNSTREAM(MODIFIER||2415||966|MICAL3|protein_coding|CODING|ENST00000400561||1),DOWNSTREAM(MODIFIER||2422||966|MICAL3|protein_coding|CODING|ENST00000444520||2),DOWNSTREAM(MODIFIER||2425||966|MICAL3|protein_coding|CODING|ENST00000444520||1)', 'GT:DP:PL', ['1/2:81:281,5,9,58,0,115,338'])
    >>> next(n)
    ('TTCC', 'T', 'EFF=CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE28E|196|MICAL3|protein_coding|CODING|ENST00000498573|1|1|WARNING_TRANSCRIPT_NO_START_CODON),CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE4E|188|MICAL3|protein_coding|CODING|ENST00000578984|1|1|WARNING_TRANSCRIPT_INCOMPLETE),CODON_CHANGE_PLUS_CODON_DELETION(MODERATE||gaggaa/gaa|EE998E|2002|MICAL3|protein_coding|CODING|ENST00000441493|21|1),DOWNSTREAM(MODIFIER||2415||966|MICAL3|protein_coding|CODING|ENST00000400561||1),DOWNSTREAM(MODIFIER||2425||966|MICAL3|protein_coding|CODING|ENST00000444520||1);ORIG_ALLELES=TTCC/T,TTCCTCC;ORIG_ALLELE_i=1', 'GT:PL', ['1/.:281,5,9'])

    >>> next(n)
    ('TTCC', 'TTCCTCC', 'EFF=CODON_INSERTION(MODERATE||gag/gaGGAg|E28EE|196|MICAL3|protein_coding|CODING|ENST00000498573|1|2|WARNING_TRANSCRIPT_NO_START_CODON),CODON_INSERTION(MODERATE||gag/gaGGAg|E4EE|188|MICAL3|protein_coding|CODING|ENST00000578984|1|2|WARNING_TRANSCRIPT_INCOMPLETE),CODON_INSERTION(MODERATE||gag/gaGGAg|E998EE|2002|MICAL3|protein_coding|CODING|ENST00000441493|21|2),DOWNSTREAM(MODIFIER||2412||966|MICAL3|protein_coding|CODING|ENST00000400561||2),DOWNSTREAM(MODIFIER||2422||966|MICAL3|protein_coding|CODING|ENST00000444520||2);ORIG_ALLELES=TTCC/T,TTCCTCC;ORIG_ALLELE_i=2', 'GT:PL', ['./1:281,58,115'])

    >>> next(n)
    Traceback (most recent call last):
    ...
    StopIteration


    """
    if len(alts) == 1:
        yield ref, alts[0], infostr(info, ()), frmt, gts
        raise StopIteration

    if not isinstance(info, dict):
        info = OrderedDict((kv[0], (kv[1] if len(kv) > 1 else None)) for kv in (x.split('=') for x in info.split(';')))

    fmts = frmt.split(':')
    gts = [OrderedDict(zip(fmts, x.split(':'))) for x in gts]

    for i, g in enumerate(gts):
        if 'GT' in g:
            gts[i]['__sep'] = "|" if "|" in g['GT'] else "/"
            gts[i]['GT'] = g['GT'].split(gts[i]['__sep'])

    for i, alt in enumerate(alts, start=1):
        if any('GT' in g for g in gts):
            fixed_genos = fix_genos(gts, i)
        else:
            fixed_genos = gts
        if fixed_genos:
            fields = [f for f in fixed_genos[0].keys() if f != "__sep"] # ordered
        else:
            fields = []
        info_sub = fix_info(info, alt, alts)
        # now we temporarily replace info[field] with ann

        info_sub['ORIG_ALLELES'] = '%s/%s' % (ref, ",".join(alts))
        info_sub['ORIG_ALLELE_i'] = str(i)
        ret = ref, alt, infostr(info_sub), ":".join(fields), fmt_genos(fixed_genos)
        yield ret


def leftalign(chrom, pos, ref, alt, fa, max_shift=1000):        
    seq = fa[chrom][max(0, pos - max_shift - 1):pos + len(ref) - 1]
    assert seq.endswith(ref), (chrom, pos, ref, alt, seq[-10:])
    return _leftalign(pos, ref, alt, seq)[:3]


def _leftalign(pos, ref, alt, seq):
    """
    
    simple implementation from the vt paper:
    # actual variant is 2-base CA insertion.
    Last argument indicates whether we ran out of sequence and therefore did not
    finish left-aligning before running out of sequence. (False is bad).

    >>> _leftalign(123, 'CAC', 'C', 'GGGCACACAC')
    (118, 'GCA', 'G', True)

    # run out of sequence!
    >>> _leftalign(123, 'CAC', 'C', 'CACACAC')
    (119, 'CAC', 'C', False)

    >>> _leftalign(123, 'CCA', 'CAA', 'ACCCCCCA')
    (123, 'CC', 'CA', True)

    # have to left-trim after left-align
    >>> normalize(*_leftalign(123, 'CCA', 'CAA', 'ACCCCCCA')[:3], left_only=True)
    (124, 'C', 'A')

    >>> _leftalign(123, 'C', 'A', 'ACCCCCC')
    (123, 'C', 'A', True)

    """
    assert seq.endswith(ref)
    assert ref != alt
    seq = seq[:-len(ref)]
    ref, alt = list(ref), list(alt)
    j = 0

    quit = False
    while j < len(seq) and not quit:
        quit = True

        if ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
            quit = False

        if len(ref) == 0 or len(alt) == 0:
            j += 1
            ref = [seq[-j]] + ref
            alt = [seq[-j]] + alt
            quit = False

    return pos - j, "".join(ref), "".join(alt), quit

def normalize(pos, ref, alt, left_only=False):
    """simplify a ref/alt a la vt normalize so that ref=CC, alt=CCT becomes
    ref=C, alt=CT. this helps in annotating variants.

    This code relies on the observation by Eric Minikel that many annotation
    misses can be addressed by removing common suffix and prefixes.
    (http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/)

    >>> normalize(123, 'T', 'C')
    (123, 'T', 'C')

    >>> normalize(123, 'CC', 'CCT')
    (124, 'C', 'CT')

    >>> normalize(123, 'TCCCCT', 'CCCCT')
    (123, 'TC', 'C')

    >>> normalize(123, 'TCCCCTA', 'CCCCT')
    (123, 'TCCCCTA', 'CCCCT')

    >>> normalize(123, 'TCCCCTA', 'CCCCTA')
    (123, 'TC', 'C')

    >>> normalize(123, 'AAATCCCCTA', 'AAACCCCTA')
    (125, 'AT', 'A')

    >>> normalize(123, 'CAAATCCCCTAG', 'AAACCCCTA')
    (123, 'CAAATCCCCTAG', 'AAACCCCTA')
    """
    if len(ref) == len(alt) == 1:
        return pos, ref, alt

    # logic for trimming from either end is the same so we just reverse the
    # string to trim the suffix (i == 0) and correct position when doing prefix
    # (i == 1). To support alleles that have already been right-trimmed from
    # _leftalign, we allow the left_only argument to just do prefix-trimming.
    if left_only:
        sides = (1,)
    else:
        sides = (0, 1)
    for i in sides:
        if i == 0: # suffix so flip
            ref, alt = ref[::-1], alt[::-1]

        n, min_len = 0, min(len(ref), len(alt))
        while n + 1 < min_len and alt[n] == ref[n]:
            n += 1

        alt, ref = alt[n:], ref[n:]
        if i == 0: # flip back
            ref, alt = ref[::-1], alt[::-1]
        else: # add to position since we stripped the prefix
            pos += n

    return pos, ref, alt


def leftnorm(chrom, pos, ref, alt, fa=None):
    """
    this is the normalization function that should be used.
    if no fasta is present, then it just normalizes. Otherwise
    it left-aligns and then normalizes.
    """
    if fa is None:
        return normalize(pos, ref, alt)

    return normalize(*leftalign(chrom, pos, ref, alt, fa), left_only=True)
    

def main(fh, fa=None, buffer_size=5000):

    # heap so that variants remain in sorted order even after normalizing.
    import heapq

    if fa:
        from pyfaidx import Fasta
        fa = Fasta(fa, as_raw=True, read_ahead=40000)

    def gen_variants(fh):
        for line in fh:
            fields = line.rstrip().split("\t")
            pos, ref, alts, info = (int(fields[1]), fields[3], fields[4].split(","), fields[7]) 
            if len(fields) > 8:
                frmt, gts = fields[8], fields[9:]
            else:
                frmt, gts = '', []

            if len(alts) > 1:
                for sref, salt, sinfo, sfmt, sgts in varsplit(ref, alts, info, frmt, gts):

                    fields[1], fields[3], fields[4] = leftnorm(fields[0], pos, sref, salt, fa)
                    fields[7] = sinfo
                    if len(fields) > 8:
                        fields[8], fields[9:] = sfmt, sgts
                    nfields = fields[:]
                    nfields.insert(2, len(fields[3]))
                    yield nfields
            else:
                fields[1], fields[3], fields[4] = leftnorm(fields[0], pos, ref, alts[0], fa)
                fields.insert(2, len(fields[3]))
                yield fields
    def _print(ov):
        ov[1] = str(ov[1])
        ov.pop(2)
        print "\t".join(ov)
    
    var_buffer = []
    last_chrom = None

    for line in fh:
        # just output header
        if line[0] == "#":
            print line.rstrip()
        if line.startswith("#CHROM\t"): break

    for v in gen_variants(fh):
        # side step chromosome ordering problems.
        if v[0] != last_chrom:
            while var_buffer:
                _print(heapq.heappop(var_buffer))

        if len(var_buffer) == buffer_size:
            ov = heapq.heappushpop(var_buffer, v)
            _print(ov)
        else:
            heapq.heappush(var_buffer, v)
        last_chrom = v[0]

    while var_buffer:
        _print(heapq.heappop(var_buffer))

if __name__ == "__main__":
    import sys, os
    import doctest
    sys.stderr.write(str(doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS | doctest.REPORT_ONLY_FIRST_FAILURE,
            verbose=0)) + "\n")
    import gzip

    if len(sys.argv) >= 2 and (os.path.isfile(sys.argv[1]) or sys.argv[1] == "-"):
        if sys.argv[1] == "-": fh = sys.stdin
        elif sys.argv[1].endswith(".gz"): fh = gzip.open(sys.argv[1])
        else: fh = open(sys.argv[1])

        ref = None
        if len(sys.argv) == 3 and os.path.isfile(sys.argv[2]):
            ref = sys.argv[2]

        sys.exit(main(fh, ref))

