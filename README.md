vcftidy
=======

bring some order to the *fruit salad* that is VCF.

About
-----

VCF is the standard for representing variants, but different aligners use different conventions for key pieces of information, e.g. alt depths in the sample fields. vcftidy aims to improve this by

 1. putting ref and alt read depths into the AD field for each genotype (a la GATK) and setting Number=A for the header.
 2. splitting multiple alts
 3. normalizing (trimming and left-aligning) the variants.

Where 2 and 3 should greatly reduce false negatives due to incorrect annotations. If you have a common error in a VCF, please open an issue so that we can address it in `vcftidy`.

Use
---

  $ python vcftidy.py $VCF $REFERENCE\_FASTA > $TIDY\_VCF

See Also
--------

	[vt](https://github.com/atks/vt) and associated [paper](http://t.co/J48Irh4wEe) does a nice job *decomposing* and *normalizing* variants.
