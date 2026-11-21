from pysam import VariantFile as VCF
from pysam import AlignmentFile as BAM
import sys
def countUniqueBarcodes(barcodes):
    if len(barcodes) == 0:
        return 0
    bc_starting_list = list(barcodes.keys())
    for bc in bc_starting_list:
        if not barcodes.get(bc.split("+")[1] + "+" + bc.split("+")[0]):
            barcodes[bc.split("+")[1] + "+" + bc.split("+")[0]] = 1
    bcs = sorted(list(barcodes.keys()))
    take_ind = [0]
    for nn in range(1, len(bcs)):
        string_diff = [a == b for a, b in zip(list(bcs[nn]), list(bcs[nn - 1]))].count(
            False
        )
        if string_diff > 1:
            take_ind.append(nn)
    bcs_uniq = [bcs[_] for _ in take_ind]
    count = 0
    bc_dict = dict()
    for bc in bcs_uniq:
        if not bc_dict.get(bc.split("+")[1] + "+" + bc.split("+")[0]):
            bc_dict[bc] = 1
            count += 1
    return count

def extractDepthSnv(bam, chrom, pos, ref, alt, params):
    has_barcode = False
    for seq in bam.fetch(chrom, pos - 1, pos):
        if len(seq.query_name.split("_")[-1].split("+")) == 2:
            has_barcode = True
        break
    altAlleleCount = 0
    refAlleleCount = 0
    otherAlleleCount = 0
    for pileupcolumn in bam.pileup(chrom, pos - 1, pos, min_base_quality=30):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if (
                    pileupread.is_refskip
                    or pileupread.alignment.is_secondary
                    or pileupread.alignment.is_supplementary
                    or pileupread.alignment.is_duplicate
                    or pileupread.alignment.mapping_quality <= params["mapq"]
                ):
                    continue
                if pileupread.is_del:
                    otherAlleleCount += 1
                elif (
                    pileupread.alignment.query_sequence[pileupread.query_position]
                    == alt
                ):
                    altAlleleCount += 1
                elif (
                    pileupread.alignment.query_sequence[pileupread.query_position]
                    == ref
                ):
                    refAlleleCount += 1
                else:
                    otherAlleleCount += 1
    depth = refAlleleCount + altAlleleCount + otherAlleleCount
    return altAlleleCount, refAlleleCount, depth
def extractDepthIndel(bam, chrom, pos, ref, alt, params):
    indel_size = len(alt) - len(ref)
    altAlleleCount = 0
    refAlleleCount = 0
    otherAllelCount = 0
    for pileupcolumn in bam.pileup(chrom, pos - 1, pos, min_base_quality=30):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if (
                    pileupread.is_del
                    or pileupread.is_refskip
                    or pileupread.alignment.is_secondary
                    or pileupread.alignment.is_supplementary
                    or pileupread.alignment.is_duplicate
                    or pileupread.alignment.mapping_quality <= params["mapq"]
                ):
                    continue
                if pileupread.indel == indel_size:
                    altAlleleCount += 1
                elif pileupread.indel == 0:
                    refAlleleCount += 1
                else:
                    otherAllelCount += 1
    depth = refAlleleCount + altAlleleCount + otherAllelCount
    return altAlleleCount, refAlleleCount, depth

normalbam = BAM(sys.argv[2],'rb')
samplebam = BAM(sys.argv[3],'rb')
vcf = VCF(sys.argv[1],'r')
vcf_out = VCF(sys.argv[4],'w',header=vcf.header)
params = {"mapq":30}
last_chrom = "0"
for rec in vcf.fetch():
    if rec.chrom != last_chrom:
        print(rec)
        last_chrom = rec.chrom
    if len(rec.alts) >=2: continue
    if len(rec.ref) == 1 and len(rec.alts[0]) == 1:
        na,nr,ndp = extractDepthSnv(normalbam,rec.chrom,rec.pos,rec.ref,rec.alts[0],params)
        ta,tr,tdp = extractDepthSnv(samplebam,rec.chrom,rec.pos,rec.ref,rec.alts[0],params)
    else:
        na,nr,ndp = extractDepthIndel(normalbam,rec.chrom,rec.pos,rec.ref,rec.alts[0],params)
        ta,tr,tdp = extractDepthIndel(samplebam,rec.chrom,rec.pos,rec.ref,rec.alts[0],params)
    if ndp >= 10 and na == 0 and ta >= 1:
        vcf_out.write(rec)
        

