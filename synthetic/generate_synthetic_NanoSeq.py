import sys
from pysam import VariantFile as VCF
from pysam import AlignmentFile as BAM
import random
import copy
from gzip import open as gzopen
import statistics

mutdonor_fq1 = sys.argv[1]
mutdonor_fq2 = sys.argv[2]
seqdonor_fq1 = sys.argv[3]
seqdonor_fq2 = sys.argv[4]
mutdonor_bam = sys.argv[5]
seqdonor_bam = sys.argv[6]
muts = sys.argv[7]
out_prefix = sys.argv[8]
truth_vcf = sys.argv[9]
out_fq1 = sys.argv[10]
out_fq2 = sys.argv[11]
mut_num = int(sys.argv[12])
normal_bam = sys.argv[13]
vcf = VCF(muts,'rb')
vcf.header.add_meta('INFO', items=[('ID',"RF"), ('Number',1), ('Type','Integer'), ('Description','read family size')])
vcf.header.add_meta('INFO', items=[('ID',"RT"), ('Number',1), ('Type','Integer'), ('Description','read family size top strand')])
vcf.header.add_meta('INFO', items=[('ID',"RB"), ('Number',1), ('Type','Integer'), ('Description','read family size bottom strand')])
out_vcf = VCF(truth_vcf,'w',header=vcf.header)
add_readnames = set()
remove_readnames = set()
mbam = BAM(mutdonor_bam,'rb')
sbam = BAM(seqdonor_bam,'rb')
nbam = BAM(normal_bam,'rb')
rec_num = 0
variant_list = list()
for rec in vcf.fetch():
    variant_list.append(rec)
while rec_num < mut_num:
    #rec_num += 1
    if rec_num%1000==0: print(rec_num)
    rng = random.randint(0,len(variant_list))
    rec = variant_list.pop(rng)
    if len(rec.alts) >1: continue
    alt = rec.alts[0]
    chrom = rec.chrom
    pos = rec.pos
    ref = rec.ref
    alt_read_barcodes = {}
    rf_size = 0
    rf_top = 0
    rf_bot = 0
    if len(ref) == 1 and len(alt) == 1:
        for pileupcolumn in mbam.pileup(
            chrom, pos - 1, pos, min_base_quality=1, flag_filter=2816
        ):
            if pileupcolumn.pos == pos - 1:
                break_flag = 0
                for pileupread in pileupcolumn.pileups:
                    if (
                        pileupread.is_refskip
                        or pileupread.alignment.is_secondary
                        or pileupread.alignment.is_supplementary
                        or pileupread.alignment.mapping_quality <= 30
                    ):
                        continue
                    barcodes = pileupread.alignment.query_name.split("_")[-1]
                    if (not pileupread.is_del and 
                        pileupread.alignment.query_sequence[pileupread.query_position]
                        == alt):
                        if not alt_read_barcodes.get(barcodes):
                            alt_read_barcodes[barcodes] = [pileupread.alignment]
                        else:
                            alt_read_barcodes[barcodes].append(pileupread.alignment)
                alt_read_barcodes_combined = dict()
                alt_read_barcodes_duplex = dict()
                for bc in alt_read_barcodes.keys():
                    bc1 = bc.split("+")[0]
                    bc2 = bc.split("+")[1]
                    if not alt_read_barcodes_combined.get(bc2+"+"+bc1):
                        alt_read_barcodes_combined[bc] = alt_read_barcodes[bc]
                    else:
                        alt_read_barcodes_duplex[bc] = alt_read_barcodes_combined[bc2+"+"+bc1] + alt_read_barcodes[bc]
                if len(alt_read_barcodes_duplex) == 0: continue
                break_flag = 0
                for rf in alt_read_barcodes_duplex.values():
                    for seq in rf:
                        if seq.query_name.split("_")[0]  in add_readnames:
                            break_flag = 1
                            break
                if not break_flag:
                    random_read_family_bc = random.choice(list(alt_read_barcodes_duplex.keys()))
                    random_read_family = alt_read_barcodes_duplex[random_read_family_bc]
                    random_read_family_starts = statistics.mode([_.reference_start for _ in random_read_family])
                    random_rf_bc1 = random_read_family_bc.split("+")[0]
                    random_rf_bc2 = random_read_family_bc.split("+")[1]
                    random_rf_names = list()
                    for pileupread in pileupcolumn.pileups:
                        if (
                            pileupread.is_del
                            or pileupread.is_refskip
                            or pileupread.alignment.is_secondary
                            or pileupread.alignment.is_supplementary
                            or pileupread.alignment.mapping_quality <= 30
                        ):
                            continue
                        barcodes = pileupread.alignment.query_name.split("_")[-1]
                        barcodes_set = set()
                        if barcodes == random_read_family_bc and pileupread.alignment.reference_start == random_read_family_starts:
                            random_rf_names.append(pileupread.alignment.query_name.split("_")[0])
                            rf_top += 1
                        elif barcodes == random_rf_bc2 + "+" + random_rf_bc1 and pileupread.alignment.reference_start == random_read_family_starts:
                            random_rf_names.append(pileupread.alignment.query_name.split("_")[0])
                            rf_bot += 1

                    rf_size = rf_top + rf_bot
                    add_readnames.update(set(random_rf_names))
            
    else:
        indel_size = len(alt) - len(ref)
        continue_flag = 0
        for pileupcolumn in nbam.pileup(
            chrom, pos - 1, pos, min_base_quality=30, flag_filter=2816
        ):
            if pileupcolumn.pos == pos - 1:
                for pileupread in pileupcolumn.pileups:
                    if pileupread.indel != 0:
                        continue_flag = 1
                        break
            if continue_flag:
                break
        if continue_flag: continue
        for pileupcolumn in mbam.pileup(
            chrom, pos - 1, pos, min_base_quality=1, flag_filter=2816
        ):
            if pileupcolumn.pos == pos - 1:
                break_flag = 0
                for pileupread in pileupcolumn.pileups:
                    if (
                        pileupread.is_del
                        or pileupread.is_refskip
                        or pileupread.alignment.is_secondary
                        or pileupread.alignment.is_supplementary
                        or pileupread.alignment.mapping_quality <= 30
                    ):
                        continue
                    if pileupread.indel != indel_size and pileupread.indel != 0:
                        continue_flag = 1
                        break
                    barcodes = pileupread.alignment.query_name.split("_")[-1]
                    if pileupread.indel == indel_size:
                        if not alt_read_barcodes.get(barcodes):
                            alt_read_barcodes[barcodes] = [pileupread.alignment]
                        else:
                            alt_read_barcodes[barcodes].append(pileupread.alignment)
                    else:
                        continue
                if continue_flag:
                    break
                alt_read_barcodes_combined = dict()
                alt_read_barcodes_duplex = dict()
                for bc in alt_read_barcodes.keys():
                    bc1 = bc.split("+")[0]
                    bc2 = bc.split("+")[1]
                    if not alt_read_barcodes_combined.get(bc2+"+"+bc1):
                        alt_read_barcodes_combined[bc] = alt_read_barcodes[bc]
                    else:
                        alt_read_barcodes_duplex[bc] = alt_read_barcodes_combined[bc2+"+"+bc1] + alt_read_barcodes[bc]
                if len(alt_read_barcodes_duplex) == 0: continue
                break_flag = 0
                for rf in alt_read_barcodes_duplex.values():
                    for seq in rf:
                        if seq.query_name.split("_")[0] in add_readnames:
                            break_flag = 1
                            break
                if not break_flag:
                    random_read_family_bc = random.choice(list(alt_read_barcodes_duplex.keys()))
                    random_read_family = alt_read_barcodes_duplex[random_read_family_bc]
                    random_read_family_starts = statistics.mode([_.reference_start for _ in random_read_family])
                    random_rf_bc1 = random_read_family_bc.split("+")[0]
                    random_rf_bc2 = random_read_family_bc.split("+")[1]
                    random_rf_names = list()
                    for pileupread in pileupcolumn.pileups:
                        if (
                            pileupread.is_del
                            or pileupread.is_refskip
                            or pileupread.alignment.is_secondary
                            or pileupread.alignment.is_supplementary
                            or pileupread.alignment.mapping_quality <= 30
                        ):
                            continue
                        barcodes = pileupread.alignment.query_name.split("_")[-1]
                        if barcodes == random_read_family_bc and pileupread.alignment.reference_start == random_read_family_starts:
                            random_rf_names.append(pileupread.alignment.query_name.split("_")[0])
                            rf_top += 1
                        elif barcodes == random_rf_bc2 + "+" + random_rf_bc1 and pileupread.alignment.reference_start == random_read_family_starts:
                            random_rf_names.append(pileupread.alignment.query_name.split("_")[0])
                            rf_bot += 1
                    rf_size = rf_top + rf_bot

                    add_readnames.update(set(random_rf_names))
        if continue_flag: continue
    if rf_size == 0: continue
    for pileupcolumn in sbam.pileup(
                    chrom, pos - 1, pos, min_base_quality=1, flag_filter=2816
                ):
        if pileupcolumn.pos == pos - 1:
            readnames = [pileupread.alignment.query_name.split("_")[0] for pileupread in pileupcolumn.pileups]
            if len(readnames) >= rf_size:
                random_read_names = random.sample(readnames,k=rf_size)
                rec.info.__setitem__('RF', rf_size)
                rec.info.__setitem__('RT', rf_top)
                rec.info.__setitem__('RB', rf_bot)
                out_vcf.write(rec)
                remove_readnames.update(set(random_read_names))
                rec_num += 1
add_reads = 0
remove_reads = 0
print(1)
if seqdonor_fq1[-3:] == ".gz":
    sfq1 = gzopen(seqdonor_fq1, "rt")
    sfq2 = gzopen(seqdonor_fq2, "rt")
else:
    sfq1 = open(seqdonor_fq1)
    sfq2 = open(seqdonor_fq1)
if mutdonor_fq1[-3:] == ".gz":
    mfq1 = gzopen(mutdonor_fq1, "rt")
    mfq2 = gzopen(mutdonor_fq2, "rt")
else:
    mfq1 = open(mutdonor_fq1)
    mfq2 = open(mutdonor_fq1)
fq1Out = open(out_fq1,"w")
fq2Out = open(out_fq2,"w")
lineIndex = 0
for line1, line2 in zip(mfq1, mfq2):
    if lineIndex == 0:
        name1 = line1
        name2 = line2
        lineIndex = 1
    elif lineIndex == 1:
        seq1 = line1
        seq2 = line2
        lineIndex = 2
    elif lineIndex == 2:
        lineIndex = 3
    else:
        qual1 = line1
        qual2 = line2
        lineIndex = 0
        if name1.strip("\n").strip("@") in add_readnames:
        # print(read1,read2)
            fq1Out.write(name1 + seq1 + "+\n" + qual1)
            fq2Out.write(name2 + seq2 + "+\n" + qual2)
            add_reads += 1
lineIndex = 0
for line1, line2 in zip(sfq1, sfq2):
    if lineIndex == 0:
        name1 = line1
        name2 = line2
        lineIndex = 1
    elif lineIndex == 1:
        seq1 = line1
        seq2 = line2
        lineIndex = 2
    elif lineIndex == 2:
        lineIndex = 3
    else:
        qual1 = line1
        qual2 = line2
        lineIndex = 0
        if not name1.strip("\n").strip("@") in remove_readnames:
        # print(read1,read2)
            fq1Out.write(name1 + seq1 + "+\n" + qual1)
            fq2Out.write(name2 + seq2 + "+\n" + qual2)
            remove_reads += 1


