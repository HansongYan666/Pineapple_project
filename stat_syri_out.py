import sys,os


stats = {'SYN': [0, 0, 0],
             'INV': [0, 0, 0],
             'TRANS': [0, 0, 0],
             'DUPA': [0, 0, 0],
             'DUPB': [0, 0, 0],
             'NOTALR': [0, 0, 0],
             'NOTALQ': [0, 0, 0],
             'SNP': [0, 0, 0],
             'LargeINS': [0, 0, 0],
             'LargeDEL': [0, 0, 0],
             'SmallINS': [0, 0, 0],
             'SmallDEL': [0, 0, 0],
             'CPG': [0, 0, 0],
             'CPL': [0, 0, 0],
             'HDR': [0, 0, 0],
             'TDM': [0, 0, 0]}


with open(sys.argv[1], 'r') as fin:
    for line in fin:
        line = line.strip().split('\t')
        if line[10] == 'SYN':
            stats['SYN'][0] += 1
            stats['SYN'][1] += abs(int(line[1]) - int(line[2])) + 1
            stats['SYN'][2] += abs(int(line[6]) - int(line[7])) + 1
        if line[10] == 'INV':
            stats['INV'][0] += 1
            stats['INV'][1] += abs(int(line[1]) - int(line[2])) + 1
            stats['INV'][2] += abs(int(line[6]) - int(line[7])) + 1
        if line[10] in ['TRANS', 'INVTR']:
            stats['TRANS'][0] += 1
            stats['TRANS'][1] += abs(int(line[1]) - int(line[2])) + 1
            stats['TRANS'][2] += abs(int(line[6]) - int(line[7])) + 1
        ## Discuss DUP length stat
        if line[10] in ['DUP', 'INVDP']:
            if line[11] == 'copygain':
                stats['DUPB'][0] += 1
                stats['DUPB'][2] += abs(int(line[6]) - int(line[7])) + 1
            if line[11] == 'copyloss':
                stats['DUPA'][0] += 1
                stats['DUPA'][1] += abs(int(line[1]) - int(line[2])) + 1
        if line[10] == 'NOTAL':
            try:
                stats['NOTALR'][1] += abs(int(line[1]) - int(line[2])) + 1
                stats['NOTALR'][0] += 1
            except ValueError:
                stats['NOTALQ'][2] += abs(int(line[6]) - int(line[7])) + 1
                stats['NOTALQ'][0] += 1
        if line[10] == 'SNP':
            stats['SNP'][0] += 1
            stats['SNP'][1] += abs(int(line[1]) - int(line[2])) + 1
            stats['SNP'][2] += abs(int(line[6]) - int(line[7])) + 1
        if line[10] == 'INS':
            length = abs(int(line[6]) - int(line[7]))
            if length > 50:
                stats['LargeINS'][0] += 1
                stats['LargeINS'][2] += abs(int(line[6]) - int(line[7]))
            else:
                stats['SmallINS'][0] += 1
                stats['SmallINS'][2] += abs(int(line[6]) - int(line[7]))
        if line[10] == 'DEL':
            length = abs(int(line[1]) - int(line[2]))
            if length > 50:
                stats['LargeDEL'][0] += 1
                stats['LargeDEL'][1] += abs(int(line[1]) - int(line[2]))
            else:
                stats['SmallDEL'][0] += 1
                stats['SmallDEL'][1] += abs(int(line[1]) - int(line[2]))
        # Discuss Copychange statistics
        if line[10] == 'CPG':
            stats['CPG'][0] += 1
            # stats['CPG'][1] += abs(int(line[1]) - int(line[2])) + 1
            stats['CPG'][2] += abs(int(line[6]) - int(line[7])) + 1
        if line[10] == 'CPL':
            stats['CPL'][0] += 1
            stats['CPL'][1] += abs(int(line[1]) - int(line[2])) + 1
            # stats['CPL'][2] += abs(int(line[6]) - int(line[7])) + 1
        if line[10] == 'HDR':
            stats['HDR'][0] += 1
            stats['HDR'][1] += abs(int(line[1]) - int(line[2])) + 1
            stats['HDR'][2] += abs(int(line[6]) - int(line[7])) + 1
        if line[10] == 'TDM':
            stats['TDM'][0] += 1
            stats['TDM'][1] += abs(int(line[1]) - int(line[2])) + 1
            stats['TDM'][2] += abs(int(line[6]) - int(line[7])) + 1


with open(sys.argv[2], 'w') as fout:
    fout.write('#Structural annotations\n')
    fout.write('#Variation_type\tCount\tLength_ref\tLength_qry\n')
    fout.write('{}\t{}\t{}\t{}\n'.format('Syntenic regions', stats['SYN'][0], stats['SYN'][1], stats['SYN'][2]))
    fout.write('{}\t{}\t{}\t{}\n'.format('Inversions', stats['INV'][0], stats['INV'][1], stats['INV'][2]))
    fout.write('{}\t{}\t{}\t{}\n'.format('Translocations', stats['TRANS'][0], stats['TRANS'][1], stats['TRANS'][2]))
    fout.write('{}\t{}\t{}\t{}\n'.format('Duplications (reference)', stats['DUPA'][0], stats['DUPA'][1], '-'))
    fout.write('{}\t{}\t{}\t{}\n'.format('Duplications (query)', stats['DUPB'][0], '-', stats['DUPB'][2]))
    fout.write('{}\t{}\t{}\t{}\n'.format('Not aligned (reference)', stats['NOTALR'][0], stats['NOTALR'][1], '-'))
    fout.write('{}\t{}\t{}\t{}\n'.format('Not aligned (query)', stats['NOTALQ'][0], '-', stats['NOTALQ'][2]))
    fout.write('\n\n#Sequence annotations\n')
    fout.write('#Variation_type\tCount\tLength_ref\tLength_qry\n')
    fout.write('{}\t{}\t{}\t{}\n'.format('SNPs', stats['SNP'][0], stats['SNP'][1], stats['SNP'][2]))
    fout.write('{}\t{}\t{}\t{}\n'.format('Large Insertions', stats['LargeINS'][0], '-', stats['LargeINS'][2]))
    fout.write('{}\t{}\t{}\t{}\n'.format('Large Deletions', stats['LargeDEL'][0], stats['LargeDEL'][1], '-'))
    fout.write('{}\t{}\t{}\t{}\n'.format('Small Insertions', stats['SmallINS'][0], '-', stats['SmallINS'][2]))
    fout.write('{}\t{}\t{}\t{}\n'.format('Small Deletions', stats['SmallDEL'][0], stats['SmallDEL'][1], '-'))
    fout.write('{}\t{}\t{}\t{}\n'.format('Copygains', stats['CPG'][0], '-', stats['CPG'][2]))
    fout.write('{}\t{}\t{}\t{}\n'.format('Copylosses', stats['CPL'][0], stats['CPL'][1], '-'))
    fout.write('{}\t{}\t{}\t{}\n'.format('Highly diverged', stats['HDR'][0], stats['HDR'][1], stats['HDR'][2]))
    fout.write('{}\t{}\t{}\t{}\n'.format('Tandem repeats', stats['TDM'][0], stats['TDM'][1], stats['TDM'][2]))

# def parser_vcf_info(info):
#     dic = {}
#     for i in info.strip().split(";"):
#         tmp = i.split("=")
#         dic[tmp[0]] = tmp[1]
#     return dic
#
# n = 0
# for i in open(sys.argv[1]):
#     if i.startswith('#'):
#         continue
#     #Chr01	9801	SNP18419	A	T	.	PASS	END=9801;ChrB=Chr01;StartB=13518;EndB=13518;Parent=SYN1;VarType=ShV;DupType=.
#     tmp = i.strip("\n").split("\t")
#     id = tmp[2]
#     info = tmp[7]
#     dic = parser_vcf_info(info)
#     if "Parent" in dic.keys():
#         parent = dic["Parent"]
#     else:
#         print(i.strip())
#         continue
#     # if parent != ".":
#     #     continue
#     if parent.startswith("SYN"):
#         continue
#     if tmp[4] == "<TRANS>" or tmp[4] == "<INVTR>":
#     # if "TRAN" in tmp[4]:
#         n+=1
# print(n)
# stats = {'SYN': [0, 0, 0],
#          'INV': [0, 0, 0],
#          'TRANS': [0, 0, 0],
#          'DUPA': [0, 0, 0],
#          'DUPB': [0, 0, 0],
#          'NOTALR': [0, 0, 0],
#          'NOTALQ': [0, 0, 0],
#          'SNP': [0, 0, 0],
#          'INS': [0, 0, 0],
#          'DEL': [0, 0, 0],
#          'CPG': [0, 0, 0],
#          'CPL': [0, 0, 0],
#          'HDR': [0, 0, 0],
#          'TDM': [0, 0, 0]}
#
# for line in open(sys.argv[1]):
#     if line.startswith("#"):
#         continue
#     line = line.strip().split('\t')
#     id = line[2]
#     info = line[7]
#     svtype = line[4]
#     svtype = svtype.strip(">").strip("<")
#     dic = parser_vcf_info(info)
#     if "Parent" in dic.keys():
#         parent = dic["Parent"]
#     else:
#         continue
#     #Chr01	11546	DEL18463	CTT	C	.	PASS	END=11548;ChrB=Chr01;StartB=15262;EndB=15262;Parent=SYN1;VarType=ShV;DupType=.
#     start1 = int( line[1])
#     end1 = int(dic["END"])
#     try:
#         start2 = int(dic["StartB"])
#         end2 = int(dic["EndB"])
#     except:
#         start2 = 0
#         end2 = 0
#     if svtype == 'SYN':
#         stats['SYN'][0] += 1
#         stats['SYN'][1] += abs(start1 - end1) + 1
#         stats['SYN'][2] += abs(start2 - end2) + 1
#     if svtype == 'INV':
#         stats['INV'][0] += 1
#         stats['INV'][1] += abs(start1 - end1) + 1
#         stats['INV'][2] += abs(start2 - end2) + 1
#
#     if svtype in ['TRANS', 'INVTR']:
#         stats['TRANS'][0] += 1
#         stats['TRANS'][1] += abs(start1 - end1) + 1
#         stats['TRANS'][2] += abs(start2 - end2) + 1
#
#     ## Discuss DUP length stat
#     if svtype in ['DUP', 'INVDP']:
#         duptype = dic["DupType"]
#         if duptype == 'copygain':
#             stats['DUPB'][0] += 1
#             stats['DUPB'][2] += abs(start2 - end2) + 1
#         if duptype == 'copyloss':
#             stats['DUPA'][0] += 1
#             stats['DUPA'][1] += abs(start1 - end1) + 1
#
#     if svtype == 'NOTAL':
#         try:
#             stats['NOTALR'][1] += abs(start1 - end1) + 1
#             stats['NOTALR'][0] += 1
#         except ValueError:
#             stats['NOTALQ'][2] += abs(start2 - end2) + 1
#             stats['NOTALQ'][0] += 1
#
#     if id.startswith("SNP"):
#         stats['SNP'][0] += 1
#         stats['SNP'][1] += abs(start1 - end1) + 1
#         stats['SNP'][2] += abs(start2 - end2) + 1
#
#     if id.startswith("INS"):
#         stats['INS'][0] += 1
#         length = abs(start2 - end2)
#         if length > 50:
#             stats['INS'][0] += 1
#             stats['INS'][2] += abs(start2 - end2)
#         else:
#             stats['SNP'][0] += 1
#             # stats['INS'][2] += abs(start2 - end2)
#     if id.startswith("DEL"):
#         length = abs(start2 - end2)
#         if length > 50:
#             stats['DEL'][0] += 1
#             stats['DEL'][1] += abs(start1 - end1)
#         else:
#             stats['SNP'][0] += 1
#             stats['SNP'][1] += abs(start1 - end1)
#
#     # Discuss Copychange statistics
#     if svtype == 'CPG':
#         stats['CPG'][0] += 1
#         # stats['CPG'][1] += abs(start1 - end1) + 1
#         stats['CPG'][2] += abs(start2 - end2) + 1
#     if svtype == 'CPL':
#         stats['CPL'][0] += 1
#         stats['CPL'][1] += abs(start1 - end1) + 1
#         # stats['CPL'][2] += abs(start2 - end2) + 1
#
#     if svtype == 'HDR':
#         stats['HDR'][0] += 1
#         stats['HDR'][1] += abs(start1 - end1) + 1
#         stats['HDR'][2] += abs(start2 - end2) + 1
#     if svtype == 'TDM':
#         stats['TDM'][0] += 1
#         stats['TDM'][1] += abs(start1 - end1) + 1
#         stats['TDM'][2] += abs(start2 - end2) + 1
# fout = open(sys.argv[2], 'w')
# fout.write('#Structural annotations\n')
# fout.write('#Variation_type\tCount\tLength_ref\tLength_qry\n')
# fout.write('{}\t{}\t{}\t{}\n'.format('Syntenic regions', stats['SYN'][0], stats['SYN'][1], stats['SYN'][2]))
# fout.write('{}\t{}\t{}\t{}\n'.format('Inversions', stats['INV'][0], stats['INV'][1], stats['INV'][2]))
# fout.write('{}\t{}\t{}\t{}\n'.format('Translocations', stats['TRANS'][0], stats['TRANS'][1], stats['TRANS'][2]))
# fout.write('{}\t{}\t{}\t{}\n'.format('Duplications (reference)', stats['DUPA'][0], stats['DUPA'][1], '-'))
# fout.write('{}\t{}\t{}\t{}\n'.format('Duplications (query)', stats['DUPB'][0], '-', stats['DUPB'][2]))
# fout.write('{}\t{}\t{}\t{}\n'.format('Not aligned (reference)', stats['NOTALR'][0], stats['NOTALR'][1], '-'))
# fout.write('{}\t{}\t{}\t{}\n'.format('Not aligned (query)', stats['NOTALQ'][0], '-', stats['NOTALQ'][2]))
# fout.write('\n\n#Sequence annotations\n')
# fout.write('#Variation_type\tCount\tLength_ref\tLength_qry\n')
# fout.write('{}\t{}\t{}\t{}\n'.format('SNPs', stats['SNP'][0], stats['SNP'][1], stats['SNP'][2]))
# fout.write('{}\t{}\t{}\t{}\n'.format('Insertions', stats['INS'][0], '-', stats['INS'][2]))
# fout.write('{}\t{}\t{}\t{}\n'.format('Deletions', stats['DEL'][0], stats['DEL'][1], '-'))
# fout.write('{}\t{}\t{}\t{}\n'.format('Copygains', stats['CPG'][0], '-', stats['CPG'][2]))
# fout.write('{}\t{}\t{}\t{}\n'.format('Copylosses', stats['CPL'][0], stats['CPL'][1], '-'))
# fout.write('{}\t{}\t{}\t{}\n'.format('Highly diverged', stats['HDR'][0], stats['HDR'][1], stats['HDR'][2]))
# fout.write('{}\t{}\t{}\t{}\n'.format('Tandem repeats', stats['TDM'][0], stats['TDM'][1], stats['TDM'][2]))

