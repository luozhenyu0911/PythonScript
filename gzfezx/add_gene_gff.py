import re
import sys
from collections import defaultdict
 
"""
python add_gene_gff.py Beauregard_v3.hc_gene_models.gff3 > Beauregard_v3.hc_gene_models_add_geneid.gff3
"""
# 读入文件
# 第一步
# 将每个mRNA的起始终止坐标都写入以geneid作为key的字典中
# 以便在第二步中提取最小和最大的坐标作为geneid最长的坐标
gene_dict=defaultdict(list)
with open(sys.argv[1], 'r') as f:
    gene_id = ''
    # 遍历每一行
    for line in f:
        if line.startswith('##'):
            # 忽略以##开头的注释行
            continue
        elif line.startswith('#'):
            # 对于以#开头的注释行，直接输出
            print(line.strip())
        else:
            # 对于数据行
            fields = line.strip().split('\t')
            # 如果是mRNA行
            if fields[2] == 'mRNA':
                # 获取mRNA ID
                mrna_id = re.findall(r'ID=([^;]+)', fields[8])[0]
                # 生成新的基因ID
                # gene_id = re.sub(r'\.\d+', '', mrna_id)
                gene_id = mrna_id.split('.')[0]+"."+mrna_id.split('.')[1]+"."+mrna_id.split('.')[2]
                gene_dict[gene_id].append(int(fields[3]))
                gene_dict[gene_id].append(int(fields[4]))

# 第二步
with open(sys.argv[1], 'r') as f:
    #lines = f.readlines()
    # 初始化基因ID
    gene_id = ''
    # 遍历每一行
    for line in f:
        if line.startswith('##'):
            # 忽略以##开头的注释行
            continue
        elif line.startswith('#'):
            # 对于以#开头的注释行，直接输出
            print(line.strip())
        else:
            # 对于数据行
            fields = line.strip().split('\t')
            # 如果是mRNA行
            if fields[2] == 'mRNA':
                # 获取mRNA ID
                mrna_id = re.findall(r'ID=([^;]+)', fields[8])[0]
                # 生成新的基因ID
                # gene_id = re.sub(r'\.\d+', '', mrna_id)
                gene_id = mrna_id.split('.')[0]+"."+mrna_id.split('.')[1]+"."+mrna_id.split('.')[2]
                if gene_id in gene_dict:
                    # 添加gene行
                    gene_line = [fields[0], fields[1], 'gene', str(min(gene_dict[gene_id])), str(max(gene_dict[gene_id])), fields[5], fields[6], fields[7], f'ID={gene_id};Name={gene_id}']
                    print('\t'.join(gene_line))
                    # 修改mRNA行
                    mrna_line = fields
                    mrna_line[2] = 'mRNA'
                    mrna_line[8] = f'ID={mrna_id};Name={mrna_id};Parent={gene_id};'
                    # mrna_line[8] = re.sub(r';source_id=[^;]+', '', mrna_line[8])
                    print('\t'.join(mrna_line))
                    del gene_dict[gene_id]
                else:
                    mrna_line = fields
                    mrna_line[2] = 'mRNA'
                    mrna_line[8] = f'ID={mrna_id};Name={mrna_id};Parent={gene_id};'
                    # mrna_line[8] = re.sub(r';source_id=[^;]+', '', mrna_line[8])
                    print('\t'.join(mrna_line))
            else:
                # 对于其他行直接输出
                print(line.strip())
