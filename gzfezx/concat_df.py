import pandas as pd
import glob
import argparse
example_text = '''example:
   python 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--glob_xx','-g',type=str,help="",required= True,metavar='')
parser.add_argument('--columns', '-c',type= int,help="how many columns to keep",required= True,metavar='')
parser.add_argument('--prefix', '-p',type= str,help="",required= True,metavar='')
args = parser.parse_args()

glob_xx = args.glob_xx
prefix = args.prefix
columns = args.columns

file_list = glob.glob("*"+glob_xx)
# print(len(argv))
# print(argv[4])
# 创建一个空的DataFrame作为初始值
result_df = pd.DataFrame()
num=[]
for i in range(int(columns)):
    num.append(i)
# 循环对多个文件两两合并
for filename in file_list:
    df1 = pd.read_csv(filename, sep='\t', header=None, usecols=num)
    name = filename.split('.')[0]
    df1['sample'] = name
    df1.columns = ['region', 'count', 'sample']
    df1['percentages'] = df1['count'] / df1['count'].sum()
    df1['percentages'] = df1['percentages'].apply(lambda x: "{:.1%}".format(x))
    if result_df.empty:
        result_df = df1
    else:
        result_df = pd.concat([result_df, df1], ignore_index=True)
        
result_df.to_csv(prefix, sep='\t', index=False, header=True)