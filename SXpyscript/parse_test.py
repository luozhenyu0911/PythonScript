import argparse

parser = argparse.ArgumentParser()
parser.description='According to a given number to extract reads that shorter than this number.'
parser.add_argument('-l','--list', nargs='+', help='<Required> Set flag',type=int,required=True)
args = parser.parse_args()

if args.list:
    print type(args.list[0])

