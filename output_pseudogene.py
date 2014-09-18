#!/usr/bin/python
import sys

# add frameshift information to pseudogene informaiton.
def output_pseudogene(in_file_name, frameshift_dict, out_file_name):
    in_file = open(in_file_name, 'Ur')
    out_file = open(out_file_name, 'w')    
    for line in in_file:
        items = line.rstrip().split()
        key = '_'.join([items[0], items[2], items[3]]) 
        if key in frameshift_dict:
            out_file.write(line.rstrip()+' '+frameshift_dict[key]+'\n')
        else:
           out_file.write(line.rstrip()+' 0\n')
    in_file.close()
    out_file.close()

# get frameshift information.
def get_frameshift_dict(in_file_name):
    frameshift_dict = {}
    in_file = open(in_file_name, 'Ur')
    for line in in_file:
        if line[0] == '>':
            items = line.rstrip().split()
            key = items[0][1:] 
            frameshift_num_str = items[3].split('=')[-1]
            frameshift_dict[key] = frameshift_num_str
    in_file.close()
    return frameshift_dict

def main():
    if len(sys.argv) < 4:
        sys.stderr.write('Usage: %s <pseudogene file> <hmmframe file>'
                         ' <output file>\n' % (sys.argv[0],))
        sys.exit(1)
    pseudogene_file_name = sys.argv[1]
    hmmframe_file_name = sys.argv[2]
    output_file_name = sys.argv[3]
    frameshift_dict = get_frameshift_dict(hmmframe_file_name)
    output_pseudogene(pseudogene_file_name, frameshift_dict, output_file_name)

if __name__ == '__main__':
    main()
