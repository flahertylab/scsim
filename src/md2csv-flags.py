import os, sys

def convert_to_SCs(flags):
    s_flags = []
    with open(flags, 'r') as f:
        for line in f:
            line = line[0:len(line)-1]
            if line[0] == '|':
                SCs = list(filter(None, line.split('|')))
                for sc in SCs:
                    s_flags.append([sc])
            else:
                SNVs = list(filter(None, line.split('|')))
                loc = int(SNVs.pop(0))
                for i, flag in enumerate(SNVs):
                    if flag == 'True':
                        flag = 1
                    else:
                        flag = 0
                    s_flags[i].append((loc, flag))
    return s_flags

def write_file(flags_list, REF_DIR):
    with open(os.path.join(REF_DIR, 'true_mutations_all_sc.csv'), 'w+') as f:
        f.write('sample,pos,snv_flag\n')
        for SC in flags_list:
            for flag in SC:
                if flag[0] is not 'S':
                    f.write('{},{},{}\n'.format(int(SC[0][2:]), flag[0], flag[1]))

if __name__ == '__main__':
    write_file(convert_to_SCs(sys.argv[1]), sys.argv[2])
