#!/usr/bin/python

import sys, getopt, codecs
import locale, os
import glob

def main(argv):
    if_directed = 0
    if_weighted = 0
    try:
        opts, args = getopt.getopt(argv,"hd:w:",["ifdirected=","ifrequireweight="])
    except getopt.GetoptError:
        print ('processing.py -d <0 for undirected> -w <0 for not writing weights>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('processing.py -d <0 for undirected> -w <0 for not writing weights>')
            sys.exit()
        elif opt in ("-d", "--ifdirected"):
            if_directed = arg
        elif opt in ("-w", "--ifrequireweight"):
            if_weighted = arg
    print ('Directed:"', if_directed)
    print ('Weighted "', if_weighted)
	
    cur_directory = os.getcwd()
    file_list = glob.glob(cur_directory + "\*.txt")
    print ('Processing data files in ' + cur_directory)
    
    for file_path in file_list:
        file_name = os.path.basename(file_path)  
        file_name = "w_" + file_name
        file_name_d = "d_" + file_name
        output_path_d = cur_directory + "\processed" + "\\" + file_name_d        
        output_path = cur_directory + "\processed" + "\\" + file_name
        
        output_file_d = open(output_path_d, 'w', newline='')
        output_file = open(output_path, 'w', newline='')
        
        dic = {}
        count_d = [0]
        count = [0]
        num_edges = 0
        with open(file_path) as f:
            for line in f:
                edge = line.split()
                if edge[0] != edge[1]:
                    num_edges += 1
                    if edge[0] in dic:
                        count[dic[edge[0]]] += 1
                    else:
                        dic[edge[0]] = len(dic) + 1
                        count.append(1)
                        count_d.append(0)
                    if edge[1] in dic: 
                        count[dic[edge[0]]] += 1
                        count_d[dic[edge[1]]] += 1
                    else:
                        dic[edge[1]] = len(dic) + 1
                        count.append(1)
                        count_d.append(1)
            print(str(len(dic)) + " " + str(num_edges), file = output_file_d)
            print(str(len(dic)) + " " + str(num_edges * 2), file = output_file)
        with open(file_path) as f:
            for line in f:
                edge = line.split()
                if edge[0] != edge[1]:
                    edge_output = str(dic[edge[0]]) + " " + str(dic[edge[1]]) + " "
                    edge_output_1 = str(dic[edge[1]]) + " " + str(dic[edge[0]]) + " " + str(1/float(count[dic[edge[0]]]))
                    edge_output_d = edge_output + str(1/float(count_d[dic[edge[1]]]))
                    edge_output = edge_output + str(1/float(count[dic[edge[1]]]))
                    print(edge_output_d , file = output_file_d)
                    print(edge_output , file = output_file)
                    print(edge_output_1 , file = output_file)
            output_file.close()
            output_file_d.close()
                
if __name__ == "__main__":
   main(sys.argv[1:])