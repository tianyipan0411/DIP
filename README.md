## DIP
Experiments for ICDM Paper "Dynamic Propagation Rates: New Dimension to Viral Marketing in Online Social Networks"


# Compiling
Run make in the root folder
# To create the datasets:
1. Download desired network files from snap.stanford.edu in txt format and place in the folder "Data Processing"
2. Run processing.py, which can renumber the node indices and generate directed/undirected, weighted/unweighted graphs as needed.
3. Use ./el2bin [input file] [output file] to convert txt to bin
[input file]: the path to text file in edge list format: the first line contains the number of nodes n and number of edges m, each of the next m lines describes an edge following the format: [src] [dest] [weight]. Node index starts from 1.
[output file]: the path to binary output file
    
# Solve the DIP problem
Run ./dip -i [input file] -m [model name] -t [num threads] -s [speedup rate] -sr [speedup threshold]  -cr [required threshold] -o [output file]
[input file]: path to the .bin file of the network
[model name]: either "IC" or "LT"
[num threads]: number of threads for the program
[speedup rate]: how many times the propagation speed increases. Take values in {1,2,3,4,5,6}. The values are indices. Refer to the paper for actual values.
[speedup threshold]: number of activated users to trigger the speedup. Take values in {1,2,3,4,5}. The values are indices. Refer to the paper for actual values, which varies for different networks.
[required threshold]: number of activated users as the final requirement. Take values in {1,2,3,4,5,6}. The values are indices. Refer to the paper for actual values, which varies for different networks.
[output file]: the path to the output file that records the seed set, running time and other information of interest

The results will be in the folder Result/[network name]

# Verify the seed set
Run ./verify -i [result file] -o [folder name]
[result file]: path to the output file of DIP. 
[folder name]: folder to store the output. The path will be Final/[folder name]/[file corresponding to the result file]
   


    
