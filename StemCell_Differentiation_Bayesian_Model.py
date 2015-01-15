#dwest25, David West, Final Assignment
#Bayesian Modeling of Stem Cell Differentiation
import math

def main():
    infile_1 = open ("ipsc1.dat", "r")
    infile_2 = open ("ipsc2.dat", "r")
    outfile_1 = open("ipsc_results.dat", 'w')
    outfile_2 = open("ipsc_results2.dat", 'w')
    file_list = infile_1.read().split('\n') #creates a list of lines from file
    file_list2 = infile_2.read().split('\n') #creates a list of lines from file
    one_node_network(file_list)

    multi_node_network(file_list, 3, outfile_1)
    multi_node_network(file_list2, 3, outfile_2)
    

def one_node_network(f_list):
    print('Network with one parent node \n')
    file_list = f_list


    #loops through all Xn's (proteins)
    for Xn in range (6):
        off = count(Xn*2, '0', file_list)
        on = count(Xn*2, '1', file_list)
        n_values = [[off[0], off[1]], [on[0], on[1]]] #list of values with n0j and n1j values for X1
        log2P = compute_log2P(n_values)

        #table header
        print('\n\tX =', Xn+1, '\n\n**** j=0 ****')
        print('n0:', off[0],
              '\nn1j:', off[1],
              '\nP(N=0):', off[2],
              '\nP(N=1):', off[3])

        print('\n**** j=1 ****')
        print('n0:', on[0],
              '\nn1j:', on[1],
              '\nP(N=0):', on[2],
              '\nP(N=1):', on[3])

        dev = round(math.fabs(1373.52+log2P), 3) #rounded to 3 digits for simpler print
        perc_dev = round(dev/1373.52*100, 3)
        print('\nlog2(P) = ', log2P)
        print('Deviation = ', dev)
        print('Percent Deviation = ', perc_dev,'%')
        print('_______________________________')

#counts n0j and n1j for any single regulatory protein Xn in a state j
def count (Xn, j, file_list):
    n0j = 0
    n1j = 0
    for i in range (len(file_list)-1):
        x_value = file_list[i][Xn] #Xn values
        N = file_list[i][12]
        if (x_value == j):
            if (N == '0'):
                n0j += 1
            else:
                n1j += 1


    P_0 = round(n0j/(n0j+n1j), 3) #computes P(N=0) and rounds to 3 digits (for simpler print)
    P_1 = round(n1j/(n0j+n1j), 3) #computes P(N=1)
        
    return n0j, n1j, P_0, P_1

def multi_node_network(f_list,n, oFile):
    print("\n\nMulti-Node Network")
    file_list = f_list
    networks = []
    network_info = [] #info about network in corresponding index of networks[]. Format (n0j, n1j, occurences)
    probs0 = [] #list of all P(N=0) values
    probs1 = [] #list of all P(N=1) values
    data = []

    for i in range (len(file_list)-1):
        cur_network = file_list[i][:n*2-1]
        cur_N = file_list[i][12]
        if (cur_network in networks):
            index = networks.index(cur_network)
            network_info[index][2] +=1 #increments occurences
            if (cur_N == '0'):
                network_info[index][0] +=1 #increments n0j
            elif(cur_N == '1'):
                network_info[index][1] +=1 #increments n1j
                
        else:
            networks.append(cur_network)
            if (cur_N == '0'):
                network_info.append([1,0, 1]) #note: occurences is not vital for the problem except as a verification
            elif(cur_N == '1'):
                network_info.append([0,1, 1])

    print(compute_log2P(network_info))

    for i in range(len(networks)):
        P_0 = round(network_info[i][0]/(network_info[i][0]+network_info[i][1]), 3)
        P_1 = round(network_info[i][1]/(network_info[i][0]+network_info[i][1]), 3)
        probs0.append(P_0)
        probs1.append(P_1)
        data.append([networks[i],network_info[i][0], network_info[i][1], P_0, P_1])
        print(data[i], file = oFile)


   
#Sterling's approzimation
def s_approx (n):
    ans = (n+.5)*math.log(n,2)-n+.5*math.log(2*math.pi,2)
    return ans

#takes a list in format (n0j, n1j) as parameter
def compute_log2P(data):
    total = 0
    
    for i in range (len(data)):
        n0j = data[i][0]
        n1j = data[i][1]
        #uses approximation if n0j or n1j is large
        if (n0j>50 or n1j >50):
            num = s_approx(n0j) + s_approx(n1j)
            denom = s_approx(n0j+n1j)
            partial_sum =  num - denom
        else:
            partial_sum = math.log(math.factorial(n0j)*math.factorial(n1j)/math.factorial(n0j+n0j),2)
        total += partial_sum
    return total-math.log(len(data),2)*55
main()
