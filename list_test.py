def create_list(n):
    network_list = []
    net = n
    for i in range(6):
        for j in range (6-i):
            curnet = net[0:j]+ (j)*[1] + net[j:]
            network_list.append(curnet)
    return network_list
        
print(create_list([0,0,0,0,0,0]))
