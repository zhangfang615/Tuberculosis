eventlist_file = file("./output_new/eventlist0.26_100_0.01_0.0018_0.05_0.75.txt")
line = eventlist_file.readline().strip()
cluster_lists = []
resistant_set = set()
cluster_dic = {}
# cluster_lists.append(list("a"))
# cluster_lists.append(list("b"))
# cluster_lists[0].append("c")
i = 0
while line:
    fields = line.split(" ")
    if fields[0] == "1" or fields[0] == "3" or fields[0] == "2":
        transmit = str(fields[1])
        transmitted = str(fields[2])
        if not transmit in cluster_dic.keys():
            new = set()
            new.add(transmit)
            cluster_lists.append(new)
            cluster_dic[transmit] = len(cluster_lists)-1
        c_list = cluster_dic[transmit]
        if not transmitted in cluster_dic.keys():
            cluster_lists[c_list].add(transmitted)
            cluster_dic[transmitted] = c_list
        else:
            d_list = cluster_dic[transmitted]
            if not c_list == d_list:
                cluster_lists[d_list].remove(transmitted)
                cluster_lists[c_list].add(transmitted)
                cluster_dic[transmitted] = c_list
    elif fields[0] == "5":
        removed = str(fields[1])
        if removed in cluster_dic.keys():
            cluster_lists[cluster_dic[removed]].remove(removed)
            cluster_dic.pop(removed)
    elif fields[0] == '4':
        resistant = str(fields[1])
        if resistant in cluster_dic.keys():
            if cluster_dic[resistant] not in resistant_set:
                new = set()
                cluster_lists[cluster_dic[resistant]].remove(resistant)
                new.add(resistant)
                cluster_lists.append(new)
                cluster_dic[resistant] = len(cluster_lists)-1
                resistant_set.add(len(cluster_lists)-1)
        else:
            new = set()
            new.add(resistant)
            cluster_lists.append(new)
            cluster_dic[resistant] = len(cluster_lists) - 1
            resistant_set.add(len(cluster_lists) - 1)
    line = eventlist_file.readline().strip()
    i += 1
    if i % 10000 == 0:
        print i

output_unresi = file("./output_new/0.26_100_0.01_0.018_0.05_0.75_unresi.txt","w")
output_resi = file("./output_new/0.26_100_0.01_0.018_0.05_0.75_resi.txt", "w")
i = 0
while i < len(cluster_lists):
    if i in resistant_set and not len(cluster_lists[i]) == 0 :
        resi = "".join(s+" " for s in cluster_lists[i])
        output_resi.writelines(resi+"\n")
    elif i not in resistant_set and not len(cluster_lists[i]) == 0 :
        unresi = "".join(s + " " for s in cluster_lists[i])
        output_unresi.writelines(unresi + "\n")
    i += 1
