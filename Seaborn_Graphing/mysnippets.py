
def GenerateLists(stringArray1, stringArray2, label1num, label2repeatnum):
        ## Generate Tuples for the zip function to make dynamic analysis for specific things easier
        
        mylist1 = []
        mylist2 = []
        for string in range(0, len(stringArray1)-1):
            for iter in range(0,label1num-1):
                mylist1.append(string)
        for iter2 in range(0, (len(mylist1)-1)/ label2repeatnum):
            for iter2 in range(0, label2repeatnum-1):
                test = test