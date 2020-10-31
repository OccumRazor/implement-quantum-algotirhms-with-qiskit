def Dec2Bi(num):
# Decimal to binary.
    res=list(bin(num)[2:])
    return [int(res[i]) for i in range(len(res))]
def Bi2Dec(p):
# Binary to decimal.
    len12=len(p)
    res=0
    for i in range(len12):
        k=p[i]*2**(len12-i-1)
        res=k+res
    return res
def returnN(a,b):
    if len(Dec2Bi(a))>len(Dec2Bi(b)):
        return b,a,len(Dec2Bi(a))
    else:
        return a,b,len(Dec2Bi(b))
def Reverse(lst):
    return [ele for ele in reversed(lst)]
def vectorForm(a):# In this place, input a binary form of var a.
    L=len(a)
    temp=[]
    for i in range(L):
        if a[i]==0:
            temp.append([1,0])
        else:
            temp.append([0,1])
    optSequence=temp[L-1]
    for i in range(L-2,-1,-1):
        mid=[]
        for k in range(2):
            if temp[i][k]==0:
                for k in range(len(optSequence)):
                    mid.append(0)
            else:
                mid+=optSequence
        optSequence=mid
    #circ.initialize(vectorForm(Dec2Bi(a)),qr0)
    #circ.initialize(vectorForm([0]+Dec2Bi(b)),qr1)
    # qiskit initialization instruction, for reference, not compatiable
    # with my code.
    return optSequence

def randomPick(data):
    # The data should be a dictionary, in this case, the qiskit counts.
    keys=list(data.keys())
    seed=[]
    for i in range(len(keys)):
        seed.extend(data[keys[i]]*[int(keys[i])])
    return choice(seed)

def printCounts(res,string=1):
    keys=list(res.keys())
    for i in range(len(keys)):
        print(keys[i]+':'+str(res[keys[i]]))
    if string!=1:
        s=sum([res[keys[i]] for i in range(len(keys))])
        try:
            index=keys.index(string)
            print('\n'+str(len(keys))+' different states.')
            print('\nExpected state:')
            print(keys[index]+':'+str(res[keys[index]])+'/'+str(s))
        except ValueError:
            print('\nExpected state does not exit!')

def iptCheck(ipt):
    if isinstance(ipt,(int,list)):
        try:
            num=len(ipt)
            for i in range(num):
                if not isinstance(ipt[i],int):
                   raise TypeError('all list elements should be int')
        except TypeError:
            num=1
    else:
        raise TypeError('the target state must be either a int or a list')
    return num

def strMaster(lst):
    # TODO: extract each string in the upper function and compare them to do
    # general jobs.
    return 0
