from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit,execute,Aer
from qiskit.providers.aer import QasmSimulator
from qiskit.quantum_info.operators import Operator
from math import log,cos,sin,sqrt,pi,exp
import matplotlib.pyplot as plt
from numpy import kron,matmul,transpose,conjugate,zeros,trace,complex128,array,inf,linspace,abs
from time import time
from scipy.linalg import expm
from random import random
from scipy.fft import fft
from scipy import integrate
import xlsxwriter
from tqdm import tqdm
H=[[1/sqrt(2),1/sqrt(2)],[1/sqrt(2),-1/sqrt(2)]]
X=[[0,1],[1,0]]
ID=[[1,0],[0,1]]
DM0=[[1,0],[0,0]]
DM1=[[0,0],[0,1]]
Y=[[0,-1j],[1j,0]]
Z=[[1,0],[0,-1]]
CX=[[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]]
ket=[[[1],[0]],[[0],[1]],ID]
CH=[[1,0,0,0],[0,1/sqrt(2),0,1/sqrt(2)],[0,0,1,0],[0,1/sqrt(2),0,-1/sqrt(2)]]

def RX(t):
    return [[cos(t),-1j*sin(t)],[-1j*sin(t),cos(t)]]
def twoQubitState(t):# Parametrized two qubit state
    # twoQubitState(t)|00>=cost|00>-isint|11>
    #CH=kron(DM0,ID)+kron(DM1,H)
    #return [[cos(t)],[0],[0],[-1j*sin(t)]]
    #return matmul(CH,kron(RX(t),ID))
    return matmul(CX,kron(RX(t),ID))

def WState(t):# Parametrized three qubit W State
    # cos(theta)|100>+sin(theta)/sqrt(2)|010>+sin(theta)/sqrt(2)|001>
    XNCNC=kron(ID,kron(DM1,DM1))+kron(ID,kron(DM0,DM1))+kron(ID,kron(DM1,DM0))+kron(X,kron(DM0,DM0))
    U=matmul(kron(kron(DM0,ID)+kron(DM1,Y),ID),kron(RX(t),kron(ID,ID)))
    U=matmul(kron(ID,kron(DM0,ID)+kron(DM1,H)),U)
    U=matmul(kron(kron(ID,DM0)+kron(X,DM1),ID),U)
    U=matmul(kron(ID,kron(ID,DM0)+kron(X,DM1)),U)
    U=matmul(XNCNC,U)
    return U
def M(p,q):# function for nQubitWState
    return [[sqrt(p)/sqrt(p+q),sqrt(q)/sqrt(p+q)],[-sqrt(q)/sqrt(p+q),sqrt(p)/sqrt(p+q)]]
def spaceExpansion(mat,n,idx):# function for nQubitWState
    if idx[0]>0:
        for _ in range(idx[0]):
            mat=kron(ID,mat)
    if idx[-1]<n-1:
        for _ in range(idx[-1]+1,n):
            mat=kron(mat,ID)
    return mat
def nQubitWState(n):
    m=M(1,n-1)
    m=spaceExpansion(m,n,[0])
    for i in range(1,n-1):
        mat=kron(DM0,ID)+kron(DM1,M(1,n-1-i))
        mat=spaceExpansion(mat,n,[i-1,i])
        m=matmul(mat,m)
    for i in range(n-1,0,-1):
        mat=spaceExpansion(CX,n,[i-1,i])
        m=matmul(mat,m)
    m=matmul(spaceExpansion(X,n,[0]),m)
    return m
def nQubitGHZState(n):# of course n must be at least 2, otherwise what do you want?
    U=matmul(CX,kron(H,ID))
    newCX=kron(ID,CX)
    while n>2:
        U=kron(U,ID)
        U=matmul(newCX,U)
        newCX=kron(ID,newCX)
        n-=1
    return U

def D2B(n,N):
    # n current number
    # N maximum
    # to make sure the output have the same length
    l=int(log(N,2))
    opt=[]
    for i in range(l):
        if n <2**(l-i-1):opt.append(0)
        else:
            opt.append(1)
            n-=2**(l-i-1)
    return opt
def B2D(li):
    l=len(li)
    res=0
    for i in range(l):res+=li[i]*2**(l-i-1)
    return res
def D2Q(num,length):# convert decimal number to quaternary number
    res=[0 for _ in range(length)]
    for i in range(length-1,-1,-1):
        if num>=4**i:
            res[length-i-1]=int((num-num%(4**i))/4**i)
            num=num%(4**i)
    return res
def mulKron(ops):
    res=ops[0]
    for i in range(1,len(ops)):res=kron(res,ops[i])
    return res
def listStringSum(li):
    text=li[0]
    for i in range(1,len(li)):text+=li[i]
    return text
def PauliDecomposition(Ham):
    n=int(log(len(Ham),2))
    set=[ID,X,Y,Z]
    dic={0:'I',1:'X',2:'Y',3:'Z'}
    nQubitPauli=[]
    remaindOps=[]
    coes=[]
    for i in range(4**n):
        idx=D2Q(i,n)
        ops=mulKron([set[idx[i]] for i in range(n)])
        temp=trace(matmul(ops,Ham))/(2**n)
        if temp!=0:
            coes.append(temp)
            nQubitPauli.append([dic[idx[i]] for i in range(n)])
            remaindOps.append(ops)
    print(f'number of nonzero term(s): {len(coes)}')
    for i in range(len(coes)):print(f'operator: {listStringSum(nQubitPauli[i])}, coefficient: {coes[i]}')

def quad(li):# hamiltonian companian function to avoid shallow copy
    res=[]# there are specific deep copy, from copy import deepcopy
    for j in range(len(li)):
        for _ in range(4):res.append([li[j][0],li[j][1]])
    return res
def hamiltonian(num_qubit,idle):
    # this Hamiltonian calculates Tr(rho^2_gamma) 
    # the variable num_qubit is actually half the number of qubit required
    loc=idle
    if isinstance(idle,int):loc=[idle]
    dim=2**(2*num_qubit)
    Har=[[0 for _ in range(dim)] for _ in range(dim)]
    # the list jk denotes the jk part of |ijk ijk> initially,
    # and finally it turns to |ijk i'jk>
    jk=[[D2B(i,2**(num_qubit-len(loc)))*2,D2B(j,2**(num_qubit-len(loc)))*2] for i in range(2**(num_qubit-1)) for j in range(2**(num_qubit-1))]
    for i in range(len(loc)):
        jk=quad(jk)
        for j in range(len(jk)):
            c=j%4
            jk[j][0]=jk[j][0][:loc[i]]+[int(c/2)]+jk[j][0][loc[i]:loc[i]+num_qubit-1]+[c%2]+jk[j][0][loc[i]+num_qubit-1:]
            jk[j][1]=jk[j][1][:loc[i]]+[int(c/2)]+jk[j][1][loc[i]:loc[i]+num_qubit-1]+[c%2]+jk[j][1][loc[i]+num_qubit-1:]
    for i in range(len(jk)):Har[B2D(jk[i][0])][B2D(jk[i][1])]=1
    return Har

def densityMatrix(state):
    dm=[]
    for i in range(len(state)):
        dm.append([])
        for j in range(len(state)):
            dm[i].append(complex(state[i]*conjugate(state[j])))
    return dm
def partialTrace(state,idle):
    # idle here refers to those traced qubits.
    n=int(log(len(state),2))
    traceList=[i for i in range(n) if i not in idle]
    li=[]
    for i in range(n):
        if i in traceList:li.append(0)
        else:li.append(1)
    N=2**sum(li) # number of basis
    rho=densityMatrix(state)
    res=complex128(zeros([2**len(traceList),2**len(traceList)]))
    for i in range(N):
        li=D2B(i,N)
        for j in traceList:
            li.insert(j,2)
        basis=ket[li[-1]]
        for k in range(n-2,-1,-1):
            basis=kron(ket[li[k]],basis)
        res+=matmul(transpose(basis),matmul(rho,basis))   
    return res

def swapTestRes(ipt,rep):
    trrho2=2*ipt['0']/rep-1 # the try sentence isn't required here, p0>=0.5
    return sqrt(2*(1-trrho2))

def swapTest(u0,idle,rep=10**6):
    n=int(log(len(u0),2))
    u0=Operator(u0)
    ctrl=QuantumRegister(1)
    d0Reg=QuantumRegister(n)
    d1Reg=QuantumRegister(n)
    c2Reg=ClassicalRegister(1)
    qc=QuantumCircuit(ctrl,d0Reg,d1Reg,c2Reg)
    qc.append(u0,[d0Reg[i] for i in range(n)])
    qc.append(u0,[d1Reg[i] for i in range(n)])
    qc.h(ctrl)
    for i in range(n):
        if i not in idle:
            qc.cswap(ctrl,d0Reg[i],d1Reg[i])
    qc.h(ctrl)
    qc.barrier()
    qc.measure(ctrl,c2Reg)
    simulator=Aer.get_backend('qasm_simulator')
    res=execute(qc,simulator,shots=rep).result().get_counts()
    return swapTestRes(res,rep)

def WStateConcurrence(N,k):
    # the GME concurrence has a simple formula,
    # it only depends on the number of qubits N
    # and the number of idle qubits k
    return sqrt(2*(1-(N**2-2*N*k+2*k**2)/N/N))

def completes(ele,otherPart,N):# to check whether a term should be added to idle
    for i in range(len(otherPart)):
        if len(ele)+len(otherPart[i])==N:
            temp=ele+otherPart[i]
            count=0
            for j in range(N):
                if j in temp:count+=1
            if count>=N:return False
        if len(ele)==len(otherPart[i]):
            count1=0
            for j in range(len(ele)):
                if ele[j] in otherPart[i]:count1+=1
            if count1==len(ele):return False
    return True
def genPartition(N):# for n qubits, 2^(n-1)-1 terms for n>2
    n=int(N/2)
    li=[[i] for i in range(N)]
    idles=[[i] for i in range(N)]
    for i in range(1,n+1):
        temp=[[j] for j in range(N)]
        for j in range(i-1):
            for k in range(len(temp)):
                if len(temp[k])==j+1:
                    for l in range(N):
                        if li[l][0] not in temp[k]:temp.append(temp[k]+li[l])
        for j in range(N,len(temp)):
            if completes(temp[j],idles,N) and temp[j] not in idles: idles.append(temp[j])
    return idles

def expectation(state,ham,time,rep):
    eiht=expm(time*1j*array(ham))
    dm=densityMatrix(kron(state,state))
    return trace(matmul(eiht,dm))*(1+(1-2*random())/sqrt(rep))
def h(x):
    if x**2>=1:return 0
    else:return exp(-1/(1-x**2))
def window(x):
    if x>=0 and x<1:return 1
    else:return 0
def fUnint(xp):
    return 5*h(2*xp-point)*window(xp)#/4.05
def F(x,num_sample):
    coes=[]
    xs=[x-2*(num_sample-m)/num_sample for m in range(num_sample)]+[x+2*m/num_sample for m in range(num_sample)]
    #xs=[-x/2+x*i/num_sample for i in range(len(num_sample))]
    #print(f'xs:{xs}')
    global point
    for i in range(len(xs)):
        point=xs[i]
        coes.append(integrate.quad(fUnint,-inf,inf)[0])
    coes=fft(coes)
    return coes
def timeSeries(theta,num_qubit=3,num_sample=50,rep=10**5):
    exp_TS=[]
    coes=F(0.75,num_sample)
    for i in range(len(theta)):
        if num_qubit==3:unitary=WState(theta[i])
        else:unitary=twoQubitState(theta[i])
        state=matmul(unitary,[[1]]+[[0] for _ in range(2**num_qubit-1)])
        exp=[]
        for j in range(-num_sample,num_sample):
            exp.append(coes[j]*expectation(state,array(hamiltonian(num_qubit,0)),j,rep))
        exp_TS.append(sum(exp))
    #norm=exp_TS[0]
    norm=max(exp_TS)
    for i in range(len(exp_TS)):exp_TS[i]/=norm#sqrt(2*(1-exp_TS[i]/norm))
    return exp_TS

def varyingTheta(num_points=100,num_qubit=3,rep=10**4):
    #theta=[i*pi/(num_points-1) for i in range(num_points-1)]+[pi]
    theta=linspace(0,pi/2,num_points)
    idle=genPartition(num_qubit)
    GME_Classical=[]
    GME_Quantum=[]
    for i in range(num_points):
        if num_qubit==3:unitary=WState(theta[i])
        else:unitary=twoQubitState(theta[i])
        state=matmul(unitary,[[1]]+[[0] for _ in range(2**num_qubit-1)])
        #print(state)
        iGME=[]
        qGME=[]
        for j in range(len(idle)):
            pdm=partialTrace(state,idle[j])
            #print(pdm)
            temp=2*(1-trace(matmul(pdm,pdm)))
            #print(temp)
            #temp=trace(matmul(pdm,pdm))
            if abs(temp.imag)>1e-7:print('imaginary part might exist')
            iGME.append(sqrt(abs(temp)))
            swap_test_res=swapTest(unitary,idle[j],rep)
            qGME.append(swap_test_res)
        GME_Classical.append(min(iGME))
        GME_Quantum.append(min(qGME))
        #GME_Har.append(min(hGME))
    x=[theta[i]/pi for i in range(num_points)]
    plt.plot(x,GME_Classical,'r',label='partial trace')
    plt.plot(x,GME_Quantum,'b-.',label='swap test')
    GME_TS=timeSeries(theta,num_qubit,10**2,rep)
    for i in range(len(GME_TS)):GME_TS[i]=sqrt(2*(1-GME_TS[i]))
    plt.plot(x,GME_TS,'g*',label='time series')
    plt.xlabel(r'$\theta/\pi$')
    plt.ylabel('GME Concurrence')
    #plt.ylabel('purity')
    plt.legend(loc='best')
    plt.savefig('with time series.png')
    plt.savefig('with time series.eps',format='eps')
    plt.show()

varyingTheta(101,2,10**4)

def varyingQubit(rep=10**4):
    start=2
    end=5
    GME_Classical=[]
    GME_Quantum=[]
    for num_qubit in range(start,end):
        print(f'{num_qubit}-qubit GME concurrence')
        idle=genPartition(num_qubit)
        unitary=nQubitWState(num_qubit)
        state=matmul(unitary,[[1]]+[[0] for _ in range(2**num_qubit-1)])
        iGME=[]
        qGME=[]
        t0=0
        print(len(idle))
        for j in tqdm(range(len(idle))):
            pdm=partialTrace(state,idle[j])
            temp=2*(1-trace(matmul(pdm,pdm)))
            if abs(temp.imag)>1e-5:
                print('imag part goes high')
            iGME.append(sqrt(temp.real))
            t1=time()
            qGME.append(swapTest(unitary,idle[j],rep))
            t2=time()
            t0+=t2-t1
        print(f'{num_qubit}-qubit GME concurrence requires time: {t0}')
        GME_Classical.append(min(iGME))
        GME_Quantum.append(min(qGME))
    x=[i for i in range(start,end)]
    plt.plot(x,GME_Classical,'r',label='partial trace')
    plt.plot(x,GME_Quantum,'b-.',label='swap test')
    plt.xlabel('number of qubits')
    plt.ylabel('GME Concurrence')
    plt.legend(loc='best')
    plt.show()
    #plt.savefig('n qubit w state.png')
    #plt.savefig('n qubit w state.eps',format='eps')

#varyingQubit(10**4)

def unnecessaryBipartitions(rep=10**4):
    # for W state and GHZ state, they are highly symmetrical, it is reasonable to 
    # assume that a lot of bipartitions give the same result.
    start=2
    end=13
    y1=[]
    y2=[]
    for num_qubit in range(start,end):
        idle=genPartition(num_qubit)
        unitary=nQubitWState(num_qubit)
        state=matmul(unitary,[[1]]+[[0] for _ in range(2**num_qubit-1)])
        for j in range(len(idle)-1):
            if len(idle[j])!=len(idle[j+1]):
                print(f'idle qubits:{idle[j]}')
                pdm=partialTrace(state,idle[j])
                temp=2*(1-trace(matmul(pdm,pdm)))
                y1.append(temp)
                if abs(temp.imag)>1e-5:
                    print('imag part goes high')
                print(f'partial trace result: {sqrt(temp.real)}')
                print(f'formula result: {WStateConcurrence(num_qubit,len(idle[j]))}')
                y2.append(WStateConcurrence(num_qubit,len(idle[j])))
    x=[i for i in range(len(y1))]
    plt.plot(x,y1,'r')
    plt.plot(x,y2,'bo')
    plt.show()

#unnecessaryBipartitions()

def varyingQubitWithFewerBipartition(rep=10**4):
    start=2
    end=16
    GME_Classical=[]
    GME_Quantum=[]
    idle=[0]
    t0=time()
    for num_qubit in range(start,end):
        unitary=nQubitWState(num_qubit)
        state=matmul(unitary,[[1]]+[[0] for _ in range(2**num_qubit-1)])
        pdm=partialTrace(state,idle)
        temp=2*(1-trace(matmul(pdm,pdm)))
        if abs(temp.imag)>1e-5:
            print('imag part goes high')
        GME_Classical.append(sqrt(temp.real))
        print(f'current precise result: {GME_Classical[-1]}')
        t2=time()
        GME_Quantum.append(swapTest(unitary,idle,rep))
        t3=time()
        print(f'current C-Swap result: {GME_Quantum[-1]}')
        print(f'time spent for {2*num_qubit} swap test: {t3-t2}s')
    t1=time()
    print(f'time spent: {t1-t0}s')
    print(f'current precise result: {GME_Classical}')
    print(f'current approximate result: {GME_Quantum}')
    x=[i for i in range(start,end)]
    plt.plot(x,GME_Classical,'r',label='partial trace')
    plt.plot(x,GME_Quantum,'b-.',label='swap test')
    plt.xlabel('number of qubits')
    plt.ylabel('GME Concurrence')
    plt.legend(loc='best')
    plt.savefig('n qubit w state with one bipartition.png')
    plt.savefig('n qubit w state with one bipartition.eps',format='eps')
    plt.show()

#varyingQubitWithFewerBipartition(10**5)
