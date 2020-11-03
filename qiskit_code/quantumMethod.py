from qiskit import *
from qiskit.circuit.library.standard_gates import SwapGate,CU1Gate,XGate,U1Gate
from math import pi,sqrt
from qiskit.quantum_info.operators import Operator
import numpy as np

def ini(circ,qr,ipt):
    # Input binary form, and append [0] ahead for qr1 block.
    for i in range(len(ipt)):
        if ipt[len(ipt)-i-1]:
            circ.x(qr[i])
    return 0
def diffusion(n): # Matrix representation of the diffusion transformation.
    N=2**n# for Grover search.
    return [N*[2/N] for i in range(N)]-np.identity(N)
def phaseFlip(reg,theta):
    # reg is a one qubit register (input a qubit is also OK).
    if reg.__len__()!=1:
        raise TypeError('The input quantum register contains more than one qubit.')
    phaseCirc=QuantumCircuit(reg,name='p\nh\na\ns\ne\n')
    phaseCirc.append(U1Gate(theta),reg)
    phaseCirc.append(XGate(),reg)
    phaseCirc.append(U1Gate(theta),reg)
    phaseCirc.append(XGate(),reg)
    return phaseCirc.to_instruction()
def CPhaseFlip(qReg,reg,theta):# If all ancilla qubits equal one, flip the phase of querry REG.
    # In this place, reg is a one qubit quantum register, and qReg is a n qubits quantum register.
    phaseCirc=QuantumCircuit(qReg,reg,name='p\nh\na\ns\ne\n')
    num=qReg.__len__()
    IN=[qReg[i] for i in range(num)]+[reg[0]]
    CU1Gate=U1Gate(theta).control(num)
    CXGate=XGate().control(num)
    phaseCirc.append(CU1Gate,IN)
    phaseCirc.append(CXGate,IN)
    phaseCirc.append(CU1Gate,IN)
    phaseCirc.append(CXGate,IN)
    return phaseCirc.to_instruction()
def amulitpdeAmplification(query,criteria,ancilla,n):# for Grover search.
    AACirc=QuantumCircuit(query,criteria,ancilla,name='A\nA\n')
    AACirc.h(query)
    from qiskit_code.Grover import check
    CHECK=check(query,criteria,ancilla,n)
    lst=range(n)# This looks rather awkward, I may (but not likely) try to change this later.
    AC=[ancilla[i] for i in lst]
    QUERY=[query[i] for i in lst]
    CHECKIN=QUERY+[criteria[i] for i in lst]+AC
    nCP=CPhaseFlip(ancilla,QuantumRegister(1),pi)
    nCX=XGate().control(n)# Generate a controlled-X gate with n controls.
    D=Operator(diffusion(n))
    for i in range(int(pi*sqrt(2**n)/8+1)):#[1] Iteration times.
        AACirc.append(CHECK,CHECKIN)
        AACirc.append(nCP,AC+[query[0]])
        AACirc.append(CHECK,CHECKIN)
        AACirc.append(D,query)
# [1]M. Boyer, G. Brassard, P. Hoyer & A. Tapp, Tight bounds on quantum searching,
#Proceedings, PhysComp 1996
    return AACirc.to_instruction()
def bigCSwap(c,t0,t1,l):# Controlled Swap of two qubit blocks.
    BCSC=QuantumCircuit(c,t0,t1,name='b\ni\ng\nC\nS\nw\na\np\n')
    for i in range(l):
        BCSC.cswap(c,t0[i],t1[i])
    return BCSC.to_instruction()

def bigSwap(t0,t1,l):# Swap of two qubit blocks.
    BSC=QuantumCircuit(t0,t1,name='b\ni\ng\nS\nw\na\np\n')
    for i in range(l):
        BSC.swap(t0[i],t1[i])
    return BSC.to_instruction()

def c_cx(c0,c1,target,n):# CCX where the second control and the target are two blocks
    c_cxC=QuantumCircuit(c0,c1,target,name='C\no\nn\nt\nr\no\nl\nl\ne\nd\n-\nC\nX\n')
    for i in range(n):
        c_cxC.ccx(c0,c1[i],target[i])
    return c_cxC.to_instruction()

def bigCCSwap(c0,c1,reg1,reg2,l):# Controlled-Controlled Swap of two qubit blocks.
    bigCCSwapC=QuantumCircuit(c0,c1,reg1,reg2,name='C\nC\nS\nw\na\np\n')
    ccswap=SwapGate().control(2)
    for i in range(l):
        bigCCSwapC.append(ccswap,[c0[0],c1[0],reg1[i],reg2[i]])
    return bigCCSwapC.to_instruction()

carry_q=QuantumRegister(4)
carryC=QuantumCircuit(carry_q,name='c\na\nr\nr\ny\n')
carryC.ccx(carry_q[1],carry_q[2],carry_q[3])
carryC.cx(carry_q[1],carry_q[2])
carryC.ccx(carry_q[0],carry_q[2],carry_q[3])
CARRY=carryC.to_instruction()

sum_q=QuantumRegister(3)
sumC=QuantumCircuit(sum_q,name='s\nu\nm')
sumC.cx(sum_q[1],sum_q[2])
sumC.cx(sum_q[0],sum_q[2])
SUM=sumC.to_instruction()

def add(q0,q1,q2,l):
# A quantum plain adder, as the main part of the oracle.
# Vedral, V., Barenco, A. and Ekert, A., 1996.
# Quantum networks for elementary arithmetic operations. Physical Review A, 54(1), p.147.
    add_circ=QuantumCircuit(q0,q1,q2,name='a\nd\nd')
    for i in range(l-1):
        add_circ.append(CARRY,[q2[i],q0[i],q1[i],q2[i+1]])
    add_circ.append(CARRY,[q2[l-1],q0[l-1],q1[l-1],q1[l]])
    add_circ.cx(q0[l-1],q1[l-1])
    add_circ.append(SUM,[q2[l-1],q0[l-1],q1[l-1]])
    RCARRY=CARRY.reverse_ops()#inverse()
    for i in range(l-2,-1,-1):
        add_circ.append(RCARRY,[q2[i],q0[i],q1[i],q2[i+1]])
        add_circ.append(SUM,[q2[i],q0[i],q1[i]])
    return add_circ.to_instruction()

def sub(q0,q1,q2,l):
    RCARRY=CARRY.reverse_ops()
    sub_circ=QuantumCircuit(q0,q1,q2,name='s\nu\nb')
    for i in range(l):
        sub_circ.append(SUM,[q2[i],q1[i],q0[i]])
        if i==l-1:
            sub_circ.cx(q0[i],q1[i])
        sub_circ.append(CARRY,[q2[i],q1[i],q0[i],q2[i+1]])
    for i in range(l-2,-1,-1):
        sub_circ.append(RCARRY,[q2[i],q1[i],q0[i],q2[i+1]])
    sub_circ.x(q2[l])
    sub_circ.swap(q0,q1)
    return sub_circ.to_instruction()

def adderMod(qr0,qr1,ac,Nr,swap_ac,t,l,ADD,SUB):
    # 0<=a,b<N
    AMC=QuantumCircuit(qr0,qr1,ac,Nr,swap_ac,t,name='a\nd\nd\ne\nr\nM\no\nd\n')
    BigCSwap=bigCSwap(t,qr0,swap_ac,l)
    BigSwap=bigSwap(qr0,Nr,l)
    lst=range(l)
    ADDIN=[qr0[i] for i in lst]+[qr1[i] for i in lst]+[ac[i] for i in lst]
    BigSwapIN=[qr0[i] for i in lst]+[Nr[i] for i in lst]
    BigCSwapIN=[t[0]]+[qr0[i] for i in lst]+[swap_ac[i] for i in lst]
    AMC.append(ADD,ADDIN)
    AMC.append(BigSwap,BigSwapIN)
    AMC.append(SUB,ADDIN)
    AMC.x(qr1[l-1])
    AMC.cx(qr1[l-1],t)
    AMC.x(qr1[l-1])
    AMC.append(BigCSwap,BigCSwapIN)
    AMC.append(ADD,ADDIN)
    AMC.append(BigCSwap,BigCSwapIN)
    AMC.append(BigSwap,BigSwapIN)
    AMC.append(SUB,ADDIN)
    AMC.cx(qr1[l-1],t)
    AMC.append(ADD,ADDIN)
    return AMC.to_instruction()

def c_mtpMOD(circ,qr0,qr1,qr2,ac,Nr,swap_ac,t,cReg,xReg,l,n):
    ADD=add(qr0,qr1,ac,l)
    SUB=sub(qr0,qr1,ac,l)
    AddMOD=adderMod(qr0,qr1,ac,Nr,swap_ac,t,l,ADD,SUB)
    iAddMOD=AddMOD.reverse_ops()
    BigCCSwap=bigCCSwap(cReg,t,qr0,swap_ac,l)
    CCX=c_cx(cReg,xReg,qr1,n)
    lst=range(l)
    AddMODIN=[qr0[i] for i in lst]+[qr1[i] for i in lst]+[ac[i] for i in lst]
    AddMODIN+=[Nr[i] for i in lst]+[swap_ac[i] for i in lst]+[t[0]]
    for i in range(n):
        BigCCSwapIN=[cReg,xReg[i]]+[qr0[i] for i in lst]+[qr2[i] for i in lst]
        circ.append(BigCCSwap,BigCCSwapIN)
        circ.append(AddMOD,AddMODIN)
        circ.append(BigCCSwap,BigCCSwapIN)
    circ.x(cReg)
    CCXIN=[cReg[0]]+[xReg[i] for i in range(n)]+[qr1[i] for i in range(l)]
    circ.append(CCX,CCXIN)
    circ.x(cReg)
    return 0

def expMod(qr0,qr1,ac,circ,N,l):
    BigSwap=bigSwap(xReg,qr1,l)
    lst=range(l)
    BigSwapIN=[xReg[i] for i in lst]+[qr1[i] for i in lst]
    C_MtpMOD=c_mtpMOD(qr0,qr1,qr2,ac,Nr,swap_ac,t,cReg,xReg,l,n)
    iC_MtpMOD=C_MtpMOD.reverse_ops()
    for i in range(m):
        MtpMODIN
        circ.append(C_MtpMOD,MtpMODIN)
        circ.append(BigSwap,BigSwapIN)
        circ.append(C_MtpMOD,MtpMODIN)
    return None

def qft(qReg):
# Michael Nielsen and Isaac Chuang (2000). Quantum Computation and Quantum
# Information. Cambridge: Cambridge University Press. ISBN 0-521-63503-9.
# OCLC 174527496.      P219, section 5.1 The quantum Fourier transform
# https://qiskit.org/documentation/stubs/qiskit.circuit.library.QFT.html
    qft_circ=QuantumCircuit(qReg,name='Q\nF\nT\n')
    num=qReg.__len__()
    for i in range(num-1,-1,-1):
        qft_circ.h(qReg[i])
        for j in range(i):
            qft_circ.append(CU1Gate(pi/2**(i-j)),[qReg[i],qReg[j]])
    # Reverse the qubit order
    for i in range(int(num/2)):# int(0.5)=0, so odd/even does not matters
        qft_circ.swap(qReg[i],qReg[num-1-i])
    return qft_circ.to_instruction()

