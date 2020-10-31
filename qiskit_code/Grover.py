from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit,Aer,execute
from qiskit.circuit.library.standard_gates import XGate,PhaseGate,U1Gate
from qiskit.providers.aer import QasmSimulator
import matplotlib.pyplot as plt
from math import pi,sqrt
import numpy as np
from qiskit_code.quantumMethod import ini,diffusion,amulitpdeAmplification
from qiskit_code.classicalMethod import Dec2Bi,printCounts

def check(query,criteria,ancilla,n): # Check entangled state of two reg are in the
    # same state quantum mechanically.
    # This plays the role of the oracle, although this is rather trivial.
    checkCirc=QuantumCircuit(query,criteria,ancilla,name='c\nh\ne\nc\nk\n')
    checkCirc.ccx(query,criteria,ancilla)
    checkCirc.x(query)
    checkCirc.x(criteria)
    checkCirc.ccx(query,criteria,ancilla)
    checkCirc.x(query)
    checkCirc.x(criteria)
    return checkCirc.to_instruction()
def phaseFlip(qReg,reg,n):# If all ancilla qubits equal one, flip the phase of querry REG.
    phaseCirc=QuantumCircuit(qReg,reg,name='p\nh\na\ns\ne\n')
    IN=[qReg[i] for i in range(n)]+[reg[0]]
    CU1Gate=U1Gate(pi).control(n)
    CXGate=XGate().control(n)
    phaseCirc.append(CU1Gate,IN)
    phaseCirc.append(CXGate,IN)
    phaseCirc.append(CU1Gate,IN)
    phaseCirc.append(CXGate,IN)
    return phaseCirc.to_instruction()

def Grover_main(state):
# Grover, Lov K. "A fast quantum mechanical algorithm for database search."
# Proceedings of the twenty-eighth annual ACM symposium on Theory of computing. 1996.
    n=len(Dec2Bi(state))
    query=QuantumRegister(n)
    criteria=QuantumRegister(n)
    ancilla=QuantumRegister(n)
    cReg=ClassicalRegister(n)
    cReg2=ClassicalRegister(n)
    circ=QuantumCircuit(query,criteria,ancilla,cReg)
    ini(circ,criteria,Dec2Bi(state))
    lst=range(n)
    AC=[ancilla[i] for i in lst]
    QUERY=[query[i] for i in lst]
    IN=QUERY+[criteria[i] for i in lst]+AC
    AA=amulitpdeAmplification(query,criteria,ancilla,n)
    circ.append(AA,IN)
    circ.measure(query,cReg)
    simulator=Aer.get_backend('qasm_simulator')
    res=execute(circ,simulator).result().get_counts()
    printCounts(res,bin(state)[2:])

def draw(N,x,state):# The latter two function draws the process of amplitude amplification.
    plt.xlim(-1,N)
    plt.ylim(-1.2,1.2)
    plt.text(-0.8,1.1,'amplitude of each state')
    plt.bar(x,state)
    plt.pause(0.5)
    plt.clf()
def classicalGrover(j):
    # input 'j' is the number that you want to find out.
    n=6# Number of qubits
    N=2**n
    D=-np.identity(N)+2/N# Thediffusion matrix.
    state=[1/sqrt(N) for i in range(N)]# Initial state.
    x=np.arange(N)# To draw the plot.
    draw(N,x,state)
    for i in range(int(pi*sqrt(N)/8+1)):
        state[j]*=-1# Controlled phase flip, of course classical.
        draw(N,x,state)
        state=np.dot(D,state)
        draw(N,x,state)    
    plt.xlim(-1.2,N)
    plt.ylim(-1.2,1.2)
    plt.text(-0.8,1.1,'amplitude of each state')
    plt.text(-0.8,0.95,'k^2: '+str(state[j]*state[j]))
    plt.bar(x,state)
# In this program, the input of the Grover search sholud be a single
# integer. If its binary form has n bits, then the searching range will
# be [0,2^n-1].
