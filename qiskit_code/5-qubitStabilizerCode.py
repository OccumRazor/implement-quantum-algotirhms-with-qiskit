from qiskit import QuantumRegister,ClassicalRegister,QuantumCircuit,Aer,execute
from qiskit.providers.aer import QasmSimulator
from qiskit.circuit.library.standard_gates import CU1Gate
from numpy import pi
from qiskit_code.classicalMethod import printCounts,Dec2Bi,modifyExpValue
from qiskit_code.quantumMethod import ini
from qiskit.aqua.operators import StateFn,I
from numpy import real
def physicalQubits(ipt):
    qr=QuantumRegister(5)
    circ=QuantumCircuit(qr)
    if ipt==1:
        circ.x(qr[0])
    # controlled phase flip - if the input state is |1>,
    # then flip the global phase by pi
    CU1=CU1Gate(pi)
    circ.append(CU1,[qr[0],qr[1]])
    circ.cx(qr[0],qr[1])
    circ.append(CU1,[qr[0],qr[1]])
    circ.cx(qr[0],qr[1])
    
    circ.h(qr[4])
    circ.s(qr[4])
    # g1
    circ.cz(qr[4],qr[3])
    circ.cz(qr[4],qr[1])
    circ.cy(qr[4],qr[0])
    
    circ.h(qr[3])
    #g2
    circ.cz(qr[3],qr[2])
    circ.cz(qr[3],qr[1])
    circ.cx(qr[3],qr[0])

    circ.h(qr[2])
    #g3
    circ.cz(qr[2],qr[4])
    circ.cz(qr[2],qr[3])
    circ.cx(qr[2],qr[0])

    circ.h(qr[1])
    circ.s(qr[1])
    #g4
    circ.cz(qr[1],qr[4])
    circ.cz(qr[1],qr[2])
    circ.cy(qr[1],qr[0])
    return circ.to_gate()

def checkPhases():
    operator=I.tensorpower(5)
    for i in range(32):
        qr=QuantumRegister(5)
        circ=QuantumCircuit(qr)
        ini(circ,qr,Dec2Bi(i))
        # generate the state vector of the test state
        psi=StateFn(circ)
        # generate the state vector of the 0L state
        phi1=StateFn(physicalQubits(0))
        thisNumber=bin(i)[2:]
        if len(thisNumber)<5:
            thisNumber=(5-len(thisNumber))*'0'+thisNumber
        print('expectation value for state '+thisNumber)
        print(' and the physical qubits of 0:')
        exp01=(~phi1@operator@psi).eval()
        exp1=modifyExpValue(real(exp01))
        print(exp1)
        # generate the state vector of the 1L state
        phi2=StateFn(physicalQubits(1))
        print(' and the physical qubits of 1:')
        exp02=(~phi2@operator@psi).eval()
        exp2=modifyExpValue(real(exp02))
        print(exp2)

def stabilized(ipt):
    operator=I.tensorpower(5)
    qr=QuantumRegister(5)
    circ=QuantumCircuit(qr)
    circ.append(physicalQubits(ipt),qr)
    # stabilized:
    phi=StateFn(circ)
    # g1 check
    circ.cz(qr[4],qr[3])
    circ.cz(qr[4],qr[1])
    circ.cy(qr[4],qr[0])
    psi=StateFn(circ)
    print((~phi@operator@psi).eval())
    #g2 check
    circ.cz(qr[3],qr[2])
    circ.cz(qr[3],qr[1])
    circ.cx(qr[3],qr[0])
    psi=StateFn(circ)
    print((~phi@operator@psi).eval())
    #g3 check
    circ.cz(qr[2],qr[4])
    circ.cz(qr[2],qr[3])
    circ.cx(qr[2],qr[0])
    psi=StateFn(circ)
    print((~phi@operator@psi).eval())
    #g4 check
    circ.cz(qr[1],qr[4])
    circ.cz(qr[1],qr[2])
    circ.cy(qr[1],qr[0])
    psi=StateFn(circ)
    print((~phi@operator@psi).eval())

#checkPhase()
stabilized(0)
stabilized(1)

