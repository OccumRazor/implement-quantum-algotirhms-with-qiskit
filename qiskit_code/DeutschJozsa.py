from qiskit import QuantumRegister,QuantumCircuit
from qiskit.aqua.operators import StateFn
from qiskit.aqua.operators import I
from qiskit_code.quantumMethod import add,ini
from qiskit_code.classicalMethod import Dec2Bi

def DeutschJozsa(l,method):
# Deutsch, D. and Jozsa, R., 1992. Rapid solution of problems by quantum computation.
# Proceedings of the Royal Society of London. Series A: Mathematical and Physical Sciences,
# 439(1907), pp.553-558.
    # The input 'l' is the equivalent to the 'N' in the original paper of 
    # David Deutsch and Richard Jozsa, and 'method' denotes the 'unknown'
    # function, if you input 'balanced' then it will be balanced and otherwise
    # it will be constant.
    qr0=QuantumRegister(l)
    qr1=QuantumRegister(l+1)
    # One qubit larger to carry.
    ac=QuantumRegister(l) # Ancilla.
    t0=QuantumRegister(1)
    circ=QuantumCircuit(qr0,qr1,ac,t0)
    circ.h(qr0)
    if method=='balanced':
        print('balanced oracle')
        ini(circ,qr1,Dec2Bi(2**(l-1)))
    else:
        print('constant oracle')
        ini(circ,qr1,Dec2Bi(0))
    lst=range(l)
    QIN1=[qr0[i] for i in lst]+[qr1[i] for i in range(l+1)]+[ac[i] for i in lst]
    ADD=add(qr0,qr1,ac,l)
    circ.append(ADD,QIN1)# Role of the U unitary
    circ.cx(qr1[l],t0)# Role of the U unitary
    circ.z(t0)# The S unitary.
    circ.cx(qr1[l],t0)# Role of the U unitary
    circ.append(ADD.inverse(),QIN1)# Role of the U unitary
    psi=StateFn(circ)
    phiReg0=QuantumRegister(l)
    phiReg1=QuantumRegister(l+1)
    phiReg2=QuantumRegister(l)
    t1=QuantumRegister(1)
    phiCirc=QuantumCircuit(phiReg0,phiReg1,phiReg2,t1)  
    phiCirc.h(phiReg0)
    if method=='balanced':
        ini(circ,qr1,Dec2Bi(2**(l-1)))
    else:
        ini(circ,qr1,Dec2Bi(0))
    phi=StateFn(phiCirc)
    operator=I.tensorpower(3*l+2)
    expectation_value=(~psi@operator@phi).eval()
    print(expectation_value)
#DeutschJozsa('constant')
#DeutschJozsa('balanced')
